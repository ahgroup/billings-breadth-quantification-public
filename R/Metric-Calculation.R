###
# Metrics Calculation
# Zane
# 2025-05-03
# Helper functions for calculating the vaccine response metrics
###

# A function to calcualte the metrics for a given model fit and metadata list
# Since there's so many subsample models it's way more memory efficient to
# do it this way and map rather than bind all the model fits together and
# do a list operation.
calculate_metrics_for_model <- function(model_fit, metadata, h = 0.01) {
	# Model cases and metrics to get per case
	# Logistic: SCR
	# Linear: GMT and AUC

	# Need this because of how targets passes the metadata list
	metadata <- metadata[[1]]

	# First get the model epreds
	model_epreds <- tidybayes::epred_draws(
		model_fit,
		newdata = tibble::tibble(d_norm = seq(0, 1, h)),
		ndraws = brms::ndraws(model_fit),
		re_formula = NA
	) |>
		dplyr::rename(x = d_norm)

	# Create output holder
	metrics_out <- list()

	# Skip the homologous/full models that shouldn't exist, remove this
	# once its fixed
	if (metadata$model == "full" && metadata$dataset == "homologous") {
		return(list(
			season = metadata$season,
			metric = metadata$metric,
			dataset = metadata$dataset,
			model = metadata$model,
			censoring = metadata$censoring,
			stats = metrics_out
		))
	}

	# Dispatch the correct metric calculations based on the model metadata
	if (metadata$model == "logistic") {
		# Calculate the seroconversion rate
		metrics_out$scr <- mean_scr(model_epreds)
	} else if (metadata$model == "reduced") {
		# Calculate the GMT
		metrics_out$gmt <- intercept(model_epreds)
	} else if (metadata$model == "full") {
		# Calculate our three model metrics
		metrics_out$auc <- AUC(model_epreds)
		metrics_out$pat <- prop_above_threshold(model_epreds, "log_posttiter")
		metrics_out$int <- intercept(model_epreds)
	} else {
		cli::cli_abort(c(
		paste0(
			"{.var metadata$model} should be one of {.val reduced}, {.val full}, ",
			"or {.val logistic}."
			),
		"i" = "\nYou passed {.val {metadata$model}} instead."
		))
	}

	# If this is a homologous model, add the 0 to the end of the name
	if (metadata$dataset == "homologous") {
		names(metrics_out) <- paste0(names(metrics_out), "0")
	}

	# If the subset ID is not there make a fake one
	if (is.null(metadata$subsample_id)) {metadata$subsample_id <- "00"}

	# Now format the output to be bound to dataframe later
	out <- list(
		season = metadata$season,
		metric = metadata$metric,
		dataset = metadata$dataset,
		model = metadata$model,
		censoring = metadata$censoring,
		subsample_id = metadata$subsample_id,
		stats = list(metrics_out)
	)

	return(out)
}

bind_calculated_metrics <- function(full_data_metrics, subsample_metrics) {

	# Bind all the little lists together into one dataframe
	metrics_dirty <- dplyr::bind_rows(
		dplyr::bind_rows(full_data_metrics),
		dplyr::bind_rows(subsample_metrics)
	) |>
		dplyr::mutate(subsample_id = factor(subsample_id, ordered = TRUE))

	# Unnest the stuff into a rectangular format
	metrics_unnested <-
		metrics_dirty |>
		tidyr::unnest_longer(stats) |>
		tidyr::unnest(stats)

	# Fix the logistic censoring by binding stuff
	stats_not_logistic <- dplyr::filter(metrics_unnested, censoring != "logistic")
	stats_logistic <- dplyr::filter(metrics_unnested, censoring == "logistic")

	stats_logistic_fixed <- dplyr::bind_rows(
		stats_not_logistic,
		dplyr::mutate(stats_logistic, censoring = "yes"),
		dplyr::mutate(stats_logistic, censoring = "no")
	) |>
		dplyr::mutate(
			stat_norm = minmax(stats),
			censoring = factor(censoring, c("yes", "no")),
			stats_id = factor(
				stats_id,
				levels = c("gmt", "gmt0", "scr", "scr0", "int", "pat", "auc"),
				labels = c("GMT", "GMT0", "SCR", "SCR0", "INT", "PAT", "AUC")
			)
		) |>
		dplyr::arrange(season, metric, censoring, subsample_id, stats_id)

	return(stats_logistic_fixed)

}

summary_fun <- function(x, all = TRUE, ...) {
	if (isTRUE(all)) {
		stats <- tibble::tibble(
			mean = mean(x, ...),
			sd = sd(x, ...),
			cv = sd / mean,
			q1 = quantile(x, probs = c(0.25)),
			q3 = quantile(x, probs = c(0.75)),
			midhinge = mean(q1, q3),
			iqr = IQR(x),
			qcd = iqr / (2 * midhinge)
		)
		ci <- ggdist::mean_hdci(x)
		out <- dplyr::bind_cols(stats, ci)
		return(out)
	} else if (isFALSE(all)) {
		return(ggdist::mean_hdci(x))
	} else {
		rlang::abort("'all' should be TRUE or FALSE!")
	}
}

summarize_calculated_metrics <- function(bound_metrics_samples) {
	# Summarize over the posterior samples
	stats_summary <-
		bound_metrics_samples |>
		dplyr::summarize(
			summary_fun(stats),
			.by = c(season, metric, stats_id, subsample_id, censoring)
		)

	stats_summary_normalized <-
		bound_metrics_samples |>
		dplyr::summarize(
			summary_fun(stat_norm),
			.by = c(season, metric, stats_id, subsample_id, censoring)
		)

	out <- dplyr::bind_rows(
		"raw" = stats_summary,
		"normalized" = stats_summary_normalized,
		.id = "stat_group"
	)

	return(out)
}

# ICC analysis of the subsamples ####
# First we need to create a nested df so we can dispatch ICC models in parallel
nest_stats_dataframe <- function(
		bound_metrics_samples, random_seeds, global_seed) {
	vars_select <-
		bound_metrics_samples |>
		dplyr::select(
			season, metric, stats_id, subsample_id, censoring, y = stat_norm
		) |>
		dplyr::filter(subsample_id != "00") |>
		droplevels()

	# Make an overall group
	iccs_data_with_overall <- dplyr::bind_rows(
		vars_select,
		dplyr::mutate(vars_select, season = "Overall")
	)

	# And now nest the data
	nested_df <-
		iccs_data_with_overall |>
		tidyr::nest(icc_data = -c(season, metric, stats_id, censoring))

	set.seed(global_seed)
	nested_df$seed <- sample(random_seeds, nrow(nested_df), replace = FALSE)

	return(nested_df)
}

generate_icc_control_arguments <- function(is_test_run) {
	if(isTRUE(is_test_run)) {
		cmdstan_sampling_arguments <-
			sampling_arguments <- list(
				chains = 4L,
				iter_warmup = 20,
				iter_sampling = 20,
				max_treedepth = 10,
				adapt_delta = 0.8
			)
	} else if (isFALSE(is_test_run)) {
		cmdstan_sampling_arguments <-
			sampling_arguments <- list(
				chains = 4L,
				iter_warmup = 250,
				iter_sampling = 1000,
				max_treedepth = 10,
				adapt_delta = 0.95
			)
	} else {
		cli::cli_abort(c(
			"{.var is_test_run} must be TRUE or FALSE, not {.val {is_test_run}}."
		))
	}

	return(cmdstan_sampling_arguments)
}

fit_icc_model <- function(icc_data, fitting_seed, fitting_arguments) {
	fit <- brms::brm(
		formula = y ~ 1 + (1 | subsample_id),
		family = "gaussian",
		data = icc_data,
		prior = c(
			brms::prior(normal(0.5, 0.25), class = "Intercept"),
			brms::prior(student_t(3, 0, 0.25), class = "sd"),
			brms::prior(student_t(3, 0, 0.25), class = "sigma")
		),
		seed = fitting_seed,
		chains = fitting_arguments$chains,
		cores = fitting_arguments$chains,
		warmup = fitting_arguments$iter_warmup,
		iter = fitting_arguments$iter_warmup + fitting_arguments$iter_sampling,
		algorithm = "sampling",
		control = list(
			adapt_delta = fitting_arguments$adapt_delta,
			max_treedepth = fitting_arguments$max_treedepth
		),
		backend = "cmdstanr"
	)
}

calculate_icc_from_model <- function(icc_model) {
	post <- brms::as_draws_df(
		icc_model,
		variable = c("sigma", "sd_subsample_id__Intercept")
	)

	var_subsample <- post$sd_subsample_id__Intercept ^ 2
	var_err <- post$sigma ^ 2

	ICC <- var_subsample / (var_subsample + var_err)

	return(ICC)
}

bind_iccs <- function(nested_stats_df, icc_results) {
	icc_data_nested <- nested_stats_df |>
		dplyr::select(season, metric, stats_id, censoring) |>
		tibble::add_column(ICC_res = icc_results)

	icc_data <- tidyr::unnest(icc_data_nested, ICC_res)

	return(icc_data)
}


summarize_icc_results <- function(icc_results_df) {
	out <-
		icc_results_df |>
		dplyr::summarise(
			ggdist::mean_hdci(ICC_res),
			.by = c(season, metric, stats_id, censoring)
		)

	return(out)
}


