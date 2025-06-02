###
# Code for multivariate ICC model
# Zane
# 2025-05-12
# Andreas wants posterior comparisons of "goodness" for each of the ICCs
# across curent and new metrics, so we have to fit multivariate models.
###

# We need to reformate the nested stats DF into a form usable for the
# multivariate ICC model
reformat_nested_stats_df <- function(nested_stats_df, random_seeds, global_seed) {
	new_factors <- nested_stats_df |>
		dplyr::mutate(
			group = dplyr::case_match(
				stats_id,
				"GMT" ~ "old-tot",
				"GMT0" ~ "old-mag",
				"SCR" ~ "old-brd",
				"INT" ~ "new-mag",
				"PAT" ~ "new-brd",
				"AUC" ~ "new-tot"
			)
		) |>
		tidyr::separate(col = group, into = c("set", "att"), sep = "-") |>
		dplyr::filter(stats_id != "SCR0") |>
		dplyr::mutate(
			set = factor(
				set,
				levels = c("old", "new"),
				labels = c("Current", "Novel")
			),
			att = factor(
				att,
				levels = c("mag", "brd", "tot"),
				labels = c("Magnitude", "Breadth", "Total Strength")
			)
		)

	new_nested <-
		# Now we unnest the icc_data so we can reformat it
		new_factors |>
		dplyr::select(-stats_id, -seed) |>
		tidyr::unnest(icc_data) |>
		dplyr::filter(season != "Overall") |>
		droplevels() |>
		# Create a stupid id variable so that the pivot works correctly
		# and doesn't give useless list column
		# I don't know why I had to use this combo of .by variables i just
		# messed with it until it gave the right result
		dplyr::mutate(
			row_id = dplyr::row_number(),
			.by = c(season, metric, censoring, subsample_id, att, set)
		) |>
		# Pivot to have one column for current metric result and one column
		# for novel metric result
		tidyr::pivot_wider(
			names_from = set,
			values_from = y
		) |>
		# Remove the stupid id variable
		dplyr::select(-row_id) |>
		# And renest
		tidyr::nest(icc_data = c(subsample_id, Current, Novel))

	set.seed(global_seed)
	new_nested$fitting_seed <- sample(random_seeds, size = nrow(new_nested))

	return(new_nested)
}

# Function to fit the multivariate models
fit_mv_icc_model <- function(icc_data, fitting_seed, fitting_arguments) {
	require(brms, quietly = TRUE)

	fit <- brms::brm(
		formula = brms::bf(mvbind(Current, Novel) ~ 1 + (1 |p| subsample_id)) +
			set_rescor(TRUE),
		family = "gaussian",
		data = icc_data,
		prior = c(
			brms::set_prior(
				"normal(0.5, 0.25)",
				class = "Intercept", resp = c("Current", "Novel")
			),
			brms::set_prior(
				"normal(0.5, 0.25)",
				class = "Intercept"
			),
			brms::set_prior(
				"student_t(3, 0, 0.25)",
				class = "sd", resp = c("Current", "Novel")
			),
			brms::set_prior(
				"student_t(3, 0, 0.5)",
				class = "sigma",  resp = c("Current", "Novel")
			),
			brms::set_prior(
				"lkj(2)",
				class = "cor"
			),
			brms::set_prior(
				"lkj(2)",
				class = "rescor"
			)
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

# And functions to summarize the multivariate models
calculate_mv_iccs_from_model <- function(icc_model) {
	post <- icc_model |>
		brms::as_draws_df(
			variable = c(
				"sd_subsample_id__Current_Intercept",
				"sd_subsample_id__Novel_Intercept",
				"sigma_Current",
				"sigma_Novel"
			)
		) |>
		tibble::as_tibble() |>
		dplyr::transmute(
			var_re_Current = sd_subsample_id__Current_Intercept ^ 2,
			var_re_Novel = sd_subsample_id__Novel_Intercept ^ 2,
			var_sigma_Current = sigma_Current ^ 2,
			var_sigma_Novel = sigma_Novel ^ 2
		) |>
		dplyr::mutate(
			icc_Current = var_re_Current / (var_re_Current + var_sigma_Current),
			icc_Novel = var_re_Novel / (var_re_Novel + var_sigma_Novel),
			.keep = "unused"
		) |>
		dplyr::mutate(icc_diff = icc_Current - icc_Novel)

	return(post)
}

calculate_standardized_icc_difference <- function(mv_icc_result) {
	d <- atanh(mv_icc_result$icc_Current) - atanh(mv_icc_result$icc_Novel)
	#Z <- d * sqrt(nrow(mv_icc_result) - 3)

	return(d)
}

bind_mv_iccs <- function(nested_stats_df_multivariate, mv_icc_results) {
	icc_data_nested <- nested_stats_df_multivariate |>
		dplyr::select(season, metric, censoring, att) |>
		tibble::add_column(ICC_res = mv_icc_results)

	icc_data <- tidyr::unnest(icc_data_nested, ICC_res)

	return(icc_data)
}

# Finally functions to compute bayesian contrast measurements
# similar to p value or whatever
summarize_mv_iccs <- function(mv_icc_data) {
	icc_current_summary <- mv_icc_data |>
		dplyr::summarize(
			ggdist::mean_hdci(icc_Current),
			.by = c(season, metric, censoring, att)
		) |>
		dplyr::transmute(
			season, metric, censoring, att,
			est_curr = y,
			lwr_curr = ymin,
			upr_curr = ymax
		)

	icc_novel_summary <- mv_icc_data |>
		dplyr::summarize(
			ggdist::mean_hdci(icc_Novel),
			.by = c(season, metric, censoring, att)
		) |>
		dplyr::transmute(
			season, metric, censoring, att,
			est_novel = y,
			lwr_novel = ymin,
			upr_novel = ymax
		)

	icc_diff_summary <- mv_icc_data |>
		dplyr::summarize(
			ggdist::mean_hdci(icc_diff),
			.by = c(season, metric, censoring, att)
		) |>
		dplyr::transmute(
			season, metric, censoring, att,
			est_diff = y,
			lwr_diff = ymin,
			upr_diff = ymax
		)

	rope_test <- mv_icc_data |>
		dplyr::summarize(
			bayestestR::rope(icc_diff, range = c(-0.1, 0.1), ci_method = "ETI"),
			.by = c(season, metric, censoring, att)
		) |>
		dplyr::transmute(
			season, metric, censoring, att,
			p_in_rope = ROPE_Percentage,
			p_out_of_rope = 1 - ROPE_Percentage
		)

	pd_test <- mv_icc_data |>
		dplyr::summarize(
			#bayestestR::pd(icc_diff),
			#med = median(icc_diff),
			p_positive = mean(icc_diff > 0),
			.by = c(season, metric, censoring, att)
		)

	jb <- c("season", "metric", "censoring", "att")
	all_effects <- mv_icc_data |>
		dplyr::distinct(season, metric, censoring, att) |>
		dplyr::left_join(icc_current_summary, by = jb) |>
		dplyr::left_join(icc_novel_summary, by = jb) |>
		dplyr::left_join(icc_diff_summary, by = jb) |>
		dplyr::left_join(rope_test, by = jb) |>
		dplyr::left_join(pd_test, by = jb)

	return(all_effects)
}
