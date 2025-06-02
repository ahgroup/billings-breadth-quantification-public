###
# Metric functions for model results
# Zane
# 2024-10-03
# Moved here since they were duplicated in 09 and 11
###

# Define a summary function for calculating point estimates and CIs
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

# Calculates the area under each individual antibody landscape
# For full models, this is the AUC of the summary landscape
# For reduced models, this is the same as the intercept
AUC <- function(preds) {
	AUC <-
		preds |>
		dplyr::group_by(.draw) |>
		dplyr::summarize(
			out = pracma::trapz(x, .epred),
			.groups = "drop"
		) |>
		dplyr::pull(out)
	return(AUC)
}

# Get the intercept of each individual antibody landscape
# For full models, this is the intercept of the regression line
# For reduced models, this is (also the intercept) the estimated mean of the
# outcome across all strains, so the mean TI and post-titer (for the
# respective outcome) for the homologous strain only.
intercept <- function(preds) {
	intercept <-
		preds |>
		dplyr::ungroup() |>
		dplyr::filter(x == 0) |>
		dplyr::select(.draw, intercept = .epred) |>
		dplyr::pull(intercept)

	return(intercept)
}

# Find the proportion of antibody landscape over 40 (3 in log units)
# This is the breadth metric in the "new" metric set.
prop_above_threshold <- function(preds, outcome_name) {
	threshold <- dplyr::if_else(outcome_name == "log_posttiter", 3, 2)
	prop_above_threshold <-
		preds |>
		dplyr::ungroup() |>
		dplyr::mutate(above_threshold = as.numeric(.epred >= threshold)) |>
		dplyr::summarize(
			out = mean(above_threshold),
			.by = ".draw"
		) |>
		dplyr::pull(out)

	return(prop_above_threshold)
}

# Get the mean seroconversion rate across all of the strains.
mean_scr <- function(preds) {
	mean_scr <-
		preds |>
		dplyr::ungroup() |>
		dplyr::summarize(
			out = mean(.epred),
			.by = ".draw"
		) |>
		dplyr::pull(out)

	return(mean_scr)
}

# Need to fix old metrics
# Using current definitions of seroconversion and seroprotection, we don't
# need to correct them for censoring because censoring can only be a lower value
# so they would have seroconverted only.
# old metrics to get: mean rise, SCR, SPR only mean rise should change if we
# take censoring. Need to get mean rise over all strains and mean rise for
# only homologous strain.
# Magnitude: old is homologous GMT, new is linear regression intercept
# Breadth: old is SCR/SPR, new is prop_above_threshold
# Total Strength: old is GMT across strains, new is AUC

# Compute all of the metrics for a given data set
metric_set <- function(
		data, reduced, full, homologous, seroconversion, outcome_name, ...
) {
	#full <- metrics_input$full

	# Now calculate the list of metrics by invoking all of the functions.
	metrics_list <- tibble::tibble(
		AUC = AUC(full),
		intercept = intercept(full),
		prop_above_threshold = prop_above_threshold(full, outcome_name),
		mean_titer_homologous = intercept(homologous),
		mean_scr = mean_scr(seroconversion),
		mean_titer_all_strains = intercept(reduced)
	) |>
		dplyr::mutate(
			.draw = 1:dplyr::n(),
			.before = dplyr::everything()
		)

	# bind into one DF and clean up
	metrics <-
		metrics_list |>
		tidyr::pivot_longer(
			cols = -.draw,
			names_to = "metric",
			values_to = "metric_value"
		) |>
		dplyr::mutate(
			metric = factor(
				metric,
				levels = c(
					"intercept", "prop_above_threshold", "AUC",
					"mean_titer_homologous", "mean_scr", "mean_titer_all_strains"
				),
				labels = c(
					"Magnitude|New", "Breadth|New", "Total Strength|New",
					"Magnitude|Old", "Breadth|Old", "Total Strength|Old"
				)
			)
		) |>
		tidyr::separate_wider_delim(
			cols = metric,
			delim = "|",
			names = c("metric_name", "metric_set")
		)

	return(metrics)
}

# END OF FILE ###
