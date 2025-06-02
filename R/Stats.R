###
# Functions for calculating various statistics
# Zane
# 2025-05-05
###

# Breadth summary statistics ####
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

# Math / stats helpers ####
## coefficient of variation
covar <- function(x, correct = NULL, na.rm = FALSE, ...) {

	cv <- sd(x, na.rm) / mean(x, na.rm = na.rm, ...)

	if (is.null(correct)) {
		return(cv)
	} else {
		correct_san <- tolower(correct)[[1]]
		n <- length(x)
	}

	if (correct_san == "n") {
		cv <- (1 + 1/(4*n)) * cv
	} else if (correct_san == "l") {
		cv <- sqrt(exp(var(log(x))) - 1)
	} else {
		stop("'correct' argument should be '(n)ormal' or '(l)og-normal.")
	}

	return(cv)
}

## quartile coefficient of dispersion
qcd <- function(x, ...) {
	q <- quantile(x, probs = c(0.25, 0.75), ...)
	q1 <- q[[1]]
	q3 <- q[[2]]
	return((q3 - q1) / (q3 + q1))
}

## gini coefficient
gini <- function(x, na.rm = FALSE) {
	n <- length(x)
	xbar <- mean(x, na.rm = na.rm)
	x_rank <- sort(x)

	out <- 2 / (n^2 * xbar) * sum(1:n * (x_rank - xbar))
	return(out)
}
