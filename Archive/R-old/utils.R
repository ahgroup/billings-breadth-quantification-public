###
# Utility functions
# Zane Billings
# 2024-08-12
# Small functions that are used in multiple places
###

# Function to save a dataset to csv and rds at the same time
write_csv_and_rds <- function(x, file) {
	readr::write_rds(x, paste0(file, ".Rds"), compress = "xz")
	readr::write_csv(x, paste0(file, ".csv"))

	invisible()
}



# Function to calculate pairwise temporal distances.
temporal_dist <- function(x) {
	# Use the short names for the actual names of the matrix.
	n <- hgp::replace_strain_names(x, from = "analysis", to = "short")

	# Set up an empty matrix to hold results
	res <- matrix(
		nrow = length(x),
		ncol = length(x)
	)

	# First get the years of isolation -- these are always the last 4 characters
	# in the name of the strain.
	x_years <- x |>
		as.character() |>
		stringr::str_sub(-4, -1) |>
		as.numeric()

	# Calculate the lower triangle of the distance matrix for all unique
	# combinations of two strains
	for (i in 2:nrow(res)) {
		for (j in 1:(i-1)) {
			res[[i, j]] <- abs(x_years[[i]] - x_years[[j]])
		}
	}

	# Set the diagonal to zero
	diag(res) <- 0

	# Quickly fill in the upper triangle
	out <- res |> as.dist() |> as.matrix()

	# Set the row and column names
	rownames(out) <- n
	colnames(out) <- n

	return(out)
}

# Function to compute distances based on cartographic maps.
racmaps_map_to_distances <- function(map) {
	coords <- Racmacs::agCoords(map)
	strains <- rownames(coords)
	x <- coords[, 1]
	y <- coords[, 2]

	res <- matrix(
		nrow = length(strains),
		ncol = length(strains),
		dimnames = list(strains, strains)
	)

	# Calculate the lower triangle of the distance matrix
	for (i in 2:nrow(res)) {
		for (j in 1:(i-1)) {
			x_part <- (x[[j]] - x[[i]])^2
			y_part <- (y[[j]] - y[[i]])^2
			dist <- sqrt(x_part + y_part)
			res[[i, j]] <- dist
		}
	}

	# We know all the diagonal entries are 0 (distance of a point to itself) so
	# set those without doing the calculation.
	diag(res) <- 0

	# Exploit R's dist class to automatically fill in the upper triangle
	out <- res |> as.dist() |> as.matrix()

	return(out)
}

# Convert a distance matrix into ggplot-approved format
tidy_dist_mat <- function(d) {
	out <- d |>
		tibble::as_tibble(rownames = "strain1") |>
		tidyr::pivot_longer(
			cols = -strain1,
			names_to = "strain2",
			values_to = "d"
		) |>
		# Order variable factors
		dplyr::mutate(
			strain1 = forcats::fct_inorder(strain1),
			strain2 = forcats::fct_inorder(strain2) |> forcats::fct_rev()
		)

	return(out)
}

# Next we'll write a function that fits this model to a given data set.
# This takes in the data and the model, along with a list of set sampling
# arguments that is (in this project) universal to all called Bayesian models.
# Any arguments that are needed for a specific model but not all models in the
# project can be passed in the ... argument.
fit_model_to_dataset <- function(data, model, sampling_arguments, ...) {
	dots <- list(...)
	args <- c("data" = list(data), sampling_arguments, dots)
	fit <- do.call(model$sample, args)
	return(fit)
}

# Reconstruct a brms model from each fit to get access to brms behind the
# scene stuff. Takes a stanfit (rstan-type) object, and the formula and
# data arguments that would be used to create a brms model, and basically
# tells brms that the fits are already instead of that stanfit object.
# Returns a brms object with samples in it.
cmdstan_to_brms <- function(stanfit, brms_formula, brms_data) {
	brms_holder <- brms::brm(
		brms_formula,
		data = brms_data,
		empty = TRUE,
		backend = "cmdstanr"
	)
	brms_holder$fit <- stanfit
	brms_holder <- brms::rename_pars(brms_holder)
	return(brms_holder)
}

# Function for min-max scaling of a vector
minmax <- function(x, ...) {
	return(
		(x - min(x, ...)) / (max(x, ...) - min(x, ...))
	)
}

# Function for coefficient of variation
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

# Function for the quartile coefficient of dispersion
qcd <- function(x, ...) {
	q <- quantile(x, probs = c(0.25, 0.75), ...)
	q1 <- q[[1]]
	q3 <- q[[2]]
	return((q3 - q1) / (q3 + q1))
}

# Function for gini coefficient
gini <- function(x, na.rm = FALSE) {
	n <- length(x)
	xbar <- mean(x, na.rm = na.rm)
	x_rank <- sort(x)

	out <- 2 / (n^2 * xbar) * sum(1:n * (x_rank - xbar))
	return(out)
}

# END OF FILE ####
