###
# Simulation functions
# Zane
# 2025-01-15
# Functions for generative simulations of HAI data
###

# First we'll write a function that can simulate titers as a linear function
# of antigenic distance, which we scale from zero to one. Then, we need to
# write a function that handles generating data from multiple labs and
# subsampling distances from a "universe" of antigenic distance values.
one_titer_sim <- function(N = 1e4L, mean = 3, sd = 2) {
	sim <-
		tibble::tibble(
			# Assume log(titer) is drawn from a normal distribution
			raw_log_titer = rnorm(N, mean, sd),
			# If we observe a titer with log(titer) < 1 (LOD), mark it as 0
			trunc_log_titer = ifelse(raw_log_titer >= 1, raw_log_titer, 0),
			# The assay is dilution based, so we only observe the floor of each
			# value.
			rounded_titer = floor(trunc_log_titer),
			# Now final observed titer is equal to this transformation.
			sim_titer = 5 * 2 ^ rounded_titer
		) |>
		dplyr::arrange(raw_log_titer)

	return(sim)
}

sim_with_distance <- function(
		N = 1e2L, distances = seq(0, 1, 0.1),
		intercept = 4, slope = -3, sd = 2
) {
	sim_means <-
		tibble::tibble(
			d = distances,
			mu = intercept + slope * distances
		) |>
		dplyr::mutate(
			sim = purrr::map(mu, \(x) one_titer_sim(N, mean = x, sd = sd))
		) |>
		tidyr::unnest(sim)

	return(sim_means)
}

one_titer_sim_lm <- function(
		N = 1e4L, distances, seed = 859947L,
		intercept, slope_dist, slope_pre, dist_pre_interaction, post_sd = 2,
		pre_mean, pre_sd
) {
	set.seed(seed)
	sim <-
		tidyr::expand_grid(
			id = seq(1, N),
			d = distances
		) |>
		dplyr::mutate(
			pre_titer = rnorm(dplyr::n(), pre_mean, pre_sd),
			mu = intercept + slope_dist * d + slope_pre * pre_titer +
				dist_pre_interaction * d * pre_titer,
			# Assume log(titer) is drawn from a normal distribution
			raw_log_titer = rnorm(dplyr::n(), mu, post_sd),
			# If we observe a titer with log(titer) < 1 (LOD), mark it as 0
			trunc_log_titer = ifelse(raw_log_titer >= 1, raw_log_titer, 0),
			# The assay is dilution based, so we only observe the floor of each
			# value.
			rounded_titer = floor(trunc_log_titer),
			# Now final observed titer is equal to this transformation.
			sim_post_titer = 5 * 2 ^ rounded_titer,
			# Do the same thing to the pre titer
			trunc_pre_titer = ifelse(pre_titer >= 1, pre_titer, 0),
			rounded_pre_titer = floor(trunc_pre_titer),
			sim_pre_titer = 5 * 2 ^ rounded_pre_titer
		) |>
		dplyr::arrange(id, d)

	return(sim)
}

# hierarchical_sim_with_distance <- function(
# 		N = 1e2L, distances = seq(0, 1, 0.1), seed = 859947L,
# 		intercept_mu = 4, intercept_sigma = 1,
# 		slope_dist_mu = -3, slope_dist_sigma = 1,
# 		dist_pre_interaction = 0.25,
# 		slope_pre_mu = 1, slope_pre_sigma = 0.5,
# 		post_sd = 2
# ) {
# 	sim <- purrr::map(
# 		distances,
# 		\(x) one_titer_sim_lm(
# 			N,
# 			intercept = rnorm(N, intercept_mu, intercept_sigma),
# 			slope_dist = rnorm(N, slope_dist_mu, slope_dist_sigma),
# 			distances = x,
# 			post_sd = post_sd,
# 			seed = seed
# 		)
# 	) |>
# 		dplyr::bind_rows() |>
# 		dplyr::mutate(d = as.numeric(d), .before = dplyr::everything()) |>
# 		dplyr::arrange(d, id)
#
# 	return(sim)
# }

sim_labs <- function(
		individuals_per_lab = 100L,
		strains_per_lab = 9L,
		n_labs = 10L,
		strain_universe = seq(0.02, 1, 0.02),
		sim_fun = hierarchical_sim_with_distance,
		global_seed = 101L,
		strain_count_method = "literal",
		...
) {
	sim_fun_args <- rlang::list2(...)

	if (strain_count_method == "pois") {
		set.seed(global_seed)
		strains_per_lab_vec <- rpois(n_labs, strains_per_lab)
	} else if (strain_count_method == "literal") {
		strains_per_lab_vec <- rep(strains_per_lab, times = n_labs)
	}

	lab_panels <- purrr::map2(
		1:n_labs, strains_per_lab_vec,
		\(idx, n_strains) sample(strain_universe, n_strains, replace = FALSE) |>
			c(0) |>
			sort()
	)

	# Generate a seed for each lab -- this ensures we get reproducible results
	# but each lab is not the same
	set.seed(global_seed)
	lab_seeds <- sample(100000L:999999L, size = n_labs)

	lab_sim <- purrr::map2(
		lab_panels, lab_seeds,
		# Have to use inject to pass in the spliced sim_fun_args list
		\(this_panel, this_seed) rlang::inject(sim_fun(
			N = individuals_per_lab,
			distances = this_panel,
			seed = this_seed,
			!!!sim_fun_args
		))
	)

	lab_names <- paste0("Lab ", pad_numbers(1:n_labs))

	out <-
		lab_sim |>
		rlang::set_names(lab_names) |>
		dplyr::bind_rows(.id = "lab")

	return(out)
}
