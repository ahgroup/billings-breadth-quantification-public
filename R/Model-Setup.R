####
# brms model formula and prior definitions
# Zane Billings
# 2025-04-28
####

suppressPackageStartupMessages({
	library(cmdstanr)
	library(brms)
	library(rstan)
})

generate_brms_control_arguments <- function(is_test_run) {
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
				max_treedepth = 12,
				adapt_delta = 0.99
			)
	} else {
		cli::cli_abort(c(
			"{.var is_test_run} must be TRUE or FALSE, not {.val {is_test_run}}."
		))
	}

	return(cmdstan_sampling_arguments)
}

generate_brms_model_info <- function() {
	# Set up the model formulas
	brms_model_formulas <- list(
		"reduced, no censoring" = brms::brmsformula(
			y ~ 1,
			family = "gaussian"
		),
		"reduced, with censoring" = brms::brmsformula(
			y | cens(c, y2) ~ 1,
			family = "gaussian"
		),
		"full, no censoring" = brms::brmsformula(
			y ~ 1 + d_norm + (1 + d_norm | subject_id),
			family = "gaussian"
		),
		"full, with censoring" = brms::brmsformula(
			y | cens(c, y2) ~ 1 + d_norm + (1 + d_norm | subject_id),
			family = "gaussian"
		),
		"logistic" = brms::brmsformula(
			seroconversion ~ 1,
			family = brms::bernoulli(link = "logit")
		)
	)

	# Set up the families, all gaussian except logistic
	brms_model_families <- list(
		"reduced, no censoring" = brms::brmsfamily("gaussian"),
		"reduced, with censoring" = brms::brmsfamily("gaussian"),
		"full, no censoring" = brms::brmsfamily("gaussian"),
		"full, with censoring" = brms::brmsfamily("gaussian"),
		"logistic" = brms::brmsfamily("bernoulli", "logit")
	)

	# Set up the priors -- brms will only let us pass priors that are accepted by
	# a model with no extras so we have to specify which ones we need at each step.
	brms_model_priors <- list(
		"reduced, no censoring" = c(
			# Prior for intercepts
			brms::prior(
				"student_t(3, 0, 3)", class = "Intercept"
			),
			# Prior for residual variance
			brms::prior(
				"student_t(3, 0, 1)", class = "sigma", lb = 0
			)
		),
		"reduced, with censoring" = c(
			# Prior for intercepts
			brms::prior(
				"student_t(3, 0, 3)", class = "Intercept"
			),
			# Prior for residual variance
			brms::prior(
				"student_t(3, 0, 1)", class = "sigma", lb = 0
			)
		),
		"full, no censoring" = c(
			# Prior for intercepts
			brms::prior(
				"student_t(3, 0, 3)", class = "Intercept"
			),
			# Prior for slopes
			brms::prior(
				"student_t(3, 0, 1)", class = "b"
			),
			# Prior for variances, we need two. One for covariances and one for global
			# variance parameter
			brms::prior(
				"student_t(3, 0, 1)", class = "sigma", lb = 0
			),
			brms::prior(
				"student_t(3, 0, 1)", class = "sd", lb = 0
			),
			# Priors for cholesky factors
			brms::prior(
				"lkj_corr_cholesky(2)", class = "L"
			)
		),
		"full, with censoring" = c(
			# Prior for intercepts
			brms::prior(
				"student_t(3, 0, 3)", class = "Intercept"
			),
			# Prior for slopes
			brms::prior(
				"student_t(3, 0, 1)", class = "b"
			),
			# Prior for variances, we need two. One for covariances and one for global
			# variance parameter
			brms::prior(
				"student_t(3, 0, 1)", class = "sigma", lb = 0
			),
			brms::prior(
				"student_t(3, 0, 1)", class = "sd", lb = 0
			),
			# Priors for cholesky factors
			brms::prior(
				"lkj_corr_cholesky(2)", class = "L"
			)
		),
		"logistic" = c(
			# Prior for intercepts
			brms::prior(
				"student_t(3, 0, 3)", class = "Intercept"
			)
		)
	)

	brms_model_info <- tibble::tibble(
		model = c('reduced', 'reduced', 'full', 'full', 'logistic'),
		censoring = c('no', 'yes', 'no', 'yes', 'logistic'),
		formula = brms_model_formulas,
		family = brms_model_families,
		priors = brms_model_priors
	)

	return(brms_model_info)
}

nest_model_data <- function(model_data) {
	# First make a copy that has only the homologous data
	homologous_model_data <- model_data |>
		dplyr::filter(as.character(vaccine_strain) == as.character(assay_strain)) |>
		dplyr::mutate(homologous_only = TRUE)

	# Add the homologous_only label to the other one
	all_model_data <- model_data |>
		dplyr::mutate(homologous_only = FALSE)

	# create nested data based on the subsets
	model_strata <- dplyr::bind_rows(
		all_model_data,
		homologous_model_data
	) |>
		tidyr::nest(dat = -c(season, metric, homologous_only)) |>
		dplyr::mutate(
			dataset = factor(
				homologous_only,
				levels = c(TRUE, FALSE),
				labels = c("homologous", "full")
			)
		) |>
		dplyr::select(-homologous_only)

	return(model_strata)
}

generate_model_metadata <- function(
		brms_model_info, nested_model_data,
		random_seeds, global_seed
	) {
	# Cross with the other model info to get a data frame with one row per
	# model we need to fit.
	# Should have nrow(augmented_model_data) * number models * k_{dose group}
	# rows after crossing
	all_model_combs <- tidyr::crossing(
		nested_model_data,
		tidyr::nesting(brms_model_info)
	) |>
		# We don't need to run the full model on homologous data
		dplyr::filter(!(model == "full" & dataset == "homologous"))

	# Randomly chose a global seed
	# Got this number between 1 and 10k on random dot org, which is also
	# where the random seeds came from
	# there are 10k random seeds so we can randomly sample from the list of
	# them
	set.seed(random_seeds[global_seed])

	# Add seeds and clean up
	model_metadata <-
		all_model_combs |>
		dplyr::mutate(
			# Add the seed column
			fitting_seed = sample(
				random_seeds,
				size = dplyr::n()
			),
			model = factor(model),
			censoring = factor(censoring)
		) |>
		dplyr::arrange(model, censoring, dataset, season, metric) |>
		dplyr::relocate(dat, .after = censoring)

	return(model_metadata)
}

listify_df <- function(df) {
	out <- purrr::pmap(df, \(...) list(...))

	return(out)
}

model_dispatch <- function(
		md_input, fitting_arguments
	) {
	metadata_list <- md_input[[1]]
	model_fit <- brms::brm(
		formula = metadata_list$formula,
		data = metadata_list$dat,
		family = metadata_list$family,
		prior = metadata_list$priors,
		seed = metadata_list$fitting_seed,
		sample_prior = "yes",
		chains = fitting_arguments$chains,
		cores = fitting_arguments$chains,
		warmup = fitting_arguments$iter_warmup,
		iter = fitting_arguments$iter_warmup + fitting_arguments$iter_sampling,
		algorithm = "sampling",
		control = list(
			adapt_delta = fitting_arguments$adapt_delta,
			max_treedepth = fitting_arguments$max_treedepth
		),
		backend = "cmdstanr",
		save_pars = brms::save_pars(all = TRUE),
		silent = 0,
		refresh = 0
	)

	return(model_fit)
}

generate_brms_model_filenames <- function(
		model_metadata, save_dir
	) {
	ensure_directory_exists(save_dir)

	s <- substr(model_metadata$season, 1, 4) |>
		as.integer()
	s2 <- paste0(s, "-", s+1)

	cens <- dplyr::case_match(
		model_metadata$censoring,
		"yes" ~ "model_wcc",
		"no" ~ "model_ncc",
		"logistic" ~ "model"
	)

	mn <- paste(
		s2, model_metadata$metric, model_metadata$model,
		cens, model_metadata$dataset, sep = "_"
	)

	n <- nrow(model_metadata)

	i <-
		seq_len(n) |>
		stringr::str_pad(width = nchar(n), side = "left", pad = "0")

	mf <- paste0("m", i, "_", mn, ".qs2")

	fn <- here::here(save_dir, mf)

	return(fn)
}

nest_homologous_data <- function(model_data) {
	homologous_model_data <-
		model_data |>
		dplyr::filter(as.character(assay_strain) == as.character(vaccine_strain)) |>
		dplyr::mutate(
			assay_strain = forcats::fct_drop(assay_strain),
			vaccine_strain = forcats::fct_drop(vaccine_strain)
		)

	model_strata <- homologous_model_data |>
		nest_model_data()

	return(model_strata)
}
