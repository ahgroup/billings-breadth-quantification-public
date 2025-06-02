# generate stancode with brms

library(brms)
library(cmdstanr)

model_path <- here::here("stan", "simple-distance-model.stan")

brms_model_formula <- brms::brmsformula(
	y | cens(c, y2) ~ 1 + x + (1 + x | id),
	family = "gaussian",
	decomp = "QR"
)

brms_priors <- c(
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
)

cmdstan_sampling_arguments <-
	sampling_arguments <- list(
		#seed = 6735136L, # from random dot org
		chains = 4L,
		iter_warmup = 10,
		iter_sampling = 10,
		max_treedepth = 10,
		adapt_delta = 0.8
	)

generate_stan_code <- function(formula, priors, data,
															 file = here::here("stan", "linear-model.stan")) {
		model_fit <- brms::brm(
			formula = formula,
			data = data,
			prior = priors,
			backend = 'cmdstanr',
			empty = TRUE
		)

		model_fit |>
			brms::stancode() |>
			writeLines(file)

		invisible()
}




