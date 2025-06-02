###
# Generate Stan code from processed data
# Zane Billings
# 2024-08-20
###
# generate stancode with brms

suppressPackageStartupMessages({
	library(brms)
	library(cmdstanr)
})

model_path <- here::here("stan")







readr::write_rds(
	brms_model_info,
	here::here("results", "data", "brms-model-info.Rds")
)

cmdstan_sampling_arguments <-
	sampling_arguments <- list(
		chains = 4L,
		iter_warmup = 100,
		iter_sampling = 250,
		max_treedepth = 10,
		adapt_delta = 0.8
	)

generate_stan_code <- function(model_info, data, pth_base) {

	empty_brms_specs <- purrr::pmap(
		model_info,
		\(formula, priors, ...) brms::brm(
			formula = formula,
			data = data,
			prior = priors,
			backend = 'cmdstanr',
			empty = TRUE
		)
	)

	file_paths <- paste0(
		pth_base, "/", model_info$file_name,
		".stan"
	)

	catch <- purrr::walk2(
		empty_brms_specs, file_paths,
		\(spec, f) spec |>
			brms::stancode() |>
			writeLines(f)
	)

	invisible(empty_brms_specs)
}




