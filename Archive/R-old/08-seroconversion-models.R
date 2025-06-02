###
# Models for estimating seroconversion rate
# Zane Billings
# 2024-08-29
# Even though the seroconversion rate doesn't have a censoring correction,
# we still want to estimate it in a Bayesian way to ensure we have a sampling
# distribution from which we can create credible intervals. So to do that we
# can fit a brms model with a bernoulli family and a logit link, with only
# a global intercept and no parameters. Then the estimated intercept will be
# an estimate of the seroconversion rate across all of the strains included
# in a given model and we'll have HMC samples we can use to construct CrIs
# of any functions of interest.
###

# Setup ####
## Package dependencies ####
box::use(
	here,
	readr,
	dplyr,
	purrr,
	forcats,
	tidyr,
	brms,
	cmdstanr,
	qs
)

source(here::here("R", "functions", "utils.R"))
source(here::here("R", "functions", "hai-helpers.R"))
source(here::here("R", "functions", "logistic-model-generation.R"))
source(here::here("R", "functions", "model-runner.R"))

# Maximum number of brms draws
BRMS_DRAWS <- 1000L

## Load data ####
# Data for full sample results
full_model_data <-
	here::here("results", "data", "analysis-data-formatted.qs") |>
	qs::qread()

# Data for subsample results
subsample_model_data <-
	here::here("results", "data", "subsample-model-input.qs") |>
	qs::qread() |>
	dplyr::select(-brms_data)

# Format the data with the seroconversion outcomes
full_model_data_nested <-
	full_model_data |>
	# Remove the posttiter and titerincrease outcome related outcomes
	dplyr::select(-c(outcome_name, outcome_value, c, y, y2)) |>
	# Remove the rows that were duplicated for each outcome
	dplyr::distinct() |>
	tidyr::nest(data = -c(strain_type, season, method, vaccine_strain)) |>
	# Add a unique ID for each group/row
	dplyr::mutate(
		group_id = dplyr::row_number() |> pad_numbers(),
		.before = dplyr::everything()
	) |>
	dplyr::mutate(
		# Now select only the columns we'll pass to the model -- this is the data
		# we would give to brms::brm() as the data argument if we were using
		# brms directly.
		brms_data = purrr::map(
			data,
			\(d) d |>
				dplyr::select(x = norm_d, y = seroconversion)
		)
	)

# Do the same thing for the subsample models
subsample_model_data_nested <-
	subsample_model_data |>
	dplyr::filter(outcome_name == "log_posttiter") |>
	dplyr::select(-outcome_name) |>
	dplyr::mutate(
		data = purrr::map(
			data,
			\(d) dplyr::select(d, -c(outcome_value, c, y, y2))
		),
		brms_data = purrr::map(
			data,
			\(d) dplyr::select(d, x = norm_d, y = seroconversion)
		)
	)

## Generate the logistic model stan code ####
# Generate the stan code for cmdstan based on a dummy brms model
generate_stan_code(
	brms_model_info,
	data = full_model_data_nested$brms_data[[1]],
	pth_base = model_path
)
stan_model_files <- paste0(
	model_path, "/", brms_model_info$file_name, ".stan"
)

# Run the models ####
# Random seeds for parallel running from random dot org
seed_vec <- c(210220L)

# Run the models -- see the model running script for more of an explanation.
full_data_model_output <-
	run_models(
		model_data = full_model_data_nested,
		model_info = brms_model_info,
		model_files = stan_model_files,
		random_seed = seed_vec,
		dump_file = here::here(
			"results", "large-files", "dump"
		),
		time_files_dir = here::here(
			"results", "times", "seroconversion-full-data-models"
		),
		stan_csv_pth_base = here::here(
			"results", "large-files", "stanfits", "seroconversion-full-data-models"
		),
		brms_file_directory = here::here(
			"results", "large-files", "brms-fits", "seroconversion-full-data-models"
		)
	)

# Save the model data results
full_data_model_processed <-
	model_result_processing(
		full_data_model_output,
		full_model_data_nested,
		here::here(
			"results", "large-files", "model-results",
			"seroconversion-full-data-model-results"
		)
	)

# Run the subsample models ####
# L'Ecuyer-CMRG random seeds for parallel running from random dot org
seed_vec <- c(567105L)

# Run the models and get the predictions
subsample_model_output <-
	run_models(
		model_data = subsample_model_data_nested,
		model_info = brms_model_info,
		model_files = stan_model_files,
		random_seed = seed_vec,
		dump_file = here::here(
			"results", "large-files", "dump"
		),
		time_files_dir = here::here(
			"results", "times", "seroconversion-subsample-models"
		),
		stan_csv_pth_base = here::here(
			"results", "large-files", "stanfits", "seroconversion-subsample-models"
		),
		brms_file_directory = here::here(
			"results", "large-files", "brms-fits", "seroconversion-subsample-models"
		)
	)

# Save the subsample model results
subsample_model_processed <-
	model_result_processing(
		subsample_model_output,
		subsample_model_data_nested,
		here::here("results", "large-files", "model-results",
							 "seroconversion-subsample-model-results")
	)

# END OF FILE ####
