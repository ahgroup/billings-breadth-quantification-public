###
# Homologous strain intercept-only models
# Zane Billings
# 2024-08-26
# In order to calculate the mean rise corrected for censoring, we need to
# fit the censoring-corrected intercept-only model just to the homologous
# strain data -- the intercept parameter of this model estimates the GMT
# after correcting for censoring.
###

# Setup ###
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
source(here::here("R", "functions", "stan-code-generation.R"))
source(here::here("R", "functions", "model-runner.R"))

# Maximum number of brms draws
BRMS_DRAWS <- 1000L

# Data cleaning ####
## Full data set cleaning ####
cohort_data <-
	here::here("results", "data", "full-data-model-input.qs") |>
	qs::qread()

# Filter each of the nested data frames to only contain the homologous data
cohort_data_homologous_only <-
	cohort_data |>
	dplyr::mutate(
		# First get a list of which data points will be included
		homologous_obs = purrr::map2(
			data, vaccine_strain,
			\(d, v) which(as.character(v) == as.character(d$strain_name))
		),
		# Now filter the three dataset columns to only contain those
		# observations
		dplyr::across(
			c(data, brms_data),
			\(col) purrr::map2(
				col, homologous_obs,
				\(d, ind) d[ind, ]
			)
		)
	)

qs::qsave(
	cohort_data_homologous_only,
	here::here("results", "data", "homologous-full-data-model-input.qs")
)

## Subsample data set cleaning ####
subsample_data <-
	here::here("results", "data", "subsample-model-input.qs") |>
	qs::qread()

# Do the same thing for the subsample data
subsample_data_homologous_only <-
	subsample_data |>
	dplyr::mutate(
		# First get a list of which data points will be included
		homologous_obs = purrr::map2(
			data, vaccine_strain,
			\(d, v) which(as.character(v) == as.character(d$strain_name))
		),
		# Now filter the three dataset columns to only contain those
		# observations
		dplyr::across(
			c(data, brms_data),
			\(col) purrr::map2(
				col, homologous_obs,
				\(d, ind) d[ind, ]
			)
		)
	)

qs::qsave(
	subsample_data_homologous_only,
	here::here("results", "data", "homologous-subsample-model-input.qs")
)

# General model setup ####
# Now we need to load the stan file directories -- we only wanted to fit the
# reduced models though.
# Make a list of the model file names.
brms_model_info <-
	brms_model_info |>
	dplyr::filter(startsWith(file_name, "reduced"))
stan_model_files <- paste0(
	model_path, "/", brms_model_info$file_name, ".stan"
)

# Run the full data models ####
# Random seeds for parallel running from random dot org
seed_vec <- c(994810L, 345326L)

# Run the models -- see the model running script for more of an explanation.
full_data_model_output <-
	run_models(
		model_data = cohort_data_homologous_only,
		model_info = brms_model_info,
		model_files = stan_model_files,
		random_seed = seed_vec,
		dump_file = here::here(
			"results", "large-files", "dump"
		),
		time_files_dir = here::here(
			"results", "times", "homologous-full-data-models"
		),
		stan_csv_pth_base = here::here(
			"results", "large-files", "stanfits", "homologous-full-data-models"
		),
		brms_file_directory = here::here(
			"results", "large-files", "brms-fits", "homologous-full-data-models"
		)
	)

# Save the model data results
full_data_model_processed <-
	model_result_processing(
		full_data_model_output,
		cohort_data_homologous_only,
		here::here("results", "large-files", "homologous-full-data-model-results")
	)

# Run the subsample models ####
# L'Ecuyer-CMRG random seeds for parallel running from random dot org
seed_vec <- c(075017L, 354547L)

# Run the models and get the predictions
subsample_model_output <-
	run_models(
		model_data = subsample_data_homologous_only,
		model_info = brms_model_info,
		model_files = stan_model_files,
		random_seed = seed_vec,
		dump_file = here::here(
			"results", "large-files", "dump"
		),
		time_files_dir = here::here(
			"results", "times", "homologous-subsample-models"
		),
		stan_csv_pth_base = here::here(
			"results", "large-files", "stanfits", "homologous-subsample-models"
		),
		brms_file_directory = here::here(
			"results", "large-files", "brms-fits", "homologous-subsample-models"
		)
	)

# Save the subsample model results
subsample_model_processed <-
	model_result_processing(
		subsample_model_output,
		subsample_data_homologous_only,
		here::here("results", "large-files", "model-results",
							 "homologous-subsample-model-results")
	)

# END OF FILE ####
