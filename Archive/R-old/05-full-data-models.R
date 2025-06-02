###
# Titer vs. distance models on the entire dataset
# Zane Billings
# 2024-08-13
# Fit Bayesian linear regression models to titer vs. antigenic distance.
# We'll use the entire dataset, but stratify by subtype and season. So for
# each season in the data, we will fit two separate models, one for H1 and one
# for H3.
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
source(here::here("R", "functions", "stan-code-generation.R"))
source(here::here("R", "functions", "model-runner.R"))

# Maximum number of brms draws
BRMS_DRAWS <- 1000L

## Load data ####
# First load the joined data.
joined_data <-
	here::here("data", "processed", "joined-data.Rds") |>
	readr::read_rds()

# Data preprocessing ####
# First, we want to normalize our antigenic distance data. This will make
# intepreting the data and comparing the models easier. Since we'll fit
# models separately for each combination of season and subtype, we'll
# normalize the antigenic distance data separately within each of those strata.
model_data <-
	joined_data |>
	# Select the variables we'll use for modeling -- we won't use any of the
	# demographic data so there's no point dragging it around.
	dplyr::select(
		row_id, unique_id, subject_id, season, study,
		strain_name, strain_type, seroconversion,
		pretiter, posttiter, fold_change, vaccine_strain, method, d
	) |>
	# Now calculated the normalized distance
	dplyr::group_by(season, strain_type, method) |>
	dplyr::mutate(norm_d = d / max(d)) |>
	# We also need to create a numeric subject ID to pass to Stan, so let's do
	# that now. There's a few ways to do this but the easiest way is to turn
	# it into a factor, which is actually an integer under the hood that we
	# can grab, so each level gets assigned a unique integer value (in order
	# of the factor levels, we call the forcats function to put them in the
	# order the subjects appear in the data).
	dplyr::mutate(
		numeric_id = as.integer(forcats::fct_inorder(as.factor(subject_id))),
		.after = subject_id
	) |>
	dplyr::ungroup()

# We want to fit a separate model to each combination of subtype and season.
# We also want to fit models that use each of the different antigenic distance
# measurements, and it's possible that we would like to fit models to both
# the post-vaccination titer and the fold change as outcomes.
# First we'll deal with the two outcome variables, because we want to pivot
# them to long form. We also need to log transform them, so we'll go ahead
# and do that.
# TODO add the HAI functions to hgp
model_data_long_outcome <-
	model_data |>
	# First put the outcomes on the log scale, using normal log transform
	# for fold change, and the HAI log transform log2(x/5) for posttiter.
	dplyr::mutate(
		log_posttiter = hai_to_log_scale(posttiter),
		log_pretiter = hai_to_log_scale(pretiter),
		titer_increase = log2(fold_change),
		.keep = "unused"
	) |>
	# Now pivot the two outcomes to long format
	tidyr::pivot_longer(
		c(log_posttiter, titer_increase),
		names_to = "outcome_name",
		values_to = "outcome_value"
	)

# Now do the hai formatting -- because of the way I wrote this function it
# is a little bit annoying.
model_data_hai_formatted <-
	model_data_long_outcome |>
	# First we create two nested data frames, one for post-vac titer outcomes
	# and one for titer increase outcomes.
	tidyr::nest(data = -outcome_name) |>
	# Now add a column for the HAI transformed data, we'll just replace the
	# existing data column since we don't lose anything. This takes care of the
	# correct formatting based on the outcome name.
	dplyr::mutate(
		data =
			purrr::map2(
				outcome_name, data,
				\(outcome_name, data)
				format_hai_data(
					data,
					post_titer = "outcome_value",
					pre_titer = "log_pretiter",
					log_scale = TRUE,
					log_out = TRUE,
					increase = dplyr::if_else(
						outcome_name == "titer_increase", TRUE, FALSE
					)
				)
			)
	) |>
	# Now remove the nesting so we get a flat data frame again
	tidyr::unnest(data) |>
	# Some of the titer increases are NA and will contribute nothing to the model
	# so we drop those ones.
	tidyr::drop_na()

# Save this so we can use it as the starting point for later analyses
qs::qsave(
	model_data_hai_formatted,
	here::here("results", "data", "analysis-data-formatted.qs")
)

# Now we nest the data because that makes it a lot easier to fit a model
# for each stratum. Once the data is nested we need to prep it to be passed
# to brms. That mostly means we need to create a new dataset that only has
# exactly the columns the model will see, and we need to format the
# censored outcome using the special brms censoring format.
nested_model_data <-
	model_data_hai_formatted |>
	dplyr::select(-seroconversion) |>
	tidyr::nest(data = -c(strain_type, season, method, outcome_name,
												vaccine_strain)) |>
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
				dplyr::select(x = norm_d, c, y, y2, id = numeric_id)
		)
	)

qs::qsave(
	nested_model_data,
	here::here("results", "data", "full-data-model-input.qs")
)

# Setup the bayesian model ####
# Generate the stan code for cmdstan based on a dummy brms model
generate_stan_code(
	brms_model_info,
	data = nested_model_data$brms_data[[1]],
	pth_base = model_path
)

# Fit the models ####
# Make a list of the model file names.
stan_model_files <- paste0(model_path, "/", brms_model_info$file_name, ".stan")

# We need four different random seeds, one for each of our sequential
# parallel runs. doFuture uses the L'Ecuyer-CMRG parallel random seed algorithm
# to ensure reproducibility. These seeds were generated on random.org.
seed_vec <- c(
	617637L,
	488629L,
	896670L,
	403975L
)

# Run the models -- see the model running script for more of an explanation.
model_output <-
	run_models(
		model_data = nested_model_data,
		model_info = brms_model_info,
		model_files = stan_model_files,
		random_seed = seed_vec,
		dump_file = here::here("results", "large-files", "dump"),
		time_files_dir = here::here("results", "times", "full-data-models"),
		stan_csv_pth_base = here::here(
			"results", "large-files", "stanfits", "full-data-models"
		),
		brms_file_directory = here::here(
			"results", "large-files", "brms-fits", "full-data-models"
		)
	)

# Model postprocessing ####
# Now we're out of parallel mode. So we want to finish up our postprocessing,
# save the final prediction results.
model_processed <-
	model_result_processing(
		model_output,
		nested_model_data,
		here::here("results","large-files", "model-results",
							 "full-data-model-results")
	)

# END OF FILE ####
