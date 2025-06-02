###
# Panel Subsampling and Models
# Zane Billings
# 2024-08-15
# In order to analyze the robustness of the metric sets, we will create several
# simulated studies, which are subsampled from our overall cohort data set.
# Then, we can calculate our metric set on each of the labs so we can
# study how the metrics vary over the simulated studies.
###

# Setup ####
## Packages and dependencies ####
box::use(
	here,
	qs,
	dplyr,
	purrr,
	readr,
	tidyr,
	brms,
	cmdstanr
)

# Other code dependencies
source(here::here("R", "functions", "utils.R"))
source(here::here("R", "functions", "hai-helpers.R"))
source(here::here("R", "functions", "stan-code-generation.R"))
source(here::here("R", "functions", "model-runner.R"))

# Set a random seed from random dot org since we're doing subsampling
set.seed(6190982)

# Maximum number of brms draws
BRMS_DRAWS <- 1000L

## Data Loading ####
cohort_data <-
	here::here("results", "data", "analysis-data-formatted.qs") |>
	qs::qread()

make_strain_subset <- function(d, homologous, N, K) {
	# Get a list of all the strains in the dataset
	strain_list <- d |>
		# Remove the homologous strain so we don't subsample it
		dplyr::filter(strain_name != homologous) |>
		# Get only the strain name column
		dplyr::pull(strain_name) |>
		# Convert to character to prevent factor weirdness
		as.character() |>
		unique()

	# Now create N subsamples of K strains each
	strain_samples <-
		purrr::map(
			1:N,
			\(n) c(
				homologous,
				sample(strain_list, size = K, replace = FALSE)
			)
		) |>
		do.call(what = rbind)

	return(strain_samples)
}

make_individual_subset <- function(d, N, K) {
	# Get a list of all the individuals in the dataset
	id_list <- d |>
		dplyr::pull(subject_id) |>
		unique()

	# Now create N subsamples of K strains each
	id_samples <-
		purrr::map(
			1:N,
			\(n) sample(id_list, size = K, replace = FALSE)
		) |>
		do.call(what = rbind)

	return(id_samples)
}

# Subsampling ####
# First we'll set up the parameters to use for the subsampling.
subsampling_parms <- list(
	# Total number of simulated studies to create
	num_studies = 10L,
	# Total number of individuals to include in each study
	num_indiv_each = 100L,
	# Total number of heterologous strains to include in each study
	num_strains_each = 9L
)

# Now we actually get the subsamples. We first nest the data using the same
# grouping as last time (recall that each group will have a model fitted
# separately -- in this case, one model per group per subsample).
subsample_data <-
	cohort_data |>
	tidyr::nest(
		data = -c(strain_type, season, method, outcome_name, vaccine_strain)
	) |>
	dplyr::mutate(
		# Number of strains used (excluding homologous) in each season
		strain_counts = purrr::map_int(
			data,
			\(d) dplyr::n_distinct(d$strain_name) - 1
		)
	) |>
	# We only want to do the subsampling analysis for seasons which have at
	# least the number of strains we're interested in -- many of the later seasons
	# had a much smaller heterologous panel and thus we don't need to use them.
	dplyr::filter(
		strain_counts > subsampling_parms$num_strains_each
	) |>
	# Now we'll use the sampling functions we wrote before to get random samples
	# of the strains and of the individuals from each of the subset groups.
	dplyr::mutate(
		# Make the strain subsamples
		strains = purrr::map2(
			vaccine_strain, data,
			\(vaccine_strain, d) d |>
				make_strain_subset(
					as.character(vaccine_strain),
					subsampling_parms$num_studies,
					subsampling_parms$num_strains_each
				)
		),
		# Make the individual subsamples
		individuals = purrr::map(
			data,
			\(d) d |>
				make_individual_subset(
					subsampling_parms$num_studies,
					subsampling_parms$num_indiv_each
				)
		),
		# Reconstruct the actual dataset for each subsample -- this looks
		# complicated, but the outer pmap() ensures that we loop through each row
		# of the dataset we're in, and we need pmap() because we need to map over
		# three arguments in parallel -- the data, the strains to select, and the
		# individuals to select.
		# The inner map() then loops through the number of subsamples we want to
		# construct, which is defined in the subsampling_parms list. This is the
		# same as the length of each list-element of the strains and subject_id
		# columns. So we go through each row of the dataset, and filter the
		# current data to only contain the individuals and subjects specified in
		# each of the random samples we already got.
		subsample = purrr::pmap(
			list(data = data, strains = strains, individuals = individuals),
			\(data, strains, individuals, ...) purrr::map(
				1:subsampling_parms$num_studies,
				\(i) data |>
					dplyr::filter(
						(strain_name %in% strains[i, ]) &
							(subject_id %in% individuals[i, ])
					)
			)
		)
	)

# Now the dataset has a column where each row is a 10-element list (one
# element per subsample for the given season). We can expand those into 10
# rows of a column where each element is a data frame and all the other
# info is repeated, which is preferable for us because we can run that data
# frame through our pipeline.
# We can also throw away all the columns we don't need.
# Then we do the brms formatting from the previous models.
subsample_data_expanded <-
	subsample_data |>
	tidyr::unnest(subsample) |>
	dplyr::mutate(
		subsample_id = rep(
			1:(subsampling_parms$num_strains_each + 1),
			times = dplyr::n() / (subsampling_parms$num_strains_each + 1)
		) |>
			pad_numbers()
	) |>
	dplyr::select(
		season, strain_type, vaccine_strain, method, outcome_name,
		subsample_id, data = subsample
	) |>
	# Now that we have a dataframe where each row contains a nested tibble
	# which we want to pass to a model as the data object, we need to format
	# the data in the way brms expects, since we used brms to generate our
	# stan code.
	dplyr::mutate(
		# Now select only the columns we'll pass to the model -- this is the data
		# we would give to brms::brm() as the data argument if we were using
		# brms directly.
		brms_data = purrr::map(
			data,
			\(d) d |>
				dplyr::select(x = norm_d, c, y, y2, id = numeric_id) |>
				# Remove the NAs from the titer increase datasets
				tidyr::drop_na()
		)
	)

qs::qsave(
	subsample_data_expanded,
	here::here("results", "data", "subsample-model-input.qs")
)

# Model Fitting ####
# Now that the data are cleaned up, we need to fit our models.
# See the model_runner.R file for more details.

# Make a list of the model file names.
stan_model_files <- paste0(model_path, "/", brms_model_info$file_name, ".stan")

# We need four different random seeds, one for each of our sequential
# parallel runs. doFuture uses the L'Ecuyer-CMRG parallel random seed algorithm
# to ensure reproducibility. These seeds were generated on random.org.
seed_vec <- c(
	863688L,
	205021L,
	351237L,
	156225L
)

# Run the models
model_output <-
	run_models(
		model_data = subsample_data_expanded,
		model_info = brms_model_info,
		model_files = stan_model_files,
		random_seed = seed_vec,
		dump_file = here::here("results", "large-files", "dump"),
		time_files_dir = here::here("results", "times", "subsample-models"),
		stan_csv_pth_base = here::here(
			"results", "large-files", "stanfits", "subsample-models"
		),
		brms_file_directory = here::here(
			"results", "large-files", "brms-fits", "subsample-models"
		)
	)

# Model Postprocessing ####
# Now we're out of parallel mode. So we want to finish up our postprocessing,
# save the final prediction results.
# again see model_runner.R for details.
model_processed <-
	model_result_processing(
		model_output,
		subsample_data_expanded,
		here::here("results", "large-files", "model-results",
							 "subsample-model-results")
	)

# Save the subsample IDs vector since it doesn't get saved
qs::qsave(
	subsample_data_expanded$subsample_id,
	here::here("results", "data", "subsample-IDs.qs")
)

# END OF FILE ####
