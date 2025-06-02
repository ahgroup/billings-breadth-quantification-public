###
# Subsample generation
# Zane
# 2025-05-02
# Using the main dataset, generate the subsamples that we need
# We want m subsamples per cohort that each have n individuals subsampled from
# our size N_c cohort and k different heterologous strains plus the
# homologous strain.
###

# Helper function that gets a subset of the strains for a given dataset
make_strain_subset <- function(
		dataset,
		homologous_strain_name,
		num_subsamples,
		num_strains_per_subsample
	) {
	# Get a list of all the strains in the dataset
	strain_list <- dataset |>
		# Remove the homologous strain so we don't subsample it
		dplyr::filter(assay_strain != homologous_strain_name) |>
		# Get only the strain name column
		dplyr::pull(assay_strain) |>
		# Convert to character to prevent factor weirdness
		as.character() |>
		unique()

	# Now create N subsamples of K strains each
	# Add the homologous strain to each subsample as well
	strain_samples <-
		purrr::map(
			1:num_subsamples,
			\(n) c(
				homologous_strain_name,
				sample(strain_list, size = num_strains_per_subsample, replace = FALSE)
			)
		) |>
		do.call(what = rbind)

	return(strain_samples)
}

# Helper function that gets a subset of the individuals for a given dataset
make_individual_subset <- function(
		dataset,
		num_subsamples,
		num_ids_per_subsample
	) {
	# Get a list of all the individuals in the dataset
	id_list <- dataset |>
		dplyr::pull(subject_id) |>
		unique()

	# Now create N subsamples of K individuals each
	id_samples <-
		purrr::map(
			1:num_subsamples,
			\(n) sample(id_list, size = num_ids_per_subsample, replace = FALSE)
		) |>
		do.call(what = rbind)

	return(id_samples)
}

generate_subsampling_parms_list <- function(
		subsamples_per_cohort,
		individuals_per_subsample,
		heterologous_strains_per_subsample
	) {
	subsampling_parms <- list(
		# Total number of simulated studies to create
		num_studies = as.integer(subsamples_per_cohort),
		# Total number of individuals to include in each study
		num_indiv_each = as.integer(individuals_per_subsample),
		# Total number of heterologous strains to include in each study
		num_strains_each = as.integer(heterologous_strains_per_subsample)
	)

	return(subsampling_parms)
}

make_subsample_dataset <- function(
	model_data,
	subsampling_parms,
	subsampling_seed
	) {
	# Set the prng seed for the subsampling
	set.seed(subsampling_seed)

	# Now we actually get the subsamples. We first nest the data using the same
	# grouping as last time (recall that each group will have a model fitted
	# separately -- in this case, one model per group per subsample).
	subsample_data <-
		model_data |>
		tidyr::nest(
			data = -c(season, metric, vaccine_strain)
		) |>
		# Get a count of how many heterologous strains were able to be sampled
		# each season
		dplyr::mutate(
			# Number of strains used (excluding homologous) in each season
			strain_counts = purrr::map_int(
				data,
				\(d) dplyr::n_distinct(d$assay_strain) - 1
			)
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
							(assay_strain %in% strains[i, ]) &
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
	subsample_data_expanded <-
		subsample_data |>
		tidyr::unnest(subsample) |>
		dplyr::mutate(
			subsample_id = rep(
				1:(subsampling_parms$num_studies),
				times = dplyr::n() / (subsampling_parms$num_studies)
			) |>
				zero_pad_numbers() |>
				factor(ordered = TRUE)
		) |>
		# We can also throw away all the columns we don't need.
		# We aren't actually throwing away info because it's all encoded in the data.
		dplyr::select(
			season, vaccine_strain, metric, strain_counts,
			subsample_id, dat = subsample
		)

	# Finally we need to make the dataset column so we can get the homologous
	# only models as well.

	# First make a copy that has only the homologous data
	# It's faster to unnest, filter, renest than it is to map over the dfs.
	col_names <- names(subsample_data_expanded)
	col_names <- col_names[!col_names == "dat"]

	homologous_model_data <-
		subsample_data_expanded |>
		tidyr::unnest(dat) |>
		dplyr::filter(
			as.character(vaccine_strain) == as.character(assay_strain)
		) |>
		tidyr::nest(dat = !dplyr::all_of(col_names)) |>
		dplyr::mutate(homologous_only = TRUE)

	# Add the homologous_only label to the other one
	all_model_data <- subsample_data_expanded |>
		dplyr::mutate(homologous_only = FALSE)

	# create nested data based on the subsets
	out <- dplyr::bind_rows(
		all_model_data,
		homologous_model_data
	) |>
		dplyr::mutate(
			dataset = factor(
				homologous_only,
				levels = c(TRUE, FALSE),
				labels = c("homologous", "full")
			),
			.before = strain_counts
		) |>
		dplyr::select(-homologous_only)

	return(out)
}
