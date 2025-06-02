###
# Antigenic Distance Calculation
# Zane Billings
# 2024-08-13
# Calculates three different antigenic distance measurements for all of the
# strains we have available.
# 1. temporal: absolute different in isolation years between two strains.
# 2. p-epitope: a sequence based distance which is defined as the maximum
# proportion of different sites across 5 different epitope regions.
# 3. cartography: uses MDS to construct a cartographic map which positions the
# strains in 2D space (in antigenic units). From the map we can extract
# Euclidean distances between all strains.
###

# Setup ####

## Packages and dependencies ####
box::use(
	here,
	readr,
	dplyr,
	tidyr,
	Racmacs[...]
)

source(here::here("R", "functions", "utils.R"))
source(here::here("R", "functions", "p-epitope-calculator.R"))

## Global settings ####
# Number of cores that can be used in parallel
N_CORES <- parallelly::availableCores(constraints = "connections", omit = 2L)
# Number of times to rerun the cartography optimization routine
N_OPT_RUNS <- 5L
N_OPT_EACH <- 10L

## Data loading ####
# We need to load the virus data
sequence_data <-
	here::here("data", "processed", "sequence-data.Rds") |>
	readr::read_rds()

# And the cohort data
cohort_data <-
	here::here("data", "processed", "analysis-data.Rds") |>
	readr::read_rds()

# And we'll also give an easier to access name to the strain names data
strain_names <-
	sequence_data |>
	# Remove the accession_no column that we don't need for analysis
	dplyr::select(-accession_no) |>
	# Add a subtype column that we can use later, and also add a column for
	# the longer strain name.
	dplyr::mutate(
		analysis_name = hgp::replace_strain_names(
			strain_name, from = "short", to = "analysis"
		),
		subtype = hgp::replace_strain_names(
			strain_name, from = "short", to = "subtype"
		),
		.after = strain_name
	) |>
	# Remove the overall strain, we don't need that here
	dplyr::filter(subtype != "") |>
	dplyr::mutate(dplyr::across(dplyr::where(is.factor), forcats::fct_drop))

## Create a data frame to hold our distance measurements ####
# Note that we need separate holders for the H1 strains and for the H3 strains
# So first we need to make data frames that contain all combinations of two
# strains from the input data.
distance_data_nested <-
	strain_names |>
	tidyr::nest(sequence_data = -subtype) |>
	# Add the cohort data as a column
	tibble::add_column(
		cohort_data = cohort_data |>
			dplyr::group_by(strain_type) |>
			dplyr::group_split() |>
			as.list()
	)

# Temporal distance ####
# The temporal distance is the easiest to calculate. This distance is defined as
# the absolute difference between the years of isolation between two strains.
# Now we calculate the pairwise distances using the temporal method. Note that
# the temporal_dist() function is from the utils.R script.
# It extracts the isolation year from the strain names and then computes the
# absolute difference.
# Now for H1 and H3 separately, we want to calculate the pairwise distances.
# Outer map iterates over all subtypes (right now just H1 and H3)
distance_data_nested$d_temporal <-
	purrr::map(
		distance_data_nested$sequence_data,
		# Inner map iterates over all combos
		\(d) d |>
			dplyr::pull(analysis_name) |>
			as.character() |>
			temporal_dist()
	)

# p-Epitope distance ####
# Before we can compute the p-epitope distance, we have to align the
# sequences correctly, which involves a few steps.
preprocess_sequences_for_pepi <- function(sequence_data, subtypes) {
	# Standardize subtypes
	subtype <- subtypes |> as.character() |> toupper() |> substr(1, 2)

	# Pull out raw sequences
	seqs_raw <- sequence_data$protein_sequence

	# Get the end of the signal peptide depending on what subtype we're on
	signal_peptipe_end <-
		dplyr::case_match(
			subtype,
			"H1" ~ 18,
			"H3" ~ 17,
			.default = NA_real_
		)

	# If there are any unregonized subtypes, stop running the code!
	if (any(is.na(signal_peptipe_end))) {
		rlang::abort("Unrecognized subtype!")
	}

	# Remove the signal peptide. For the one sequence that is not full length,
	# we know we can skip this. Otherwise remove the signal peptide.
	seqs_processed <-
		dplyr::if_else(
			sequence_data$strain_name == "MI/85",
			seqs_raw,
			stringr::str_sub(
				seqs_raw,
				signal_peptipe_end,
				-1
			)
		)

	# For H1, we also need to add a gap to a specific site if the length is
	# 548 after removing the signal peptide.
	if (subtype == "H1") {
		seqs_processed <- dplyr::if_else(
			nchar(seqs_processed) == 548,
			# This regular expression adds a dash at spot 130.
			gsub('^(.{129})(.*)$', '\\1-\\2', seqs_processed),
			seqs_processed
		)
	}

	# Now give the sequences the correct names.
	names(seqs_processed) <- sequence_data$strain_name

	return(seqs_processed)
}

distance_data_nested$aligned_seqs <-
	purrr::map2(
		distance_data_nested$sequence_data,
		distance_data_nested$subtype,
		preprocess_sequences_for_pepi
	)

distance_data_nested$d_pepitope <-
	purrr::map2(
		distance_data_nested$aligned_seqs,
		distance_data_nested$subtype,
		dist_pepi
	)

# Cartographic distance ####
## Titer matrix creation ####
# First we need to process the cohort data to make the titer table.
# For this project, we'll just use all of the data that we have, but in the
# future there are a lot of sensitivity tests and analyses we want to do.
distance_data_nested$titer_matrix <-
	purrr::map(
		distance_data_nested$cohort_data,
		\(d) d |>
			# First pick the variables we'll need, which are just the serum IDs
			# (individual measurement IDs), antigen IDs (which strain the assay
			# was against) and finally the post-vaccination titer measurement.
			dplyr::select(unique_id, posttiter, strain_name) |>
			# We'll treat the subjects differently for each person-year of data, treating
			# individuals as independent -- MDS does not accomodate a correlation
			# structure for measurements by individuals.
			# First we need to pivot, so we have the person-year ID as rows and the
			# strains as columns.
			tidyr::pivot_wider(names_from = strain_name, values_from = posttiter) |>
			# Racmacs cannot accurately place n-D coordinates for individuals
			# who have less than n + 1 non-censored non-missing measurements. So
			# we need to remove individuals who have fewer. Because this is a rowwise
			# operation it is a bit annoying in dplyr.
			# This removes individuals who do not have at least 3 non-missing
			# observations which are 10 or higher.
			dplyr::rowwise() |>
			dplyr::filter(
				sum(tidyr::replace_na(dplyr::c_across(-unique_id), 0) >= 10) >= 3
			) |>
			dplyr::ungroup() |>
			# Now, recode missing and censored values into Racmacs format.
			dplyr::mutate(
				# First we have to turn all the numbers into character vectors to format
				# the data how Racmacs wants it.
				dplyr::across(
					dplyr::everything(),
					as.character
				),
				# Now replace the missing values with a "*", which is what Racmacs wants.
				dplyr::across(
					dplyr::everything(),
					\(x) tidyr::replace_na(x, "*")
				),
				# Now replace the censored values <LOD with "<10".
				dplyr::across(
					dplyr::everything(),
					\(x) dplyr::if_else(x == "5", "<10", x)
				)
			) |>
			# Convert the person-year IDs to an actual rowname instead of a column
			# It has to be quoted because for some reason this is the only tidyverse
			# function set without data masking?
			tibble::column_to_rownames("unique_id") |>
			# Convert to matrix to improve speed
			as.matrix() |>
			# And transpose since Racmacs wants individuals ("sera") as columns
			t()
	)

## Map optimization ####
# Now with the underconstrained sera removed we can fit the maps and run the
# optimization algorithm. Note we shouldn't use dimensional annealing because we
# would have underconstrained sera in higher dimensions.
# We'll run each map through N_OPT_EACH rounds of optimization, and we'll
# restart the optimization process N_OPT_RUNS times, which can help us rule
# out bad starting conditions. Then from the results we can choose the best
# overall map.
make_optimized_maps <- function(titer_matrix) {
	optimized_maps <- purrr::map(
		1:N_OPT_RUNS,
		\(idx) {
			# We'll get a warning about unstable results, but this is because we have
			# a lot of sera. So suppress that warning.
			suppressWarnings(
				Racmacs::make.acmap(
					titer_table = titer_matrix,
					ag_names = rownames(titer_matrix),
					sr_names = colnames(titer_matrix),
					number_of_dimensions = 2L,
					number_of_optimizations = N_OPT_EACH,
					minimum_column_basis = "none",
					options = RacOptimizer.options(num_cores = N_CORES),
					verbose = FALSE
				)
			)
		},
		.progress = "Iterating map optimization"
	)

	return(optimized_maps)
}
# Since this is the most time consuming step, we'll also record how long it
# took on our machine.
{
	start <- Sys.time()
	all_maps <-
		purrr::map(
			distance_data_nested$titer_matrix,
			make_optimized_maps
		)
	end <- Sys.time()
	duration <- difftime(end, start)
}
# Save result to file since it will be large
purrr::walk2(
	all_maps,
	paste0(here::here(
		"Results", "Large-Files"), "/", distance_data_nested$subtype,
		"-map-optimizations.Rds"
	),
	\(x, f) readr::write_rds(f, file = f, compress = "xz")
)
readr::write_rds(
	duration,
	here::here("Results", "Times", "map-optimization-time.Rds")
)

## Get best maps ####
# Now we want to get only the lowest stress maps. Getting the best maps from
# each of the "inner" optimization runs is easy cause there's a function
# for it, but we have to extract and get the minimum stress for the "outer"
# optimization runs.
# This function takes in a list of acmaps, chooses the best optimization of
# each one, and computes the stress of the best map.
get_best_map_from_optim_run <- function(map_list) {
	# Get the best optimization (lowest stress) from the list
	best_map_per_iter <- purrr::map(map_list, Racmacs::keepBestOptimization)

	# Now compute the stress of each of the best acmaps
	stresses <- purrr::map(best_map_per_iter, Racmacs::mapStress)

	# Return the best map with the lowest stress -- the best of the best
	out <- best_map_per_iter[[which.min(stresses)]]
	return(out)
}

# Get the best overall acmap for each subtype
best_acmaps <-
	purrr::map(
		all_maps,
		get_best_map_from_optim_run
	)
# And write to file
readr::write_rds(
	best_acmaps,
	here::here("Results", "Data", "Best-Maps.Rds")
)

## Get distance matrices from acmaps ####
distance_data_nested$d_cartographic <-
	purrr::map(
		best_acmaps,
		racmaps_map_to_distances
	)

# Save distances to file ####
# Finally, we only want to export the three distance metrics and the subtype.
# So we need to clean the data up a bit.
tidy_distance_data <-
	distance_data_nested |>
	dplyr::select(subtype, dplyr::starts_with("d_")) |>
	# Put all of the three distance matrices into tidy data format instead of
	# the matrix format they are currently in
	# tidy_dist_mat() is a function from utils.R that does just that so we
	# apply it to each of our columns
	dplyr::mutate(
		dplyr::across(
			dplyr::starts_with("d_"),
			\(x) purrr::map(x, tidy_dist_mat)
		)
	) |>
	# First we want to pivot our three distance columns longer, so we have
	# one column of method names and one column of distance matrices
	tidyr::pivot_longer(
		cols = dplyr::starts_with("d_"),
		names_to = "method",
		names_prefix = "d_",
		values_to = "distance_matrix"
	) |>
	# Now we can unnest the data, meaning we expand the list column back into
	# regular columns, and we are done!
	tidyr::unnest(
		cols = distance_matrix
	)

# Now we write the distance data to file!
write_csv_and_rds(
	tidy_distance_data,
	here::here("Results", "Data", "distance-data")
)

# END OF FILE ####
