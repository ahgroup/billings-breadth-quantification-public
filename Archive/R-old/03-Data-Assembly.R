###
# Data Assembly
# 2024-08-12
# Zane Billings
# Combines the imported cohort data with the generated antigenic distance data.
###

# Setup ####
box::use(
	readr,
	here,
	dplyr
)

source(here::here("R", "functions", "utils.R"))

# Load the cohort data
cohort_data <- here::here("data", "processed", "analysis-data.Rds") |>
	readr::read_rds()

# Load the distance data
dist_data <- here::here("results", "data", "distance-data.Rds") |>
	readr::read_rds()

# Data joining ####
# For each assay measurement in the cohort data, we want to add the distance
# measurement corresponding to the distance between the vaccine strain and
# the assay strain. That means each one record in the cohort data should have
# EXACTLY 3 matches in the distance data. Most rows in the distance data
# won't get used, and that's OK. However, some rows in the distance
# data will be used multiple times, making this a many-to-many merge.
# What we need to check at the end to make sure it worked right is that
# the number of rows in the joined data is exactly 3x as much as the
# original cohort data.

# The first thing we need to do is a bit of data manipulation to the cohort
# data, because right now the vaccine strain is stored in two separate columns
# in the cohort data.
# After we solve that issue we can do the join.
joined_data <-
	cohort_data |>
	# First create a new variable for the vaccine name that matches the current
	# strain type of the assay that a given record represents.
	dplyr::mutate(
		vaccine_strain = dplyr::if_else(
			strain_type == "H1N1",
			h1n1_vaccine_strain,
			h3n2_vaccine_strain
		)
	) |>
	# Now join to the distance data.
	dplyr::left_join(
		dist_data,
		by = c(
			"strain_type" = "subtype",
			"vaccine_strain" = "strain1",
			"strain_name" = "strain2"
		),
		relationship = "many-to-many"
	)

if (nrow(joined_data) != 3 * nrow(cohort_data)) {
	rlang::abort('Something about the join has gone wrong!')
}

# Write to file.
write_csv_and_rds(
	joined_data,
	here::here("data", "processed", "joined-data")
)


# END OF FILE ####
