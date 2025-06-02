###
# Data import
# Zane Billings
# 2024-08-12
# Loads the following data sources from file:
#  - UGAFluVac cohort data
#  - Influenza strain accession numbers and names
# Then, we clean and process each dataset as needed,
# and save it in Rds and csv formats.
###

# Setup ####
box::use(
	here,
	readr,
	readxl,
	dplyr,
	hgp,
	forcats,
	purrr,
	tidyr,
	janitor
)

source(here::here("R", "functions", "utils.R"))

# UGAFluVac data cleaning ####

## Analysis dataset ####
# Load in the data -- the input file is clean_data.Rds, which was downloaded
# from the UGAFluVac repo on 2024-08-11.
raw_data <-
	here::here("data", "raw", "clean_data.Rds") |>
	readr::read_rds()

# Subset the data to get the subcohort used for our study, and only keep the
# variables we will actually need to use.
analysis_data <-
	raw_data |>
	# Inclusion/exclusion criteria:
	dplyr::filter(
		# 1. only individuals who got intramuscular vaccine will be used
		dose != "Intradermal",
		# 2. only use individuals who enrolled in the study before 2018, since they
		# have way more heterologous assays.
		season <= "2017 - 2018",
		# 3. only use flu A data
		strain_type %in% c("H1N1", "H3N2")
	) |>
	# Select the variables we will use in the analysis -- we do not need any
	# demographic data for our analytic method, but we need a few things for
	# the demographic information table.
	dplyr::select(
		row_id, unique_id = uniq_id, subject_id, season, study, id, dose,
		age, birth_year, sex, race = race_factor, ethnicity,
		h1n1_vaccine_strain, h3n2_vaccine_strain, strain_name, strain_type,
		pretiter, posttiter, fold_change, seroconversion
	) |>
	# Final tweaks
	dplyr::mutate(
		# Over all of the factor variables, drop the levels that have zero records
		dplyr::across(dplyr::where(is.factor), forcats::fct_drop),
		# Convert all of the strain names to short name format and order correctly
		dplyr::across(
			c(h1n1_vaccine_strain, h3n2_vaccine_strain, strain_name),
			hgp::replace_strain_names
		)
	)

# Get the row_ids that are included in the analysis dataset
analysis_data_rows_used <- dplyr::pull(analysis_data, row_id)

## Metadata reporting dataset ####
# Make the dataset with everything needed for the human subjects demographics
# metadata
reporting_data <- raw_data |>
	# First, make sure we only include individuals who were included in the
	# analysis subset
	dplyr::filter(row_id %in% analysis_data_rows_used) |>
	# Select the variables needed for the CIVICs human subjects demographics
	# data standard construction
	dplyr::select(
		subject_id, season, sex, gender, age, birth_year, study, ethnicity,
		race = race_civics_standard, dose
	) |>
	dplyr::distinct()

# Import virus accession and sequence data ####
sequence_data_raw <-
	here::here("data", "raw", "strain-sequences.xlsx") |>
	readxl::read_excel()

sequence_data <-
	sequence_data_raw |>
	# clean up the column names
	janitor::clean_names() |>
	# Select only the relevant (protein) columns
	dplyr::select(
		strain_name, protein_sequence, protein_sequence_source,
		protein_sequence_accession_number
	) |>
	# Combine the accession number columns
	tidyr::unite(
		col = "accession_no",
		c(protein_sequence_source, protein_sequence_accession_number),
		sep = ": "
	) |>
	# Cleanup
	dplyr::mutate(
		# Convert strain name to short format
		strain_name = hgp::replace_strain_names(strain_name, from = "full"),
		# Remove any newlines from the sequences
		protein_sequence = stringr::str_remove_all(protein_sequence, "\\s"),
		# Remove any excess whitespace from all columns
		dplyr::across(
			dplyr::where(is.character),
			stringr::str_squish
		),
		# Add a column with the sequence length
		sequence_length = nchar(protein_sequence),
		# Add an indicator showing that the MI/85 strain is not full length
		full_length = ifelse(strain_name == "MI/85", FALSE, TRUE)
	)

# Write the datasets to file ####
purrr::walk2(
	list(analysis_data, reporting_data, sequence_data),
	list(
		here::here("data", "processed", "analysis-data"),
		here::here("data", "processed", "reporting-data"),
		here::here("data", "processed", "sequence-data")
	),
	\(x, f) write_csv_and_rds(x, f)
)

# END OF FILE ####
