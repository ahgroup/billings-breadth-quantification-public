###
# Data Cleaning
# Zane Billings
# 2025-04-28
# New (targetized) version of data cleaning script that imports the UGAFluVac
# cohort data along with the antigenic distance data from the standalone repo.
# All of the distance calculation code that was in this repo was moved to the
# common distance calculation repo.
###

prep_cohort_data <- function(raw_cohort_data) {
	# First step is basic filtering with our study's inclusion/exclusion criteria
	dat_filtered <-
		raw_cohort_data |>
		dplyr::filter(
			# Only look at SD individuals here
			dose == "SD",
			# Remove individuals with infinite titerincrease (1 person who has a
			# pretiter of 0, probably a data entry error)
			is.finite(titerincrease),
			# Only use data from H1N1 in this analysis.
			strain_type == "H1N1",
			# Only the useful seasons
			season <= "2017 - 2018"
		) |>
		# Select the variables we'll need
		dplyr::select(
			row_id, unique_id = uniq_id, subject_id, season, year, study, dose,
			age, birth_year, sex, race = race_collapsed,
			vaccine_strain = h1n1_vaccine_strain, assay_strain = strain_name,
			pretiter, posttiter, fold_change, seroconversion, seroprotection,
			# We need these for the CIVICs reporting part
			gender, bmi, ethnicity, race_civics_standard
		) |>
		# Final tweaks
		dplyr::mutate(
			# there are currently two people who have multiple different race_collapsed
			# variables so they need to be standardized, we'll do that manually.
			race = dplyr::if_else(
				subject_id %in% c("0045", "0072"),
				"Hispanic or Latino",
				race
			) |>
				forcats::fct_infreq() |>
				forcats::fct_na_value_to_level("Unknown"),
			# Fix the typo where the race says "American American"
			race = forcats::fct_recode(
				race, "Black or African American" = "Black or American American"
			),
			# Fix the weird Shangdong typo
			assay_strain = replace(
				as.character(assay_strain),
				assay_strain == "H3N2-Shangdong-1993",
				"H3N2-Shandong-1993"
			),
			# Now do the strain name replacement
			dplyr::across(c(vaccine_strain, assay_strain), hgp::replace_strain_names),
			# Log outcome variables
			log_pretiter = hgp::hai_to_log_scale(pretiter),
			log_posttiter = hgp::hai_to_log_scale(posttiter),
			titer_increase = log2(fold_change)
		)

	return(dat_filtered)
}

join_antigenic_distance_to_cohort_data <- function(
		cleaned_cohort_data, raw_antigenic_distance_data
) {
	# Nest the antigenic distance data to make joining easier -- this allows us
	# to avoid making sure a many-to-many merge is working correctly.
	nested_agdist_data <-
		raw_antigenic_distance_data |>
		dplyr::mutate(
			# Clean up the metric variable
			# This version has more, uncomment if we need the other ones and then
			# comment the below part.
			# metric = factor(
			# 	as.character(metric),
			# 	levels = c(
			# 		"cartographic", "dominant_pepitope", "p_all_epitope", "grantham",
			# 		"flu_sub", "hamming", "absolute_temporal"
			# 	),
			# 	labels = c(
			# 		"Cartographic", "p-Epitope", "p-All-Epitope", "Grantham", "FLU Sub",
			# 		"Hamming", "Temporal"
			# 	)
			metric = factor(
				as.character(metric),
				levels = c(
					"cartographic", "dominant_pepitope", "absolute_temporal"
				),
				labels = c(
					"Cartographic", "p-Epitope", "Temporal"
				)
			)
		) |>
		# And drop the metrics we don't want
		tidyr::drop_na(metric) |>
		tidyr::nest(antigenic_distances = c(metric, d)) |>
		dplyr::select(-type_subtype)

	# Now join the antigenic distances
	dat_joined <-
		cleaned_cohort_data |>
		# Now join the antigenic distances
		dplyr::left_join(
			nested_agdist_data,
			by = c("vaccine_strain" = "Strain1", "assay_strain" = "Strain2"),
			relationship = "many-to-one"
		) |>
		# Now unnest the antigenic distance data
		tidyr::unnest(antigenic_distances)

	return(dat_joined)
}

# Now we need to make a few additional changes before we can pass the
# data to our models.
create_model_data <- function(
		cleaned_cohort_data, raw_antigenic_distance_data
	) {
	# First join the antigenic distances to the cohort data
	joined_data <- join_antigenic_distance_to_cohort_data(
		cleaned_cohort_data,
		raw_antigenic_distance_data
	)

	# Now do some model transformations and cleanup
	dat_model <-
		joined_data |>
		dplyr::mutate(
			# Remove the extra levels from the vaccine strain variable
			vaccine_strain = forcats::fct_drop(vaccine_strain),
			# We need to make versions of the time and birth year variables that are
			# closer to being scale-free than the current versions -- models that
			# have those numbers in the thousands have worse conditioning problems that
			# can lead to numerical issues in an otherwise fine model. So we'll scale
			# the year variable by subtracting 2013, the first year, and we'll scale
			# the birth_year variable with minmax transformation.
			year_c = year - 2013,
			birth_year_c = minmax(birth_year),
			# We'll also minmax scale the age.
			age_c = minmax(age),
			# Finally we'll create indicator variables for the categorical covariates
			# that will go in the model
			# Indicator variable for sex, 0 is Male, 1 is female.
			sex_i = as.numeric(sex == "Female"),
			# Indicator variable race/ethnicity, 0 is non-Hispanic White, 1 is Other
			race_i = as.numeric(race != "White"),
			# Short version of the season for plotting
			season_short = gsub(" - 20([0-9]{2})", "/\\1", as.character(season)) |>
				factor(ordered = TRUE)
		) |>
		hgp::format_hai_data(post_titer = "log_posttiter") |>
		# We also want to minmax scale the antigenic distances, but we want to
		# do it separately by metric.
		dplyr::group_by(metric) |>
		dplyr::mutate(d_norm = minmax(d)) |>
		dplyr::ungroup() |>
		droplevels.data.frame()

	return(dat_model)
}

# This function creates the formatted NIH CIVICs reporting data that we have
# to submit to the CIVICs hub
create_civics_reporting_data <- function(prepped_cohort_data) {
	dat_used <- prepped_cohort_data |>
		dplyr::rename(id = subject_id) |>
		dplyr::distinct(
			season, study, age, birth_year, gender, sex, bmi, dose, ethnicity,
			race_civics_standard, year, id
		)

	# Get the minimum and maximum age in the study
	dat_ages <- dat_used |>
		dplyr::group_by(id, study) |>
		dplyr::mutate(
			Min_Age = min(age),
			Max_Age = max(age)
		) |>
		dplyr::ungroup()

	dat_reporting <-
		dat_used |>
		dplyr::transmute(
			Study_Code = "TBD",
			Subject_ID = id,
			Cohort_ID = paste(study, year, dose, sep = "_"),
			Sex_Assigned_at_Birth = ifelse(
				is.na(sex), "Unknown", as.character(sex)
			),
			Gender = gender,
			Min_Age = dat_ages$Min_Age,
			Max_Age = dat_ages$Max_Age,
			Subject_Age_Unit = "Years",
			Birth_Year = birth_year,
			Subject_Age_Event = "Age at enrollment",
			Subject_Phenotype = "Not Collected",
			Subject_Location = ifelse(
				study == "UGA", "GA", study
			),
			Ethnicity = ifelse(
				is.na(ethnicity),
				"Unknown",
				ethnicity
			),
			Race = race_civics_standard,
			Subject_Description = "Not Provided"
		) |>
		dplyr::mutate(dplyr::across(dplyr::where(is.character), as.factor))

	return(dat_reporting)
}

# Save the final model data to file
write_model_data_to_file <- function(model_data, filename_base) {
	filenames <- paste0(filename_base, c(".qs2", ".csv"))
	qs2::qs_save(
		model_data,
		file = filenames[[1]]
	)
	readr::write_csv(
		model_data,
		file = filenames[[2]]
	)

	return(filenames)
}
