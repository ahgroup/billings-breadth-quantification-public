###
# Data Cleaning For Model Fitting
# Zane Billings
# 2022-02-21
# Get
###

# SETUP ####

# Declare dependencies
box::use(
	readr,
	tidyr,
	dplyr
)

# Data cleaning ####

# Read in the cleaned data set
clean_data <-
	readr::read_rds(
		here::here("data", "processed", "distance_data.rds")
	)

# Basic processing to clean up some of Amanda's data
dat <-
	clean_data |>
	# Remove unnecessary text from the vaccine_type column
	dplyr::mutate(
		vaccine_type = stringr::str_remove(vaccine_type, "_vaccine_fullname$") |>
			stringr::str_to_upper() |>
			factor()
	) |>
	dplyr::filter(
		# We only care about the rows where the strains are comparable -- it does
		# not make sense to compare H3N2 responses to the H1N1 vaccine.
		vaccine_type == strain_type,
		# Get only the Ag distance methods that Amanda used.
		method %in% c("cart_2d_post", "p_epi", "year"),
		# We decided to only use seasons from before 2019 where there is a wide
		# historical panel to choose from
		season < 2019
	) |>
	# Drop any unused factor levels
	dplyr::mutate(
		dplyr::across(tidyselect:::where(is.factor), forcats::fct_drop)
	) |>
	# These are the raw scale variables, but we already have the log2 scale
	# so we don't need these
	dplyr::select(-pretiter, -postiter) |>
	# Pivot the data so we can filter by the outcome of interest.
	tidyr::pivot_longer(
		cols = c(postvactiter, titerincrease),
		names_to = "outcome",
		values_to = "y"
	) |>
	# Normalize the antigenic distance measurements within vaccine group
	dplyr::group_by(vaccine_fullname, method) |>
	dplyr::mutate(
		norm_dist = distance / max(distance)
	) |>
	dplyr::ungroup() |>
	# Recode the dose to be an integer since Stan cannot handle categorical
	# variables. Note that SD = 1 and HD = 2!
	dplyr::mutate(
		dose = dose |>
			as.character() |>
			factor(levels = c("SD", "HD")) |>
			as.integer()
	)

# In order to fit multiple models, we can create nested data frames. We'll fit
# separate models for each vaccine strain/outcome/distance method combination.
# The "model level" variables should be outside of the nested data and the
# "individual level" and "group level" variables should be inside of the
# nested data.
dat2 <-
	dat |>
	# Remove variables that we don't really need -- they don't provide any
	# extra information so we don't need them in the modeling data.
	# We'll also do some renaming here to make my life easier.
	dplyr::select(
		-strain_type, -strain_short, -vac_short, -vaccine_type,
		vaccine_strain = vaccine_fullname,
		assay_strain = strains_fullname
	) |>
	# Create a numeric ID variable -- this has to be GROUPWISE PER MODEL or else
	# Stan won't know what to do, because it uses the variable as a literal
	# index.
	dplyr::group_by(vaccine_strain, method, outcome) |>
	dplyr::mutate(
		id = uniq_id |>
			factor() |>
			forcats::fct_inorder() |>
			as.integer()
	) |>
	dplyr::ungroup()

# Next we'll make a lookup table of ID conversions before we get rid of the
# original uniq_id column, so we can join back if we need to.
id_lut <-
	dat2 |>
	dplyr::select(vaccine_strain, method, outcome, id, uniq_id) |>
	dplyr::distinct()

readr::write_rds(
	id_lut,
	here::here("Data", "Processed", "Modeling-Data", "id-lut.Rds")
)

# Now we can get rid of the uniq_id and nest the data.
dat_nested <-
	dat2 |>
	dplyr::select(-uniq_id) |>
	tidyr::nest(
		dat = c(id, season, assay_strain, y, norm_dist, distance, prevactiter,
						age, dose)
	)

# Eventually the data will need to be turned into a list to pass to the
# Stan models, but it is easier to do further processing operations on
# the data frame, so that will be done ad-hoc before modeling.

readr::write_rds(
	dat_nested,
	here::here("Data", "Processed", "Modeling-Data", "model-data.Rds")
)
