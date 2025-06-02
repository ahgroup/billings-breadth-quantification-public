###
# Data Sampler for Methods
# Zane
# 2024-01-17
# Get a randomized subsample of data from Ted's cohort to use for
# the simple methods example
###

# Load full data and get only one season
dat <-
	here::here("data", "processed", "joined-data.Rds") |>
	readr::read_rds() |>
	dplyr::filter(season == "2016 - 2017")

# Get a random subsample
set.seed(2351)
subjects <- sample(dat$subject_id, size = 100L, replace = FALSE)
dat_use <-
	dat |>
	dplyr::filter(
		subject_id %in% subjects,
		strain_type == "H3N2",
		method == "cartographic"
	) |>
	dplyr::transmute(
		id = subject_id |>
			forcats::fct_inorder() |>
			as.integer() |>
			pad_numbers(),
		assay_strain = strain_name,
		pretiter, posttiter,
		distance = d
	)

# Save to file
qs2::qd_save(
	dat_use,
	here::here("products", "simple-methods", "methods-ex-dat.qs2")
)
