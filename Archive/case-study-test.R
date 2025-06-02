###
# Test analysis of BQ methods on HD/SD
# Zane
# 2025-05-17
####

cli::cli_h1("Starting setup")

source(here::here("R", "Data-Cleaning.R"))
source(here::here("R", "Utils.R"))
source(here::here("R", "Stats.R"))
library(brms)
library(cmdstanr)

file_location <- here::here("Archive", "test-case-study-fits.Rds")

dat <- targets::tar_read("raw_UGAFluVac_data")
raw_antigenic_distance_data <- targets::tar_read("raw_antigenic_distance_data")

make_comparison_data <- function(raw_cohort_data, raw_antigenic_distance_data) {
	cleaned_cohort_data <-
		dat |>
		dplyr::filter(
			dose %in% c("SD", "HD"),
			is.finite(titerincrease),
			strain_type == "H1N1",
			season <= "2017 - 2018",
			age >= 65
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

	# First join the antigenic distances to the cohort data
	joined_data <- join_antigenic_distance_to_cohort_data(
		cleaned_cohort_data,
		raw_antigenic_distance_data
	)

	# Now do some model transformations and cleanup
	comparison_data <-
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

	nested_data <-
		comparison_data |>
		dplyr::bind_rows(dplyr::mutate(comparison_data, season = "Overall")) |>
		dplyr::select(season, y, c, y2, d_norm, subject_id, row_id, dose, metric) |>
		tidyr::nest(this_data = -c(season, dose, metric)) |>
		dplyr::arrange(season, dose)

	return(nested_data)
}

fit_model <- function(model_data, ...) {
	model_fit <- brms::brm(
		formula = y | cens(c, y2) ~ 1 + d_norm + (1 + d_norm | subject_id),
		family = "gaussian",
		data = model_data,
		prior = c(
			brms::prior(
				"student_t(3, 0, 3)", class = "Intercept"
			),
			# Prior for slopes
			brms::prior(
				"student_t(3, 0, 1)", class = "b"
			),
			# Prior for variances, we need two. One for covariances and one for global
			# variance parameter
			brms::prior(
				"student_t(3, 0, 1)", class = "sigma", lb = 0
			),
			brms::prior(
				"student_t(3, 0, 1)", class = "sd", lb = 0
			),
			# Priors for cholesky factors
			brms::prior(
				"lkj_corr_cholesky(2)", class = "L"
			)
		),
		seed = 7317L,
		#sample_prior = "yes",
		chains = fitting_arguments$chains,
		cores = fitting_arguments$chains,
		warmup = fitting_arguments$iter_warmup,
		iter = fitting_arguments$iter_warmup + fitting_arguments$iter_sampling,
		algorithm = "sampling",
		control = list(
			adapt_delta = fitting_arguments$adapt_delta,
			max_treedepth = fitting_arguments$max_treedepth
		),
		backend = "cmdstanr"
	)

	return(model_fit)
}

if (file.exists(file_location)) {
	cli::cli_h1("Reading fitted models from file")
	nested_data <- make_comparison_data(dat, raw_antigenic_distance_data)
	res <- readr::read_rds(file_location)
} else {
	nested_data <- make_comparison_data(dat, raw_antigenic_distance_data)

	fitting_arguments <-
		sampling_arguments <- list(
			chains = 12L,
			iter_warmup = 250L,
			iter_sampling = 250L,
			max_treedepth = 12L,
			adapt_delta = 0.99
		)

	res <- list()
	n <- nrow(nested_data)
	cli::cli_h1("Starting model fitting")

	for (i in 1:n) {
		this_season <- nested_data$season[[i]]
		this_dose <- nested_data$dose[[i]]
		cli::cli_h2(paste(
			"Fitting model", i, "of", n, "for: season", this_season, "and dose",
			this_dose
		))
		season_data <- nested_data$this_data[[i]]
		season_model_fit <- fit_model(season_data)
		res[[i]] <- season_model_fit
		rm(season_model_fit, this_season, season_data)
		gc(verbose = FALSE)
	}

	readr::write_rds(
		res,
		here::here("Archive", "test-case-study-fits.Rds")
	)
}

cli::cli_h1("Computing preds and stats")

calculate_new_metrics <- function(epreds, ...) {
	auc <- AUC(epreds) |> ggdist::mean_hdci()
	pat <- prop_above_threshold(epreds, "log_posttiter") |> ggdist::mean_hdci()
	int <- intercept(epreds) |> ggdist::mean_hdci()

	out <- dplyr::bind_rows(
		auc = auc,
		pat = pat,
		int = int,
		.id = "stat"
	)

	return(out)
}

stuff <-
	nested_data |>
	# Calculate the predictions from each model
	dplyr::mutate(
		model_fit = res,
		model_epreds = purrr::map(
			model_fit,
			\(x) tidybayes::epred_draws(
				x,
				newdata = tibble::tibble(d_norm = seq(0, 1, 0.01)),
				ndraws = brms::ndraws(x),
				re_formula = NA
			) |>
				dplyr::rename(x = d_norm)
		)
	) |>
	# Get the metrics
	dplyr::mutate(stats = purrr::map(model_epreds, calculate_new_metrics))

summarized_stuff <-
	stuff |>
	dplyr::select(season, dose, metric, stats) |>
	tidyr::unnest(stats)

tbl_data <-
	summarized_stuff |>
	dplyr::transmute(
		season, dose, metric, stat,
		est = paste0(
			formatC(y, format = "f", digits = 2),
			" (",
			formatC(ymin, format = "f", digits = 2),
			", ",
			formatC(ymax, format = "f", digits = 2),
			")"
		)
	) |>
	tidyr::pivot_wider(names_from = stat, values_from = est)


summarized_stuff

case_study_plot <-
	summarized_stuff |>
		dplyr::filter(metric == "Cartographic") |>
	dplyr::mutate(
		season = gsub("([0-9]{4})\\s*-\\s*[0-9]{2}([0-9]{2})", "\\1/\\2", season),
		season = factor(season, ordered = TRUE),
		stat = factor(
			stat,
			levels = c("int", "pat", "auc"),
			labels = c(
				"Magnitude\n(Intercept)",
				"Breadth\n(Prop. above threshold)",
				"Total Strength\n(AUC)"
			)
		)
	) |>
	ggplot2::ggplot() +
	ggplot2::aes(
		x = season,
		y = y,
		ymin = ymin,
		ymax = ymax,
		color = dose,
		shape = dose
	) +
	ggplot2::geom_pointrange(
		position = ggplot2::position_dodge(width = 0.5)
	) +
	ggplot2::facet_wrap(~stat, scales = "free_y") +
	ggplot2::labs(
		x = "Season",
		y = "Statistic",
		color = "Fluzone dose",
		shape = "Fluzone dose"
	) +
	ggplot2::scale_color_viridis_d(begin = 0.2, end = 0.8) +
	hgp::theme_ms() +
	ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))

ggplot2::ggsave(
	filename = here::here("Archive", "test-case-study-plot.png"),
	width = 16,
	height = 8
)

