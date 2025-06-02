###
# Descriptive analyses
# Zane Billings
# 2025-04-28
# Demographics tables and other summary statistics and plots needed to
# understand the dataset.
###

make_assay_counts_table <- function(model_data, file_path) {
	distinct_data <-
		model_data |>
		dplyr::distinct(subject_id, season_short, study, assay_strain)

	distinct_data_w_overall <-
		dplyr::bind_rows(
			distinct_data,
			dplyr::mutate(distinct_data, season_short = "Overall"),
			dplyr::mutate(distinct_data, study = "Overall"),
			dplyr::mutate(distinct_data, season_short = "Overall", study = "Overall")
		)

	tbl_data <-
		distinct_data_w_overall |>
		dplyr::summarize(
			n = dplyr::n(),
			.by = c(season_short, study)
		) |>
		dplyr::mutate(
			p = n / max(n) * 100,
			p_f = formatC(p, format = "f", digits = 0),
			#n_a = stringr::str_pad(n, 5, side = "left", pad = " "),
			n_f = paste0(n, "\n(", p_f, "%)")
		) |>
		dplyr::transmute(season_short, study, n = n_f) |>
		tidyr::pivot_wider(
			names_from = season_short,
			values_from = n,
			names_glue = "Season_{season_short}"
		) |>
		dplyr::mutate(
			dplyr::across(
				.cols = dplyr::starts_with("Season"),
				.fns = \(x) tidyr::replace_na(as.character(x), "----")
			),
			# dplyr::across(
			# 	.cols = c(`Season_2013/14`, `Season_2014/15`),
			# 	.fns = \(x) ifelse(study == "FL" & dose == "HD", 0, x)
			# ),
			dplyr::across(
				.cols = dplyr::starts_with("Season"),
				.fns = \(x) stringr::str_pad(x, 4, "left", " ")
			),
			dplyr::across(
				.cols = dplyr::starts_with("Season"),
				.fns = \(x) stringr::str_replace(x, "----", " — ")
			)
		)

	assay_tbl <-
		tbl_data |>
		dplyr::rename(
			"Study Site" = study
		) |>
		flextable::flextable() |>
		flextable::separate_header() |>
		flextable::merge_v(j = 1) |>
		flextable::hline(i = 1, j = 2:7, part = "header") |>
		flextable::align(j = 2:7, align = "center") |>
		flextable::fix_border_issues() |>
		fit_flextable_to_page()

	save_file_as_qs2(
		assay_tbl,
		file_path
	)

	return(file_path)
}

make_person_counts_table <- function(model_data, file_path) {
	distinct_data <-
		model_data |>
		dplyr::distinct(subject_id, season_short, study)

	distinct_data_w_overall <-
		dplyr::bind_rows(
			distinct_data,
			dplyr::mutate(distinct_data, season_short = "Overall"),
			dplyr::mutate(distinct_data, study = "Overall"),
			dplyr::mutate(distinct_data, season_short = "Overall", study = "Overall")
		)

	tbl_data <-
		distinct_data_w_overall |>
		dplyr::summarize(
			n = dplyr::n(),
			.by = c(season_short, study)
		) |>
		dplyr::mutate(
			p = n / max(n) * 100,
			p_f = formatC(p, format = "f", digits = 0),
			#n_a = stringr::str_pad(n, 5, side = "left", pad = " "),
			n_f = paste0(n, "\n(", p_f, "%)")
		) |>
		dplyr::transmute(season_short, study, n = n_f) |>
		tidyr::pivot_wider(
			names_from = season_short,
			values_from = n,
			names_glue = "Season_{season_short}"
		) |>
		dplyr::mutate(
			dplyr::across(
				.cols = dplyr::starts_with("Season"),
				.fns = \(x) tidyr::replace_na(as.character(x), "----")
			),
			# dplyr::across(
			# 	.cols = c(`Season_2013/14`, `Season_2014/15`),
			# 	.fns = \(x) ifelse(study == "FL" & dose == "HD", 0, x)
			# ),
			dplyr::across(
				.cols = dplyr::starts_with("Season"),
				.fns = \(x) stringr::str_pad(x, 4, "left", " ")
			),
			dplyr::across(
				.cols = dplyr::starts_with("Season"),
				.fns = \(x) stringr::str_replace(x, "----", " — ")
			)
		)

	# Side quest to get the number of strains in each season
	strain_counts <-
		model_data |>
		dplyr::distinct(season, assay_strain) |>
		dplyr::count(season) |>
		tidyr::pivot_wider(names_from = season, values_from = n) |>
		do.call(what = c)

	total_unique_strains <- dplyr::n_distinct(model_data$assay_strain)

	# Add the strains count row to the end of the tbl_data
	new_row <- rlang::set_names(
		as.character(c("Num. strains", strain_counts, total_unique_strains)),
		colnames(tbl_data)
	)

	person_tbl <-
		dplyr::bind_rows(tbl_data, new_row) |>
		dplyr::rename(
			"Study Site" = study
		) |>
		flextable::flextable() |>
		flextable::separate_header() |>
		flextable::merge_v(j = 1) |>
		flextable::hline(i = 1, j = 2:7, part = "header") |>
		flextable::align(j = 2:7, align = "center") |>
		flextable::fix_border_issues() |>
		fit_flextable_to_page()

	save_file_as_qs2(
		person_tbl,
		file_path
	)

	return(file_path)
}

get_unique_individual_count <- function(model_data, file_path) {
	out <-
		dplyr::bind_rows(
			model_data,
			dplyr::mutate(model_data, study = "Overall")
		) |>
		dplyr::distinct(study, subject_id) |>
		dplyr::count(study) |>
		tibble::deframe()

	save_file_as_qs2(
		out,
		file_path
	)

	here::here(file_path)
}

# Figure showing raw pre and post titers for all strains
make_titer_plot <- function(model_data, file_path) {
	set.seed(370)
	titer_data <-
		model_data |>
		dplyr::distinct(
			assay_strain, vaccine_strain, unique_id, log_pretiter, log_posttiter
		) |>
		tidyr::pivot_longer(
			cols = c(log_pretiter, log_posttiter),
			names_to = "time",
			values_to = "titer"
		) |>
		dplyr::mutate(
			time = factor(
				time,
				levels = c("log_pretiter", "log_posttiter"),
				labels = c("Pre-vaccination", "Post-vaccination")
			),
			jx = as.integer(assay_strain) + rnorm(dplyr::n(), 0, 0.12),
			jy = hgp::hai_to_natural_scale(titer + rnorm(dplyr::n(), 0, 0.12))
		)

	labs_s <- droplevels(model_data$assay_strain) |> levels()
	titer_plot <-
		titer_data |>
			ggplot2::ggplot() +
		ggplot2::aes(
			x = jx,
			y = jy
		) +
		ggplot2::geom_point(
			size = 0.5,
			alpha = 0.075,
			#position = ggplot2::position_jitter(0.4, 0.25, 370)
		) +
		ggplot2::facet_grid(vaccine_strain ~ time) +
		ggplot2::scale_x_continuous(
			breaks = 1:length(labs_s),
			labels = labs_s
		) +
		ggplot2::scale_y_continuous(
			trans = "log2",
			breaks = hgp::hai_to_natural_scale(seq(0, 10, 2))
		) +
		ggplot2::labs(
			x = "Assay Strain",
			y = "Reciprocal HAI titer"
		) +
		hgp::theme_ms() +
		ggplot2::theme(
			legend.position = "bottom",
			legend.key.width = ggplot2::unit(0.5, "in"),
			axis.text.x = ggplot2::element_text(
				size = 14, hjust = 1, vjust = 1, angle = 45
			),
			axis.text.y = ggplot2::element_text(size = 18),
			strip.text = ggplot2::element_text(size = 28),
			axis.title = ggplot2::element_text(size = 28),
			legend.title = ggplot2::element_text(size = 24),
			legend.text = ggplot2::element_text(size = 18)
		)

	ensure_directory_exists(file_path)
	ggplot2::ggsave(
		plot = titer_plot,
		filename = file_path,
		width = 13,
		height = 13
	)

	return(file_path)
}

# Demographics table ####


# Strain names table ####
make_strain_name_table <- function(model_data, file_name) {
	s <- unique(model_data$assay_strain)

	tbl_data <- tibble::tibble(
		"Short Name" = s,
		"Full Name" = hgp::replace_strain_names(s, from = "short", to = "full")
	)

	out_tbl <-
		tbl_data |>
		dplyr::arrange(`Short Name`) |>
		dplyr::distinct() |>
		flextable::flextable() |>
		flextable::merge_v(j = 1) |>
		flextable::valign(j = 1, valign = "top") |>
		flextable::fontsize(size = 10, part = "all") |>
		flextable::padding(padding.top = 0.5, padding.bottom = 0.5) |>
		fit_flextable_to_page()

	save_file_as_qs2(
		out_tbl,
		file_name
	)

	return(file_name)
}

# Vaccine strains for each year table ####
make_vaccine_strain_table <- function(model_data, file_name) {
	tbl_data <-
		model_data |>
		dplyr::distinct(season_short, vaccine_strain) |>
		dplyr::rename(
			Season = season_short,
			Vaccine = vaccine_strain
		)

	tbl_out <-
		flextable::flextable(tbl_data) |>
		fit_flextable_to_page()

	save_file_as_qs2(
		tbl_out,
		file_name
	)

	return(file_name)
}

# Demographics table ####
make_demographics_table <- function(model_data, file_name) {
	set_gtsummary_theme()

	demo_data <-
		model_data |>
		# Select the demographics columns we need
		dplyr::distinct(
			# Cohort grouping
			subject_id, season = season_short,
			# Info for table
			dose, sex, race, age, birth_year, study
		)

	tbl_out <-
		gtsummary::tbl_summary(
			data = demo_data,
			by = season,
			include = -subject_id,
			label = list(
				dose = "Dose",
				sex = "Sex assigned at birth",
				race = "Race/Ethnicity",
				age = "Age (years)",
				birth_year = "Birth year",
				study = "Study site"
			),
			statistic = list(
				gtsummary::all_continuous() ~ "{median} ({min} - {max})"
			),
			digits = list(
				gtsummary::all_continuous() ~ 0,
				gtsummary::all_categorical() ~ 0
			)
		) |>
		gtsummary::add_overall(last = TRUE) |>
		gtsummary::as_flex_table() |>
		fit_flextable_to_page()

	save_file_as_qs2(
		tbl_out,
		file_name
	)

	return(file_name)
}
