###
# Data Summary
# Zane Billings
# 2024-08-12
# Calculates descriptive and summary statistics for the combined dataset.
###

# Setup ####

## Declare dependencies ####
box::use(
	gtsummary
)

## gtsummary and ggplot setup ####
gtsummary::theme_gtsummary_language("en", big.mark = "")

## Data loading ####
# Load cohort data
cohort_data <-
	here::here("data", "processed", "analysis-data.Rds") |>
	readr::read_rds()

# Load reporting data
reporting_data <-
	here::here("data", "processed", "reporting-data.Rds") |>
	readr::read_rds()

# Basic counts ####
## Total number of assays in each study with number of strains ####
n_assays_by_study <-
	cohort_data |>
	dplyr::distinct(id, season, study, strain_name, dose) |>
	dplyr::summarize(
		n = dplyr::n(),
		n_strains = dplyr::n_distinct(strain_name),
		.by = c(study, dose)
	)

n_assays_overall <-
	cohort_data |>
	dplyr::distinct(id, season, study, strain_name, dose) |>
	dplyr::summarize(
		study = "Overall",
		n = dplyr::n(),
		n_strains = dplyr::n_distinct(strain_name),
		.by = dose
	)

n_assays_study_with_overall <-
	dplyr::bind_rows(n_assays_by_study, n_assays_overall)

qs::qsave(
	n_assays_study_with_overall,
	here::here("results", "tables", "n-assays-by-study.qs")
)

## Total number of person years in each study ####
n_person_years_by_study <-
	cohort_data |>
	dplyr::distinct(id, season, study, dose) |>
	dplyr::summarize(
		n = dplyr::n(),
		.by = c(study, dose)
	)

n_person_years_by_dose <-
	cohort_data |>
	dplyr::distinct(id, season, study, dose) |>
	dplyr::summarize(
		study = "Overall",
		n = dplyr::n(),
		.by = dose
	)

n_person_years_overall <-
	cohort_data |>
	dplyr::distinct(id, season, study, dose) |>
	dplyr::summarize(
		study = "Overall",
		dose = "Overall",
		n = dplyr::n()
	)

n_person_years_study_with_overall <-
	dplyr::bind_rows(
		n_person_years_by_study,
		n_person_years_by_dose,
		n_person_years_overall
	)

qs::qsave(
	n_person_years_study_with_overall,
	here::here("results", "tables", "n-person-years-by-study.qs")
)

## Total number of individuals in each study ####
n_individuals_by_study <-
	cohort_data |>
	dplyr::distinct(id, study) |>
	dplyr::summarize(
		n = dplyr::n(),
		.by = study
	)

n_individuals_overall <-
	cohort_data |>
	dplyr::distinct(id, study) |>
	dplyr::summarize(
		study = "Overall",
		n = dplyr::n()
	)

n_individuals_study_with_overall <-
	dplyr::bind_rows(n_individuals_by_study, n_individuals_overall)

qs::qsave(
	n_individuals_study_with_overall,
	here::here("results", "tables", "n-individuals-by-study.qs")
)

# Demographic information table ####
## Table with all person-years by strain/study ####
table1_person_years_study <-
	reporting_data |>
	dplyr::mutate(
		season = season |>
			as.character() |>
			gsub(pattern = " - ", replacement = "/") |>
			gsub(pattern = "20", replacement = ""),
		race = dplyr::if_else(
			startsWith(race, "OTH"),
			"Other",
			race
		) |>
			factor(levels = c(
				"White", "Black or African American", "Asian",
				"American Indian or Alaska Native",
				"Native Hawaiian or Other Pacific Islander", "Other",
				"Not Specified"
			)),
		dplyr::across(dplyr::where(is.factor), forcats::fct_drop)
	) |>
	gtsummary::tbl_strata(
		strata = study,
		.tbl_fun = \(x) x |>
			gtsummary::tbl_summary(
				by = season,
				include = c("sex", "race", "ethnicity", "age", "birth_year"),
				label = list(
					"sex" ~ "Sex assigned at birth",
					"race" ~ "Race",
					"ethnicity" ~ "Ethnicity",
					"age" ~ "Age",
					"birth_year" ~ "Birth year"
				),
				statistic = list(
					gtsummary::all_continuous() ~ "{median} ({min}, {max})",
					gtsummary::all_categorical() ~ "{n} ({p}%)"
				),
				missing = "no"
			) |>
			gtsummary::add_overall(last = TRUE),
		.header = "**{strata}**, N = {n}"
	)

table1_person_years_overall <-
	reporting_data |>
	dplyr::mutate(
		study = "Overall",
		season = season |>
			as.character() |>
			gsub(pattern = " - ", replacement = "/") |>
			gsub(pattern = "20", replacement = ""),
		race = dplyr::if_else(
			startsWith(race, "OTH"),
			"Other",
			race
		) |>
			factor(levels = c(
				"White", "Black or African American", "Asian",
				"American Indian or Alaska Native",
				"Native Hawaiian or Other Pacific Islander", "Other",
				"Not Specified"
			)),
		dplyr::across(dplyr::where(is.factor), forcats::fct_drop)
	) |>
	gtsummary::tbl_summary(
		by = season,
		include = c("sex", "race", "ethnicity", "age", "birth_year"),
		label = list(
			"sex" ~ "Sex assigned at birth",
			"race" ~ "Race",
			"ethnicity" ~ "Ethnicity",
			"age" ~ "Age",
			"birth_year" ~ "Birth year"
		),
		statistic = list(
			gtsummary::all_continuous() ~ "{median} ({min}, {max})",
			gtsummary::all_categorical() ~ "{n} ({p}%)"
		)
	) |>
	gtsummary::add_overall(last = TRUE) |>
	gtsummary::modify_spanning_header(
		gtsummary::all_stat_cols() ~ '**All Studies**, N = {N}'
	)

gtsummary::tbl_merge(
	list(
		table1_person_years_study,
		table1_person_years_overall
	),
	tab_spanner = FALSE
)

## Table 1 individuals table ####
table1_person_years_study <-
	reporting_data |>
	dplyr::mutate(
		season = season |>
			as.character() |>
			gsub(pattern = " - ", replacement = "/") |>
			gsub(pattern = "20", replacement = "") |>
			factor(ordered = TRUE),
		race = dplyr::if_else(
			startsWith(race, "OTH"),
			"Other",
			race
		) |>
			factor(levels = c(
				"White", "Black or African American", "Asian",
				"American Indian or Alaska Native",
				"Native Hawaiian or Other Pacific Islander", "Other",
				"Not Specified"
			)),
		dplyr::across(dplyr::where(is.factor), forcats::fct_drop)
	) |>
	gtsummary::tbl_strata(
		strata = study,
		.tbl_fun = \(x) x |>
			gtsummary::tbl_summary(
				by = season,
				include = c("sex", "race", "ethnicity", "age", "birth_year"),
				label = list(
					"sex" ~ "Sex assigned at birth",
					"race" ~ "Race",
					"ethnicity" ~ "Ethnicity",
					"age" ~ "Age",
					"birth_year" ~ "Birth year"
				),
				statistic = list(
					gtsummary::all_continuous() ~ "{median} ({min}, {max})",
					gtsummary::all_categorical() ~ "{n} ({p}%)"
				),
				missing = "no"
			) |>
			gtsummary::add_overall(last = TRUE),
		.header = "**{strata}**, N = {n}"
	)

table1_person_years_overall <-
	reporting_data |>
	dplyr::mutate(
		study = "Overall",
		season = season |>
			as.character() |>
			gsub(pattern = " - ", replacement = "/") |>
			gsub(pattern = "20", replacement = "") |>
			factor(ordered = TRUE),
		race = dplyr::if_else(
			startsWith(race, "OTH"),
			"Other",
			race
		) |>
			factor(levels = c(
				"White", "Black or African American", "Asian",
				"American Indian or Alaska Native",
				"Native Hawaiian or Other Pacific Islander", "Other",
				"Not Specified"
			)),
		dplyr::across(dplyr::where(is.factor), forcats::fct_drop)
	) |>
	gtsummary::tbl_summary(
		by = season,
		include = c("sex", "race", "ethnicity", "age", "birth_year"),
		label = list(
			"sex" ~ "Sex assigned at birth",
			"race" ~ "Race",
			"ethnicity" ~ "Ethnicity",
			"age" ~ "Age",
			"birth_year" ~ "Birth year"
		),
		statistic = list(
			gtsummary::all_continuous() ~ "{median} ({min}, {max})",
			gtsummary::all_categorical() ~ "{n} ({p}%)"
		)
	) |>
	gtsummary::add_overall(last = TRUE) |>
	gtsummary::modify_spanning_header(
		gtsummary::all_stat_cols() ~ '**All Studies**, N = {N}'
	)

gtsummary::tbl_merge(
	list(
		table1_person_years_study,
		table1_person_years_overall
	),
	tab_spanner = FALSE
) |>
	gtsummary::modify_caption(paste0(
		"Demographic information for all person-years included in the study. ",
		"Column headers, e.g. 13/14, refer to the influenza season, and spanning ",
		"headers, e.g. FL, refer to the study site. Note that the PA and FL ",
		"study sites ceased recruitment after the 2016/17 influenza season, which ",
		"is when the UGA study site began recruitment."
	)) |>
	gtsummary::as_flex_table() |>
	flextable::fontsize(size = 8, part = "all") |>
	flextable::padding(padding = 0, part = "all") |>
	flextable::bold(
		i = c(1, 4, 12, 16, 17),
		j = 1,
		part = "body"
	)

# Non gtsummary version to see if I can do what I want
table1_reporting_data_studies <-
	reporting_data |>
	dplyr::mutate(
		season = season |>
			as.character() |>
			gsub(pattern = " - ", replacement = "/") |>
			gsub(pattern = "20", replacement = "") |>
			factor(ordered = TRUE),
		race = dplyr::if_else(
			startsWith(race, "OTH"),
			"Other",
			race
		) |>
			factor(levels = c(
				"White", "Black or African American", "Asian",
				"American Indian or Alaska Native",
				"Native Hawaiian or Other Pacific Islander", "Other",
				"Not Specified"
			)),
		dplyr::across(dplyr::where(is.factor), forcats::fct_drop)
	)

table1_reporting_data <-
	dplyr::bind_rows(
		table1_reporting_data_studies,
		table1_reporting_data_studies |>
			dplyr::mutate(study = "All Studies")
	) |>
	dplyr::mutate(
		study = factor(study, levels = c("FL", "PA", "UGA", "All Studies"))
	)

fit_flextable_to_page <- function(ft, pgwidth = 6.5){

	ft_out <- flextable::autofit(ft)

	ft_out <- flextable::width(
		ft_out,
		width = (dim(ft_out)$widths*pgwidth) /
			(flextable::flextable_dim(ft_out)$widths)
	)
	return(ft_out)
}

table1_person_year_function <- function(d, cols_to_delete) {
	suppressMessages({
	out <- d |>
		gtsummary::tbl_strata(
			strata = study,
			.tbl_fun = \(x) x |>
				gtsummary::tbl_summary(
					by = season,
					include = c("sex", "race", "ethnicity", "age", "birth_year"),
					label = list(
						"sex" ~ "Sex assigned at birth",
						"race" ~ "Race",
						"ethnicity" ~ "Ethnicity",
						"age" ~ "Age",
						"birth_year" ~ "Birth year"
					),
					statistic = list(
						gtsummary::all_continuous() ~ "{median} ({min}, {max})",
						gtsummary::all_categorical() ~ "{n} ({p}%)"
					),
					missing = "no"
				) |>
				gtsummary::add_overall(last = TRUE),
			.header = "**{strata}**, N = {n}"
		) |>
		gtsummary::as_flex_table() |>
		flextable::fontsize(size = 8, part = "all") |>
		flextable::padding(padding = 0, part = "all") |>
		flextable::bold(
			i = c(1, 4, 12, 16, 17),
			j = 1,
			part = "body"
		)
	})

	if(length(cols_to_delete) >= 1) {
		out <- out |> flextable::compose(
			i = 1:17,
			j = cols_to_delete,
			value = flextable::as_paragraph(""),
			part = "body"
		)
	}

	out <- out |>
		flextable::autofit() |>
		fit_flextable_to_page()

	return(out)
}


table1_person_year_list <-
	table1_reporting_data |>
	dplyr::group_split(study) |>
	purrr::map2(
		list(6, 6, 2:4, numeric(0)),
		\(x, y) table1_person_year_function(x, y)
	)

# normalize <- function(x, ...) {
# 	return(x / sum(x, ...))
# }
#
# tbl2_header <- table1_py_bot$header$dataset[1, ] |>
# 	as.character()
#
# pars <- flextable::as_paragraph(
# 	stringr::str_extract(tbl2_header, "\\*\\*(.*)\\*\\*", group = 1) |>
# 		tidyr::replace_na(replace = " ") |>
# 	flextable::as_b(),
# 	c("", rep(", ", 12)),
# 	stringr::str_extract(tbl2_header, ", (.*)$", group = 1) |>
# 		tidyr::replace_na(replace = " ")
# )
#
# table1_py_top_hbot <-
# 	table1_py_top |>
# 	add_body_row(
# 		values = pars,
# 		top = FALSE
# 	)
#
# do.call(
# 	table1_py_bot$body$dataset |> as.list(),
# 	what = \(...) flextable::add_body(table1_py_top_hbot, ..., top = FALSE)
# ) |>
# 	flextable::merge_h(18) |>
# 	flextable::border_remove() |>
# 	flextable::hline(part = "header", i = 1:2)
#
#
# rbind(table1_py_top, table1_py_bot)


geo_mean <- function(x, ...) {
	return(exp(mean(log(x, ...))))
}

geo_sd <- function(x, ...) {
	return(exp(sd(log(x, ...))))
}

cohort_data |>
	dplyr::mutate(study_season = paste0(study, ": ", season)) |>
	gtsummary::tbl_continuous(
		variable = posttiter,
		by = study,
		include = c(season),
		statistic = gtsummary::everything() ~ "{geo_mean} ({geo_sd})"
	)

# Vaccine strains table ####

# Historical panel table ####

# END OF FILE ####
