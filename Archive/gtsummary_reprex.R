# Gtsummary reprex

# First make the data
dat <-
	tidyr::nesting(
		A = c(1, 2),
		B = c(1, 2),
		C = c(2, 3)
	) |>
	tidyr::pivot_longer(
		dplyr::everything(),
		names_to = "study_site",
		values_to = "year"
	) |>
	tidyr::expand_grid(id = 1:3) |>
	dplyr::mutate(
		y = rnorm(dplyr::n()),
		study_site = factor(study_site),
		year = factor(year, ordered = TRUE),
		id = as.character(id)
	)

gtsummary::tbl_strata(
	data = dat,
	strata = study_site,
	.tbl_fun = \(x) x |>
		gtsummary::tbl_summary(
			by = year,
			include = c(y),
			type = list(y ~ "continuous")
		)
)
