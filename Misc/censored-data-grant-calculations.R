###
# Censored data grant calculations
# Zane & Andreas
# 2024-09-27
# We want to show how censored data can affect important calculations. So for
# each vaccine, we'll calculate the homologous and overall GMT and average
# titer increase. We'll do this naively (using 5 as the LoD), and with the
# censoring correction.
###

# Declare dependencies
library(ggplot2)
library(patchwork)

# Load data
dat <-
	here::here("data", "processed", "joined-data.Rds") |>
	readr::read_rds() |>
	dplyr::mutate(
		season_first = substr(as.character(season), 1, 4),
		season_label = paste0(study, ": ", season_first)
	)

do_the_analysis <- function(
		.data = dat,
		outcome, lod = 5,
		strata = c("study", "season", "strain_type", "vaccine_strain"),
		upper_limit = 5
	) {

	dat_local <-
		.data |>
		dplyr::mutate(
			# Mark which pre and post values are censored before we mess with the
			# lod coding
			log_post_cens = !(posttiter == 5),
			log_pre_cens = !(pretiter == 5),
			# Update posttiter lod coding
			posttiter = ifelse(
				posttiter == 5, lod, posttiter
			),
			log_post = log2(posttiter / lod),
			# Update pretiter lod coding
			pretiter = ifelse(
				pretiter == 5, lod, pretiter
			),
			log_post = log2(pretiter / lod),
			# Recalculate naive TI after checking censoring
			fold_change = posttiter / pretiter,
			titerincrease = log2(fold_change)
		) |>
		# Compute Yang's censoring bounds for TI
		dplyr::mutate(
			ti_lwr = dplyr::case_when(
				log_pre_cens == 0 & log_post_cens == 0 ~ -Inf,
				log_pre_cens == 0 & log_post_cens == 1 ~ log2(posttiter / 10),
				log_pre_cens == 1 & log_post_cens == 0 ~ -Inf,
				log_pre_cens == 1 & log_post_cens == 1 ~ titerincrease
			),
			ti_upr = dplyr::case_when(
				log_pre_cens == 0 & log_post_cens == 0 ~ Inf,
				log_pre_cens == 0 & log_post_cens == 1 ~ Inf,
				log_pre_cens == 1 & log_post_cens == 0 ~ log2(10 / pretiter),
				log_pre_cens == 1 & log_post_cens == 1 ~ titerincrease
			),
			ti_cens = dplyr::case_when(
				log_pre_cens == 0 & log_post_cens == 0 ~ 3,
				log_pre_cens == 0 & log_post_cens == 1 ~ 0,
				log_pre_cens == 1 & log_post_cens == 0 ~ 2,
				log_pre_cens == 1 & log_post_cens == 1 ~ 1
			)
		)

	# Nest the data so we can easily fit a model per each stratum
	dat_nested <-
		dat_local |>
		dplyr::filter(
			as.character(strain_name) == as.character(vaccine_strain)
		) |>
		tidyr::nest(data = -tidyselect::all_of(strata))

	# Set some scaling factors and formulas correctly based on the outcome
	if (outcome == "post") {
		naive_fmla <- as.formula(
			log_post ~ 1
		)
		surv_fmla <- as.formula(
			survival::Surv(
				ifelse(log_post_cens == 0, 10, log_post),
				log_post_cens, type = "left") ~ 1
		)
		#cutoff <- 40
		mult <- 5
		#outcome_label <- "Average postvaccination titer"
	} else if (outcome == "ti") {
		naive_fmla <- as.formula(
			titerincrease ~ 1
		)
		surv_fmla <- as.formula(
			survival::Surv(time = ti_lwr, time2 = ti_upr, type = "interval2") ~ 1
		)
		#cutoff <- 4
		mult <- 1
		#outcome_label <- "Vaccine-induced antibody titer increase (fold change)"
	} else {stop("wrong outcome")}

	# Fit the two models
	dat_models <-
		dat_nested |>
		dplyr::mutate(
			naive = purrr::map(
				data,
				\(d) lm(naive_fmla, data = d)
			),
			cens = purrr::map(
				data,
				\(d) survival::survreg(
					surv_fmla,
					dist = "gaussian",
					data = d
				)
			)
		)

	# Get predictions -- note that the newdata value doesn't matter but does
	# have to exist.
	dat_preds <-
		dat_models |>
		dplyr::mutate(
			naive_est = purrr::map(
				naive,
				\(m) broom::augment(
					m, newdata = data.frame(x = 1),
					se_fit = TRUE
				)
			),
			cens_est = purrr::map(
				cens,
				\(m) broom::augment(
					m, newdata = data.frame(x = 1),
					type.predict = "response",
					se.fit = TRUE
				)
			)
		) |>
		dplyr::select(-naive, -cens) |>
		tidyr::pivot_longer(c(naive_est, cens_est)) |>
		tidyr::unnest(value) |>
		dplyr::mutate(
			est = .fitted,
			# Calculate wald-type CI
			lwr = .fitted - 1.96 * .se.fit,
			upr = .fitted + 1.96 * .se.fit,
			# And back transform
			dplyr::across(
				c(est, lwr, upr),
				\(x) mult * 2 ^ x
			),
			# Clean up method
			name = factor(
				name,
				levels = c("naive_est", "cens_est"),
				labels = c("naive", "corrected")
			)
		) |>
		dplyr::select(-data, -x, -.se.fit, -.fitted)

	return(dat_preds)
}

make_plot <- function(
		analysis_results, outcome_label, lower_limit, upper_limit, title, cutoff,
		mult = 1
	) {
	plt <-
		ggplot(analysis_results) +
		aes(
			x = est, xmin = lwr, xmax = upr,
			y = season_label, color = name, shape = name, linetype = name
		) +
		geom_vline(xintercept = cutoff) +
		geom_pointrange(position = position_dodge(width = 0.75)) +
		hgp::theme_ms(text_size_axis_text = 14) +
		#facet_wrap(~strain_type) +
		scale_x_continuous(
			breaks = scales::breaks_log(base = 2, n = 6),
			trans = "log2"
		) +
		coord_cartesian(
			xlim = mult * 2 ^ c(lower_limit, upper_limit),
		) +
		scale_color_brewer(palette = "Dark2") +
		labs(
			x = outcome_label,
			y = NULL,
			color = "Model",
			shape = "Model",
			linetype = "Model",
			title = title
		)

	return(plt)
}

# Commented out because this was us trying different things
# post_5 <- do_the_analysis(dat, "post")
# post_0 <- do_the_analysis(dat, "post", lod = 0.01)
# post_10 <- do_the_analysis(dat, "post", lod = 10)
#
#
# ti_5 <- do_the_analysis(dat, "ti", lod = 5)
# ti_0 <- do_the_analysis(dat, "ti", lod = 0.01)
# ti_10 <- do_the_analysis(dat, "ti", lod = 10)

# list(ti_5, ti_0, ti_10) |>
# 	purrr::map()
# 	patchwork::plot_layout(axes = "collect", guides = "collect") &
# 	theme(legend.position = "bottom")

# ti_1 <- do_the_analysis(dat,"ti", lod = 1)
# ti_01 <- do_the_analysis(dat,"ti", lod = 0.01)
# ti_10 <- do_the_analysis(dat,"ti", lod = 1e-10)

# res <-
# 	purrr::map(
# 		c(10, 5, 1),
# 		\(x) do_the_analysis(dat_pafl, "ti", lod = x)
# 	)

# res_plot <-
# 	res |>
# 	purrr::map(\(r) purrr::pluck(r, "plot")) |>
# 	purrr::reduce(`/`) +
	# patchwork::plot_layout(axes = "collect", guides = "collect") &
	# theme(legend.position = "bottom")

# res_plot


dat_pafl_h1 <-
	dat |>
	dplyr::filter(study != "UGA", strain_type == "H1N1") |>
	dplyr::mutate(
		season_label = factor(
			season_label,
			label = paste0(c("PA: ", "FL: "), rep(seq(2013, 2016), each = 2))
		)
	)


dat_pafl_h3 <-
	dat |>
	dplyr::filter(study != "UGA", strain_type == "H3N2") |>
	dplyr::mutate(
		season_label = factor(
			season_label,
			label = paste0(c("PA: ", "FL: "), rep(seq(2013, 2016), each = 2))
		)
	)

pafl_h3_ti_lod10 <- do_the_analysis(
	dat_pafl_h3,
	"ti",
	lod = 10,
	strata = c("season_label", "strain_type", "vaccine_strain")
)

pafl_h3_ti_lod00 <- do_the_analysis(
	dat_pafl_h3,
	"ti",
	lod = 1,
	strata = c("season_label", "strain_type", "vaccine_strain")
)

pafl_h1_ti_lod10 <- do_the_analysis(
	dat_pafl_h1,
	"ti",
	lod = 10,
	strata = c("season_label", "strain_type", "vaccine_strain")
)

pafl_h1_ti_lod00 <- do_the_analysis(
	dat_pafl_h1,
	"ti",
	lod = 1,
	strata = c("season_label", "strain_type", "vaccine_strain")
)

{
	plot1_h3 <- make_plot(
		pafl_h3_ti_lod10,
		"Average fold change post-vaccination",
		0.9, 3.4,
		"Censored values set to LoD",
		cutoff = 4
	)

	plot2_h3 <- make_plot(
		pafl_h3_ti_lod00,
		"Average fold change post-vaccination",
		0.9, 3.4,
		"Censored values set to 0 (after log transformation)",
		cutoff = 4
	)

	combined_plot_h3 <-
		(plot1_h3 / plot2_h3)	+
		patchwork::plot_layout(axes = "collect", guides = "collect") &
		theme(legend.position = "bottom")
}

ggsave(
	filename = here::here("misc", "censored-data-grant-example-h3.png"),
	plot = combined_plot_h3,
	width = 8,
	height = 6,
	dpi = 300,
	units = "in"
)

{
	plot1_h1 <- make_plot(
		pafl_h1_ti_lod10,
		"Average fold change post-vaccination",
		0.4, 2.1,
		"Censored values set to LoD",
		cutoff = 4
	)

	plot2_h1 <- make_plot(
		pafl_h1_ti_lod00,
		"Average fold change post-vaccination",
		0.4, 2.1,
		"Censored values set to 0 (after log transformation)",
		cutoff = 4
	)

	combined_plot_h1 <-
		(plot1_h1 / plot2_h1)	+
		patchwork::plot_layout(axes = "collect", guides = "collect") &
		theme(legend.position = "bottom")
}

ggsave(
	filename = here::here("misc", "censored-data-grant-example-h1.png"),
	plot = combined_plot_h1,
	width = 8,
	height = 6,
	dpi = 300,
	units = "in"
)

h1_res <-
	dplyr::bind_rows(
		"lod" = pafl_h1_ti_lod10,
		"zero" = pafl_h1_ti_lod00,
		.id = "id"
	) |>
	dplyr::filter(!(id == "zero" & name == "corrected")) |>
	dplyr::mutate(
		model = factor(
			dplyr::case_when(
				id == "lod" & name == "naive" ~ "lod",
				id == "lod" & name == "corrected" ~ "corrected",
				id == "zero" & name == "naive" ~ "zero"
			),
		levels = c("lod", "zero", "corrected"),
		labels = c(
			"Censored values set to LoD",
			"Censored values set to zero",
			"Corrected for censoring"
		)
	)
)

alltogether_h1 <-
	ggplot(h1_res) +
	aes(
		x = est, xmin = lwr, xmax = upr,
		y = season_label, color = model, shape = model, linetype = model
	) +
	geom_vline(xintercept = 4) +
	geom_pointrange(
		position = position_dodge(width = 0.75),
		linewidth = 1.25, fatten = 7
	) +
	geom_point(
		position = position_dodge(width = 0.75),
		size = 2, color = "white"
	) +
	hgp::theme_ms(
		base_size = 14,
		text_size_axis_text = 18,
		text_size_axis_title = 20,
		text_size_legend = 28,
		text_size_legend_title = 22
	) +
	scale_x_continuous(
		breaks = scales::breaks_log(base = 2, n = 6),
		trans = "log2"
	) +
	coord_cartesian(
		xlim =  2 ^ c(0.4, 2.1),
	) +
	scale_color_brewer(palette = "Dark2") +
	labs(
		x = "Average fold change post-vaccination",
		y = NULL,
		color = "Model",
		shape = "Model",
		linetype = "Model"
	) +
	guides(
		color = guide_legend(ncol = 1)
	)

ggsave(
	filename = here::here("misc", "censored-data-grant-example-h1-combined.png"),
	plot = alltogether_h1,
	width = 8,
	height = 6,
	dpi = 300,
	units = "in"
)

h3_res <-
	dplyr::bind_rows(
		"lod" = pafl_h3_ti_lod10,
		"zero" = pafl_h3_ti_lod00,
		.id = "id"
	) |>
	dplyr::filter(!(id == "zero" & name == "corrected")) |>
	dplyr::mutate(
		model = factor(
			dplyr::case_when(
				id == "lod" & name == "naive" ~ "lod",
				id == "lod" & name == "corrected" ~ "corrected",
				id == "zero" & name == "naive" ~ "zero"
			),
			levels = c("lod", "zero", "corrected"),
			labels = c(
				"Censored values set to LoD",
				"Censored values set to zero",
				"Corrected for censoring"
			)
		)
	)

alltogether_h3 <-
	ggplot(h3_res) +
	aes(
		x = est, xmin = lwr, xmax = upr,
		y = season_label, color = model, shape = model, linetype = model
	) +
	geom_vline(xintercept = 4) +
	geom_pointrange(
		position = position_dodge(width = 0.75),
		linewidth = 1.25, fatten = 7
	) +
	geom_point(
		position = position_dodge(width = 0.75),
		size = 2, color = "white"
	) +
	hgp::theme_ms(
		base_size = 14,
		text_size_axis_text = 18,
		text_size_axis_title = 20,
		text_size_legend = 28,
		text_size_legend_title = 22
	) +
	#facet_wrap(~strain_type) +
	scale_x_continuous(
		breaks = scales::breaks_log(base = 2, n = 6),
		trans = "log2"
	) +
	coord_cartesian(
		xlim =  2 ^ c(0.9, 3.4),
	) +
	scale_color_brewer(palette = "Dark2") +
	labs(
		x = "Average fold change post-vaccination",
		y = NULL,
		color = "Model",
		shape = "Model",
		linetype = "Model"
	) +
	guides(
		color = guide_legend(ncol = 1)
	)

ggsave(
	filename = here::here("misc", "censored-data-grant-example-h3-combined.png"),
	plot = alltogether_h3,
	width = 8,
	height = 6,
	dpi = 300,
	units = "in"
)



# Combined versions
