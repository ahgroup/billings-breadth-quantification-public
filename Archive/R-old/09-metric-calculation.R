###
# Metrics Calculation
# Zane Billings
# 2024-08-21
# Calculating breadth metrics from the different models
###

# Setup ####
box::use(
	here,
	qs,
	dplyr,
	ggdist,
	pracma
)

source(here::here("R", "functions", "utils.R"))

# Data processing ####
# Function for loading all of the data we need since there is A LOT
create_metrics_data <- function() {
	# Load the full model fits
	full_model_pth <- here::here(
		"results","large-files", "model-results", "full-data-model-results"
	)
	full_model_data <- qs::qread(here::here(full_model_pth, "data.qs"))
	full_model_groups <- qs::qread(here::here(full_model_pth, "groups.qs"))
	full_model_preds <- qs::qread(here::here(full_model_pth, "preds.qs"))

	# Homologous full model fits
	h_full_model_pth <- here::here(
		"results","large-files", "model-results",
		"homologous-full-data-model-results"
	)
	h_full_model_preds <- qs::qread(here::here(h_full_model_pth, "preds.qs"))

	# Seroconversion full model fits
	s_full_model_pth <- here::here(
		"results", "large-files", "model-results",
		"seroconversion-full-data-model-results"
	)
	s_full_model_groups <- qs::qread(here::here(s_full_model_pth, "groups.qs"))
	s_full_model_preds <- qs::qread(here::here(s_full_model_pth, "preds.qs"))
	# Because the seroconversion model doesn't care about censoring or which model
	# outcome we used, there are less rows than in the other data frames. So we
	# have to do a bit of cleanup and then join this one to the others.
	seroconversion_data_full_model <-
		s_full_model_groups |>
		dplyr::select(-vaccine_strain, -group_id, -model) |>
		dplyr::mutate(
			seroconversion = s_full_model_preds,
			subsample_id = "full data"
		)

	# Clean up the full model data
	full_model_metrics_input <-
		full_model_groups |>
		dplyr::mutate(
			censoring_correction = ifelse(
				grepl("with", model), "censoring correction", "naive model"
			),
			model_structure = ifelse(
				grepl("full", model), "full", "reduced"
			),
			subsample_id = "full data",
			.keep = "unused"
		) |>
		tibble::add_column(
			data = full_model_data,
			preds = full_model_preds
		) |>
		tidyr::pivot_wider(names_from = model_structure, values_from = preds) |>
		tibble::add_column(
			homologous = h_full_model_preds
		) |>
		# We only care about these seasons.
		dplyr::filter(season < "2018 - 2019") |>
		# Drop extra values from season
		dplyr::mutate(season = forcats::fct_drop(season)) |>
		# Join with the seroconversion data
		dplyr::left_join(
			seroconversion_data_full_model,
			by = c('season', 'strain_type', 'method', 'subsample_id')
		)

	# Load the subset fits
	subset_model_pth <- here::here(
		"results", "large-files", "model-results", "subsample-model-results"
	)
	subsample_ids <-
		here::here("results", "data", "subsample-IDs.qs") |>
		qs::qread()
	subset_model_data <- qs::qread(here::here(subset_model_pth, "data.qs"))
	subset_model_groups <- qs::qread(here::here(subset_model_pth, "groups.qs"))
	subset_model_preds <- qs::qread(here::here(subset_model_pth, "preds.qs"))

	# Homologous subset fits
	h_subset_model_pth <- here::here(
		"results", "large-files", "model-results",
		"homologous-subsample-model-results"
	)
	h_subset_model_preds <- qs::qread(here::here(h_subset_model_pth, "preds.qs"))

	# Seroconversion fits
	s_subsample_pth <- here::here(
		"results", "large-files", "model-results",
		"seroconversion-subsample-model-results"
	)
	s_subsample_groups <- qs::qread(here::here(s_subsample_pth, "groups.qs"))
	s_subsample_preds <- qs::qread(here::here(s_subsample_pth, "preds.qs"))
	# Because the seroconversion model doesn't care about censoring or which model
	# outcome we used, there are less rows than in the other data frames. So we
	# have to do a bit of cleanup and then join this one to the others.
	seroconversion_data_subsample <-
		s_subsample_groups |>
		dplyr::select(-vaccine_strain, -model) |>
		dplyr::mutate(
			seroconversion = s_subsample_preds,
			subsample_id = subsample_ids[1:dplyr::n()]
		)

	# Clean up the subsample model data
	subsample_metrics_input <-
		subset_model_groups |>
		dplyr::mutate(
			censoring_correction = ifelse(
				grepl("with", model), "censoring correction", "naive model"
			),
			model_structure = ifelse(
				grepl("full", model), "full", "reduced"
			),
			subsample_id = rep(subsample_ids, times = dplyr::n_distinct(model)),
			.keep = "unused"
		) |>
		tibble::add_column(
			data = subset_model_data,
			preds = subset_model_preds
		)	|>
		tidyr::pivot_wider(names_from = model_structure, values_from = preds) |>
		tibble::add_column(
			homologous = h_subset_model_preds
		) |>
		# We only care about these seasons.
		dplyr::filter(season < "2018 - 2019") |>
		# Drop extra values from season
		dplyr::mutate(season = forcats::fct_drop(season)) |>
		# Join with the seroconversion data
		dplyr::left_join(
			seroconversion_data_subsample,
			by = c('season', 'strain_type', 'method', 'subsample_id')
		)

	# Put the two things together so we only need to run the metric suite once
	out <- dplyr::bind_rows(
		full_model_metrics_input,
		subsample_metrics_input
	)

	return(out)

}

# Run the function to create the data set -- this allows for all of the
# intermediate objects to be deleted and garbage collected so we don't have
# them clogging up the RAM or whatever.
metrics_input_data <- create_metrics_data()

# Metric suite calculation ####
metric_df <-
	purrr::pmap(
		metrics_input_data,
		metric_set,
		.progress = "Calculating metrics…"
	)

metrics_output <-
	metrics_input_data |>
	dplyr::select(-data, -full, -reduced, -homologous, -seroconversion) |>
	tibble::add_column(metrics_data = metric_df) |>
	tidyr::unnest(metrics_data)

qs::qsave(
	metrics_output,
	here::here("results", "large-files", "metric-results.qs")
)

test_plt <-
	metrics_output |>
	dplyr::filter(
		season == "2016 - 2017",
		subsample_id != "full data",
		outcome_name == "log_posttiter",
		method == "pepitope",
		strain_type == "H1N1"
	) |>
	dplyr::group_by(censoring_correction, metric_name, metric_set) |>
	dplyr::mutate(metric_value = minmax(metric_value)) |>
	dplyr::ungroup() |>
	dplyr::mutate(
		cens_mod = paste0(
			metric_set, " - ",
			ifelse(censoring_correction == "naive model", "Naive", "Corrected")
		),
		cens_mod = factor(
			cens_mod,
			levels = c("Old - Naive", "New - Naive", "Old - Corrected",
								 "New - Corrected")
		)
	) |>
	tidyr::nest(
		plt_data = c(subsample_id, .draw, cens_mod, metric_value),
		.by = c(season, strain_type, method, outcome_name, metric_name)
	) |>
	dplyr::mutate(
		plt = purrr::map2(
			plt_data, metric_name,
			\(d, t) d |>
				ggplot2::ggplot() +
				ggplot2::aes(
					x = metric_value,
					y = forcats::fct_rev(subsample_id)
				) +
				ggplot2::geom_boxplot(
					outlier.shape = NA,
					color = "navy"
				) +
				ggplot2::geom_point(
					position = ggplot2::position_jitter(),
					alpha = 0.03,
					size = 0.75
				) +
				ggplot2::facet_wrap(
					ggplot2::vars(cens_mod),
					ncol = 4,
					scales = "free"
				) +
				hgp::theme_ms() +
				ggplot2::labs(
					x = "Min-max scaled metric value",
					y = "Stimulated lab",
					title = paste0("Metric: ", t)
				)
		)
	)

library(patchwork)
metrics_plot_draft <-
	purrr::reduce(test_plt$plt, `+`) +
	patchwork::plot_layout(ncol = 1, axes = "collect") &
	ggplot2::theme_minimal(base_size = 14) &
	ggplot2::theme(
		plot.title = ggplot2::element_text(size = 16),
		strip.text = ggplot2::element_text(size = 15),
		axis.text = ggplot2::element_text(color = "black"),
		plot.background = ggplot2::element_rect(fill = "white", color = "white")
	)

ggplot2::ggsave(
	here::here("results", "figures", "real-data-metric-summary.png"),
	plot = metrics_plot_draft,
	width = 16,
	height = 9,
	dpi = 150,
	units = "in"
)

# icc_test <-
# 	metrics_output |>
# 	dplyr::filter(subsample_id != "full data") |>
# 	tidyr::pivot_wider(
# 		names_from = subsample_id,
# 		values_from = metric_value,
# 		names_prefix = "Subsample_"
# 	) |>
# 	tidyr::nest(
# 		subsample_matrix = dplyr::contains("Subsample_"),
# 		.by = c(season, strain_type, method, outcome_name, censoring_correction,
# 						metric_name, metric_set)
# 	)
#
# icc_test2 <-
# 	purrr::map(
# 		icc_test$subsample_matrix[1:6],
# 		psych::ICC,
# 		lmer = TRUE,
# 		missing = FALSE,
# 		.progress = TRUE
# 	)
#
# icc_test3 <-
# 	purrr::map(
# 		icc_test$subsample_matrix,
# 		\(d) irr::icc(
# 			ratings = d,
# 			model = "twoway",
# 			type = "agreement",
# 			unit = "single"
# 		),
# 		.progress = TRUE
# 	)
#
# icc_test4 <-
# 	purrr::map(
# 		icc_test3,
# 		\(l) tibble::tibble(
# 			est = l$value,
# 			lwr = l$lbound,
# 			upr = l$ubound
# 		),
# 		.progress = TRUE
# 	) |>
# 	dplyr::bind_rows()
#
# icc_data <-
# 	dplyr::bind_cols(
# 		icc_test,
# 		icc_test4
# 	)
#
# icc_data |>
# 	dplyr::filter(
# 		season == "2016 - 2017",
# 		strain_type == "H1N1",
# 		outcome_name == "log_posttiter",
# 		method == "pepitope"
# 	) |>
# 	ggplot2::ggplot() +
# 	ggplot2::aes(x = est, xmin = lwr, xmax = upr, y = metric_name,
# 							 color = metric_set) +
# 	ggplot2::geom_pointrange(
# 		position = ggplot2::position_dodge(width = 0.1)
# 	) +
# 	ggplot2::facet_wrap(~censoring_correction) +
# 	hgp::theme_ms()

# END OF FILE ####
