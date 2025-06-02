###
# Metrics Calculation on Simulation Data
# Zane Billings
# 2024-10-03
# Calculating breadth metrics from the different models on the simulated data
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
source(here::here("R", "functions", "breadth-metrics.R"))

# Data processing ####
# Function for loading all of the data we need since there is A LOT
create_sim_metrics_data <- function() {
	# Run gc on exit since these things are big -- can help performance
	on.exit({gc(verbose = FALSE)}, add = TRUE)

	# Load the heterologous model fits
	het_model_pth <- here::here(
		"results","large-files", "model-results", "simulation-models",
		"het-lin"
	)
	het_model_data <- qs::qread(here::here(het_model_pth, "data.qs"))
	het_model_groups <- qs::qread(here::here(het_model_pth, "groups.qs"))
	het_model_preds <- qs::qread(here::here(het_model_pth, "preds.qs"))

	# Homologous model fits
	hom_model_pth <- here::here(
		"results","large-files", "model-results", "simulation-models",
		"hom-lin"
	)
	hom_model_preds <- qs::qread(here::here(hom_model_pth, "preds.qs"))

	# Logistic model fits
	log_model_pth <- here::here(
		"results", "large-files", "model-results", "simulation-models",
		"logistic"
	)
	log_model_groups <- qs::qread(here::here(log_model_pth, "groups.qs"))
	log_model_preds <- qs::qread(here::here(log_model_pth, "preds.qs"))
	# Because the seroconversion model doesn't care about censoring or which model
	# outcome we used, there are less rows than in the other data frames. So we
	# have to do a bit of cleanup and then join this one to the others.
	seroconversion_data_full_model <-
		log_model_groups |>
		dplyr::select(-model) |>
		dplyr::mutate(seroconversion = log_model_preds)

	# Clean up the full model data
	metrics_input <-
		het_model_groups |>
		dplyr::mutate(
			censoring_correction = ifelse(
				grepl("with", model), "censoring correction", "naive model"
			),
			model_structure = ifelse(
				grepl("full", model), "full", "reduced"
			),
			.keep = "unused"
		) |>
		tibble::add_column(
			data = het_model_data,
			preds = het_model_preds
		) |>
		tidyr::pivot_wider(names_from = model_structure, values_from = preds) |>
		tibble::add_column(
			homologous = hom_model_preds,
			outcome_name = "log_posttiter"
		) |>
		# Join with the seroconversion data
		dplyr::left_join(
			seroconversion_data_full_model,
			by = c('lab')
		)

	return(metrics_input)
}

# Run the function to create the data set -- this allows for all of the
# intermediate objects to be deleted and garbage collected so we don't have
# them clogging up the RAM or whatever.
metrics_input_data <- create_sim_metrics_data()

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
	here::here("results", "large-files", "metric-results-simulation.qs")
)

test_plt <-
	metrics_output |>
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
	dplyr::select(-outcome_name) |>
	tidyr::nest(
		plt_data = c(lab, .draw, cens_mod, metric_value),
		.by = c(metric_name)
	) |>
	dplyr::mutate(
		plt = purrr::map2(
			plt_data, metric_name,
			\(d, t) d |>
				ggplot2::ggplot() +
				ggplot2::aes(
					x = metric_value,
					y = forcats::fct_rev(lab)
				) +
				ggplot2::geom_boxplot(
					outlier.shape = NA,
					color = "navy"
				) +
				ggplot2::geom_point(
					color = "firebrick3",
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
	here::here("results", "figures", "sim-data-metric-summary.png"),
	plot = metrics_plot_draft,
	width = 16,
	height = 9,
	dpi = 150,
	units = "in"
)

metrics_output |>
	dplyr::filter(
		metric_name == "Total Strength",
		censoring_correction == "censoring correction"
	) |>
	dplyr::summarize(
		covar(metric_value),
		.by = c(outcome_name, metric_set)
	)


# END OF FILE ####
