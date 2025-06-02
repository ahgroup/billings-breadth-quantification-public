###
# Simulation Study
# Zane Billings
# 2024-09-03
# Simulate antibody landscapes from a regression model with censoring and
# see how our methods vs the old methods perform across different amounts of
# LoD variables.
###

# Setup ####
box::use(
	cli,
	dplyr,
	here,
	purrr,
	rlang,
	tibble,
	tidyr,
	Hmisc
)

source(here::here("R", "functions", "hai-helpers.R"))
source(here::here("R", "functions", "utils.R"))
source(here::here("R", "functions", "model-runner.R"))
source(here::here("R", "functions", "breadth-metrics.R"))
source(here::here("R", "functions", "simulation-functions.R"))

random_seed <- 4820296L
set.seed(random_seed)

# Simulation ####


# simulated_latent_data <- sim_labs(
# 	individuals_per_lab = 100L,
# 	strains_per_lab = 9L,
# 	n_labs = 10L,
# 	strain_universe = seq(0.02, 1, 0.02),
# 	sim_fun = hierarchical_sim_with_distance,
# 	global_seed = random_seed,
# 	intercept_mu = 4, intercept_sigma = 3,
# 	slope_mu = -10, slope_sigma = 4,
# 	sd = 3
# )
simulated_latent_data <- sim_labs(
	individuals_per_lab = 100L,
	strains_per_lab = 9L,
	n_labs = 10L,
	strain_universe = seq(0.02, 1, 0.02),
	sim_fun = one_titer_sim_lm,
	global_seed = random_seed,
	intercept = 4,
	slope_dist = -6,
	slope_pre = -0.25,
	dist_pre_interaction = 0,
	post_sd = 3,
	pre_mean = 2,
	pre_sd = 1
)
simulated_observed_data <-
	simulated_latent_data |>
	dplyr::select(lab, id, d, sim_pre_titer, sim_post_titer) |>
	# Calculate clinical endpoints
	dplyr::mutate(
		titer_increase = log2(sim_post_titer / sim_pre_titer),
		seroprotection = as.integer(sim_post_titer >= 40),
		seroconversion = as.integer(seroprotection & (titer_increase >= 2))
	)

# Repeat the sim for multiple slopes to get different censoring amounts
slope_vals <- seq(-15, 0, 0.25)
sim_data <-
	slope_vals |>
	purrr::map(
		\(x) {
			sim_labs(
				individuals_per_lab = 100L,
				strains_per_lab = 9L,
				n_labs = 10L,
				strain_universe = seq(0.02, 1, 0.02),
				sim_fun = one_titer_sim_lm,
				global_seed = random_seed,
				intercept = 4,
				slope_dist = x,
				slope_pre = -0.25,
				dist_pre_interaction = 0,
				post_sd = 3,
				pre_mean = 2,
				pre_sd = 1
			)
		}
	) |>
	rlang::set_names(slope_vals) |>
	dplyr::bind_rows(.id = "slope_dist") |>
	dplyr::mutate(slope_dist = as.numeric(slope_dist))

pct_censored <-
	sim_data |>
	dplyr::group_by(slope_dist, lab) |>
	dplyr::summarise(
		pct_censored = mean(sim_post_titer == 5),
		.groups = "drop_last"
	) |>
	dplyr::summarise(
		ggplot2::mean_cl_normal(pct_censored),
		.groups = "drop"
	)

target_slopes <-
	tibble::tibble(
		target = seq(0.2, 0.7, 0.1),
		slope = 	purrr::map_dbl(
			target,
			\(target) pct_censored$slope_dist[which.min(abs(pct_censored$y - target))]
		)
	)

ggplot2::ggplot(pct_censored) +
	ggplot2::aes(x = slope_dist, y = y) +
	ggplot2::geom_hline(
		data = target_slopes,
		ggplot2::aes(yintercept = target),
		lty = 1, color = "darkgray"
	) +
	ggplot2::geom_vline(
		data = target_slopes,
		ggplot2::aes(xintercept = slope),
		lty = 1, color = "darkgray"
	) +
	ggplot2::geom_line() +
	ggplot2::geom_point(size = 3) +
	ggplot2::geom_point(
		data = target_slopes,
		ggplot2::aes(x = slope, y = target),
		color = "red", size = 3
	) +
	ggplot2::scale_x_continuous(
		breaks = seq(-15, 0, 1),
		minor_breaks = NULL,
		limits = c(-15, 0)
	) +
	ggplot2::scale_y_continuous(
		breaks = seq(0.0, 1.0, 0.1),
		limits = c(0.1, 0.8)
	) +
	ggplot2::labs(
		x = "Antigenic distance slope term",
		y = "Percentage of censored data points"
	) +
	hgp::theme_ms()

ggplot2::ggplot(simulated_observed_data) +
	ggplot2::aes(x = d, y = sim_post_titer)  +
	ggplot2::geom_line(
		ggplot2::aes(group = id),
		position = ggplot2::position_jitter(width = 0, height = 0.1),
		alpha = 0.25
	) +
	ggplot2::geom_count(
		shape = 21,
		fill = "white",
		stroke = 1.5
	) +
	ggplot2::scale_y_continuous(
		trans = "log2",
		breaks = hai_to_natural_scale(seq(0, 12, 2)),
		labels = hai_to_natural_scale(seq(0, 12, 2))
	) +
	ggplot2::facet_wrap(
		ggplot2::vars(lab),
		ncol = 2
	) +
	ggplot2::labs(
		x = "Simulated antigenic distance",
		y = "Simulated HAI titer"
	) +
	hgp::theme_ms()

## Data processing for models ####
# We basically need to repeat the steps from the real data analysis.
# The censoring step can be done first, everything else we can do for specific
# models.
model_data_hai_formatted <-
	sim_data |>
	# Get only certain slopes so we don't fit all a billion useless models
	dplyr::filter(slope_dist %in% target_slopes$slope) |>
	# Get only the needed columns
	dplyr::select(slope_dist, lab, id, d, sim_pre_titer, sim_post_titer) |>
	# Calculate clinical endpoints
	dplyr::mutate(
		titer_increase = log2(sim_post_titer / sim_pre_titer),
		seroprotection = as.integer(sim_post_titer >= 40),
		seroconversion = as.integer(seroprotection & (titer_increase >= 2))
	) |>
	# put the HAI data into the correct censoring format
	format_hai_data(
		post_titer = "sim_post_titer",
		log_scale = FALSE,
		log_out = TRUE,
		increase = FALSE
	)

# Model fitting ####
# Next we need to fit all of our models to the simulation. That includes:
# 1) the original four models to the complete dataset
# 2) the two models for heterologous-only data
# 3) the logistic regression models.
# Because of how the code is written it's easiest to do those in three
# separate calls.

## Heterologous data linear models ####
cli::cli_alert_info("Starting heterologous data models.")
# First we get the data nested correctly
# Next do the nesting steps
heterologous_linear_model_data <-
	model_data_hai_formatted |>
	tidyr::nest(data = -c(lab)) |>
	dplyr::mutate(
		# Now select only the columns we'll pass to the model -- this is the data
		# we would give to brms::brm() as the data argument if we were using
		# brms directly.
		brms_data = purrr::map(
			data,
			\(d) d |>
				dplyr::select(x = d, c, y, y2, id = id)
		)
	)

# Now set up the Stan files
source(here::here("R", "functions", "stan-code-generation.R"))
stan_model_files <- paste0(model_path, "/", brms_model_info$file_name, ".stan")

# We need four different random seeds, one for each of our sequential
# parallel runs. doFuture uses the L'Ecuyer-CMRG parallel random seed algorithm
# to ensure reproducibility. These seeds were generated on random.org.
seed_vec <- c(
	975571L,
	458898L,
	170276L,
	872859L
)

# Run the models -- see the model running script for more of an explanation.
het_lin_model_output <-
	run_models(
		model_data = heterologous_linear_model_data,
		model_info = brms_model_info,
		model_files = stan_model_files,
		random_seed = seed_vec,
		dump_file = here::here("results", "large-files", "dump"),
		time_files_dir = here::here(
			"results", "times", "simulation-models", "het-lin"
		),
		stan_csv_pth_base = here::here(
			"results", "large-files", "stanfits", "simulation-models", "het-lin"
		),
		brms_file_directory = here::here(
			"results", "large-files", "brms-fits", "simulation-models", "het-lin"
		)
	)

## Homologous data linear models
# Process the results
het_lin_processed <-
	model_result_processing(
		het_lin_model_output,
		heterologous_linear_model_data,
		here::here(
			"results","large-files", "model-results", "simulation-models", "het-lin"
		)
	)

gc()

## Homologous data-only linear models ####
cli::cli_alert_info("Starting homologous data models.")
# First we get the data nested correctly
# Next do the nesting steps
homologous_linear_model_data <-
	heterologous_linear_model_data |>
	dplyr::mutate(
		# First get a list of which data points will be included
		homologous_obs = purrr::map(
			data,
			\(df) which(df$d == 0)
		),
		# Now filter the three dataset columns to only contain those
		# observations
		dplyr::across(
			c(data, brms_data),
			\(col) purrr::map2(
				col, homologous_obs,
				\(d, ind) d[ind, ]
			)
		)
	)

# Now set up the Stan files
brms_model_info <-
	brms_model_info |>
	dplyr::filter(startsWith(file_name, "reduced"))
stan_model_files <- paste0(
	model_path, "/", brms_model_info$file_name, ".stan"
)

# We need four different random seeds, one for each of our sequential
# parallel runs. doFuture uses the L'Ecuyer-CMRG parallel random seed algorithm
# to ensure reproducibility. These seeds were generated on random.org.
seed_vec <- c(
	216346L,
	686363L
)

# Run the models -- see the model running script for more of an explanation.
hom_lin_model_output <-
	run_models(
		model_data = homologous_linear_model_data,
		model_info = brms_model_info,
		model_files = stan_model_files,
		random_seed = seed_vec,
		dump_file = here::here("results", "large-files", "dump"),
		time_files_dir = here::here(
			"results", "times", "simulation-models", "hom-lin"
		),
		stan_csv_pth_base = here::here(
			"results", "large-files", "stanfits", "simulation-models", "hom-lin"
		),
		brms_file_directory = here::here(
			"results", "large-files", "brms-fits", "simulation-models", "hom-lin"
		)
	)

# Process the results
hom_lin_processed <-
	model_result_processing(
		hom_lin_model_output,
		homologous_linear_model_data,
		here::here(
			"results","large-files", "model-results", "simulation-models", "hom-lin"
		)
	)

gc()

## Logistic regression models ####
cli::cli_alert_info("Starting logistic models.")
# First we get the data nested correctly
# Next do the nesting steps
logistic_model_data <-
	model_data_hai_formatted |>
	tidyr::nest(data = -c(lab)) |>
	dplyr::mutate(
		# Now select only the columns we'll pass to the model -- this is the data
		# we would give to brms::brm() as the data argument if we were using
		# brms directly.
		brms_data = purrr::map(
			data,
			\(d) d |>
				dplyr::select(x = d, y = seroconversion)
		)
	)

# Now set up the Stan files
source(here::here("R", "functions", "logistic-model-generation.R"))
stan_model_files <- paste0(
	model_path, "/", brms_model_info$file_name, ".stan"
)

# We need four different random seeds, one for each of our sequential
# parallel runs. doFuture uses the L'Ecuyer-CMRG parallel random seed algorithm
# to ensure reproducibility. These seeds were generated on random.org.
seed_vec <- c(
	738477L
)

# Run the models -- see the model running script for more of an explanation.
logistic_model_output <-
	run_models(
		model_data = logistic_model_data,
		model_info = brms_model_info,
		model_files = stan_model_files,
		random_seed = seed_vec,
		dump_file = here::here("results", "large-files", "dump"),
		time_files_dir = here::here(
			"results", "times", "simulation-models", "logistic"
		),
		stan_csv_pth_base = here::here(
			"results", "large-files", "stanfits", "simulation-models", "logistic"
		),
		brms_file_directory = here::here(
			"results", "large-files", "brms-fits", "simulation-models", "logistic"
		)
	)

# Process the results
logistic_processed <-
	model_result_processing(
		logistic_model_output,
		logistic_model_data,
		here::here(
			"results","large-files", "model-results", "simulation-models", "logistic"
		)
	)

gc()

# END OF FILE ####
