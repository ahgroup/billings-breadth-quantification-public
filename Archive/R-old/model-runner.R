###
# Model Running Code
# Zane Billings
# 2024-08-22
# Function that runs all four of the main models on a list of input data
# frames.
####

# Declare dependencies
box::use(
	brms,
	cmdstanr,
	crayon,
	dplyr,
	furrr,
	future,
	here,
	progressr,
	purrr,
	qs,
	readr,
	rlang,
	stringr,
	tibble,
	tidybayes,
	tidyr
)

# Set the handler for progress bars
# You can edit this if you want them to look different or make noise
# I like cli because it gives accurate ETA estimates
progressr::handlers(progressr::handler_txtprogressbar(clear = FALSE))

# Now we'll run the fitting routine once for every row of nested_model_data.
# Each model will run 4 NUTS chains in parallel, since the computer has more
# than 4 cores we can run multiple model jobs in parallel using a nested
# parallelization scheme.
# We need to do this separately for each of the four models we want to fit, so
# we'll do those sequentially.

## Model fitting loop ####
# Now we'll loop through each of the models. Doing this as a loop because it's
# a bit complicated and unnecessary to make it a function for a map call.
# We put the entire loop inside of the function to avoid internal loop
# variables being assigned to the global scope.
run_models <- function(
		model_data,
		model_info,
		model_files,
		random_seed,
		dump_file = here::here("results", "large-files", "dump"),
		time_files_dir,
		stan_csv_pth_base,
		brms_file_directory,
		parallel = TRUE,
		models_to_run = seq_along(model_files),
		max_brms_draws = NULL
	) {

	# Check that BRMS_DRAWS global exists
	if (is.null(max_brms_draws)) {
		if (exists("BRMS_DRAWS")) {
			cli::cli_alert_info("Getting number of draws from BRMS_DRAWS global.")
			max_brms_draws <- BRMS_DRAWS
		} else {
			cli::cli_alert_info("Setting number of brms draws to default, 1000.")
			max_brms_draws <- 1000L
		}
	} else if (!((is.numeric(max_brms_draws)) & (length(max_brms_draws) == 1))) {
		cli::cli_alert_warning("brms_draws not found.")
		cli::cli_alert_info("Setting number of brms draws to default, 1000.")
		max_brms_draws <- 1000L
	}

	oa_st <- Sys.time()
	paste0("Starting model setup!") |>
		cli::cli_h1()
	paste0("Current time: ", format(oa_st, "%Y-%m-%d %H:%M")) |>
		cli::cli_alert_info()
	paste0("Starting setup steps.") |>
		cli::cli_alert_info()

	# Create an output list to hold the final stuff
	interp_list <- vector(mode = "list", length = length(stan_model_files))
	preds_list <- vector(mode = "list", length = length(stan_model_files))

	# Create a temporary location for the parallel process to dump the
	# CSV files created by cmdstanr
	cli::cli_alert_info(
		"Dump file path: {.path {dump_file}}"
	)
	cli::cli_alert_danger("Stop now if you don't want that file emptied!")
	Sys.sleep(5)
	unlink(dump_file, recursive = TRUE) # CLEARS THE DUMP FILE!!!
	dir.create(dump_file, showWarnings = FALSE, recursive = TRUE)

	# Get total number of models to process
	n <- length(model_files)

	# Width to pad model names
	chr_width <- nrow(model_data) |> nchar()

	# Make sure directories exist
	dir.create(time_files_dir, showWarnings = FALSE, recursive = TRUE)

	# Setup parallel processing
	if(isTRUE(parallel)) {
		source(here::here("R", "functions", "parallel-plan.R"))
	} else {
		rlang::inform(
			c("i" = "parallel arg is not TRUE, running models sequentially.")
		)
		suppressPackageStartupMessages({
			library(future)
			library(doFuture)
		})
		future::plan("sequential")
	}

	# This part runs the model fitting loop -- one iteration for each of the
	# four models.
	for (i in seq_along(model_files)) {
		if(!(i %in% models_to_run)) next

		current_model_file <- model_files[[i]]
		cli::cli_h2(c(
			crayon::bold(crayon::green(paste("Starting model", i, "of", n)))
		))
		cli::cli_alert_info("Model file: {.path {current_model_file}}")

		# Now get the current brms formula and model name
		current_bf <- model_info$formula[[i]]
		model_name <- model_info$file_name[[i]]

		# Now compile the current model and set it up for cmdstanr
		cli::cli_alert_info("Compiling stan model: ")
		current_model <- cmdstanr::cmdstan_model(current_model_file)

		# Create a list of data input files, a model will be fit to each one. We can
		# use the convenience function make_standata() from brms to create correctly
		# formatted data for cmdstanr.
		input_list <- purrr::map(
			model_data$brms_data,
			\(d) brms::make_standata(current_bf, data = d)
		)

		# Create list of file names for saving results
		model_number <-
			stringr::str_pad(
				1:length(input_list),
				width = chr_width, side = "left", pad = "0"
			)

		# Now fit the models, saving the CSV files to the holding location
		cli::cli_alert_info("Starting model fitting to each dataset.")
		start <- Sys.time()
		progressr::with_progress({
			p <- progressr::progressor(
				along = input_list
			)

			fitted_models <-
				foreach(
					dataset = input_list
				) %dofuture% {
					p()
					fit_model_to_dataset(
						dataset,
						current_model,
						cmdstan_sampling_arguments,
						parallel_chains = future::nbrOfWorkers(),
						refresh = 0,
						diagnostics = NULL,
						show_messages = FALSE,
						show_exceptions = FALSE,
						output_dir = dump_file
					)
				} %seed% {
					random_seed[[i]]
				} %packages% {
					'cmdstanr'
				}
		})
		end <- Sys.time()
		duration <- difftime(end, start)

		# Save the time stamp to file
		readr::write_rds(
			duration,
			here::here(time_files_dir,
								 paste0(model_name, "-time.Rds"))
		)

		# Model postprocessing
		# Now we have the annoying task of processing the models. First we'll save
		# the cmdstan CSV for each one.
		# With the fitted models, we need to get the useful stuff out
		# of them, which is namely the predictions. Doing that with a cmdstan model
		# is a bit annoying, it's actually easier to convert it into a brms object
		# (since we used brms to generate the stan code anyways and didn't edit the
		# likelihood part of the stan code) and use the postprocessing stuff that is
		# built into brms.

		# Save all of the cmdstan csvs to the known location -- it is really
		# annoying to have to do this in two steps, but we do AFAIK because we
		# don't get granular control over the names in dopar.
		# Note that this step is typically very fast.
		cli::cli_alert_info("Saving cmdstan CSVs to file.")
		model_csv_dirs <- paste0(
			here::here(stan_csv_pth_base, model_name),
			"/model", model_number
		)
		progressr::with_progress({
			p <- progressr::progressor(
				along = fitted_models
			)

			csv_out <-
				foreach(
					model_dir = model_csv_dirs,
					model = fitted_models
				) %dofuture% {
					p()
					dir.create(model_dir, showWarnings = FALSE, recursive = TRUE)
					suppressMessages(
						model$save_output_files(
							dir = model_dir,
							basename = "model",
							timestamp = FALSE,
							random = FALSE
						)
					)
				} %seed% {
					random_seed[[i]]
				} %packages% {
					'cmdstanr'
				}
		})

		# Now we have to use the brms function to load those CSVs that we just wrote
		# to a known location. This actually loads them as rstan 'stanfit' objects
		# instead of cmdstanr objects.
		cli::cli_alert_info("Loading cmdstan CSVs into stanfit objects.")
		progressr::with_progress({
			p <- progressr::progressor(
				along = fitted_models

			)

			stanfit_objects <-
				foreach(
					dir = model_csv_dirs,
					model = fitted_models
				) %dofuture% {
					p()
					list.files(dir, full.names = TRUE) |>
						brms::read_csv_as_stanfit(
							variables = model$metadata()$stan_variables,
							model = current_model,
							exclude = "",
							algorithm = "sampling"
						)
				} %seed% {
					random_seed[[i]]
				} %packages% {
					'cmdstanr'
				}
		})

		# Now those stanfit objects contain all the fits, so we can build empty
		# brms objects using the same formula and data, and shove the stanfit object
		# directly into where brms stores the samples.
		cli::cli_alert_info("Converting stanfit objects to brms objects.")
		progressr::with_progress({
			p <- progressr::progressor(
				along = stanfit_objects
			)

			brms_import <-
				foreach(
					samples = stanfit_objects,
					data = model_data$brms_data
				) %dofuture% {
					p()
					cmdstan_to_brms(
						stanfit = samples,
						brms_formula = current_bf,
						brms_data = data
					)
				} %seed% {
					random_seed[[i]]
				} %packages% {
					'cmdstanr'
				}
		})

		# Save the brms models to file in case we need them again
		cli::cli_alert_info("Saving brms objects to file.")
		cur_brms_dir <- here::here(brms_file_directory, model_name)
		dir.create(cur_brms_dir, showWarnings = FALSE, recursive = TRUE)
		brms_file_names <- paste0(
			cur_brms_dir,
			"/model", model_number
		)
		progressr::with_progress({
			p <- progressr::progressor(
				along = brms_import
			)

			brms_filesave <-
				foreach(
					model = brms_import,
					model_dir = brms_file_names
				) %dofuture% {
					p()
					readr::write_rds(
						x = model,
						file = paste0(model_dir, ".qs"),
						compress = "xz"
					)
				} %seed% {
					random_seed[[i]]
				}
		})


		# Now we can get the predictions for all of our models.
		# First, get the predictions on the interpolated data.
		# That means we need to create the fake data that we want to predict on.
		# This is for plotting, so we want to interpolate.
		# If we're using the full model with antigenic distance in it, we want
		# a nice grid of interpolated distance values to make a plot.
		#preds_newdata <- tibble::tibble(x = seq(0, 1, 0.01))

		# Now actually get the predictions
		rlang::inform(c("i" = "Getting predictions on interpolated data."))
		progressr::with_progress({
			p <- progressr::progressor(
				along = input_list,
				label = "Fitting models",
			)

			brms_interp <-
				foreach(
					m = brms_import
				) %dofuture% {
					p()
					tidybayes::epred_draws(
						m, newdata = tibble::tibble(x = seq(0, 1, 0.01)),
						ndraws = min(max_brms_draws, brms::ndraws(m)),
						re_formula = NA
					)
				} %seed% {
					random_seed[[i]]
				}
		})

		# Next, get the predictions on the actual data.
		# rlang::inform(c("i" = "Getting predictions on actual data."))
		# brms_preds <- purrr::map2(
		# 	brms_import, model_data$brms_data,
		# 	\(m, d) {
		# 		tidybayes::epred_draws(
		# 			m, newdata = d,
		# 			ndraws = min(BRMS_DRAWS, brms::ndraws(m)),
		# 			re_formula = NA
		# 		)
		# 	}
		# )

		# End the parallel processing and close the cluster workers
		future::plan(sequential)

		# Save them into the lists for storage
		interp_list[[i]] <- brms_interp
		#preds_list[[i]] <- brms_preds

		rlang::inform(c(
			"v" = crayon::green(paste("Finished model", i, "of", n, "\n"))
		))

		# end of this huge loop
	}

	# Return the two prediction lists that we need to process further
	out <- list(
		"interpolated" = interp_list,
		"actual" = preds_list
	)

	cli::cli_h1("Successfully completed all model runs!")
	oa_et <- Sys.time()
	paste0("Current time: ", format(oa_et, "%Y-%m-%d %H:%M")) |>
		cli::cli_alert_info()
	readr::write_rds(difftime(oa_et, oa_st), here::here(time_files_dir, "overall.Rds"))
	return(out)
}

model_result_processing <- function(
	model_output,
	data,
	storage_dir
) {
	rlang::inform(c("!" = "Starting postprocessing."))
	model_preds <-
		# First, we need to add the predictions (on both the interpolated and
		# real data) and add them to the nested model data (the input data). We need
		# to do this once for each of the four models we fit. We'll also add a column
		# for the name, so we can tell which model the predictions came from.
		purrr::pmap(
			# pmap() requires a list as input, and we name the variables because we can
			# use the names to make sure we're referring to the correct variable in our
			# mapped function. We need to pass in the two different data lists and the
			# name of the model, all of these things are length (number of models).
			list(
				"interp" = model_output$interpolated,
				#"preds" = model_output$actual,
				"name" = brms_model_info$file_name
			),
			# In the mapped function, we add the interp and preds lists as columns to
			# the nested model data, so the predictions go with the data that the model
			# was fit to. We also add the ID column for the model name, which is (at
			# the time this function sees it) a length one character vector that gets
			# recycled.
			\(interp, preds, name) {
				data |>
					dplyr::select(-brms_data) |>
					tibble::add_column(
						interp_preds = interp,
						#data_preds = preds,
						model = name
					)
			}
		) |>
		# Now we have one data frame for each of the models that we fit, and since
		# they have a column identifying the model we can flatten the structure into
		# one tibble which is easier to use and more efficient to store.
		dplyr::bind_rows()



	# Next we have to format the data for saving. Make sure the destination
	# file exists.
	dir.create(
		storage_dir,
		showWarnings = FALSE,
		recursive = TRUE
	)

	rlang::inform(c("i" = "Formatting data."))
	# We want to do the expansion and saving for the data, interp, and preds
	# columns so we'll use purrr.
	model_preds_list <-
		list(
			model_preds |>
				dplyr::select(-data, -interp_preds),
			model_preds$data,
			model_preds$interp_preds
		)

	rlang::inform(c("i" = "Saving files to disk."))
	# Save the expanded results to disk, using qs "quick serialization" format
	purrr::walk2(
		paste0(
			storage_dir,
			"/", c("groups", "data", "preds"), ".qs"
		),
		model_preds_list,
		\(pth, x) qs::qsave(x, pth)
	)

	rlang::inform(c("v" = "Finished postprocessing successfully."))
	# Invisibility return the expanded preds, since the side effect is more
	# important than the return.
	invisible(model_preds_list)
}

# END OF FILE ####
