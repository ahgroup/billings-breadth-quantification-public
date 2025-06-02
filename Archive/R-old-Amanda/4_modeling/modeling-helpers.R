###
# Modeling Helper Functions
# Zane
# 2023-02-21
# This script contains functions that are used to fit the models.
###

#' Get completed version of filename to save a model
#'
#' @param home the top-level directory. dir should be a subdirectory of home.
#' This defaults to the Rproj location so you don't necessarily have to
#' type the whole here::here(fasjdhflkasjdflkajsdfkl) stuff out each time.
#' @param dir the current-level directory where model files will be saved
#' @param idx index of the current model
#' @param ext the file extension, defaults to .Rds
#'
#' @return A string with the concatenated filename.
#' @export
#'
get_fn <- function(home = here::here(), dir, idx, ext = ".Rds") {
	paste0(
		home, "/", dir, "/Model",
		stringr::str_pad(idx, width = 2, side = "left", pad = "0"), ext
	)
}

# TODO update this function so that it saves the diagnostics to an actual
# data file, not a random text file that is hard to deal with.
#' Save rstan chain diagnostics to file
#'
#' @param fit a fitted `rstan` model
#' @param fn filename for where to save the diagnostics rules
#'
#' @return The fitted model object, invisibily
#' @export
#'
get_rstan_diagnostics <- function(fit, fn) {
	sink(fn)
	sink(stdout(), type = "message")
	rstan::check_hmc_diagnostics(fit)
	sink()
	sink(stderr(), type = "message")
	invisible(fit)
}

#' Fit a cmdstanr model to every dataframe in a list.
#'
#' @param dat A list of data frames.
#' @param model The cmdstan model to fit.
#' @param dir The directory where model results should be saved (note that
#' you MUST save files to a directory, the model results are far too large
#' to be held in memory simultaneously on a normal machine.)
#' @param rstan_out Should model files be saved as rstan objects? If set to
#' anything other than TRUE, results will be saved as cmdstanr objects.
#' The result of the functions I wrote typically require the result to be
#' saved as an rstan object, similar to how Richard McElreath's rethinking
#' package works.
#'
#' @return The cmdstan model (invisibly). All other effects of this function
#' are through side effects.
#' @export
#'
fit_all_models <- function(dat, model, dir, rstan_out = TRUE, verbose = TRUE,
													 ctrl = NULL, ...) {
	starttime <- Sys.time()

	# Function for sample from the model at each iteration of the loop
	sample_from_model <- function(d, id, n, model, verbose = TRUE, ctrl = NULL) {
		if(isTRUE(verbose)) {
			paste0("Starting model ", id, " of ", n, "!\n") |>
				crayon::white() |>
				crayon::bgBlue() |>
				cat()
		}

		if(is.null(ctrl)) {
			ctrl <- list(
				seed = 370,
				threads_per_chain = 2,
				chains = 8,
				parallel_chains = 8,
				iter_warmup = 1250,
				iter_sampling = 1250,
				adapt_delta = 0.95
			)
		}

		sample_args <- c(data = list(d), ctrl)
		out <- do.call(model$sample, sample_args)

		return(out)
	}

	# Easy error handling for the model, if it runs great, if there is an
	# error in the fit it will return "uh-oh!" without crashing the
	# entire loop
	possibly_sample_from_model <- purrr::possibly(sample_from_model, "uh-oh!")

	# Need this constant for message printing
	n_mods <- nrow(dat)

	purrr::iwalk(
		dat[["dat"]],
		\(d, idx) {
			# Get the file name to save at
			fn <- paste0(
				here::here("Results", "_Out"), "/", dir, "/Model",
				stringr::str_pad(idx, side = "left", width = "2", pad = "0"), ".Rds"
			)

			# Invoke the model fitting routine
			mod <- possibly_sample_from_model(d, idx, n_mods, model, verbose, ctrl)

			# Fit models with cmdstan, but then read the model back as an Rstan object
			# So I can use the rstan::extract function that McElreath uses to get the
			# output in the same format as rethinking.
			if (isTRUE(rstan_out)) {
				mod <- rstan::read_stan_csv(mod$output_files())
			}

			# Save the results to disk
			readr::write_rds(mod, fn, compress = "gz")

			# Cleanup this iteration of the loop to save ram
			rm(mod)
			invisible(gc())
		}
	)
	stoptime <- Sys.time()

	message("Fitting took ", difftime(stoptime, starttime) |> format())

	invisible(model)
}

#' Run a batch of models. Parallel in the same sense as pmin/pmax.
#'
#' @param dat_list A list of data frames to fit on.
#' @param model_list A list of models to fit.
#' @param dir_list A list of directories where files should be saved.
#' @param ... Other arguments to pass to `fit_all_models()`.
#'
#' All arguments in this list should be ENCLOSED LISTS. E.g. you should pass
#' a list of length 1 containing a list of 84 data frames instead of just a'
#' list of length 84 containing 84 data frames. Then the model in the first
#' slot of the model list will be fit on EACH dataframe in the first list.
#'
#' @return
#' @export
#'
#' @examples
run_all_model_batches <- function(dat_list, model_list, dir_list, ctrl, ...) {
	start_time <- Sys.time()
	# Get lengths
	ld <- length(dat_list)
	lm <- length(model_list)
	lf <- length(dir_list)

	# dir_list MUST BE the same length as the maximum length of either
	# dat_list or model_list!
	if (lf != max(ld, lm)) {
		stop("You need to specify one directory per dat/model combo.")
	}

	# If dat_list and model_list are not the same length and one of them
	# is not length one, need to stop.
	if(ld != lm) {
		if(!(ld == 1) & !(lm == 1)) {
			stop("Must specify the same number of data sets and models ",
					 "or one of them should be length one.\n",
					 "If you wanted to cross data/models you should do that ",
					 "before passing to this function.")
		}
		}

	# Pwalk over the three lists, fitting one model per combo.
	purrr::pwalk(
		.l = list(d = dat_list, m = model_list, f = dir_list,
							idx = seq_along(dir_list)),
		.f = \(d, m, f, idx) {
			paste0("Starting batch ", idx, " of ", lf, "!\n") |>
				crayon::bgRed() |>
				cat()

			fit_all_models(
				dat = d,
				model = m,
				dir = f,
				ctrl = ctrl,
				...
			)
		}
	)

	stop_time <- Sys.time()
	message("Batch fitting took ", difftime(stop_time, start_time) |> format(),
					"total.")

	invisible(model_list)

}

#' Compile a list of cmdstanr models
#'
#' @param model_list A list of cmdstanr model objects.
#' @param compile_args Arguments to pass to `x$compile()` from `cmdstanr`.
#'
#' @return The list, invisibily.
#' @export
#'
compile_model_list <- function(model_list, compile_args = NULL) {
	purrr::walk2(
		model_list, seq_along(model_list),
		\(x, idx) {
			paste0("Compiling model: ", names(model_list)[[idx]], " (", idx, " of ",
						 length(model_list), ")\n") |>
				crayon::bgMagenta() |>
				cat()

			if (is.null(compile_args)) {
				compile_args <- list(
					cpp_options = list(stan_threads = TRUE),
					quiet = TRUE,
					pedantic = TRUE
				)
			}

			do.call(x$compile, compile_args)
		}
	)

	invisible(model_list)
}

# END OF FILE ####
