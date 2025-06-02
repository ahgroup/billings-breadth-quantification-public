###
# Overall simple linear models
# Zane
# 2022-02-21
# In this script, I'll fit simple linear regression models for the
# relationship between the outcomes and the antigenic distance that
# do not take any other covariates into account. For this script, we'll
# consider three different models: complete pooling, no pooling, and partial
# pooling.
###

# SETUP ####

box::use(
	readr,
	tidyr,
	dplyr,
	cmdstanr[...],
	rstan[...]
)

source(here::here("R", "4_modeling", "modeling-helpers.R"))

# Data cleaning ####
model_data <- readr::read_rds(
	here::here("Data", "Processed", "Modeling-Data", "model-data.Rds")
)

dat_stan <-
	model_data |>
	dplyr::mutate(
		dat = purrr::map(
			dat,
			\(x) x |>
				dplyr::select(y, id, norm_dist) |>
				as.list()
		),
		dat = purrr::map(
			dat,
			\(x) c(x, N = length(x$id), k = max(x$id))
		)
	)

# For these models we need to do a bit more model transformation.

# Model setup ####
# See either the Stan code or the writeup for an explanation of the model.
model_list <- list(
	"correlated-parameters" =
		cmdstanr::cmdstan_model(
			stan_file = here::here("Stan", "partial-pooling-correlated-parms.stan"),
			exe_file  = here::here("Stan", "partial-pooling-correlated-parms.exe"),
			compile = FALSE
		),
	"partial-pooling" =
		cmdstanr::cmdstan_model(
			stan_file = here::here("Stan", "partial-pooling-across-seasons.stan"),
			exe_file  = here::here("Stan", "partial-pooling-across-seasons.exe"),
			compile = FALSE
		),
	"no-pooling" =
		cmdstanr::cmdstan_model(
			stan_file = here::here("Stan", "no-pooling-across-seasons.stan"),
			exe_file  = here::here("Stan", "no-pooling-across-seasons.exe"),
			compile = FALSE
		),
	"complete-pooling" =
		cmdstanr::cmdstan_model(
			stan_file = here::here("Stan", "complete-pooling-across-seasons.stan"),
			exe_file  = here::here("Stan", "complete-pooling-across-seasons.exe"),
			compile = FALSE
		)
)

# Unfortunately because rstan is object-oriented, this function is a
# mutator. :(
compile_model_list(model_list)

# Model fitting ####

CTRL <-
	list(
		seed = 370,
		threads_per_chain = 2,
		chains = 8,
		parallel_chains = 8,
		iter_warmup = 700,
		iter_sampling = 700,
		adapt_delta = 0.8
	)

# test_loo <-
# 	purrr::map(
# 		test_fit_l,
# 		\(x) {
			# ll <- x$draws("log_lik")
			# r_eff <- loo::relative_eff(exp(ll), cores = 8)
			# psis <- loo::loo(ll, r_eff = r_eff, cores = 8)
			# return(psis)
# 		}
# 	)
#
# test_loo_compare <- loo::loo_compare(test_loo)

run_all_model_batches(
	dat_list = list(dat_stan[1:6, ]),
	model_list = model_list,
	dir_list = as.list(c(
		"SimpleModels/correlated-parameters",
		"SimpleModels/partial-pooling",
		"SimpleModels/no-pooling",
		"SimpleModels/complete-pooling"
	)),
	TRUE, TRUE, ctrl = CTRL
)

# Next test PSIS for each of the rows. We want to see if the correlated
# parameters model is actually better in all cases.
ncores <- max(1L, parallel::detectCores() - 2L, na.rm = TRUE)
fn <-		list.files(
	here::here(
		"Results", "_Out", "SimpleModels"
	),
	full.names = TRUE,
	recursive = TRUE
)
test_psis_all <-
	purrr::imap(
		fn,
		\(x, idx) {
			paste0("Starting model ", idx, " of ", length(fn), ".\n") |>
				crayon::white() |>
				cat()

			suppressWarnings(rm(m, ll, r_eff, psis))
			invisible(gc())

			m <- readr::read_rds(x)
			message("ℹ", " Read file")
			ll <- loo::extract_log_lik(m, merge_chains = FALSE)
			message("ℹ", " Got log likelihood")
			r_eff <- loo::relative_eff(exp(ll), cores = 16)
			message("ℹ", " Got relative effective n")
			psis <- loo::loo(ll, r_eff = r_eff, cores = 16)
			message("ℹ", " Got PSIS")
			return(psis)
		},
		.progress = "Importance sampling"
	)

test_psis_res <-
	purrr::map(
	test_psis_all,
	\(x) x$estimates |>
		as.data.frame() |>
		tibble::rownames_to_column("metric") |>
		tibble::tibble()
) |>
	rlang::set_names(
		fn |>
			stringr::str_remove("^.*/SimpleModels/") |>
			stringr::str_remove("/Model.*")
	) |>
	dplyr::bind_rows(.id  = "model") |>
	dplyr::filter(metric == "looic") |>
	dplyr::select(-metric) |>
	dplyr::bind_cols(
		replicate(4, dat_stan[1:6, ], simplify = FALSE) |> dplyr::bind_rows()
	)

test_psis_res |>
	ggplot(aes(x = model, y = Estimate, ymin = Estimate - 2 * SE, ymax = Estimate + 2 * SE)) +
	geom_pointrange() +
	facet_grid(method ~ outcome) +
	zlib::theme_ms()

readr::write_rds(test_psis_all, here::here("psis-test-simple-models-first-6.Rds"))


# Model predictions processing ####

# END OF FILE ####
