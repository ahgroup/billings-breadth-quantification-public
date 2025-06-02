###
# Plan for parallelism via future
# Zane Billings
# 2024-08-20
# Sets up a plan for parallel computing using the R future package based on
# a safe number of detected cores.
###

rlang::inform("Setting up plan for future")
suppressPackageStartupMessages({
	library(future)
	library(doFuture)
})

# Get the number of available cores with some restrictions so the computer
# is still usable while the jobs are running
N_cores <- parallelly::availableCores(constraints = "connections", omit = 2L)

# Number of cores to use for each level of the nested scheme.
# From the future documentation, we know that the total number of cores that
# will be used is
# 1 + N_outer * (1 + N_inner).
# For fixed N_inner = 4, we want to solve
# N_cores >= 1 + N_outer * (5) giving us
# N_outer <= floor((N_cores - 1) / 5)
# But we want to make sure anyone running this code on their computer doesn't
# auto-trigger an error, so use the most cores we can up to 4.
N_inner <- min(N_cores, 4)
N_outer <- max(floor((N_cores - 1) / 5), 1)

# Check that the total number of cores we will use is no more than the number of
# cores we said were available
check_req_cores <- function(required_cores, total_cores) {
	rlang::inform(c("i" = paste0(
		"Total cores required: ", required_cores, "/", total_cores
	)))
	if (required_cores > total_cores) {
		rlang::warn(c("x" = "Tried to use more cores than we got!"))
		rlang::inform(
			c("i" = "Falling back to sequential plan.")
		)
		future::plan(sequential)
	} else {
		rlang::inform(c("v" = paste0("Good to go!")))
	}
}

# Set up the parallel processing plan based on the total number of cores.
if (N_cores <= 4) {
	# If there's less than 4 cores, do everything sequentially.
	req_cores <- 1
	check_req_cores(req_cores, N_cores)
	rlang::inform(
		c("i" = "Using completely sequential plan due to number of cores.")
	)
	future::plan(sequential)
} else if (N_cores <= 8) {
	# If there's only 8 cores, use only the inner parallelization plan so each
	# model still runs chains in parallel, but models are executed sequentially.
	req_cores <- 6
	check_req_cores(req_cores, N_cores)
	rlang::inform(
		c("i" = "Running models sequentially with four parallel chains.")
	)
	future::plan(
		list(
			sequential,
			tweak(multisession, workers = N_inner)
		)
	)
} else {
	# If we got more than 8 cores available we'll run our nested parallel
	# plan where multicore jobs are also run in parallel.
	req_cores <- 1 + N_outer * (1 + N_inner)
	check_req_cores(req_cores, N_cores)
	rlang::inform(c("i" = paste0(
		"Using nested parallel structure with 4 parallel chains per model and ",
		N_outer,
		" parallel modeling jobs."
	)))
	future::plan(
		list(
			tweak(multisession, workers = N_outer),
			tweak(multisession, workers = N_inner)
		)
	)
}

# END OF FILE ####
