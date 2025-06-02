# furrr test for cluster replacement
# Zane
# 2024-08-07

library(future)
library(doFuture)

N_cores <- parallelly::availableCores(constraints = "connections", omit = 1L)

input_list <- list.files(
	here::here("data", "cluster", "input"),
	full.names = TRUE
)

# Number of cores to use
# Outer has length(input_list) number of workers, so uses length(input_list) + 1
# cores total.
N_outer <- length(input_list)
# From the future documentation, we know that the total number of cores that
# will be used is
# 1 + N_outer * (1 + N_inner).
# Solving N_cores <= that, we get that for a fixed N_outer,
# N_inner <= N_cores / N_outer - 1, rounded down to get an integer.
N_inner <- N_cores %/% N_outer - 1
# Check that the total number of cores we will use is no more than the number of
# cores we said were available
req_cores <- 1 + N_outer * (1 + N_inner)
if (req_cores > N_cores) {
	rlang::abort("Tried to use more cores than we got!")
} else {
	rlang::inform(paste0("Total cores required: ", req_cores, "/", N_cores))
}

future::plan(
	list(
		tweak(multisession, workers = N_outer),
		tweak(multisession, workers = N_inner)
	)
)



dimension_test <- function(file, cores_each) {
	library(base)
	library(Racmacs)

	g <- readRDS(file)

	map_g <- Racmacs::acmap(
		titer_table = g[,-1],
		ag_names = g[[1]],
		sr_names = colnames(g)[-1]
	)

	#Dimension Testing
	d <- Racmacs::dimensionTestMap(
		map = map_g,
		options = list(num_cores = cores_each)
	)

	return(d)
}

start <- Sys.time()
res <-
	foreach(
		file = input_list
	) %dofuture% {
		dimension_test(file, nbrOfWorkers())
	} %seed% {8618654}
end <- Sys.time()


# tic()
# test <- furrr::future_map(
# 	input_list,
# 	dimension_test,
# 	cores_each = N_c,
# 	.options = furrr::furrr_options(seed = 3750000),
# 	.progress = TRUE
# )
# toc(log = TRUE)
names(res) <- gsub(
	"(D:/Research/Current/Breadth-Quant/data/cluster/input/|.rds)",
	"",
	input_list
)
out <- dplyr::bind_rows(res, .id = "map") |> tibble::as_tibble()

saveRDS(out, here::here("furrr-test.Rds"))
(elapsed <- difftime(end, start, units = "mins"))
saveRDS(elapsed, here::here("furrr-time.Rds"))
