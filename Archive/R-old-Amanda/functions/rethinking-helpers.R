###
# Functions for dealing with rethinking models
# Zane
# 2023-02-09
###


# Function for processing the posterior samples
get_preds <- function(model, sim_dat, n_samps = Inf, f = "link", ...) {
	# This is a memory intensive operation, so force the garbage collector
	# to run. Lots of people say this doesn't do anything but I don't
	# believe them.
	invisible(gc())

	# Get the predictions
	out <- list()
	if (f == "link") {
		out$samples <- rethinking::link(model, data = sim_dat, n = n_samps, ...)
	} else if (f == "sim") {
		out$samples <- rethinking::sim(model, data = sim_dat, n = n_samps, ...)
	} else {
		stop("f must be 'link' or 'sim'.")
	}

	# Get the mean and PI
	out$means <- colMeans(out$samples)
	out$PIs <- apply(out$samples, 2, rethinking::PI)

	return(out)
}

batch_process_preds <- function(batch_size, model, sim_dat, n_samps, f = "link",
																...) {

	# Get size of data to simulated on
	N <- nrow(sim_dat)

	# bp stands for "batching points" -- the indices of the dataframe that define
	# the series of bathes. We have to add N to the end to make sure we grab the
	# last batch which is allowed to be smaller.
	bp <- c(seq(1, N, batch_size), N)

	# So now bp[i, i + 1] will get the indices for the next batch.
	# That is sim_dat[bp[i, i+1], ] is the i-th batch of data which contains
	# batch_size number of rows.
}

# Fcn to put this stuff into tibble
clean_list <- function(ls, x = NULL) {
	out <-
		tibble::tibble(
			est = ls$means,
			lwr = ls$PIs[1, ],
			upr = ls$PIs[2, ]
		)
	if (!is.null(x)) {
		out <- tibble::add_column(out, x)
	}
	return(out)
}

# Conditional effects -- population level predictions that conditional on
# See rethinking pg 428
get_marginals <- function(x, post, type = 1, avar = "a", bvar = "b",
													inv_link = identity, ...) {
	# Global grand mean CI
	if (type == 1) {
		a_bar <- rowMeans(post[[avar]])
		b_bar <- rowMeans(post[[bvar]])
	}

	# New typical individual CI
	if (type == 3) {
		a_bar <- apply(post[[avar]], 1, sample, size = 1)
		b_bar <- apply(post[[bvar]], 1, sample, size = 1)
	}

	# Container for output
	out <- list()

	# Link calculation with simulated x values
	out$samples <- sapply(x, \(x) inv_link(a_bar + b_bar * x))

	# Mean and CI calculation
	out$means <- colMeans(out$samples)
	out$PIs <- apply(out$samples, 2, \(x) rethinking::PI(x, ...))

	return(out)
}




# END OF FILE ####
