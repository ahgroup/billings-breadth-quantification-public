# .Rprofile -- this will be loaded before every new session.
# Load renv
source("renv/activate.R")

tv <- function() {
	targets::tar_visnetwork(T)
}

td <- function() {
	targets::tar_make(use_crew = FALSE, callr_function = NULL, as_job = FALSE)
}

tm <- function() {
	targets::tar_make()
}

tc <- function() {
	targets::tar_make(use_crew = FALSE)
}

# END OF FILE
