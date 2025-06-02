# CEnsoring format tester

source(here::here("censoring-format.R"))
test <- tibble::tribble(
	~pretiter, ~posttiter,
	5, 20,
	20, 5,
	5, 5,
	20, 20,
	10, 40,
	40, 10
)

# raw brms format
# Expected answers
# c = c('interval', 'left', 'left', 'interval', 'interval', 'interval')
#   = c(2, -1, -1, 2, 2, 2)
# y = c(20, 10, 10, 20, 40, 10)
# y2 = c(40, 10, 10, 40, 80, 20)
format_hai_data(
	test,
	post_titer = "posttiter",
	pre_titer = NULL,
	log_scale = FALSE,
	increase = FALSE
)

# Expected same as above but on log scale
# y = c(2, 1, 1, 2, 3, 1)
# y2 = c(3, 1, 1, 3, 4, 2)
format_hai_data(
	test,
	post_titer = "posttiter",
	pre_titer = NULL,
	log_scale = FALSE,
	increase = FALSE,
	log_out = TRUE
)

# test for titer increase
# Expected fold changes for the test data are in these intervals:
# (0, 4), (1/4, Inf), (0, Inf), (1/2, 2), (1/4, 4), (1/4, 4)
# so the expected answers are:
# c = c('left', 'right',
# y = c(0, 1/4, 0,)
# y2 = c()
format_hai_data(
	test,
	post_titer = "posttiter",
	pre_titer = "pretiter",
	log_scale = FALSE,
	increase = TRUE,
	log_out = FALSE
)

format_hai_data(
	test,
	post_titer = "posttiter",
	pre_titer = "pretiter",
	log_scale = FALSE,
	increase = TRUE,
	log_out = TRUE
)
