# testing if slope gives an error when there's only one predictor value
set.seed(791737080L)
# create a data frame with one predictor value
test_data <-
	tibble::tibble(
		x = rep(1, 100)
	) |>
	dplyr::mutate(
		y = 3 + 10 * x + rnorm(dplyr::n(), 0, 2)
	)

test_lm <- lm(y ~ x, data = test_data)
summary(test_lm)

test_rstanarm <- rstanarm::stan_glm(y ~ x, data = test_data)
summary(test_rstanarm)

test_brm <- brms::brm(y ~ x, data = test_data)
summary(test_brm)

# so we can see that if there is only one x value the slope cannot be
# estimated at all, which is what I thought but thought I might need to convince
# Andreas.
