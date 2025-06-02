# Simulate two data series from the same linear model
# but in the second one the x-values are closer together
dat <- tibble::tibble(
	x1 = c(0, 0.2, 0.4, 0.6, 0.8, 1.0),
	y1 = 6 - 5 * x1,
	x2 = c(0, 0.1, 0.2, 0.3, 0.4, 0.7),
	y2 = 6 - 5 * x2
)

# GMT is now distorted
mean(dat$y1)
mean(dat$y2)

# But lm estimates are the same
lm(y1 ~ x1, data = dat)
lm(y2 ~ x2, data = dat)

# Do it again with a bit of noise
set.seed(100)
dat2 <- dat |>
	dplyr::mutate(
		z1 = rnorm(dplyr::n()) + y1,
		z2 = rnorm(dplyr::n()) + y2
	)

# GMT is still distorted
mean(dat2$z1)
mean(dat2$z2)

# But lm estimates are very similar
lm(z1 ~ x1, data = dat2)
lm(z2 ~ x2, data = dat2)
