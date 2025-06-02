library(reprex)
library(brms)

brms::brm(
	formula = bf(
		mpg ~ 1,
		decomp = 'QR'
	),
	data = mtcars
)
