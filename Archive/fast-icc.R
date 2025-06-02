# Frequentist Icc model with bootstrap ICCs
# not "fast" but way faster than bayesian
calculate_fast_icc_from_model <- function(icc_data) {
	icc_model <- lme4::lmer(
		y ~ 1 + (1 | subsample_id),
		data = icc_data[[1]] # Have to have the [[1]] because of targets?
		# probably because idk what I am doing....
	)
	icc_res <- performance::icc(icc_model, ci = 0.95)

	out <- setNames(
		icc_res$ICC_adjusted,
		c("est", "lwr", "upr")
	) |>
		as.list() |>
		tibble::as_tibble()

	return(out)
}
