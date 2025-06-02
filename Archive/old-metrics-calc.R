# Calculate metrics ####
# Next we need to calculate the breadth metrics on the full data models.
# First, calculate the old metrics: mean TI and seroconversion rate.

# And the new metrics
metrics_calc <-
	full_model_preds |>
	dplyr::select(-data) |>
	dplyr::mutate(
		interp_preds = purrr::map(
			interp_preds,
			\(d) d |> dplyr::ungroup() |> dplyr::group_by(.draw)
		),
		data_preds = purrr::map(
			data_preds,
			\(d) d |> dplyr::ungroup() |> dplyr::group_by(.draw)
		),
		#AUC Calculation -- no weighting
		AUC = purrr::map(
			interp_preds,
			# Calculate the pointwise AUC using the trapezoidal method
			\(d) d |>
				dplyr::summarize(
					out = pracma::trapz(x, .epred),
					.groups = "drop"
				) |>
				dplyr::summarize(ggdist::mean_hdci(out))
		),
		# Intercept of line
		intercept = purrr::map(
			interp_preds,
			\(d) d |>
				dplyr::group_split() |>
				purrr::map_dbl(\(d) d$.epred[d$x == 0]) |>
				ggdist::mean_hdci()
		),
		# percent of values above threshold (3)
		p40 = purrr::map(
			interp_preds,
			\(d) d |>
				dplyr::group_split() |>
				purrr::map_dbl(\(d) mean(d$.epred >= 3)) |>
				ggdist::mean_hdci()
		),
		# TODO FINISH THIS AFTER WE FIX TITER INCREASE FOR OUTCOME
		# seroconversion rate
		# scr = purrr::map_dbl(
		# 	data_preds,
		# 	\(d) d |>
		# 		dplyr::mutate(
		# 			seroconversion =
		# 		)
		# 		dplyr::group_split() |>
		# 		purrr::map_dbl(\(d) mean(d$.epred >= 3)) |>
		# 		ggdist::mean_hdci()
		# ),
		# # geometric mean titer
		# gmt = purrr::map_dbl(
		# 	lm_preds,
		# 	\(d) geo_mean(d$y_nat)
		# ),
		# # homologous only GMT
		# gmt_hom = purrr::map_dbl(
		# 	lm_preds,
		# 	\(d) d |> dplyr::filter(d == 0) |> dplyr::pull(y_nat) |> geo_mean()
		# )
	) |>
	dplyr::select(-interp_preds, -data_preds)

# Now we have to collect all that crap together somehow since I couldn't
# think of a better way to do this.
metrics_tidy <-
	metrics_calc |>
	dplyr::rowwise() |>
	dplyr::mutate(
		metrics = dplyr::c_across(c(AUC, intercept, p40)) |>
			dplyr::bind_rows() |>
			dplyr::mutate(metric = c("AUC", "intercept", "p40")) |>
			list(),
		.keep = "unused"
	) |>
	dplyr::ungroup() |>
	tidyr::unnest(metrics)

write_csv_and_rds(
	metrics_tidy,
	here::here("results", "data", "full-model-metrics")
)
