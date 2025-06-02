###
# Plots for CEIRR 2023
# Zane
# 2023-08-23
###

library(ggplot2)
ggplot2::theme_set(
	ggplot2::theme_bw() +
		ggplot2::theme(
			plot.background = ggplot2::element_rect(fill = "white", color = "white"),
			axis.text = ggplot2::element_text(size = 12, color = "black"),
			axis.title = ggplot2::element_text(size = 22),
			plot.subtitle = ggplot2::element_text(
				size = 16, hjust = 0, margin = ggplot2::margin(b = 2)
			),
			plot.title = ggplot2::element_text(
				size = 19, hjust = 0, margin = ggplot2::margin(b = 2, l = 2)
			),
			plot.caption = ggplot2::element_text(size = 14),
			strip.text = ggplot2::element_text(
				size = 16, hjust = 0.5, margin = ggplot2::margin(b = 2, t = 2)
			),
			panel.spacing = ggplot2::unit(2, "lines"),
			legend.position = "bottom",
			legend.text = ggplot2::element_text(size = 16, color = "black"),
			legend.title = ggplot2::element_text(size = 18, color = "black"),
			plot.margin = ggplot2::margin(t = 6, r = 6, b = 6, l = 6)
		)
)

# The regression plot with the distribution on it ####
set.seed(100)
ds = c(0, 0.17, 0.22, 0.30, 0.36, 0.41, 0.45, 0.52, 0.61, 0.76, 0.83)
res <- matrix(nrow = 100, ncol = length(ds))
for (i in 1:nrow(res)) {
	res[i, ] <- rnorm(length(ds), 5 - 5 * ds, 1)
}

res[, 1:10] |>
	`colnames<-`(ds[1:10]) |>
	tibble::as_tibble() |>
	tidyr::pivot_longer(cols = dplyr::everything()) |>
	ggplot() +
	aes(x = name, y = value) +
	geom_point(size = 3, alpha = 0.25, color = "blue",
						 shape = 21, fill = "white", stroke = 2) +
	theme_void() +
	theme(plot.background = ggplot2::element_rect(fill = NA, color = NA))


ggsave(here::here("Products/CEIRR2023-Presentation/Figures/Rplot001.png"))

# Example histograms ####

png(here::here("Products/CEIRR2023-Presentation/Figures/hist-continuous.png"))
hist(res[, 11], breaks = seq(-3, 4, 0.25), main = NULL, xlab = NULL, ylab = NULL, cex.axis = 2)
dev.off()

png(here::here("Products/CEIRR2023-Presentation/Figures/hist-floored.png"))
barplot(
	table(floor(res[, 11])),
	main = NULL, xlab = NULL, ylab = NULL, cex.axis = 2, cex.names = 2
)
dev.off()

x <- floor(res[, 11])
z <- ifelse(x <= 0, 0, x)

png(here::here("Products/CEIRR2023-Presentation/Figures/hist-cens.png"))
barplot(
	table(z),
	main = NULL, xlab = NULL, ylab = NULL, cex.axis = 2, cex.names = 2
)
dev.off()

# Metric plots for ted's data ###########################

clean_data <-
	readr::read_rds(
		here::here("data", "processed", "distance_data.rds")
	)

dat <-
	clean_data |>
	dplyr::mutate(
		vaccine_type = stringr::str_remove(vaccine_type, "_vaccine_fullname$") |>
			stringr::str_to_upper() |>
			factor()
	) |>
	dplyr::filter(
		vaccine_type == strain_type,
		# Get only the Ag distance methods that Amanda used.
		method %in% c("cart_2d_post", "p_epi", "year"),
		vaccine_fullname == "H1N1-California-2009",
		dose == "SD"
	) |>
	dplyr::mutate(
		dplyr::across(tidyselect:::where(is.factor), forcats::fct_drop)
	) |>
	# Remove the columns I know we do not need here
	dplyr::select(
		-c(pretiter, postiter, vac_short, strain_short)
	)

# Pivot the data so we can fit pre/post models at the same time
dat_models <-
	dat |>
	tidyr::pivot_longer(
		cols = c(postvactiter, titerincrease),
		names_to = "outcome",
		values_to = "y"
	) |>
	# Turn the categorical variables into integer indexes, this is required
	# for index coding in stan
	dplyr::mutate(
		id = uniq_id |>
			factor() |>
			forcats::fct_inorder() |>
			as.integer()
	) |>
	# Normalize the antigenic distance measurements within vaccine group
	dplyr::group_by(vaccine_fullname, method) |>
	dplyr::mutate(
		norm_dist = distance / max(distance)
	) |>
	dplyr::ungroup()

# We want to fit a model for every unique combination of:
# vaccine strain; outcome; distance method.
# (Models will be fitted ACROSS seasons, not per-season.)
# We'll do this by nesting the data into strata to get a per-stratum
# data frame and then mapping the model fitting function over the
# nested data frames.
dat_nested <-
	dat_models |>
	# Select only the data we need RIGHT NOW to prevent Stan from throwing a fit,
	# it will often get angry over unused data that is in the wrong format.
	dplyr::select(
		# Variables that should go into the Stan model. These DO need to be
		# processed to be NUMERIC.
		id = uniq_id, y, x = norm_dist, p = prevactiter,
		# Stratum/nesting variables -- these do NOT need to be processed.
		vaccine_fullname, method, outcome,
		# Strain name variable to filter on later. doesn't needed to be
		# processed as it will be removed in a later step
		s = strains_fullname
	) |>
	# Cleaning up the factor variabels that have to become indices
	# it has to be done in groups to ensure there are no missing levels!
	dplyr::group_by(dplyr::across(!c(id, x, y, s))) |>
	dplyr::mutate(
		id = id |> factor() |> forcats::fct_inorder() |> as.integer()
	) |>
	dplyr::ungroup() |>
	# This line creates the stratified data frames. Analogous to base R split()
	# but organized better.
	tidyr::nest(dat = c(id, x, y, s, p)) |>
	# TODO EVENTUALLY DO THIS FOR ALL OUTCOMES
	dplyr::filter(outcome == "titerincrease") |>
	dplyr::select(-outcome)

{ # Group this bit so that the seed set always runs with the sampling bit
	set.seed(370)
	N <- 1000
	k <- 10
	strains <- unique(dat_models$strains_fullname)
	# Generate all of the combinations of 10 strains
	# This produces an array of size k by (N choose k)
	all_subsamples <- combn(strains, 10)
	# Now randomly sample column indices to ensure we don't sample the
	# same sub-panel twice.
	# array, but this makes it easier to join with the dataframe correctly.
	subsamples <-
		purrr::map(
			sample.int(ncol(all_subsamples), N),
			\(i) all_subsamples[, i]
		)

}

dat_stan <-
	# Do the crossing part as described
	tidyr::expand_grid(dat_nested, subsamples) |>
	# Filter each row and remove strain variable
	dplyr::mutate(
		dat = purrr::map2(
			dat, subsamples,
			\(d, l) d |>
				dplyr::filter(s %in% l) |>
				dplyr::select(-s)
		),
		# Convert to list for stan
		dat = purrr::map(
			dat,
			\(x) x |>
				na.omit() |>
				as.list()
		),
		# Add the variable N that we need to pass
		dat = purrr::map(
			dat,
			\(d) c(d, N = length(d$id))
		)
	)

h <- here::here("Results", "_Out", "Panel-Subsampling")

# It's easier to clean the data if do this bit separately
dat_preds <-
	dat_stan |>
	dplyr::mutate(
		# First get the model predictions
		preds = purrr::map(
			list.files(paste0(h, "/Preds"), full.names = TRUE),
			\(x) readr::read_rds(x)
		)
	)

dat_parms <-
	dat_preds |>
	dplyr::mutate(
		# Now get the slope and intercept
		parameters = 	purrr::map(
			list.files(paste0(h, "/Summary"), full.names = TRUE),
			\(f) {
				p <- readr::read_rds(f)
				parms <-
					p |>
					dplyr::filter(parameter %in% c("a", "b")) |>
					dplyr::select(parameter, est = mean) |>
					tibble::tibble()
			}
		),
		# Compute the mean outcome for comparison to AUCS
		# ymean = purrr::map(
		# 	dat,
		# 	\(d) ggplot2::mean_cl_normal(d$y, conf.int = .89) |>
		# 		rlang::set_names(c("est", "lwr", "upr")) |>
		# 		tibble::tibble()
		# )
		ymean = purrr::map_dbl(
			dat,
			\(d) mean(d$y)
		),
		p_sc = purrr::map_dbl(
			dat,
			\(d) mean(
				( (d$p == 0) & (d$y >= 3) ) |
					( (d$p >= 1) & (d$y >= 2) )
			)
		),
		p_sp = purrr::map_dbl(
			dat,
			\(d) mean((d$y - d$p) >= 3)
		),
		p40 = purrr::map_dbl(
			preds,
			\(d) mean(d$est >= 3)
		)
	)

dat_stats <-
	dat_parms |>
	tidyr::unnest(parameters) |>
	tidyr::pivot_wider(
		names_from = parameter,
		values_from = est
	) |>
	dplyr::left_join(
		dat_auc,
		by = dplyr::join_by(vaccine_fullname, method, dat, subsamples, preds)
	) |>
	# tidyr::unnest(ymean, names_sep = "_") |>
	# tidyr::pivot_longer(cols = !(vaccine_fullname:preds)) |>
	# tidyr::separate(name, into = c("name", "stat")) |>
	# tidyr::pivot_wider(names_from = stat, values_from = value) |>
	# dplyr::bind_rows(dat_auc) |>
	tidyr::pivot_longer(
		cols = c(ymean, a, b, p_sc, p_sp, p40, dplyr::starts_with("AUC"))
	) |>
	dplyr::left_join(
		# Get the short names back
		readr::read_rds(here::here("data", "processed", "virus_info.rds")),
		by = c("vaccine_fullname" = "analysis_name")
	) |>
	dplyr::mutate(
		subtype = substr(vaccine_fullname, 1, 4),
		method = factor(
			method,
			levels = c("year", "p_epi", "cart_2d_post"),
			labels = c("Year difference", "p-Epitope", "Antigenic cartography")
		),
		name = factor(
			name,
			levels = c("ymean", "p_sc", "p_sp",
								 "AUC_unweighted", "AUC_linear", "AUC_step",
								 "a", "b", "p40"),
			labels = c("Mean TI",
								 "Prop. seroconverted",
								 "Prop. seroprotection",
								 "AUC",
								 "AUC (linear)",
								 "AUC (2 AU)",
								 "Intercept",
								 "Slope",
								 "p40")
		)
	)|>
	dplyr::rename(stat = name)

quantile_df <- function(x, coverage = c(0.5, 0.75, 1)) {
	a <- (1 - coverage)/2
	tibble::tibble(
		lwr = quantile(x, a, na.rm = TRUE),
		upr = quantile(x, 1 - a, na.rm = TRUE),
		quant = coverage
	)
}

dat_summary <-
	dat_stats |>
	dplyr::group_by(method, stat) |>
	dplyr::reframe(
		mean = mean(value),
		med = median(value),
		quantile_df(value, c(0.75, 0.95, 1))
	)

library(ggplot2)
# Manual
dat_summary |>
	dplyr::filter(
		stat %in% c("Intercept", "AUC", "p40", "Prop. seroconverted", "Mean TI")
	) |>
	dplyr::arrange(-quant) |>
	ggplot2::ggplot() +
	ggplot2::aes(x = mean, y = stat, xmin = lwr, xmax = upr, color = factor(quant)) +
	ggplot2::geom_linerange(linewidth = 3, alpha = 1) +
	ggplot2::geom_point(color = "red", shape = "|", stroke = 5, size = 5) +
	ggplot2::facet_wrap(ggplot2::vars(method), ncol = 3) +
	ggplot2::labs(
		x = NULL,
		y = NULL,
		color = "ETCI width",
		title = "Vaccine: H1N1-California-2009"
	) +
	ggplot2::scale_color_grey()

ggsave(
	here::here("Products/CEIRR2023-Presentation/Figures/metrics-plot-dat.png"),
	width = 16,
	height = 9
)

# Antigenic distance measures and lm/GAM plot ####
p1data <- readr::read_rds(here::here("Andreas-Poster-Plots/p1data.Rds"))
dat_p1 <- p1data[[1]] |>
	dplyr::filter(
		dose == "SD",
	)
plt_test <- p1data[[2]] |>
	dplyr::mutate(
		o = factor(outcome,
							 levels = c("prevactiter", "postvactiter", "titerincrease"),
							 labels = c("Pre-vaccination titer",
							 					 "Post-vaccination titer",
							 					 "Titer increase")),
	)
dat_table <- p1data[[3]]

dat_labs <-
	dat_p1 |>
	dplyr::filter(
		strain_short %in% c("CA/09", "NC/99", "SC/18", "HK/14", "Shan/93", "TX/77"),
		vac_short %in% c("CA/09", "HK/14")
	) |>
	dplyr::select(
		strain_type,
		strain_short,
		norm_dist,
		distance,
		method
	) |>
	dplyr::distinct() |>
	dplyr::mutate(
		y = -0.25,
		m = factor(
			method,
			levels = c("year", "p_epi", "cart_2d_post"),
			labels = c(
				"Year difference",
				"p-Epitope sequence distance",
				"Antigenic cartography distance"
			)
		)
	)

p1 <- list()

for (i in c("year", "p_epi", "cart_2d_post")) {

	dat <-dat_p1 |>
		dplyr::filter(
			method == i, vaccine_fullname == "H1N1-California-2009",
			outcome == "postvactiter"
		)

	dat2 <-
		plt_test |>
		dplyr::filter(
			method == i, vaccine_fullname == "H1N1-California-2009",
			outcome == "postvactiter"
		)

	if( i == "p_epi") {
		sec <-
			ggh4x::help_secondary(
				data =dat ,
				primary = norm_dist,
				secondary = distance,
				name = "Original distance measure"
			)
	} else {
		sec <-
			ggh4x::help_secondary(
				data =dat ,
				primary = norm_dist,
				secondary = distance,
				name = NULL
			)
	}

	p1s <-
		dat |>
		dplyr::mutate(
			m = factor(
				method,
				levels = c("year", "p_epi", "cart_2d_post"),
				labels = c(
					"Year difference",
					"p-Epitope sequence distance",
					"Antigenic cartography distance"
				)
			)
		) |>
		ggplot2::ggplot() +
		ggplot2::aes(x = norm_dist, y = y) +
		ggplot2::geom_point(
			alpha = 0.1,
			size = 0.5,
			position = ggplot2::position_jitter(width = 0.005, height = 0.35, seed = 370)
		) +
		# ggplot2::geom_ribbon(
		# 	data = dat2,
		# 	mapping = ggplot2::aes(
		# 		x = distance, y = est, ymin = est - 0.05, ymax = est + 0.05, group = o
		# 	),
		# 	alpha = 0.25,
		# 	color = "white"
		# ) +
		# ggplot2::geom_line(
		# 	data = dat2,
		# 	mapping = ggplot2::aes(
	# 		x = distance, y = est
	# 	),
	# 	linewidth = 1.5,
	# 	alpha = 0.85,
	# 	lineend = "round",
	# 	color = "dodgerblue2"
	# ) +
	ggplot2::geom_smooth(
		method = "lm",
		formula = y~x,
		linewidth = 1.5,
		alpha = 0.85,
		lineend = "round",
		color = "#56B4E9",
		se = FALSE
	) +
		ggplot2::geom_smooth(
			method = "gam",
			formula = y ~ s(x, bs = "cs", k = 5),
			linewidth = 1.5,
			alpha = 0.85,
			lineend = "round",
			color = "#E69F00",
			linetype = "dashed",
			se = FALSE
		) +
		ggrepel::geom_label_repel(
			data = dat_labs |>
				dplyr::filter(strain_type == 'H1N1', method == i),
			ggplot2::aes(
				x = norm_dist, y = y,
				label = strain_short
			),
			inherit.aes = FALSE,
			min.segment.length = ggplot2::unit(0, 'lines')
		) +
		ggplot2::facet_wrap(ggplot2::vars(m)) +
		ggplot2::scale_x_continuous(
			limits = c(-0.05, 1.05),
			breaks = seq(0, 1, 0.25),
			expand = ggplot2::expansion(0, 0),
			sec.axis = sec
		) +
		ggplot2::scale_y_continuous(
			limits = c(-1, 10.5),
			breaks = c(-5, 0, 5, 10),
			expand = ggplot2::expansion(0, 0)
		) +
		ggplot2::labs(
			x = "",
			y = NULL,
			color = "Outcome",
			#title = "H1N1-California-2009 (n = 773)"
		) +
		ggplot2::guides(
			color = ggplot2::guide_legend(override.aes = list(alpha = 1, size = 3))
		) +
		ggplot2::theme(
			plot.margin = ggplot2::margin(t = 6, r = 8, b = 6, l = 8)
		)

	if (i == "p_epi") {
		p1s <- p1s + ggplot2::labs(
			x = "Normalized distance measure"
		)
	}

	if (i == "year") {
		p1s <- p1s + ggplot2::labs(
			y = "titer (log2 HAI)"
		)
	}

	p1 <- c(p1, list(p1s))
}

full_p1 <-
	purrr::reduce(p1, `+`) +
	patchwork::plot_layout(guides = "collect", ncol = 3)  +
	patchwork::plot_annotation(
		title = "H1N1-California-2009 (n = 773)"
	)&
	ggplot2::theme(legend.position = "bottom"); full_p1

ggsave(
	here::here("Products/CEIRR2023-Presentation/Figures/h1n1-dist-lm.png"),
	plot = full_p1,
	width = 16,
	height = 9
)

p2 <- list()

for (i in c("year", "p_epi", "cart_2d_post")) {

	dat <- dat_p1 |>
		dplyr::filter(
			method == i, vaccine_fullname == "H3N2-Hong Kong-2014",
			outcome == "postvactiter"
		)

	dat2 <-
		plt_test |>
		dplyr::filter(
			method == i, vaccine_fullname == "H3N2-Hong Kong-2014",
			outcome == "postvactiter"
		)

	if( i == "p_epi") {
		sec <-
			ggh4x::help_secondary(
				data =dat ,
				primary = norm_dist,
				secondary = distance,
				name = "Original distance measure"
			)
	} else {
		sec <-
			ggh4x::help_secondary(
				data =dat ,
				primary = norm_dist,
				secondary = distance,
				name = NULL
			)
	}

	p2s <-
		dat |>
		dplyr::mutate(
			m = factor(
				method,
				levels = c("year", "p_epi", "cart_2d_post"),
				labels = c(
					"Year difference",
					"p-Epitope sequence distance",
					"Antigenic cartography distance"
				)
			)
		) |>
		ggplot2::ggplot() +
		ggplot2::aes(x = norm_dist, y = y) +
		ggplot2::geom_point(
			alpha = 0.1,
			size = 0.5,
			position = ggplot2::position_jitter(width = 0.005, height = 0.35, seed = 370)
		) +
		ggplot2::geom_smooth(
			method = "lm",
			formula = y~x,
			linewidth = 1.5,
			alpha = 0.85,
			lineend = "round",
			color = "#56B4E9",
			se = FALSE
		) +
		ggplot2::geom_smooth(
			method = "gam",
			formula = y ~ s(x, bs = "cs", k = 5),
			linewidth = 1.5,
			alpha = 0.85,
			lineend = "round",
			color = "#E69F00",
			linetype = "dashed",
			se = FALSE
		) +
		ggrepel::geom_label_repel(
			data = dat_labs |>
				dplyr::filter(strain_type == 'H3N2', method == i),
			ggplot2::aes(
				x = norm_dist, y = y,
				label = strain_short
			),
			inherit.aes = FALSE,
			min.segment.length = ggplot2::unit(0, 'lines')
		) +
		ggplot2::facet_wrap(ggplot2::vars(m)) +
		ggplot2::scale_x_continuous(
			limits = c(-0.05, 1.05),
			breaks = seq(0, 1, 0.25),
			expand = ggplot2::expansion(0, 0),
			sec.axis = sec
		) +
		ggplot2::scale_y_continuous(
			limits = c(-1, 10.5),
			breaks = c(-5, 0, 5, 10),
			expand = ggplot2::expansion(0, 0)
		) +
		ggplot2::labs(
			x = "",
			y = NULL,
			color = "Outcome",
			#title = "H1N1-California-2009 (n = 773)"
		) +
		ggplot2::guides(
			color = ggplot2::guide_legend(override.aes = list(alpha = 1, size = 3))
		) +
		ggplot2::theme(
			plot.margin = ggplot2::margin(t = 6, r = 8, b = 6, l = 8)
		)

	if (i == "p_epi") {
		p2s <- p2s + ggplot2::labs(
			x = "Normalized distance measure"
		)
	}

	if (i == "year") {
		p2s <- p2s + ggplot2::labs(
			y = "titer (log2 HAI)"
		)
	}

	p2 <- c(p2, list(p2s))
}

full_p2 <-
	purrr::reduce(p2, `+`) +
	patchwork::plot_layout(guides = "collect", ncol = 3) +
	patchwork::plot_annotation(title = "H3N2-Hong Kong-2014 (n = 583)") &
	ggplot2::theme(legend.position = "bottom"); full_p2

ggsave(
	here::here("Products/CEIRR2023-Presentation/Figures/h3n2-dist-lm.png"),
	plot = full_p2,
	width = 16,
	height = 9
)

leg <- cowplot::get_legend(
	full_p1 + ggplot2::theme(legend.box.margin = ggplot2::margin(t = 12))
)

full_p1 <- full_p1 & ggplot2::theme(legend.position = "none")

overall_p1 <- cowplot::plot_grid(
	full_p1, full_p2, leg,
	rel_heights = c(1, 1, .1),
	ncol = 1
)

# END OF FILE ####
