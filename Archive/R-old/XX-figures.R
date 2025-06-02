###
# Figure generation
# Zane Billings
# 2024-08-14
# This script generates all of the figures in the manuscript and supplement.
###

# Setup ####

## Packages and dependencies ####
box::use(
	ggplot2[...],
	patchwork[...],
	ggrepel,
	ggh4x,
	hgp,
	readr,
	here,
	dplyr,
	tidyr
)

## Set ggplot2 theme ####
ggplot2::theme_set(
	ggplot2::theme_bw() +
		ggplot2::theme(
			plot.background = ggplot2::element_rect(fill = "white", color = "white"),
			axis.text = ggplot2::element_text(size = 16, color = "black"),
			axis.title = ggplot2::element_text(size = 18),
			plot.subtitle = ggplot2::element_text(
				size = 16, hjust = 0, margin = ggplot2::margin(b = 2)
			),
			plot.title = ggplot2::element_text(
				size = 24, hjust = 0, margin = ggplot2::margin(b = 4)
			),
			plot.caption = ggplot2::element_text(size = 14),
			strip.text = ggplot2::element_text(
				size = 16, hjust = 0.5, margin = ggplot2::margin(b = 2, t = 2)
			),
			panel.spacing = ggplot2::unit(2, "lines"),
			legend.position = "bottom",
			legend.text = ggplot2::element_text(size = 16, color = "black"),
			legend.title = ggplot2::element_text(size = 18, color = "black")
		)
)

# Seed for jittering plots, etc.
S <- 6344430L

# Titer vs. distance full sample figure ########################################
# Make a figure with the following parameters:
# x axis: antigenic distance
# y axis: post-vaccination titer or titer increase
# facet by: antigenic distance method
# Have two x-axes, one for normalized and one on top for non-normalized distance
# Show one H1 strain and one H3 strain.

# First load the datasets. We'll load the saved results from the full data
# models fitting script.
full_data_models_res <-
	here::here("results", "data", "full-model-preds.Rds") |>
	readr::read_rds()

# Now we can unpack the cohort data and the model predictions.
# TODO decide if we should pool across studies or not
full_data_models_cohort_data <-
	full_data_models_res |>
	dplyr::select(-interp_preds, -data_preds) |>
	tidyr::unnest(data) |>
	dplyr::filter(season == "2016 - 2017") |>
	# Drop unused factor levels
	dplyr::mutate(
		dplyr::across(dplyr::where(is.factor), forcats::fct_drop)
	)

full_data_models_interp <-
	full_data_models_res |>
	dplyr::select(-data, -data_preds) |>
	dplyr::filter(season == "2016 - 2017") |>
	tidyr::unnest(interp_preds) |>
	# Calculate a credible interval instead of having multiple draws
	dplyr::group_by(season, strain_type, method, outcome_name, x)

full_data_models_means <-
	full_data_models_interp |>
	dplyr::summarize(
		ggdist::mean_hdci(.epred, .width = 0.97),
		.groups = "drop"
	)

# Load the dataset with the metrics
full_data_models_metrics <-
	here::here("results", "data", "full-model-metrics.Rds") |>
	readr::read_rds() |>
	# Doing some data cleaning and plot prep
	# First filter so that we have the same records as before
	dplyr::filter(
		season == "2016 - 2017",
		outcome_name == "log_posttiter",
		metric %in% c("p40", "intercept", "AUC")
	) |>
	# Set factor levels for some variables to make them pretty print
	dplyr::mutate(
		method_factor = factor(
			method,
			levels = c("temporal",
								 "pepitope",
								 "cartographic"),
			labels = c("Year difference",
								 "p-Epitope sequence distance",
								 "Antigenic cartography distance")
		),
		metric_factor = factor(
			metric,
			levels = c("intercept", "p40", "AUC"),
			labels = c("Magnitude", "Breadth", "Total Strength")
		)
	) |>
	dplyr::arrange(strain_type, metric_factor, method_factor) |>
	dplyr::mutate(
		lab = paste0(
			metric_factor, ": ", sprintf("%.2f", y),
			" (", sprintf("%.2f", ymin), ", ", sprintf("%.2f", ymax), ")"
		)
	) |>
	dplyr::group_by(strain_type, method_factor, method) |>
	dplyr::summarise(
		lab = paste(lab, collapse = "\n"),
		x = 0.65,
		y = 9,
		.groups = "drop"
	)

# Make a dataset with labels to put on the plot
# We want to label a few of the strains to show how using the temporal
# distance changes things for H1. So we picked a few representative strains.
full_data_models_strain_labels <-
	full_data_models_cohort_data |>
	dplyr::filter(
		strain_name %in% c("CA/09", "NC/99", "SC/18", "HK/14", "Shan/93", "TX/77"),
		vaccine_strain %in% c("CA/09", "HK/14")
	) |>
	dplyr::select(
		strain_type,
		strain_name,
		norm_d,
		d,
		method
	) |>
	dplyr::distinct() |>
	dplyr::mutate(
		y = -0.25,
		m = factor(
			method,
			levels = c("temporal", "pepitope", "cartographic"),
			labels = c(
				"Year difference",
				"p-Epitope sequence distance",
				"Antigenic cartography distance"
			)
		)
	)

# Make the CA/09 2014 ggplot
full_data_models_cohort_data |>
	dplyr::filter(
		outcome_name == "log_posttiter"
	) |>
	ggplot2::ggplot() +
	ggplot2::aes(
		x = norm_d,
		y = outcome_value
	) +
	ggplot2::geom_point(
		alpha = 0.1,
		size = 0.5,
		position = ggplot2::position_jitter(width = 0.01, height = 0.1, seed = 370)
	) +
	ggplot2::facet_wrap(ggplot2::vars(method))

# Making the figure 1 is kind of complicated because ggplot rebels against me
# for nearly every part of it. But nonetheless, through complicated feats of
# engineering, I have managed to construct it like the pyramids of old.
# This is all wrapped inside of a function to avoid a bunch of nasty stuff
# being in the global environment (it is better for everyone if all the
# intermediates get garbage collected).
make_dist_panel_figure <- function(subtype) {
	# We begin our great feat of engineering by constructing a simple empty list.
	p1 <- list()

	# We want to have one facet for each of the antigenic distance methods that
	# we used, and we want each facet to have two x-axes, one with the original
	# distance measurement (top) and one with the normalized distance measurement
	# (bottom). In ggplot2, secondary axes must be the result of a one-to-one
	# transformation of the primary axis, and cannot depend on the data. The
	# ggh4x::help_secondary() function can construct this transformation for us,
	# but we have to do it separately for each facet and then glue them together
	# at the end. So we will loop through the three distance measurements and
	# create the first three plots.
	for (i in c("temporal", "pepitope", "cartographic")) {

		# First we filter the data to get only the data for the current method.
		# We also need separate plots for H1N1 and H3N2, so we filter on that as well.
		dat_cohort <- full_data_models_cohort_data |>
			dplyr::filter(
				method == i, strain_type == subtype
			)

		# Now we can get the interpolated model data. We do the same filtering
		# as we did on the cohort data.
		dat_interp <-
			full_data_models_interp |>
			dplyr::filter(
				method == i, strain_type == subtype
			)

		dat_means <-
			full_data_models_means |>
			dplyr::filter(
				method == i, strain_type == subtype
			)

		# Filter the labels data the same way
		dat_labs <-
			full_data_models_strain_labels |>
			dplyr::filter(
				method == i, strain_type == subtype
			)

		# And finally do the exact same thing for the breadth metrics.
		dat_stats <-
			full_data_models_metrics |>
			dplyr::filter(
				method == i, strain_type == subtype
			)

		# Next we need to construct the secondary axis. This is a bit contrived,
		# but we know that the facets will go in order: "temporal", then "peitope",
		# then "cartographic". We only want one overall axis label, so we give the
		# pepitope secondary axis a name and leave the other two blank.
		if(i == "pepitope") {
			secondary_axis_scale <-
				ggh4x::help_secondary(
					data = dat_cohort,
					primary = norm_d,
					secondary = d,
					name = "Original distance measure"
				)
		} else {
			secondary_axis_scale <-
				ggh4x::help_secondary(
					data = dat_cohort,
					primary = norm_d,
					secondary = d,
					name = ""
				)
		}

		# Next we make the H1N1 ggplot.
		p1s <-
			dat_cohort |>
			ggplot2::ggplot() +
			ggplot2::aes(x = norm_d, y = outcome_value) +
			# Add a line to mark the breadth threshold
			ggplot2::geom_hline(
				alpha = 1,
				linewidth = 1,
				linetype = "dashed",
				yintercept = 3,
				color = "#bbbbbb"
			) +
			# Add data points for the actual observed values
			ggplot2::geom_point(
				alpha = 0.15,
				size = 0.5,
				position = ggplot2::position_jitter(
					width = 0.01, height = 0.1, seed = S
				)
			) +
			# Add a ribbon and a likne for the model prediction with CrI
			ggplot2::geom_line(
				data = dat_interp,
				mapping = ggplot2::aes(x = x, y = .epred, group = .draw),
				alpha = 0.25,
				color = "#cbc9e2"
			) +
			# ggplot2::geom_ribbon(
			# 	data = dat_interp,
			# 	mapping = ggplot2::aes(
			# 		x = x, y = y, ymin = ymin, ymax = ymax
			# 	),
			# 	alpha = 0.75,
			# 	fill = "lightblue"
			# ) +
			ggplot2::geom_line(
				data = dat_means,
				mapping = ggplot2::aes(
					x = x, y = y
				),
				color = "#6a51a3",
				linewidth = 1,
				alpha = 0.85,
				lineend = "round"
			) +
			# Now add text labels for the selected strains to show how they move
			# around.
			ggrepel::geom_label_repel(
				data = dat_labs,
				ggplot2::aes(
					x = norm_d, y = y,
					label = strain_name
				),
				inherit.aes = FALSE,
				min.segment.length = ggplot2::unit(0, 'lines')
			) +
			# Now add a box showing the metric calculations
			ggplot2::geom_label(
				data = dat_stats,
				ggplot2::aes(
					x = x, y = y,
					label = lab
				),
				inherit.aes = FALSE,
				alpha = 0.85
			) +
			# We can use facet_wrap() to add strip text to the top of the plot
			# even though there's only one facet. This is the easiest way to
			# get text that is spaced nicely.
			ggplot2::facet_wrap(ggplot2::vars(method_factor)) +
			# Next we customize the x and y axes and ensure they are all the same.
			ggplot2::scale_x_continuous(
				limits = c(-0.05, 1.05),
				breaks = seq(0, 1, 0.25),
				expand = ggplot2::expansion(0, 0),
				sec.axis = secondary_axis_scale
			) +
			ggplot2::scale_y_continuous(
				limits = c(-0.5, 10.5),
				breaks = c(0, 5, 10),
				expand = ggplot2::expansion(0, 0)
			) +
			# Specify the axis titles
			ggplot2::labs(
				x = "",
				y = NULL,
				color = "Outcome"
			) +
			# ggplot2::guides(
			# 	color = ggplot2::guide_legend(override.aes = list(alpha = 1, size = 3))
			# ) +
			# Theme customizations in addition to the baseline in the hgp theme.
			# In particular we need to edit the margins a bit to prevent ugly overlaps.
			# We also need to draw outlines or the axes are too confusing.
			ggplot2::theme(
				plot.margin = ggplot2::margin(t = 6, r = 8, b = 6, l = 8)
			)

		# Some more annoying plot details -- if we're on p-epitope, which is in the
		# midle, add the bottom x-axis legend.
		if (i == "pepitope") {
			p1s <- p1s + ggplot2::labs(
				x = "Normalized distance measure"
			)
		}
		# If we're on temporal distance, which is the leftmost plot, add a y-axis
		# legend.
		if (i == "temporal") {
			p1s <- p1s + ggplot2::labs(
				y = "titer (log2 HAI)"
			)
		}

		p1 <- c(p1, list(p1s))
	}

	# Create a title for the plot based on the subtype we are looking at.
	title <- paste0(
		"Vaccine: ",
		# First get the name of the vaccine strain we're looking at right now.
		dat_cohort$vaccine_strain |>
			as.character() |>
			unique() |>
			hgp::replace_strain_names(from = "short", to = "full") |>
			as.character(),
		"; 2014 - 2015 season",
		" (n = ",
		# And get the number of individuals observed
		dplyr::n_distinct(dat_cohort$subject_id),
		")"
	)

	full_p1 <-
		purrr::reduce(p1, `+`) +
		patchwork::plot_layout(guides = "collect", ncol = 3)  +
		patchwork::plot_annotation(
			title = title
		) &
		ggplot2::theme(legend.position = "bottom")

	return(full_p1)

}

p1 <- make_dist_panel_figure("H1N1")
p2 <- make_dist_panel_figure("H3N2")

dist_figure <- cowplot::plot_grid(
	p1, p2,
	rel_heights = c(1, 1),
	ncol = 1
)

ggplot2::ggsave(
	here::here("results/figures/posttiter-vs-distance.png"),
	plot = dist_figure,
	width = 12,
	height = 10
)
ggplot2::ggsave(
	here::here("results/figures/posttiter-vs-distance.tiff"),
	plot = dist_figure,
	width = 12,
	height = 10
)

################################################################################
# Panel subsampling metrics figure
################################################################################
# First load the data
subsample_metrics <-
	here::here("results", "data", "subsample-metrics.Rds") |>
	readr::read_rds()
