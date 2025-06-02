###
# Extracitn results from models
# Zane
# 2025-05-03
###

# Summary landscape figure ####
# function to process the full data model predictions
get_summary_landscape <- function(full_data_model_fits, model_metadata, h = 0.01) {
	# Probably would benefit from being branched over all model fits but it's not
	# slow enough to worry about
	aug_metadata <-
		model_metadata |>
		tibble::add_column(.fit = full_data_model_fits) |>
		dplyr::select(-formula, -family, -priors, -fitting_seed) |>
		dplyr::filter(model == "full", dataset == "full") |>
		# Get the epreds
		dplyr::mutate(
			.preds = purrr::map(
				.fit,
				\(model) tidybayes::epred_draws(
					model,
					newdata = tibble::tibble(d_norm = seq(0, 1, h)),
					ndraws = brms::ndraws(model),
					re_formula = NA
				)
			)
		)

	preds <-
		aug_metadata |>
		dplyr::select(-.fit, -dat) |>
		tidyr::unnest(.preds) |>
		dplyr::group_by(season, metric, dataset, model, censoring, d_norm) |>
		dplyr::summarise(tidybayes::mean_hdci(.epred), .groups = "drop") |>
		dplyr::mutate(
			dplyr::across(
				c(y, ymin, ymax),
				hgp::hai_to_natural_scale
			),
			censoring = factor(
				as.character(censoring),
				levels = c("no", "yes"),
				labels = c("LoD set to 5", "Censoring correction")
			)
		)

	return(preds)
}

make_summary_landscape_plot <- function(
		this_season, model_data, preds, strain_text
	) {
	this_data <- dplyr::filter(model_data, season == this_season)
	this_preds <- dplyr::filter(preds, season == this_season)
	this_strains <- dplyr::filter(strain_text, season == this_season)

	out <-
		ggplot2::ggplot(this_data) +
		ggplot2::geom_point(
			mapping = ggplot2::aes(x = d_norm, y = posttiter),
			size = 0.5,
			alpha = 0.15,
			position = ggplot2::position_jitter(0.01, 0.15, 370)
		) +
		ggplot2::geom_ribbon(
			mapping = ggplot2::aes(x = d_norm, ymin = ymin, ymax = ymax),
			data = this_preds,
			alpha = 0.5,
			fill = "gray",
			color = "black",
			linewidth = 0.5,
			linetype = 3
		) +
		ggplot2::geom_line(
			mapping = ggplot2::aes(x = d_norm, y = y),
			data = this_preds,
			linetype = 2,
			linewidth = 1.5
		) +
		ggplot2::geom_segment(
			mapping = ggplot2::aes(x = d_norm, xend = d_norm, y = 5120, yend = 2560),
			data = this_strains,
			color = "black",
			linewidth = 1
		) +
		ggplot2::geom_label(
			mapping = ggplot2::aes(x = d_norm, y = 5120, label = assay_strain),
			data = this_strains,
			size = 10,
			size.unit = "pt"
		) +
		ggplot2::scale_y_continuous(
			trans = "log2",
			breaks = hgp::hai_to_natural_scale(seq(0, 10, 2)),
			limits = hgp::hai_to_natural_scale(c(-0.15, 10.15))
		) +
		ggplot2::scale_x_continuous(
			limits = c(0, 1),
			expand = ggplot2::expansion(mult = 0.05, add = 0.025)
		) +
		ggplot2::facet_grid(censoring~metric) +
		ggplot2::labs(
			x = "Normalized antigenic distance",
			y = "Post-vaccination HAI titer",
			title = this_season
		) +
		hgp::theme_ms()

	return(out)
}

make_all_summary_landscape_plots <- function(
		preds, model_metadata,
		base_name
	) {
	model_data <- model_metadata |>
		dplyr::select(-formula, -family, -priors, -fitting_seed) |>
		dplyr::filter(model == "full", dataset == "full") |>
		tidyr::unnest(dat) |>
		dplyr::mutate(
			censoring = factor(
				as.character(censoring),
				levels = c("no", "yes"),
				labels = c("LoD set to 5", "Censoring correction")
			)
		)

	strains_text <- model_data |>
		dplyr::filter(assay_strain %in% c("USSR/77", "CA/09", "SC/18")) |>
		dplyr::distinct(season, assay_strain, d_norm, metric, censoring)

	all_seasons <- unique(model_data$season)

	plot_vec <-
		purrr::map(
			all_seasons,
			\(s) make_summary_landscape_plot(s, model_data, preds, strains_text)
		)

	plot_file_names <- here::here(
		base_name,
		paste0(
			"Summary-Landscape-",
			stringr::str_sub(all_seasons, 3, 4), "-",
			stringr::str_sub(all_seasons, 10, 11),
			".png"
		)
	)

	# Doesn't need to be dynamically branched because the overhead would be slower
	# than the time it takes to just do it
	purrr::walk2(
		plot_vec, plot_file_names,
		\(p, n) ggplot2::ggsave(
			filename = n,
			plot = p,
			width = 13,
			height = 10
		)
	)

	return(plot_file_names)
}

# Metrics table ####
make_model_metrics_table <- function(
		metrics_summary_df, example_only, file_name
) {
	filtered_df <-
		metrics_summary_df |>
		dplyr::filter(
			stat_group == "raw",
			subsample_id == "00"
		) |>
		dplyr::select(season, metric, stats_id, censoring, y, ymin, ymax)

	if(isTRUE(example_only)) {
		tbl_data <-
			filtered_df |>
			dplyr::filter(season == "2016 - 2017")
	} else if (isFALSE(example_only)) {
		tbl_data <- filtered_df
	} else {
		cli::cli_abort(c(
			"{.val example_only} should be {.val TRUE} or {.val FALSE}.",
			"x" = "Not {.val {example_only}}."
		))
	}

	data_long <-
		tbl_data |>
		dplyr::mutate(
			dplyr::across(
				c(y, ymin, ymax),
				\(x) formatC(x, digits = 2, width = 4, format = "f")
			),
			ci = paste0(y, " (", ymin, ", ", ymax, ")"),
			group = dplyr::case_match(
				stats_id,
				"GMT" ~ "old-tot",
				"GMT0" ~ "old-mag",
				"SCR" ~ "old-brd",
				"INT" ~ "new-mag",
				"PAT" ~ "new-brd",
				"AUC" ~ "new-tot"
			)
		) |>
		tidyr::separate(col = group, into = c("set", "att"), sep = "-") |>
		dplyr::select(-c(y, ymin, ymax)) |>
		dplyr::filter(stats_id != "SCR0") |>
		dplyr::mutate(
			set = factor(
				set,
				levels = c("old", "new"),
				labels = c("Current", "Novel")
			),
			att = factor(
				att,
				levels = c("mag", "brd", "tot"),
				labels = c("Magnitude", "Breadth", "Total Strength")
			),
			stats_id = factor(
				stats_id,
				levels = c("GMT0", "SCR", "GMT", "INT", "PAT", "AUC"),
				labels = c(
					"Homologous GMT",
					"Seroconversion rate",
					"Overall GMT",
					"Landscape intercept",
					"Prop. landscape above threshold",
					"Area under landscape curve"
				)
			),
			censoring = factor(
				censoring,
				levels = c("yes", "no"),
				labels = c("Censoring correction", "LoD set to 5")
			)
		) |>
		dplyr::arrange(season, metric, censoring, att, set)


	ft_data <-
		data_long |>
		dplyr::mutate(
			ultragroup = paste0(att, "_", set)
		) |>
		dplyr::select(-c(att, set, stats_id)) |>
		tidyr::pivot_wider(
			values_from = ci,
			names_from = ultragroup
		)

	# Separate stuff so we can make the table narrower for word doc--
	# make the current metrics in the ag dist column
	dat_all <-
		data_long |>
		dplyr::select(-stats_id) |>
		dplyr::mutate(
			newgroup = ifelse(set == "Novel", as.character(metric), "curr"),
			newgroupf = factor(
				newgroup,
				levels = c("curr", "Cartographic", "p-Epitope", "Temporal"),
				labels = c(
					"Current",
					"Novel (Cartographic)",
					"Novel (p-Epitope)",
					"Novel (Temporal)"
				)
			)
		) |>
		dplyr::select(season, metric = newgroupf, censoring, ci, att) |>
		dplyr::arrange(season, metric, censoring, att)

	# Have to do some manually filtering of the current metrics cause they
	# have inconsequential (0.01 or 0.02) differences due to MC error
	dat_current <-
		dat_all |>
		dplyr::filter(metric == "Current") |>
		dplyr::group_by(season, metric, censoring, att) |>
		dplyr::slice_head(n = 1) |>
		dplyr::ungroup()

	dat_current_distinct <-
		dplyr::bind_rows(
			dat_current,
			dplyr::filter(dat_all, metric != "Current")
		) |>
		dplyr::arrange(season, metric, censoring, att)

	ft_data <-
		dat_current_distinct |>
		tidyr::pivot_wider(
			names_from = c(censoring, att),
			values_from = ci
		) |>
		dplyr::rename(Season = season, "Metric Set" = metric)

	if(isTRUE(example_only)) {
		ft_data2 <- ft_data |>
			dplyr::select(-Season)
		vmerge_vals <- c(1)
		hline_cols <- c(2:7)
		vline_vals <- c(1, 4)
	} else {
		ft_data2 <- ft_data
		vmerge_vals <- c(1, 2)
		hline_cols <- c(3:8)
		vline_vals <- c(2, 5)
	}

	metrics_ft <-
		ft_data2 |>
		flextable::flextable() |>
		flextable::separate_header(split = "_") |>
		flextable::merge_v(j = vmerge_vals) |>
		flextable::valign(j = vmerge_vals, valign = "top") |>
		flextable::hline(i = 1, j = hline_cols, part = "header") |>
		flextable::vline(j = vline_vals[[1]], part = "body") |>
		flextable::vline(j = vline_vals[[2]], part = "all") |>
		flextable::fontsize(size = 8, part = "body") |>
		flextable::fontsize(size = 10, part = "header") |>
		flextable::fix_border_issues() |>
		fit_flextable_to_page()

	save_file_as_qs2(
		metrics_ft,
		file_name
	)

	return(file_name)
}

# ICC on subsample plots ####
create_icc_plots <- function(bound_metrics_df, path_base) {
	require(patchwork, quietly = TRUE)
	plot_df <-
		bound_metrics_df |>
		dplyr::filter(
			stats_id != "SCR0",
			censoring == "yes",
			subsample_id != "00"
		) |>
		dplyr::select(season, metric, stats_id, subsample_id, stat_norm) |>
		dplyr::mutate(
			set = dplyr::case_match(
				stats_id,
				"GMT" ~ "Current",
				"GMT0" ~ "Current",
				"SCR" ~ "Current",
				"INT" ~ "Novel",
				"PAT" ~ "Novel",
				"AUC" ~ "Novel"
			) |> factor(levels = c("Current", "Novel")),
			att = dplyr::case_match(
				stats_id,
				"GMT" ~ "Total Strength",
				"GMT0" ~ "Magnitude",
				"SCR" ~ "Breadth",
				"INT" ~ "Magnitude",
				"PAT" ~ "Breadth",
				"AUC" ~ "Total Strength"
			) |> factor(levels = c("Magnitude", "Breadth", "Total Strength"))
		) |>
		dplyr::select(-stats_id)

	# Separate stuff so we can make the table narrower for word doc--
	# make the current metrics in the ag dist column
	dat_all <-
		plot_df |>
		dplyr::mutate(
			newgroup = ifelse(set == "Novel", as.character(metric), "curr"),
			newgroupf = factor(
				newgroup,
				levels = c("curr", "Cartographic", "p-Epitope", "Temporal"),
				labels = c(
					"Current",
					"Novel (Cartographic)",
					"Novel (p-Epitope)",
					"Novel (Temporal)"
				)
			)
		) |>
		dplyr::select(
			season, metric = newgroupf, subsample_id, att, y = stat_norm
		) |>
		dplyr::arrange(season, metric, att)

	# Have to do some manually filtering of the current metrics cause they
	# have inconsequential (0.01 or 0.02) differences due to MC error
	# Just get the first 4000 draws since we arranged everything before
	dat_current <-
		dat_all |>
		dplyr::filter(metric == "Current") |>
		dplyr::group_by(season, metric, att, subsample_id) |>
		dplyr::slice_head(n = 4000) |>
		dplyr::ungroup()

	dat_current_distinct <-
		dplyr::bind_rows(
			dat_current,
			dplyr::filter(dat_all, metric != "Current")
		) |>
		dplyr::arrange(season, metric, att) |>
		dplyr::filter(metric != "Novel (Temporal)")

	set.seed(4589)
	final_dat_plot <-
		dat_current_distinct |>
		dplyr::group_by(season, metric, subsample_id, att) |>
		dplyr::slice_sample(n = 1000) |>
		dplyr::ungroup()

	final_dat_plot_nested <-
		final_dat_plot |>
		tidyr::nest(plt_dat = -c(season, att)) |>
		dplyr::mutate(
			overall_means = purrr::map(
				plt_dat,
				\(d) d |>
					dplyr::summarize(y = mean(y), .by = metric)
			),
			plt = purrr::map2(
				plt_dat, overall_means,
				\(d, o) d |>
					ggplot2::ggplot() +
					ggplot2::aes(x = subsample_id, y = y) +
					ggplot2::geom_hline(
						ggplot2::aes(yintercept = y),
						data = o,
						color = "red",
						linetype = "dashed",
						linewidth = 0.75
					) +
					ggplot2::geom_point(alpha = 0.01) +
					ggplot2::stat_summary(
						fun = mean,
						geom = "point",
						color = "red",
						shape = 4,
						size = 1.25,
						stroke = 1.25
					) +
					ggplot2::facet_grid(~metric) +
					hgp::theme_ms() +
					ggplot2::labs(x = NULL, y = "Minmax scaled metric value")
			),
			title = paste0(season, ": ", att),
			plt = purrr::map2(
				plt, att,
				\(p, t) p + ggplot2::ggtitle(t)
			)
		)

	plots_combined <-
		final_dat_plot_nested |>
		dplyr::summarize(
			plt_group = list(
				purrr::reduce(plt, `/`) +
					patchwork::plot_layout(axes = "collect")
			),
			.by = season
		) |>
		dplyr::mutate(
			plt_group = purrr::map2(
				plt_group, season,
				\(p, s) p +
				patchwork::plot_annotation(title = paste0("Cohort: ", s)) &
				ggplot2::theme(
					axis.text.x = ggplot2::element_text(size = 10, angle = 45),
					plot.title = ggplot2::element_text(size = 32)
				)
			)
		)

	file_names <- paste0(
		here::here(path_base, "icc-plot-"),
		stringr::str_sub(plots_combined$season, 3, 4), "-",
		stringr::str_sub(plots_combined$season, 10, 11),
		".png"
	)

	purrr::walk2(
		plots_combined$plt_group,
		file_names,
		\(p, fn) {
			ggplot2::ggsave(
				fn,
				plot = p,
				width = 13,
				height = 8.5 * 2
			)
		}
	)

	return(file_names)
}

# ICC table ####
clean_up_icc_results <- function(icc_results_df) {
	icc_res_cleaned <-
		icc_results_df |>
		dplyr::rename(
			est = y,
			lwr = ymin,
			upr = ymax
		) |>
		dplyr::select(-c(.width, .point, .interval)) |>
		dplyr::filter(stats_id != "SCR0", censoring != "no", season != "Overall") |>
		dplyr::mutate(
			set = dplyr::case_match(
				stats_id,
				"GMT" ~ "Current",
				"GMT0" ~ "Current",
				"SCR" ~ "Current",
				"INT" ~ "Novel",
				"PAT" ~ "Novel",
				"AUC" ~ "Novel"
			) |> factor(levels = c("Current", "Novel")),
			att = dplyr::case_match(
				stats_id,
				"GMT" ~ "Total Strength",
				"GMT0" ~ "Magnitude",
				"SCR" ~ "Breadth",
				"INT" ~ "Magnitude",
				"PAT" ~ "Breadth",
				"AUC" ~ "Total Strength"
			) |> factor(levels = c("Magnitude", "Breadth", "Total Strength"))
		) |>
		dplyr::select(season, metric, set, att, est, lwr, upr)

	# Separate stuff so we can make the table narrower for word doc--
	# make the current metrics in the ag dist column
	dat_all <-
		icc_res_cleaned |>
		dplyr::mutate(
			newgroup = ifelse(set == "Novel", as.character(metric), "curr"),
			newgroupf = factor(
				newgroup,
				levels = c("curr", "Cartographic", "p-Epitope", "Temporal"),
				labels = c(
					"Current",
					"Novel (Cartographic)",
					"Novel (p-Epitope)",
					"Novel (Temporal)"
				)
			)
		) |>
		dplyr::select(
			season, metric = newgroupf, att, est, lwr, upr
		) |>
		dplyr::arrange(season, metric, att)

	# Have to do some manually filtering of the current metrics cause they
	# have inconsequential (0.01 or 0.02) differences due to MC error
	# Just get the first 4000 draws since we arranged everything before
	dat_current <-
		dat_all |>
		dplyr::filter(metric == "Current") |>
		dplyr::group_by(season, metric, att) |>
		dplyr::slice_head(n = 1) |>
		dplyr::ungroup()

	dat_current_distinct <-
		dplyr::bind_rows(
			dat_current,
			dplyr::filter(dat_all, metric != "Current")
		) |>
		dplyr::arrange(season, metric, att) |>
		dplyr::filter(metric != "Novel (Temporal)")

	return(dat_current_distinct)
}

make_icc_table <- function(clean_icc_results_df, example_only, file_path) {
	dat_tbl <-
		clean_icc_results_df |>
		dplyr::rename(Season = season, "Metric Set" = metric) |>
		dplyr::filter(Season != "Overall") |>
		dplyr::mutate(
			dplyr::across(
				c(est, lwr, upr),
				\(x) formatC(x, digits = 2, width = 4, format = "f")
			),
			ci = paste0(est, " (", lwr, ", ", upr, ")")
		) |>
		dplyr::select(-c(est, lwr, upr)) |>
		tidyr::pivot_wider(names_from = att, values_from = ci)

	if (isTRUE(example_only)) {
		dat_tbl <- dat_tbl |>
			dplyr::filter(Season == "2016 - 2017") |>
			dplyr::select(-Season)

		tbl_out <- 	dat_tbl |>
			flextable::flextable()
	} else if (isFALSE(example_only)) {
		tbl_out <- dat_tbl |>
			flextable::flextable() |>
			flextable::merge_v(j = 1) |>
			flextable::valign(j = 1, valign = "top") |>
			flextable::fix_border_issues()
	} else {
		cli::cli_abort("{.arg example_only} should be {.val TRUE} or {.val FALSE}.")
	}

	final_tbl <-
		tbl_out |>
		fit_flextable_to_page()

	save_file_as_qs2(
		final_tbl,
		file_path
	)

	return(file_path)
}

# Metrics over time plot ####
## ICCs over time
make_icc_over_time_plot <- function(clean_icc_results_df, file_path) {
	clean_icc_results_df |>
		ggplot2::ggplot() +
		ggplot2::aes(
			x = season,
			y = est,
			ymin = lwr,
			ymax = upr
		) +
		ggplot2::geom_ribbon(
			ggplot2::aes(group = 1),
			alpha = 0.25, color = "transparent"
		) +
		ggplot2::geom_line(
			ggplot2::aes(group = 1)
		) +
		ggplot2::geom_point() +
		ggplot2::facet_grid(metric ~ att) +
		hgp::theme_ms()

	plt_out <-
		clean_icc_results_df |>
		ggplot2::ggplot() +
		ggplot2::aes(
			x = season,
			y = est,
			ymin = lwr,
			ymax = upr,
			color = metric
		) +
		ggplot2::geom_ribbon(
			ggplot2::aes(group = metric, fill = metric),
			alpha = 0.25, color = "transparent"
		) +
		ggplot2::geom_line(
			ggplot2::aes(group = metric)
		) +
		ggplot2::geom_point() +
		ggplot2::facet_wrap(~att, scales = "free_y") +
		ggplot2::scale_color_brewer(
			palette = "Dark2",
			aesthetics = c("color", "fill")
		) +
		hgp::theme_ms() +
		ggplot2::labs(
			x = "Season",
			y = "ICC",
			color = "Metric",
			fill = "Metric"
		)

	ggplot2::ggsave(
		filename = file_path,
		plot = out_plt,
		width = 13,
		height = 7
	)
}

## Actual metrics over time
make_metrics_over_time_plot <- function(bound_metrics_df, file_path) {
		cleaned_df <-
			bound_metrics_df |>
		dplyr::filter(
			censoring == "yes",
			subsample_id == "00",
			stats_id != "SCR0",
			season != "Overall"
		) |>
		dplyr::select(season, metric, stats_id, stats) |>
		dplyr::summarize(
			ggdist::mean_hdci(stats),
			.by = c(season, metric, stats_id)
		) |>
		dplyr::mutate(
			set = dplyr::case_match(
				stats_id,
				"GMT" ~ "Current",
				"GMT0" ~ "Current",
				"SCR" ~ "Current",
				"INT" ~ "Novel",
				"PAT" ~ "Novel",
				"AUC" ~ "Novel"
			) |> factor(levels = c("Current", "Novel")),
			att = dplyr::case_match(
				stats_id,
				"GMT" ~ "Total Strength",
				"GMT0" ~ "Magnitude",
				"SCR" ~ "Breadth",
				"INT" ~ "Magnitude",
				"PAT" ~ "Breadth",
				"AUC" ~ "Total Strength"
			) |> factor(levels = c("Magnitude", "Breadth", "Total Strength")),
			season_short = gsub(
				"([0-9]{4})\\s*-\\s*[0-9]{2}([0-9]{2})", "\\1/\\2", season
			)
		) |>
			dplyr::select(
				season = season_short, metric,
				est = y, lwr = ymin, upr = ymax, set, att
			)

		# Do the stupid comparing trick thing
		# Separate stuff so we can make the table narrower for word doc--
		# make the current metrics in the ag dist column
		dat_all <-
			cleaned_df |>
			dplyr::mutate(
				newgroup = ifelse(set == "Novel", as.character(metric), "curr"),
				newgroupf = factor(
					newgroup,
					levels = c("curr", "Cartographic", "p-Epitope", "Temporal"),
					labels = c(
						"Current",
						"Novel (Cartographic)",
						"Novel (p-Epitope)",
						"Novel (Temporal)"
					)
				)
			) |>
			dplyr::select(
				season, metric = newgroupf, att, est, lwr, upr
			) |>
			dplyr::arrange(season, metric, att)

		# Have to do some manually filtering of the current metrics cause they
		# have inconsequential (0.01 or 0.02) differences due to MC error
		# Just get the first 4000 draws since we arranged everything before
		dat_current <-
			dat_all |>
			dplyr::filter(metric == "Current") |>
			dplyr::group_by(season, metric, att) |>
			dplyr::slice_head(n = 1) |>
			dplyr::ungroup()

		dat_current_distinct <-
			dplyr::bind_rows(
				dat_current,
				dplyr::filter(dat_all, metric != "Current")
			) |>
			dplyr::arrange(season, metric, att) |>
			dplyr::filter(metric != "Novel (Temporal)")

		out_plt <-
			dat_current_distinct |>
			ggplot2::ggplot() +
			ggplot2::aes(
				x = season,
				y = est,
				ymin = lwr,
				ymax = upr,
				color = metric
			) +
			ggplot2::geom_line(
				ggplot2::aes(group = metric)
			) +
			ggplot2::geom_pointrange(
				ggplot2::aes(group = metric),
				alpha = 0.8
			) +
			ggplot2::facet_wrap(~att, scales = "free_y") +
			ggplot2::scale_color_brewer(
				palette = "Dark2",
				aesthetics = c("color", "fill")
			) +
			hgp::theme_ms() +
			ggplot2::labs(
				x = "Season",
				y = "Metric value",
				color = "Metric",
				fill = "Metric"
			)

		ggplot2::ggsave(
			filename = file_path,
			plot = out_plt,
			width = 13,
			height = 7
		)
}


# Head-to-head comparison ####
make_icc_comparisons <- function(icc_results_dataframe) {
	cleaned_df <-
		icc_results_dataframe |>
		dplyr::filter(
			censoring == "yes",
			stats_id != "SCR0",
			season != "Overall"
		) |>
		dplyr::select(season, metric, stats_id, ICC_res) |>
		dplyr::mutate(
			set = dplyr::case_match(
				stats_id,
				"GMT" ~ "Current",
				"GMT0" ~ "Current",
				"SCR" ~ "Current",
				"INT" ~ "Novel",
				"PAT" ~ "Novel",
				"AUC" ~ "Novel"
			) |> factor(levels = c("Current", "Novel")),
			att = dplyr::case_match(
				stats_id,
				"GMT" ~ "Total Strength",
				"GMT0" ~ "Magnitude",
				"SCR" ~ "Breadth",
				"INT" ~ "Magnitude",
				"PAT" ~ "Breadth",
				"AUC" ~ "Total Strength"
			) |> factor(levels = c("Magnitude", "Breadth", "Total Strength")),
			season = gsub(
				"([0-9]{4})\\s*-\\s*[0-9]{2}([0-9]{2})", "\\1/\\2", season
			)
		) |>
		dplyr::select(-stats_id)

	# Do the stupid comparing trick thing
	# Separate stuff so we can make the table narrower for word doc--
	# make the current metrics in the ag dist column
	dat_all <-
		cleaned_df |>
		dplyr::mutate(
			newgroup = ifelse(set == "Novel", as.character(metric), "curr"),
			newgroupf = factor(
				newgroup,
				levels = c("curr", "Cartographic", "p-Epitope", "Temporal"),
				labels = c(
					"Current",
					"Novel (Cartographic)",
					"Novel (p-Epitope)",
					"Novel (Temporal)"
				)
			)
		) |>
		dplyr::select(
			season, metric = newgroupf, att, ICC_res
		) |>
		dplyr::arrange(season, metric, att)

	# Have to do some manually filtering of the current metrics cause they
	# have inconsequential (0.01 or 0.02) differences due to MC error
	# Just get the first 4000 draws since we arranged everything before
	dat_current <-
		dat_all |>
		dplyr::filter(metric == "Current") |>
		dplyr::group_by(season, metric, att) |>
		dplyr::slice_head(n = 4000) |>
		dplyr::ungroup()

	dat_current_distinct <-
		dplyr::bind_rows(
			dat_current,
			dplyr::filter(dat_all, metric != "Current")
		) |>
		dplyr::arrange(season, metric, att) |>
		dplyr::filter(metric != "Novel (Temporal)")

	dat_res <- dat_current_distinct |>
		dplyr::mutate(
			.draw = dplyr::row_number(),
			.by = c(season, att, metric)
		) |>
		tidyr::pivot_wider(names_from = metric, values_from = ICC_res)

	# Now that we did the data cleaning we construct the comparisons
	metrics <- dat_current_distinct$metric |> unique()
	combs <- combn(metrics, 2)


	dat_nested <- 	dat_res |>
		tidyr::nest(data = -c(season, att))

	contrast_samples <- purrr::map(
		1:ncol(combs),
		\(j) {
			met1 <- combs[[1, j]] |> as.character()
			met2 <- combs[[2, j]] |> as.character()

			contr <- purrr::map(
				dat_nested$data,
				\(d) d[[met2]] - d[[met1]]
			)

			out <-
				tibble::add_column(dat_nested, contr = contr) |>
				dplyr::select(-data)

			return(out)
		}
	)

	contrast_df <-
		tibble::tibble(
			contrast = paste0(combs[2, ], " - ", combs[1, ]),
			sample = contrast_samples
		) |>
		tidyr::unnest(sample) |>
		tidyr::unnest(contr)

	contrast_res <-
		contrast_df |>
		dplyr::summarize(
			ggdist::mean_hdci(contr),
			bayestestR::rope(contr, range = c(-0.1, 0.1)),
			.by = c(contrast, season, att)
		)

	return(contrast_res)
}
