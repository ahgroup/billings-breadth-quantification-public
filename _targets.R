# Created by use_targets().
# Follow the comments below to fill in this target script.
# Then follow the manual to check and run the pipeline:
#   https://books.ropensci.org/targets/walkthrough.html#inspect-the-pipeline

# Load packages required to define the pipeline:
suppressWarnings({suppressPackageStartupMessages({
	library(targets)
	library(tarchetypes)
	library(crew)
	library(crew.cluster)
	library(brms)
	library(qs2)
})})

# Detect whether we're on the cluster or not -- if we are, the global variable
# "SLURM_JOB_ID" will be set, if we aren't it will not be set.
# We'll use this to dispatch the correct crew controller without having to
# manually change settings at runtime.
hpc <- !is.na(Sys.getenv("SLURM_JOB_ID", unset = NA))

# Set up crew controllers
controller_hpc <- crew.cluster::crew_controller_slurm(
	name = "hpc",
	workers = 100,
	seconds_idle = 120,  # time until workers are shut down after idle
	options_cluster = crew.cluster::crew_options_slurm(
		script_lines = c(
			"module load R/4.4.2-gfbf-2024a"
			#add additional lines to the SLURM job script as necessary here
		),
		log_output = "logs/crew_bq_log_%A.out",
		log_error = "logs/crew_bq_log_%A.err",
		memory_gigabytes_per_cpu = 4,
		cpus_per_task = 4,
		time_minutes = 1200, # wall time for each worker
		partition = "batch",
		n_tasks = 1,
		verbose = TRUE
	)
)

controller_local <- crew::crew_controller_local(
	name = "local",
	workers = 10,
	options_local = crew::crew_options_local(log_directory = "logs")
	# options_metrics = crew_options_metrics(
	# 	path = "logs",
	# 	seconds_interval = 1
	# )
)

# Set target options:
tar_option_set(
	# Packages that are required for running the pipeline but not necessarily for
	# defining the targets
	packages = c("qs2"),
	# Optionally set the default storage format. qs is fast.
	format = "qs",
	# Set up parallel options via crew
	controller = crew::crew_controller_group(
		controller_hpc,
		controller_local
	),
	resources = targets::tar_resources(
		# If we detect that we're on slurm, by default use the small hpc controller
		# otherwise use the local controller
		# The large hpc controller must be manually specified
		crew = tar_resources_crew(controller = ifelse(hpc, "hpc", "local"))
	)
)

# Run the R scripts in the R/ folder with your custom functions:
targets::tar_source()

# This global tells us if we're running a test run for debugging or a real run
# for actual running.
IS_TEST_RUN <- FALSE

# Model setup that allows us to do dynamic-within-static branching for
# better version control
# Right now it's all dynamic because that is easier to code
# brms_model_info <- generate_brms_model_info()

# Replace the target list below with your own:
list(
	# Targets for raw data
	## UGAFluVac Cohort data set
	### This data is the output of the cleaning process from the repo
	### ahgroup/UGAFluVac-data/, which includes the raw excel files.
	tarchetypes::tar_file_read(
		name = raw_UGAFluVac_data,
		command = here::here("Data", "Raw", "UGAFluVac-cohort-data.Rds"),
		read = readr::read_rds(!!.x)
	),
	## Antigenic distance data set
	### Load the clean antigenic distance data which is exported from the
	### ahgroup/influenza-antigenic-distance-data repo.
	tarchetypes::tar_file_read(
		name = raw_antigenic_distance_data,
		command = here::here("Data", "Raw", "antigenic-distance-data.Rds"),
		read = readr::read_rds(!!.x)
	),
	# Cleaning data
	targets::tar_target(
		name = prepped_cohort_data,
		command = prep_cohort_data(raw_UGAFluVac_data)
	),
	targets::tar_target(
		name = model_data,
		command = create_model_data(
			prepped_cohort_data,
			raw_antigenic_distance_data
		)
	),
	targets::tar_target(
		name = reporting_data,
		command = create_civics_reporting_data(prepped_cohort_data)
	),
	# Save datasets to file
	tarchetypes::tar_file(
		name = model_data_files,
		command = write_model_data_to_file(
			model_data,
			here::here("Data", "Processed", "Model-Data")
		)
	),
	tarchetypes::tar_file(
		name = reporting_data_files,
		command = write_model_data_to_file(
			reporting_data,
			here::here("Data", "Processed", "Reporting-Data")
		)
	),
	# Descriptive analyses
	## Assay counts table
	tarchetypes::tar_file(
		name = assay_counts_table,
		command = make_assay_counts_table(
			model_data,
			here::here("Results", "Tables", "Assay-Counts-Table.qs2")
		)
	),
	tarchetypes::tar_file(
		name = person_counts_table,
		command = make_person_counts_table(
			model_data,
			here::here("Results", "Tables", "Person-Counts-Table.qs2")
		)
	),
	tarchetypes::tar_file(
		name = individual_counts,
		command = get_unique_individual_count(
			model_data,
			here::here("Results", "Tables", "Unique-Individual-Counts.qs2")
		)
	),
	## Demographics table
	tarchetypes::tar_file(
		name = demographics_table,
		command = make_demographics_table(
			model_data,
			here::here("Results", "Tables", "Demographics-Table.qs2")
		)
	),
	## Raw data figures for supplement
	tarchetypes::tar_file(
		name = titer_data_figure,
		command = make_titer_plot(
			model_data,
			here::here("Results", "Figures", "Titer-Plot.png")
		)
	),
	## Strain short/full names table
	tarchetypes::tar_file(
		name = strain_names_table,
		command = make_strain_name_table(
			model_data,
			here::here("Results", "Tables", "Strain-Names.qs2")
		)
	),
	## Vaccine strains table
	# tarchetypes::tar_file(
	# 	name = vaccine_strain_table,
	# 	command = make_vaccine_strain_table(
	# 		model_data,
	# 		here::here("Results", "Tables", "Vaccine-Strains.qs2")
	# 	)
	# ),
	# Model fitting
	## General setup
	targets::tar_target(
		name = brms_model_info,
		command = generate_brms_model_info()
	),
	targets::tar_target(
		name = brms_control_arguments,
		command = generate_brms_control_arguments(IS_TEST_RUN)
	),
	tarchetypes::tar_file_read(
		name = random_seeds,
		command = here::here("Data", "Raw", "seeds.txt"),
		read = scan(!!.x, quiet = TRUE)
	),
	targets::tar_target(
		name = nested_model_data,
		command = nest_model_data(model_data)
	),
	targets::tar_target(
		name = model_metadata,
		command = generate_model_metadata(
			brms_model_info,
			nested_model_data,
			random_seeds,
			global_seed = 8751
		)
	),
	tarchetypes::tar_file(
		name = model_metadata_file,
		command = save_file_as_qs2(
			model_metadata,
			here::here("Results", "Output", "Model-Metadata.qs2")
		)
	),
	targets::tar_target(
		name = model_metadata_list,
		command = listify_df(model_metadata)
	),
	# "Full Data" models
	# i.e. models fitted to the entire dataset for each cohort
	targets::tar_target(
		name = full_data_model_fits,
		command = model_dispatch(
			model_metadata_list,
			brms_control_arguments
		),
		pattern = map(model_metadata_list),
		iteration = "list",
		# Don't keep huge models in memory, we'll read them from file again later
		memory = "transient",
		garbage_collection = TRUE
	),
	tar_target(
		name = full_data_model_fit_file_names,
		command = generate_brms_model_filenames(
			model_metadata,
			here::here("Results", "Large-Files", "Full-Data-Models")
		)
	),
	tarchetypes::tar_file(
		name = full_data_model_fit_files,
		command = save_file_as_qs2(
			full_data_model_fits,
			full_data_model_fit_file_names
		),
		pattern = map(
			full_data_model_fits,
			full_data_model_fit_file_names
		)
	),
	# Subsampling analysis
	## Define the subsample parameters
	targets::tar_target(
		name = subsampling_parameters,
		command = generate_subsampling_parms_list(
			subsamples_per_cohort = 25,
			individuals_per_subsample = 100,
			heterologous_strains_per_subsample = 9
		)
	),
	## Now make the subsample data
	targets::tar_target(
		name = subsample_data,
		command = make_subsample_dataset(
			model_data,
			subsampling_parameters,
			# seed from random dot org
			subsampling_seed = 6133124L
		)
	),
	## Save the subsample data to a file
	tarchetypes::tar_file(
		name = subsample_data_files,
		command = write_model_data_to_file(
			subsample_data,
			here::here("Data", "Processed", "Subsample-Data")
		)
	),
	## Make the subsample model metadata
	targets::tar_target(
		name = subsample_model_metadata,
		command = generate_model_metadata(
			brms_model_info,
			subsample_data,
			random_seeds,
			global_seed = 6793
		)
	),
	tarchetypes::tar_file(
		name = subsample_model_metadata_files,
		command = save_file_as_qs2(
			subsample_model_metadata,
			here::here("Results", "Output", "Model-Metadata-Subsamples.qs2")
		)
	),
	## Listify the subsample model info
	targets::tar_target(
		name = subsample_model_metadata_list,
		command = listify_df(subsample_model_metadata)
	),
	# Subsample model fitting
	targets::tar_target(
		name = subsample_model_fits,
		command = model_dispatch(
			subsample_model_metadata_list,
			brms_control_arguments
		),
		pattern = map(subsample_model_metadata_list),
		iteration = "list",
		# Don't keep huge models in memory, we'll read them from file again later
		memory = "transient",
		garbage_collection = TRUE
	),
	tar_target(
		name = subsample_model_fit_file_names,
		command = generate_brms_model_filenames(
			subsample_model_metadata,
			here::here("Results", "Large-Files", "Subsample-Models")
		)
	),
	tarchetypes::tar_file(
		name = subsample_model_fit_files,
		command = save_file_as_qs2(
			subsample_model_fits,
			subsample_model_fit_file_names
		),
		pattern = map(
			subsample_model_fits,
			subsample_model_fit_file_names
		)
	),
	# Metrics calculation
	targets::tar_target(
		name = full_data_metrics,
		command = calculate_metrics_for_model(
			full_data_model_fits,
			model_metadata_list
		),
		pattern = map(
			full_data_model_fits,
			model_metadata_list
		),
		iteration = "list"
	),
	targets::tar_target(
		name = subsample_metrics,
		command = calculate_metrics_for_model(
			subsample_model_fits,
			subsample_model_metadata_list
		),
		pattern = map(
			subsample_model_fits,
			subsample_model_metadata_list
		),
		iteration = "list"
	),
	# Combining and cleaning the metrics
	targets::tar_target(
		name = bound_metrics_df,
		command = bind_calculated_metrics(full_data_metrics, subsample_metrics)
	),
	targets::tar_target(
		name = metrics_summary_df,
		command = summarize_calculated_metrics(bound_metrics_df)
	),
	# ICC across subsamples
	targets::tar_target(
		name = nested_stats_dataframe,
		command = nest_stats_dataframe(bound_metrics_df, random_seeds, 152L)
	),
	targets::tar_target(
		name = icc_control_arguments,
		command = generate_icc_control_arguments(IS_TEST_RUN)
	),
	targets::tar_target(
		name = icc_models,
		command = fit_icc_model(
			nested_stats_dataframe$icc_data,
			nested_stats_dataframe$seed,
			icc_control_arguments
		),
		pattern = map(nested_stats_dataframe),
		iteration = "list"
	),
	targets::tar_target(
		name = icc_results,
		command = calculate_icc_from_model(icc_models),
		pattern = map(icc_models),
		iteration = "list"
	),
	targets::tar_target(
		name = icc_results_dataframe,
		command = bind_iccs(nested_stats_dataframe, icc_results)
	),
	targets::tar_target(
		name = icc_summary_dataframe,
		command = summarize_icc_results(icc_results_dataframe)
	),
	# Model results and figures
	## Summary landscape figures
	targets::tar_target(
		name = summary_landscape_predictions,
		command = get_summary_landscape(
			full_data_model_fits,
			model_metadata
		)
	),
	tarchetypes::tar_file(
		name = summary_landscape_figure_files,
		command = make_all_summary_landscape_plots(
			summary_landscape_predictions,
			model_metadata,
			here::here("Results", "Figures")
		)
	),
	## Sample metrics table
	tarchetypes::tar_file(
		name = example_metrics_table,
		command = make_model_metrics_table(
			metrics_summary_df,
			TRUE,
			here::here("Results", "Tables", "Example-Metrics-Table.qs2")
		)
	),
	tarchetypes::tar_file(
		name = all_seasons_metrics_table,
		command = make_model_metrics_table(
			metrics_summary_df,
			FALSE,
			here::here("Results", "Tables", "All-Season-Metrics-Table.qs2")
		)
	),
	## ICC plot
	tarchetypes::tar_file(
		name = icc_plots,
		command = create_icc_plots(
			bound_metrics_df,
			here::here("Results", "Figures")
		)
	),
	# ICC Tables
	targets::tar_target(
		name = clean_icc_data_for_results,
		command = clean_up_icc_results(icc_summary_dataframe)
	),
	tarchetypes::tar_file(
		name = example_icc_table,
		command = make_icc_table(
			clean_icc_data_for_results,
			TRUE,
			here::here("Results", "Tables", "Example-ICC-Table.qs2")
		)
	),
	tarchetypes::tar_file(
		name = full_icc_table,
		command = make_icc_table(
			clean_icc_data_for_results,
			FALSE,
			here::here("Results", "Tables", "All-Seasons-ICC-Table.qs2")
		)
	),
	# Make the sequence accession table for the supplement
	tarchetypes::tar_file_read(
		name = seq_accession_data,
		command = here::here("Data", "Raw", "sequence-accession.csv"),
		read = readr::read_csv(!!.x)
	),
	tarchetypes::tar_file(
		name = seq_accession_table,
		command = make_seq_accession_table(
			seq_accession_data,
			here::here("Results", "Tables", "Seq-Accession-Table.qs2")
		)
	),
	# Contrast calculation
	targets::tar_target(
		name = contrast_dataframe,
		command = make_icc_comparisons(icc_results_dataframe)
	),
	# Make contrast tables
	# TODO
	# Plots over time
	## Metrics
	tarchetypes::tar_file(
		name = icc_over_time_plot,
		command = make_icc_over_time_plot(
			clean_icc_data_for_results,
			here::here("Results", "Figures", "ICC-Over-Time.png")
		)
	),
	## ICC
	tarchetypes::tar_file(
		name = metrics_over_time_plot,
		command = make_metrics_over_time_plot(
			bound_metrics_df,
			here::here("Results", "Figures", "Metrics-Over-Time.png")
		)
	),
	## Contrast
	## TODO
	# Endpoint creation
	## First generate the software bibliography
	tar_file(
		name = software_bibliography,
		command = generate_software_bibliography(
			here::here("Assets", "package-refs.bib")
		)
	),
	## Now render the manuscript
	tarchetypes::tar_quarto(
		name = manuscript,
		path = here::here("Products", "Manuscript", "Manuscript.qmd"),
		extra_files = c(
			here::here("Assets", "billings-breadth-quantification.bib"),
			here::here("Assets", "package-refs.bib"),
			here::here("Assets", "vancouver.csl"),
			here::here("Assets", "word-template.docx")
		),
		quiet = FALSE
	),
	## And the supplement
	tarchetypes::tar_quarto(
		name = supplement,
		path = here::here("Products", "Manuscript", "Supplement.qmd"),
		extra_files = c(
			here::here("Assets", "billings-breadth-quantification.bib"),
			here::here("Assets", "package-refs.bib"),
			here::here("Assets", "vancouver.csl"),
			here::here("Assets", "word-template.docx")
		),
		quiet = FALSE
	),
	## And the supplementary methods document
	tarchetypes::tar_quarto(
		name = supplementary_methods,
		path = here::here("Products", "Methods-Doc", "Methods.qmd"),
		extra_files = c(
			here::here("Products", "Methods-Doc", "methods-ex-dat.qs2"),
			here::here("Products", "Methods-Doc", "data-sampler.R"),
			here::here("Products", "Methods-Doc", "censoring_optimized.png"),
			here::here("Assets", "word-template.docx")
		),
		quiet = FALSE,
		cache = TRUE
	),
	## Finally render the README and we are done!
	tarchetypes::tar_quarto(
		name = README,
		path = here::here("README.qmd"),
		quiet = FALSE
	)
)
