###
# Miscellaneous Stuff
# Zane
# 2025-05-07
###

# Generate software bibliography -- this should always be right before
# the manuscript and supplement
generate_software_bibliography <- function(file_path) {
	require(softbib, quietly = TRUE)

	unlink(file_path)
	softbib::softbib(
		output = file_path
	)

	return(file_path)
}

# Make the sequence accession table
make_seq_accession_table <- function(seq_accession_data, file_path) {
	out <-
		seq_accession_data[, c(2, 4)] |>
		dplyr::rename("Strain" = "Short Name") |>
		flextable::flextable() |>
		fit_flextable_to_page()

	save_file_as_qs2(
		out,
		file_path
	)

	return(file_path)
}
