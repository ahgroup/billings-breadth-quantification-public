###
# Utility functions used in multiple scripts
# Zane Billings
# 2025-04-28
###

# Function for minmax scaling a numeric vector
minmax <- function(x, ...) {
	stopifnot(is.numeric(x))
	return((x - min(x, ...)) / (max(x, ...) - min(x, ...)))
}

# Convenience function for saving one output as a qs2 file
# for use with tar_file.
save_file_as_qs2 <- function(file_to_save, file_path) {
	file_path |>
		dirname() |>
		dir.create(showWarnings = FALSE, recursive = TRUE)
	qs2::qs_save(file_to_save, file_path)
	return(file_path)
}

# Convenience function to ensure the directory for saving a file exists
ensure_directory_exists <- function(file_path) {
	dir.create(dirname(file_path), showWarnings = FALSE, recursive = TRUE)

	invisible(file_path)
}

# Function for making flextables the correct size
fit_flextable_to_page <- function(tbl, page_width = 6.5){
	ft_out <- tbl |> flextable::autofit()

	ft_width <- (dim(ft_out)$widths * page_width) /
		(flextable::flextable_dim(ft_out)$widths)

	ft_out <- flextable::width(ft_out, width = ft_width)

	return(ft_out)
}

# Function to set the gtsummary theme
set_gtsummary_theme <- function() {
	suppressMessages({
		gtsummary::theme_gtsummary_journal(journal = "jama")
		gtsummary::theme_gtsummary_compact(font_size = 10)
		gtsummary::theme_gtsummary_language(language = "en", big.mark = "")
	})

	invisible(TRUE)
}

# Function to add zeroes to the left side of a vector of numbers
# TODO add to zlib
zero_pad_numbers <- function(x) {
	# Get the number of digits in the largest number
	max_digits <- max(nchar(x))

	# Pad the numbers
	out <- stringr::str_pad(
		x,
		width = max_digits,
		side = "left",
		pad = "0",
		use_width = TRUE
	)

	return(out)

}

head_tail <- function(dat, n = 6) {
	if (!is.numeric(n) || n >= nrow(dat)) {
		stop("'n' should be an integer and < nrow(dat)!")
	} else if (length(n) == 1) {
		head_n <- n
		tail_n <- n
	} else if (length(n) == 2) {
		head_n <- n[[1]]
		tail_n <- n[[2]]
	} else {
		stop("'n' should be length 1 or 2.")
	}

	header <- head(dat, n = head_n)
	footer <- tail(dat, n = tail_n)

	middle <- dat[0, ]
	middle[1, 1] <- "$\\ \\vdots$"

	out <- dplyr::bind_rows(header, middle, footer)
	return(out)
}
