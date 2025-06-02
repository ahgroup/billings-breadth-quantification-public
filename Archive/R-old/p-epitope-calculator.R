###
# P-epitope calculator
# Zane Billings
# 2024-05-04
# Calculate the p-epitope distance between two aligned influenza sequences
###

# First we need to create a reference data frame of p-epitope site locations
get_pepitope_sites <- function(
		subtype,
		sites = c('full', 'a', 'b', 'c', 'd', 'e', 'all_epitopes')
	) {
	
	# standardize subtype
	subtype <- stringr::str_to_lower(subtype)
	
	# The sites are different for H1 and H3. These are based on Amanda's
	# code where she referenced the actual papers.
	if (startsWith(subtype, 'h1')) {
		full <- 1:326
		
		site_a <- c(
			118, 120:122, 126:129, 132:135, 137, 139:143, 146, 147, 149, 165, 252,
			253
		)
		
		site_b <- c(
			124, 125, 152:157, 160, 162, 183:187, 189:191, 193:196
		)
		
		site_c <- c(
			34:38, 40, 41, 43:45, 269:274, 276:278, 283, 288, 292, 295, 297, 298,
			302, 303, 305:310
		)
		
		site_d <- c(
			89, 94:96, 113, 117, 163, 164, 166:174, 176, 179, 198, 200, 202, 204:216,
			222:227, 235, 237, 241, 243:245
		)
		
		site_e <- c(
			47, 48, 50, 51, 53, 54, 56:58, 66, 68:75, 78:80, 82:86, 102, 257:261,
			263, 267
		)
		
	} else if (startsWith(subtype, 'h3')) {
		full <- 1:328
		
		site_a <- c(
			122, 124, 126, 130:133, 135, 137:138, 140, 142:146, 150, 152,
			168
		)
		
		site_b <- c(
			128, 129, 155:160, 163, 165, 186:190, 192:194, 196:198
		)
		
		site_c <- c(
			44:48, 50, 51, 53, 54, 273, 275:276, 278:280, 294, 297, 299, 300,
			304:305, 307:312
		)
		
		site_d <- c(
			96, 102, 103, 117, 121, 167, 170:177, 179, 182, 201, 203, 207:209,
			212:219, 226:230, 238, 240, 242, 244, 246:248
		)
		
		site_e <- c(
			57, 59, 62, 63, 67, 75, 78, 80:83, 86:88, 91, 92, 94, 109, 260:262,
			265
		)
	} else {
		stop("'subtype' should be 'h1n1' or 'h3n2'.")
	}
	
	# Make a vector containing all epitope residues
	all <- c(site_a, site_b, site_c, site_d, site_e) |>
		unique() |>
		sort()
	
	# Make a list containing all the vectors
	res <- list(full, site_a, site_b, site_c, site_d, site_e, all)
	names(res) <- c('full', 'a', 'b', 'c', 'd', 'e', 'all_epitopes')
	
	# Now return the sites that were requested before
	out <- res[sites]
	return(out)
}

# Utility function for subsetting characters in a single string (in R language,
# getting a substring by position for a character vector of length 1)
extract_string_chars <- function(str, pos) {
	sub <- charToRaw(str)[pos]
	return(rawToChar(sub))
}

# Utility function for computing the hamming distance between two strings.
hamming_dist <- function(string1, string2, norm = TRUE) {
	if (nchar(string1) != nchar(string2)) {
		rlang::abort('The two strings must be the same length.')
	}
	
	# Split the strings into vector of single characters
	x <- stringr::str_split_1(string1, "")
	y <- stringr::str_split_1(string2, "")
	
	# Calculate the distance
	d <- sum(x != y)
	
	# Normalized if asked
	if (isTRUE(norm)) {
		return(d / nchar(string1))
	} else {
		return(d)
	}
}

# Compute the p-epitope distance between two strings
pepitope <- function(seq_1, seq_2, subtype, mode = "dominant") {
	# Get the numbers for the residues in each epitope
	p_epi_sites <- get_pepitope_sites(subtype, sites = c('a', 'b', 'c', 'd', 'e'))
	
	epi_dists <-
		purrr::map_dbl(
			p_epi_sites,
			# This function looks kind of weird because doing this calculation is
			# really annoying to make sure the Hamming dist is the same as the
			# overall Hamming dist and treats special chars the same.
			# First we subset the two strings to only contain residues from the
			# current epitope.
			\(current_sites) list(
				extract_string_chars(seq_1, current_sites),
				extract_string_chars(seq_2, current_sites)
			) |>
				do.call(what = hamming_dist)
		)
	
	# Now decide whether to return all epitopes or not
	if (mode == "dominant") {
		# Dominant p-epitope = maximum p-epitope
		p_epi <- max(epi_dists)
	} else if (mode == "average") {
		p_epi <- mean(epi_dists)
	} else {
		p_epi <- epi_dists
	}
	return(p_epi)
}

# Compute a p-epitope distance matrix for all sequences in a character vector
dist_pepi <- function(seqs, subtype) {
	# Set up an empty matrix to hold results
	res <- matrix(
		nrow = length(seqs),
		ncol = length(seqs),
		dimnames = list(names(seqs), names(seqs))
	)
	
	# Calculate the lower triangle of the distance matrix for all unique
	# combinations of two strains
	for (i in 2:nrow(res)) {
		for (j in 1:(i-1)) {
			res[[i, j]] <- pepitope(seqs[[i]], seqs[[j]], subtype)
		}
	}
	
	# Set the diagonal to zero
	diag(res) <- 0
	
	# Quickly fill in the upper triangle
	out <- res |> as.dist() |> as.matrix()
	
	return(out)
}

# test <- dist.pepi(prot_h1, 'h1')

##### END OF FILE ####
