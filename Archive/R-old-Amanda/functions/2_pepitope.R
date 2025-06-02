## Functions for 2_pepitope_calculator.Rmd

## epi() Function : Selects the characters of a specific antigenic site
## ######### Description: Selects the characters of a specific antigenic site
## from an amino acid string. Convert the string into a vector of characters,
## and select the characters that correspond to the given antigenic site.
## Re-condense the vector of characters into a string, and return that string.

epi <- function(sequence,
                antigenic_site)
{
  # Arguments "sequence" = a character string of the amino acid sequence NOT
  # starting from the Methionine start, but instead after the leader sequence
  # has been removed. Please refer to Antigenic Evergreen residues for
  # assistance.

  # "antigenic_site" = a vector of integers indicating the amino acid residues
  # of the antigenic site of interest

  vec <- seqinr::s2c(sequence)
  site <- vec[antigenic_site]
  epit <- seqinr::c2s(site)

  return(epit)

  #Value: Returns the antigenic site as a string

}


## combin() Function : Adds a column of a specific antigenic site to the overall
## dataframe #######

combin <- function(epitope_site,
                   dataframe,
                   reference_column = "ha_0",
                   ignore_starting_sequence_warning = FALSE)
{
  # "epitope_site" = epitope site name of interest.

  # "dataframe" = dataframe that will be expanded by one column that contains
  # the extracted epitope site from the full HA0 protein that is contained
  # within a column of the dataframe as string. The epitope site will be
  # returned as a string.

  # "ha_0" = The name of column that contains the full HA0 protein sequence.

  # "ignore_starting_sequence_warning" = This is the warning that the HA0
  # protein sequence does not start with a methionine. To override this warning,
  # change to TRUE. Arguments: c(TRUE, FALSE)

  # Purpose: Adds a column of a specific antigenic site to the overall
  # dataframe. Using the full amino acid sequence use the previously defined epi
  # function to select a specific antigenic site, and add that character string
  # as a column to a pre-existing dataframe

  # Pull the full amino acid sequence that does not contain the leader sequence
  colnames(dataframe)[colnames(dataframe) == reference_column] <- "ref_col"
  df <- dataframe$ref_col

  # Check to make sure that the HA0 starts with Methionine start amino acid
  #  if(ignore_starting_sequence_warning == FALSE)
  #    if(sum(toupper(str_sub(df, 1, 1)) != "D) > 0)
  #      stop("Warning: HA0 does not start with Methionine. Double check starting sequence")

  # Using the previously defined epi function, extract the select epitope residues
  xyz <- lapply(df,
                FUN = epi,
                epitope_site)
  xyz <- unlist(xyz)

  # Add the column to the previous dataframe
  df_new<- cbind(dataframe, xyz)

  # Rename the new column after the epitope site used
  name <- deparse(substitute(epitope_site))
  names(df_new)[length(names(df_new))]<- name

  # Rename the reference column back to it's original name
  colnames(df_new)[colnames(df_new) == "ref_col"] <- reference_column


  # Return the dataframe with the newly named column.
  return(df_new)
}

## string.diff.ex() Function : Determining the number of differences between
## strings of the same length######


# Description: Compare two strings (a and b) to each other to determine the
# number of differences. There is the option to not include characters (exclude)
# in the the comparison if they are unknown. There is also the option to ignore
# the case of the strings.

# Usage:

string.diff.ex <- function(a,
                           b,
                           exclude=c("x","X","?","-"),
                           ignore.case=TRUE) {

  # Arguments:
  # "a" = Input character string

  # "b" = Comparison character string for "a"

  # "exclude" = Option to exclude certain characters from comparison. These are
  # usually ambiguities in the amino acid sequence. "x", "X", and "?" are
  # commonly used characters for ambiguities

  # "ignore.case" = Option to ignore the case of the amino acids. This is set to TRUE

  if(nchar(a)!=nchar(b))
    stop("Warning: Lengths of input strings differ. Please check your input.")

  if(ignore.case==TRUE)
  {
    a<-toupper(a)
    b<-toupper(b)
  }

  diff.a<-unlist(base::strsplit(a,split=""))
  diff.b<-unlist(base::strsplit(b,split=""))
  diff.d<-rbind(diff.a,diff.b)

  for(ex.loop in 1:length(exclude))
  {
    diff.d<-diff.d[,!(diff.d[1,]==exclude[ex.loop]|diff.d[2,]==exclude[ex.loop])]
  }

  differences<-sum(diff.d[1,]!=diff.d[2,])

  return(differences)

  # Value: The number of differences between two strings will be returned.
}


## difference_matrix() Function : Determination of the number of differences
## between two strains #######

#Data should be a dataframe that is
difference_matrix <- function(data) {
  for (i in 1:length(data[,1])) {
    for (j in 1:length(data[,1])) {
      data[j,2+i] <- string.diff.ex(a = data[i,2], b = data[j,2])
    }
  }

  colnames(data)[-(1:2)] <- c(data[,1])

  return(data)
}




