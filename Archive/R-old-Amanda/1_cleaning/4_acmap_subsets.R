## "Addition of Acmap Data to Antigenic Landscape"
## "Amanda Skarlupka"
## date: "3/3/2022"

# Purpose ####

## This files creates the final dataset that contains: antigenic distances from the acmap objects, and the sequence based distances from the amino acid sequences.

# Prerequisites:

#First run:

# ../code/1_cleaning/1_cleaning_code.R
# ../code/1_cleaning/2_pepitope_calculator.Rmd
# ../code/1_cleaning/3_acmap_creation.R

library(Racmacs) #Calculating antigenic distance and creating maps
library(tidyverse) #Working with data

#Functions for working with Acmap objects
source(here::here("R/functions/antigenic-cartography.R"))

#Read in the data
#There is minor cleaning the Pre viruses have duplicate names

acmap_file_location <- "data/processed/subset_analysis/cartography"
dataset_file_location <- "data/processed/subset_analysis/datasets"
subtype <- "H1N1"

## Datasets for the subsets:
data_list <- paste0(dataset_file_location, "/", list.files(path = here::here(dataset_file_location))) #Pull all of the subset file names from the datasets folder

subset_datasets <- lapply(X = data_list, 
                          FUN = function(.) readRDS(file = here::here(.))) #Pull the actual dataframes

names(subset_datasets) <- stringr::str_match(data_list, "datasets/\\s*(.*?)\\s*.rds")[,2] #Name the subset based on the subset number

subset_datasets <- data.table::rbindlist(subset_datasets,
                                         idcol = "subset")

## Sequence Distances:
### Even for subsets the p-epitope and year distances remain the same. 
sequence_distance <- vector(mode = "list")
for (i in 1:2) {
  x <- c("h1", "h3")[i] 
sequence_distance[[x]] <- readRDS(file = here::here(paste0("data/processed/", 
                                             x, 
                                             "_seq_distance/", 
                                             x, 
                                             "_key.rds"))) %>%
  dplyr::filter(method == "p_epi") %>%
  tidyr::pivot_longer(cols = -c(analysis_name, method),
               names_to = "strains_fullname",
               values_to = "distance")
}

sequence_distance <- data.table::rbindlist(sequence_distance, 
                                           idcol = "strain_type") %>%
  mutate(strain_type = dplyr::recode(strain_type, 
                                     h1 = "H1N1", 
                                     h3 = "H3N2"))

## Distances for 2d post subset maps:

### Load the antigenic Maps 

file_list <- paste0(acmap_file_location, "/", list.files(path = here::here(acmap_file_location))) #Pull all of the subset file names from the datasets folder

foo <- function(location) {
  AcmapToDF(map = Racmacs::read.acmap(here::here(location))) %>%
    dplyr::filter(type == "antigen") %>%
    dplyr::mutate(method = "cart_2d_post",
                  subtype = stringr::str_match(location, "cartography/\\s*(.*?)\\s*_subset")[,2],
                  subset = paste0("subset_", stringr::str_match(location, "_subset_\\s*(.*?)\\s*.ace")[,2])) %>%
    dplyr::select(V1, V2, identifier, method, subtype, subset)
}

map_list <- lapply(X = file_list, foo) #Pull the actual dataframes
names(map_list) <- stringr::str_match(file_list, "cartography/\\s*(.*?)\\s*.ace")[,2] #Name the subset based on the subset number
antigen <- data.table::rbindlist(map_list)

## Calculate the Euclidean distances
antigenic_distance <- subset_datasets %>%
  dplyr::select(subset, strain_type, dplyr::ends_with("_fullname")) %>%# Identifies the vaccine strains to be used, and comparison distances
  tidyr::pivot_longer(cols = c(dplyr::ends_with("vaccine_fullname")),
                    names_to = "vaccine_strain_type",
                    values_to = "vaccine_fullname") %>%
  dplyr::mutate(vaccine_strain_type = ifelse(vaccine_strain_type == "h1n1_vaccine_fullname",
                                             "H1N1",
                                             ifelse(vaccine_strain_type == "h3n2_vaccine_fullname",
                                                    "H3N2",
                                                    "typeb"))) %>%
  dplyr::filter(vaccine_fullname != "None",
                strain_type == vaccine_strain_type) %>%
  dplyr::select(-vaccine_strain_type) %>%
  dplyr::distinct() %>%
  dplyr::right_join(antigen,
             by = c("strains_fullname" = "identifier", "strain_type" = "subtype", "subset")) %>%
  base::droplevels() %>%
  dplyr::group_by(subset, strain_type, vaccine_fullname) %>%
  dplyr::mutate(distance = sqrt((V1-V1[which(strains_fullname==paste(vaccine_fullname))])^2+(V2-V2[which(strains_fullname==paste(vaccine_fullname))])^2)) %>% #Euclidean distance from vaccine strain
  dplyr::ungroup() %>%
  dplyr::select(-c(V1:V2))


## Calculate the Sequence Distances

seq <- sequence_distance %>%
  dplyr::filter(analysis_name %in% levels(antigenic_distance$vaccine_fullname)) %>%
  dplyr::rename(vaccine_fullname = analysis_name)
  
cent <- antigenic_distance %>%
  dplyr::left_join(seq, 
                   by = c("strain_type", "strains_fullname", "vaccine_fullname")) %>%
  tidyr::pivot_longer(cols = c(distance.x, distance.y),
               values_to = "distance",
               names_to = "method") %>%
  dplyr::mutate(method = ifelse(method == "distance.x", "cart_2d_post", "p_epi")) %>%
  dplyr::select(-c(method.x, method.y)) %>%
  dplyr::distinct()

year <- antigenic_distance %>%
  dplyr::select(vaccine_fullname, strains_fullname, subset, strain_type) %>%
  dplyr::distinct() %>%
  dplyr::mutate(vac_year = as.numeric(stringr::str_sub(vaccine_fullname, start = -4)),
         test_year = as.numeric(stringr::str_sub(strains_fullname, start = -4)),
         distance = abs(test_year - vac_year),
         method = "year") %>%
  dplyr::select(-c(vac_year, test_year)) %>%
  dplyr::distinct()

methods <- cent %>%
  dplyr::full_join(year,
                   by = c("subset", "strain_type", "strains_fullname", "vaccine_fullname", "method", "distance"))

virus <- readRDS(file = here::here("data/processed/virus_info.rds"))


subset_datasets %>%
  tidyr::pivot_longer(cols = dplyr::ends_with("vaccine_fullname"),
               names_to = "vaccine_strain_type", 
               values_to = "vaccine_fullname") %>%
  dplyr::mutate(vaccine_strain_type = ifelse(vaccine_strain_type == "h1n1_vaccine_fullname",
                                             "H1N1",
                                             ifelse(vaccine_strain_type == "h3n2_vaccine_fullname",
                                                    "H3N2",
                                                    "typeb"))) %>%
  dplyr::filter(vaccine_fullname != "None",
                strain_type == vaccine_strain_type)  %>%
  dplyr::select(-vaccine_strain_type) %>%
  dplyr::left_join(methods,
                   by = c("subset", "strain_type", "strains_fullname", "vaccine_fullname")) %>%
  dplyr::distinct() %>%
  dplyr::left_join(virus[, c('analysis_name', 'short_name')], 
            by = c('vaccine_fullname' = 'analysis_name')) %>%
  dplyr::left_join(virus[, c('analysis_name', 'short_name')], 
            by = c('strains_fullname' = 'analysis_name')) %>%
  dplyr::rename(vac_short = short_name.x,
         strain_short = short_name.y) %>%
  base::droplevels() %>%
  dplyr::ungroup() %>%
  saveRDS(., file = here::here("data/processed/subset_analysis/distance/distance_data.rds"))

