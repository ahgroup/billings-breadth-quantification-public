## Data cleaning file ############

## Purpose: clean the cohort data and viral data. Saves the supplemental table for viral accession numbers 

library(tidyverse)
library(here)
library(readxl)
library(dplyr)

source(file = here::here("R/functions/1_cleaning.R"))

# Cohort data cleaning: 
# Load the original cohort data file:
# Save the cleaned data file in the processed folder:
individuals <- cleaning_all()


## Save supplemental table with accession numbers and the seasons it was used:

try <- readRDS(file = here::here("Data/processed/clean_data.rds")) %>%
  dplyr::select(season, strain_type, strains_fullname) %>% 
  dplyr::distinct() %>%
  dplyr::mutate(contain = "x")

virus_clean() %>% #Recall the original virus data to clean up the accession numbers
  dplyr::inner_join(try, by = c("analysis_name" = "strains_fullname")) %>%
  tidyr::pivot_wider(names_from = season,
                     values_from = contain) %>%
  dplyr::select(-strain_type) %>%
  base::saveRDS(file = here::here("Results/Tables/virus_accession.rds"))

## These are all of the viruses used in the 2014-2019 HAI breadth panels including the vaccine strains
virus_data <- try %>%
  dplyr::select(-c(contain, season, strain_type)) %>%
  dplyr::distinct()

# Virus data cleaning:
virus_data <- readxl::read_excel(here::here("Data/raw/accession_no.xlsx")) %>%
  dplyr::filter(analysis_name %in% virus_data$strains_fullname) # Keep only the strains used in the 2014-2019 studies

season_data <- readRDS(here::here("Data/raw/vaccines.rds")) %>%
  tidyr::pivot_longer(cols = h1n1_vaccine:h3n2_vaccine,
                      names_to = "vaccine_type",
                      values_to = "analysis_name") %>%
  dplyr::select(season, analysis_name) 

## Set the factor levels for the short names
virus <- virus_data %>%
  dplyr::select(c(subtype, full_strain_name, analysis_name, short_name, factor_order)) %>%
  dplyr::left_join(y = season_data,
                  by = "analysis_name") %>%
  dplyr::select(-season) %>%
  dplyr::distinct() %>%
  dplyr::arrange(factor_order) %>%
  dplyr::mutate(short_name = factor(factor_order, 
                             labels = short_name, 
                             levels = factor_order))

base::saveRDS(virus, file = here::here('Data/processed/virus_info.rds'))


  

## Create the subsetted datasets; Total number is based on the arg number_of_subsets

subset_data <- subset_the_data(dataframe = virus,
                number_of_subsets = 5,
                number_of_viruses = 10,
                starting_seed = 123,
                individuals_df = individuals)




## Save supplemental table with the strains used for each subset:
try <- try %>%
  dplyr::select(strain_type, strains_fullname, contain) %>%
  dplyr::distinct()

subset_data %>%
  dplyr::mutate(contain = "x") %>%
  tidyr::pivot_wider(names_from = subset,
                     values_from = contain,
                     names_prefix = "subset_") %>%
  dplyr::select(-subtype) %>%
  dplyr::right_join(try, by = c("analysis_name" = "strains_fullname")) %>%
  dplyr::select(-factor_order) %>%
  dplyr::left_join(virus[,c("analysis_name", "factor_order")], by = "analysis_name") %>%
  dplyr::select(strain_type, analysis_name, factor_order, contain, rev(starts_with("subset"))) %>%
  base::saveRDS(file = here::here("tables/supplemental/data_subsets_viruses.rds"))



