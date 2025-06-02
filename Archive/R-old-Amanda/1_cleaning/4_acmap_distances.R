## "Addition of Acmap Data to Antigenic Landscape"
## "Amanda Skarlupka"
## date: "12/29/2022"

# Purpose ####

##This files creates the final dataset that contains: antigenic distances from the acmap objects, and the sequence based distances from the amino acid sequences.

# Prerequisites:

# First run:
  
# ../code/1_cleaning/1_cleaning_code.R
# ../code/1_cleaning/2_pepitope_calculator.Rmd
# ../code/1_cleaning/3_acmap_creation.R

library(Racmacs) #Calculating antigenic distance and creating maps
library(tidyverse) #Working with data

#Functions for working with Acmap objects
source(here::here("R/functions/antigenic-cartography.R"))

#Read in the data
#There is minor cleaning the Pre viruses have duplicate names

d <- readRDS(file = here::here("data/processed/clean_data.rds"))

clean_data <- d # Keep a clean version for future reference

# H1N1 ####

subtype <- "h1"

distance <- readRDS(file = here::here(paste0("data/processed/", 
                                             subtype, 
                                             "_seq_distance/", 
                                             subtype, 
                                             "_key.rds"))) #H1N1 sequence distances

## Distances ####

#Load the different antigenic maps. These include pre and post vaccination and for 1, 2, and 3 dimensions. All maps (sd; as indicated) were created with only standard dose sera.


map_list <- vector(mode = "list")

### 1 Dimension Maps ####
map_list[["cart_1d_pre"]] <- AcmapToDF(map = read.acmap(here::here(paste0("data/processed/", 
                                                                          subtype, 
                                                                          "_maps/", 
                                                                          subtype, 
                                                                          "_pre_sd_1d.ace")))) %>%
  dplyr::filter(type == "antigen") %>%
  dplyr::mutate(method = "cart_1d_pre",
         V2 = 0,
         V3 = 0) %>%
  dplyr::select(V1, V2, V3, identifier, method)

map_list[["cart_1d_post"]] <- AcmapToDF(map = read.acmap(here::here(paste0("data/processed/", 
                                                                           subtype, 
                                                                           "_maps/", 
                                                                           subtype, 
                                                                           "_post_sd_1d.ace")))) %>%
  dplyr::filter(type == "antigen") %>%
  dplyr::mutate(method = "cart_1d_post",
         V2 = 0, 
         V3 = 0) %>%
  dplyr::select(V1, V2, V3, identifier, method)

### 2 Dimension Maps ####
map_list[["cart_2d_pre"]] <- AcmapToDF(map = read.acmap(here::here(paste0("data/processed/", 
                                                                          subtype, 
                                                                          "_maps/", 
                                                                          subtype, 
                                                                          "_pre_sd_2d.ace")))) %>%
  dplyr::filter(type == "antigen") %>%
  dplyr::mutate(method = "cart_2d_pre",
         V3 = 0) %>%
  dplyr::select(V1, V2, V3, identifier, method)

map_list[["cart_2d_post"]] <- AcmapToDF(map = read.acmap(here::here(paste0("data/processed/", 
                                                                           subtype, 
                                                                           "_maps/", 
                                                                           subtype, 
                                                                           "_post_sd_2d.ace")))) %>%
  dplyr::filter(type == "antigen") %>%
  dplyr::mutate(method = "cart_2d_post",
         V3 = 0) %>%
  dplyr::select(V1, V2, V3, identifier, method)

### 3 Dimension Maps ####
map_list[["cart_3d_pre"]] <- AcmapToDF(map = read.acmap(here::here(paste0("data/processed/", 
                                                                          subtype, 
                                                                          "_maps/", 
                                                                          subtype, 
                                                                          "_pre_sd_3d.ace")))) %>%
  dplyr::filter(type == "antigen") %>%
  dplyr::mutate(method = "cart_3d_pre") %>%
  dplyr::select(V1, V2, V3, identifier, method)

map_list[["cart_3d_post"]] <- AcmapToDF(map = read.acmap(here::here(paste0("data/processed/", 
                                                                           subtype, 
                                                                           "_maps/", 
                                                                           subtype, 
                                                                           "_post_sd_3d.ace")))) %>%
  dplyr::filter(type == "antigen") %>%
  dplyr::mutate(method = "cart_3d_post") %>%
  dplyr::select(V1, V2, V3, identifier, method)

antigen <- data.table::rbindlist(map_list)

#Warnings are generated here, they are suppressed
suppressWarnings({
abland <- d %>%
  dplyr::filter(strain_type == "H1N1", 
         !is.na(titerincrease), 
         !is.na(postvactiter)) %>%
  dplyr::select(h1n1_vaccine_fullname, strains_fullname) %>%
  dplyr::right_join(antigen, 
             by = c("strains_fullname" = "identifier")) %>%
  dplyr::filter(h1n1_vaccine_fullname != "H1N1-Guangdong Maonan-2019",
         method != "cart_3d_pre") %>%
  dplyr::group_by(h1n1_vaccine_fullname, 
           method) %>%
  dplyr::mutate(distance = sqrt((V1-V1[which(strains_fullname==paste(h1n1_vaccine_fullname))])^2+(V2-V2[which(strains_fullname==paste(h1n1_vaccine_fullname))])^2+(V3-V3[which(strains_fullname==paste(h1n1_vaccine_fullname))])^2)) %>% #Euclidean distance from vaccine strain
  dplyr::ungroup() %>%
  dplyr::select(-c(V1:V3))
})

cent <- distance %>%
  dplyr::filter(analysis_name %in% levels(abland$h1n1_vaccine_fullname)) %>%
  tidyr::pivot_longer(cols = starts_with("H1N1"), 
               names_to = "strains_fullname", 
               values_to = "distance") %>%
  dplyr::mutate(h1n1_vaccine_fullname = analysis_name) %>%
  dplyr::select(-analysis_name) %>%
  dplyr::full_join(abland,
                   by = c("method", "strains_fullname", "distance", "h1n1_vaccine_fullname")) %>%
  dplyr::select(h1n1_vaccine_fullname, strains_fullname, distance, method) %>%
  dplyr::distinct()

year <- abland %>%
  dplyr::select(h1n1_vaccine_fullname, strains_fullname) %>%
  dplyr::distinct() %>%
  dplyr::mutate(vac_year = as.numeric(stringr::str_sub(h1n1_vaccine_fullname, start = -4)),
         test_year = as.numeric(stringr::str_sub(strains_fullname, start = -4)),
         distance = abs(test_year - vac_year),
         method = "year") %>%
  dplyr::select(-c(vac_year, test_year)) %>%
  dplyr::distinct()

h1n1_methods <- cent %>%
  dplyr::full_join(year,
                   by = c("h1n1_vaccine_fullname", "strains_fullname", "distance", "method"))
  

# H3N2 ####

subtype <- "h3"

d <- clean_data %>%
  dplyr::filter(strain_type == "H3N2", !is.na(titerincrease), !is.na(postvactiter)) %>%
  dplyr::select(h3n2_vaccine_fullname, strains_fullname)

distance <- readRDS(file = here::here(paste0("data/processed/", 
                                             subtype, 
                                             "_seq_distance/", 
                                             subtype, 
                                             "_key.rds"))) #H3N2 sequence distances

map_list <- vector(mode = "list")


## Distances ####

### 1 Dimension Maps ####
map_list[["cart_1d_pre"]] <- AcmapToDF(map = read.acmap(here::here(paste0("data/processed/", 
                                                                          subtype, 
                                                                          "_maps/", 
                                                                          subtype, 
                                                                          "_pre_sd_1d.ace")))) %>%
  dplyr::filter(type == "antigen") %>%
  dplyr::mutate(method = "cart_1d_pre",
         V2 = 0,
         V3 = 0) %>%
  dplyr::select(V1, V2, V3, identifier, method)

map_list[["cart_1d_post"]] <- AcmapToDF(map = read.acmap(here::here(paste0("data/processed/", 
                                                                           subtype, 
                                                                           "_maps/", 
                                                                           subtype, 
                                                                           "_post_sd_1d.ace")))) %>%
  dplyr::filter(type == "antigen") %>%
  dplyr::mutate(method = "cart_1d_post",
         V2 = 0, 
         V3 = 0) %>%
  dplyr::select(V1, V2, V3, identifier, method)

### 2 Dimension Maps ####
map_list[["cart_2d_pre"]] <- AcmapToDF(map = read.acmap(here::here(paste0("data/processed/", 
                                                                          subtype, 
                                                                          "_maps/", 
                                                                          subtype, 
                                                                          "_pre_sd_2d.ace")))) %>%
  dplyr::filter(type == "antigen") %>%
  dplyr::mutate(method = "cart_2d_pre",
         V3 = 0) %>%
  dplyr::select(V1, V2, V3, identifier, method)

map_list[["cart_2d_post"]] <- AcmapToDF(map = read.acmap(here::here(paste0("data/processed/", 
                                                                           subtype, 
                                                                           "_maps/", 
                                                                           subtype, 
                                                                           "_post_sd_2d.ace")))) %>%
  dplyr::filter(type == "antigen") %>%
  dplyr::mutate(method = "cart_2d_post",
         V3 = 0) %>%
  dplyr::select(V1, V2, V3, identifier, method)

### 3 Dimension Maps ####
map_list[["cart_3d_pre"]] <- AcmapToDF(map = read.acmap(here::here(paste0("data/processed/", 
                                                                          subtype, 
                                                                          "_maps/", 
                                                                          subtype, 
                                                                          "_pre_sd_3d.ace")))) %>%
  dplyr::filter(type == "antigen") %>%
  dplyr::mutate(method = "cart_3d_pre") %>%
  dplyr::select(V1, V2, V3, identifier, method)

map_list[["cart_3d_post"]] <- AcmapToDF(map = read.acmap(here::here(paste0("data/processed/", 
                                                                           subtype, 
                                                                           "_maps/", 
                                                                           subtype, 
                                                                           "_post_sd_3d.ace")))) %>%
  dplyr::filter(type == "antigen") %>%
  dplyr::mutate(method = "cart_3d_post") %>%
  dplyr::select(V1, V2, V3, identifier, method)

antigen <- data.table::rbindlist(map_list)

#Suppress the warnings:
suppressWarnings({
abland <- d %>%
  dplyr::right_join(antigen, 
             by = c("strains_fullname" = "identifier")) %>%
  dplyr::group_by(h3n2_vaccine_fullname, 
           method) %>%
  dplyr::mutate(distance = sqrt((V1-V1[which(strains_fullname==paste(h3n2_vaccine_fullname))])^2 +
                           (V2-V2[which(strains_fullname==paste(h3n2_vaccine_fullname))])^2 + 
                           (V3-V3[which(strains_fullname==paste(h3n2_vaccine_fullname))])^2)) %>%
  dplyr::ungroup() %>%
  dplyr::select(-c(V1:V3))
})
cent <- distance %>%
  dplyr::filter(analysis_name %in% levels(abland$h3n2_vaccine_fullname)) %>%
  tidyr::pivot_longer(cols = starts_with("H3N2"), 
               names_to = "strains_fullname", 
               values_to = "distance") %>%
  dplyr::mutate(h3n2_vaccine_fullname = analysis_name) %>%
  dplyr::select(-analysis_name) %>%
  dplyr::full_join(abland,
                   by = c("method", "strains_fullname", "distance", "h3n2_vaccine_fullname")) %>%
  dplyr::select(h3n2_vaccine_fullname, strains_fullname, distance, method) %>%
  dplyr::distinct()

year <- abland %>%
  dplyr::select(h3n2_vaccine_fullname, strains_fullname) %>%
  dplyr::distinct() %>%
  dplyr::mutate(vac_year = as.numeric(stringr::str_sub(h3n2_vaccine_fullname, start = -4)),
         test_year = as.numeric(stringr::str_sub(strains_fullname, start = -4)),
         distance = abs(test_year - vac_year),
         method = "year") %>%
  dplyr::select(-c(vac_year, test_year)) %>%
  dplyr::distinct()


h3n2_methods <- cent %>%
  dplyr::full_join(year,
                   by = c("h3n2_vaccine_fullname", "strains_fullname", "distance", "method"))


# Saving Files ####

saveRDS(h1n1_methods,
        file = here::here("data/processed/h1_seq_distance/h1_distance_measure_key.rds"))

saveRDS(h3n2_methods,
        file = here::here("data/processed/h3_seq_distance/h3_distance_measure_key.rds"))


virus <- readRDS(file = here::here("data/processed/virus_info.rds"))


colnames(h3n2_methods)[colnames(h3n2_methods) == "h3n2_vaccine_fullname"] <- "vaccine_fullname"
colnames(h1n1_methods)[colnames(h1n1_methods) == "h1n1_vaccine_fullname"] <- 'vaccine_fullname'

all_methods <- base::rbind(h3n2_methods, 
                     h1n1_methods)

a <- clean_data %>%
  tidyr::pivot_longer(cols = ends_with("vaccine_fullname"),
               names_to = "vaccine_type", 
               values_to = "vaccine_fullname") %>%
  dplyr::filter(vaccine_fullname != "None") %>%
  dplyr::mutate(working_type = ifelse(strain_type %in% c("B-Pre", "B-Vic", "B-Yam"), 
                               "B-", 
                               ifelse(strain_type == "H1N1", 
                                      "H1", 
                                      "H3"))) %>%
  dplyr::left_join(all_methods,
                   by = c("strains_fullname", "vaccine_fullname")) %>%
  dplyr::mutate(working_type = as.character(working_type),
         vaccine_fullname = as.character(vaccine_fullname))  %>%
  dplyr::filter(str_detect(str_sub(vaccine_fullname, 
                            start = 1, 
                            end = 2), 
                    working_type)) %>%
  dplyr::select(-working_type) %>%
  dplyr::distinct() %>%
  dplyr::left_join(virus[, c('analysis_name', 'short_name')], 
            by = c('vaccine_fullname' = 'analysis_name')) %>%
  dplyr::left_join(virus[, c('analysis_name', 'short_name')], 
            by = c('strains_fullname' = 'analysis_name')) %>%
  dplyr::rename(vac_short = short_name.x,
         strain_short = short_name.y) %>%
  base::droplevels() %>%
  dplyr::ungroup()

saveRDS(a, file = here::here("data/processed/distance_data.rds"))

## Save a list of the viruses that were included in both the original data and also have associated distances:


a %>%
  dplyr::select(c(strain_type, season, titerincrease, method, distance, strain_short)) %>%
  dplyr::filter(method %in% c("cart_2d_post", "p_epi", "year"), 
         !is.na(distance),
         strain_type %in% c("H1N1", "H3N2"),
         season != 2020,
         !is.na(titerincrease)) %>%
  dplyr::select(strain_short) %>%
  base::droplevels() %>%
  base::unique() %>%
  saveRDS(file = here::here("data/processed/viruses_used.rds"))



  