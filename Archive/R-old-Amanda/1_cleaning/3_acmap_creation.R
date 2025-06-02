
# Load the packages and functions:

library(Racmacs) #Calculating antigenic distance and creating maps
library(tidyverse) #Working with data
source(here::here("R/functions/antigenic-cartography.R"))


# Create the 10 SD, post-titer Acmaps for the H1N1 and H3N2 subsets 
subset_Acmaps()

# Read and separate the the data
d <- readRDS(file = here::here("data/processed/clean_data.rds"))

d_sd <- d %>%
  dplyr::filter(dose == "SD")

# H1N1 ####
subtype <- "h1"

## Pretiter ####
### Remove Underconstrained ####
#### Entire population ####

dimension_set <- vector(mode = "list")

for (i in 1:5) {
  dimension_set[[i]] <- remove_underconstrained(d, 
                                                subtype = "H1N1", 
                                                time = "pretiter", 
                                                dimensions = i)
}

#### SD Only ####
#For only the individuals who received a standard dose

dimension_set_sd <- vector(mode = "list")

for (i in 1:5) {
  dimension_set_sd[[i]] <- remove_underconstrained(d_sd, 
                                                   subtype = "H1N1", 
                                                   time = "pretiter", 
                                                   dimensions = i)
}

### Create maps and save files ####

#This file is for running the dimension selection on. 
  
  saveRDS(dimension_set[[1]], 
          file = here::here(paste0("data/cluster/input/", 
                                   subtype, 
                                   "_pre.rds")))
  
  saveRDS(dimension_set_sd[[1]], 
          file = here::here(paste0("data/cluster/input/",
                                   subtype, 
                                   "_pre_sd.rds")))
  for (i in 1:3) {
    CreateSaveAcmap(dimension_set[[i]], 
                    dimensions = i, 
                    optimizations = 100, 
                    subtype = subtype, 
                    timing = "pre", 
                    specific_filtering = "all")
  }
  
  for (i in 1:3) {
    CreateSaveAcmap(dimension_set_sd[[i]], 
                    dimensions = i, 
                    optimizations = 100, 
                    subtype = subtype, 
                    timing = "pre", 
                    specific_filtering = "sd")
  }
  
## PostTiter ####
  
  dimension_set <- vector(mode = "list")
  
### Remove Underconstrained ####
#### Entire population ####
  #For the entire population of sera
  for (i in 1:5) {
    dimension_set[[i]] <- remove_underconstrained(d, 
                                                  subtype = "H1N1", 
                                                  time = "postiter", 
                                                  dimensions = i)
  }
  
#### SD Only ####
  #For a only the individuals who received a standard dose:
  
  dimension_set_sd <- vector(mode = "list")
  
  for (i in 1:5) {
    dimension_set_sd[[i]] <- remove_underconstrained(d_sd, 
                                                     subtype = "H1N1", 
                                                     time = "postiter", 
                                                     dimensions = i)
  }
  
### Create Maps, Save Files  ####
#This file is for running the dimension selection on. 
    
    saveRDS(dimension_set[[1]], 
            file = here::here(paste0("data/cluster/input/",
                                     subtype, 
                                     "_post.rds")))
    
    saveRDS(dimension_set_sd[[1]], 
            file = here::here(paste0("data/cluster/input/",
                                     subtype, 
                                     "_post_sd.rds")))
    
    for (i in 1:3) {
      CreateSaveAcmap(dimension_set[[i]], 
                      dimensions = i, 
                      optimizations = 100, 
                      subtype = subtype, 
                      timing = "post", 
                      specific_filtering = "all")
    }
    
    for (i in 1:3) {
      CreateSaveAcmap(dimension_set_sd[[i]], 
                      dimensions = i, 
                      optimizations = 100, 
                      subtype = subtype, 
                      timing = "post", 
                      specific_filtering = "sd")
    }
    
# H3N2 ####
    subtype <- "h3"
## Pretiter ####
### Remove Underconstrained ####
#### Entire population ####
# For the entire population of sera
    
    dimension_set <- vector(mode = "list")
    
    for (i in 1:5) {
      dimension_set[[i]] <- remove_underconstrained(d, 
                                                    subtype = "H3N2", 
                                                    time = "pretiter", 
                                                    dimensions = i)
    }
#### SD Only ####    
# For only the individuals who received a standard dose:
    
    dimension_set_sd <- vector(mode = "list")
    
    for (i in 1:5) {
      dimension_set_sd[[i]] <- remove_underconstrained(d_sd, 
                                                       subtype = "H3N2", 
                                                       time = "pretiter", 
                                                       dimensions = i)
    }
    
### Create Maps, Save Files ####
    #This file is for running the dimension selection on. 
    
    saveRDS(dimension_set[[1]], 
            file = here::here(paste0("data/cluster/input/",
                                     subtype, 
                                     "_pre.rds")))
    
    saveRDS(dimension_set_sd[[1]], 
            file = here::here(paste0("data/cluster/input/",
                                     subtype, 
                                     "_pre_sd.rds")))
    for (i in 1:3) {
      CreateSaveAcmap(dimension_set[[i]], 
                      dimensions = i, 
                      optimizations = 100, 
                      subtype = subtype, 
                      timing = "pre", 
                      specific_filtering = "all")
    }
    
    for (i in 1:3) {
      CreateSaveAcmap(dimension_set_sd[[i]], 
                      dimensions = i, 
                      optimizations = 100, 
                      subtype = subtype, 
                      timing = "pre", 
                      specific_filtering = "sd")
    }
    
## Post Titer: ####
    
    dimension_set <- vector(mode = "list")
### Remove Underconstrained ####
#### Entire population ####
#For the entire population of sera
    for (i in 1:5) {
      dimension_set[[i]] <- remove_underconstrained(d, 
                                                    subtype = "H3N2", 
                                                    time = "postiter", 
                                                    dimensions = i)
    }
#### SD Only ####
    #For a only the individuals who received a standard dose:
    
    dimension_set_sd <- vector(mode = "list")
    
    for (i in 1:5) {
      dimension_set_sd[[i]] <- remove_underconstrained(d_sd, 
                                                       subtype = "H3N2", 
                                                       time = "postiter", 
                                                       dimensions = i)
    }
    
### Create Maps, Save Files ####   
    #This file is for running the dimension selection on. 
    
    saveRDS(dimension_set[[1]], 
            file = here::here(paste0("data/cluster/input/",
                                     subtype, 
                                     "_post.rds")))
    
    saveRDS(dimension_set_sd[[1]], 
            file = here::here(paste0("data/cluster/input/",
                                     subtype, 
                                     "_post_sd.rds")))
    
    for (i in 1:3) {
      CreateSaveAcmap(dimension_set[[i]], 
                      dimensions = i, 
                      optimizations = 100, 
                      subtype = subtype, 
                      timing = "post", 
                      specific_filtering = "all")
    }
    
    for (i in 1:3) {
      CreateSaveAcmap(dimension_set_sd[[i]], 
                      dimensions = i, 
                      optimizations = 100, 
                      subtype = subtype, 
                      timing = "post", 
                      specific_filtering = "sd")
    }
    
    
# SD maps specific for each season #######
    
 #Read in the data
    
    d <- readRDS(file = here::here("data/processed/clean_data.rds"))

    # SD individuals by Season
    
    df_split <- vector(mode = "list")
    df_split <- d %>%
      dplyr::filter(dose == "SD",
             !is.na(titerincrease)) %>%
      dplyr::group_split(season, 
                  .keep = TRUE) %>%
      stats::setNames(unique(d$season))
    
## H1N1 ####
### Pre-titer: ####
    
    #Select the H1N1 only viruses and sera, run map and then deselect the underconstrained sera, and rerun map. For the 2 dimensional maps only. Cycle through seasons. 
    
    
    h1n1_vaccines <- d %>%
      dplyr::select(season, h1n1_vaccine_fullname) %>%
      dplyr::distinct()
    
    reduced <- vector(mode = "list")
    
    #For the SD receivers by season
    for (i in 1:length(df_split)) {
      x <- names(df_split[i])
      reduced[[x]] <- remove_underconstrained(df_split[[i]], 
                                              subtype = "H1N1", 
                                              time = "pretiter", 
                                              dimensions = 2)
    }
    
    
    map_list <- vector(mode = "list")
    
      
      for (i in 1:(length(reduced)-1)) { #Season 2020 has no data points for creating a map
        x <- paste0("pre_",
                    names(reduced[i]))
        map_list[[x]] <- CreateSaveAcmap(reduced[[i]], 
                                         dimensions = 2, 
                                         optimizations = 100, 
                                         subtype = subtype, 
                                         timing = "pre", 
                                         specific_filtering = "sd_season",
                                         save = FALSE)
      }
      
    
    
### Post-Titer: ####
    
    reduced <- vector(mode = "list")
    
    #For the SD receivers by season
    for (i in 1:length(df_split)) {
      x <- names(df_split[i])
      reduced[[x]] <- remove_underconstrained(df_split[[i]], 
                                              subtype = "H1N1", 
                                              time = "postiter", 
                                              dimensions = 2)
    }

      
      for (i in 1:(length(reduced)-1)) { #Season 2020 has no data points for creating a map
        x <- paste0("post_",
                    names(reduced[i]))
        map_list[[x]] <- CreateSaveAcmap(reduced[[i]], 
                                         dimensions = 2, 
                                         optimizations = 100, 
                                         subtype = subtype, 
                                         timing = "post", 
                                         specific_filtering = "sd_season",
                                         save = FALSE)
      }
      
    
    
## Save Data ####
      
      maps <- vector(mode = "list")
      
      for (i in 1:length(map_list)) {
        
        maps[[i]] <- AcmapToDF(map = map_list[[i]]) %>%
          dplyr::mutate(sera_count = sum(type == "sera")) %>%
          dplyr::filter(type == "antigen") %>%
          dplyr::mutate(method = paste0("cart_2d", ifelse(stringr::str_detect(names(map_list[i]), 
                                                              pattern = "^pre"),
                                                   "_pre",
                                                   "_post")),
                 season = stringr::str_sub(names(map_list[i]), start = -4)) %>%
          dplyr::select(V1, V2, identifier, method, season, sera_count)
        
      }
      
      antigen <- data.table::rbindlist(maps)
      df <- antigen %>%
        dplyr::mutate(season = as.numeric(season)) %>%
        dplyr::left_join(h1n1_vaccines,
                         by = "season")
      
      
      saveRDS(object = df,
              file = here::here("data/processed/h1_seq_distance/h1_seasons_distance.rds"))
    
    
    
    rm(antigen, df, h1n1_vaccines, map_list, maps, post_maps, pre_maps, reduced)
    
    
## H3N2 ####
    
### Pre-titer: ####
    
    #Select the H1N1 only viruses and sera, run map and then deselect the underconstrained sera, and rerun map. For the 2 dimensional maps only. Cycle through seasons. 
    
    subtype <- "h3"
    h3n2_vaccines <- d %>%
      dplyr::select(season, h3n2_vaccine_fullname) %>%
      dplyr::distinct()
    
    reduced <- vector(mode = "list")
    
    #For the SD receivers by season
    for (i in 1:length(df_split)) {
      x <- names(df_split[i])
      reduced[[x]] <- remove_underconstrained(df_split[[i]], 
                                              subtype = "H3N2", 
                                              time = "pretiter", 
                                              dimensions = 2)
    }
    
    
    map_list <- vector(mode = "list")
    
      
      for (i in 1:(length(reduced)-1)) { #Season 2020 has no data points for creating a map
        x <- paste0("pre_",
                    names(reduced[i]))
        map_list[[x]] <- CreateSaveAcmap(reduced[[i]], 
                                         dimensions = 2, 
                                         optimizations = 100, 
                                         subtype = subtype, 
                                         timing = "pre", 
                                         specific_filtering = "sd_season",
                                         save = FALSE)
      }
      
    
    
    
### Post-Titer: ####
    
    
    reduced <- vector(mode = "list")
    
    #For the SD receivers by season
    for (i in 1:length(df_split)) {
      x <- names(df_split[i])
      reduced[[x]] <- remove_underconstrained(df_split[[i]], 
                                              subtype = "H3N2", 
                                              time = "postiter", 
                                              dimensions = 2)
    }
    
      
      for (i in 1:(length(reduced)-1)) { #Season 2020 has no data points for creating a map
        x <- paste0("post_",
                    names(reduced[i]))
        map_list[[x]] <- CreateSaveAcmap(reduced[[i]], 
                                         dimensions = 2, 
                                         optimizations = 100, 
                                         subtype = subtype, 
                                         timing = "post", 
                                         specific_filtering = "sd_season",
                                         save = FALSE)
      }
      
    
    
    
    
## Save Data: ####
    
    maps <- vector(mode = "list")
    
    for (i in 1:length(map_list)) {
      
      maps[[i]] <- AcmapToDF(map = map_list[[i]]) %>%
        dplyr::mutate(sera_count = sum(type == "sera")) %>%
        dplyr::filter(type == "antigen") %>%
        dplyr::mutate(method = paste0("cart_2d", ifelse(stringr::str_detect(names(map_list[i]), 
                                                            pattern = "^pre"),
                                                 "_pre",
                                                 "_post")),
               season = stringr::str_sub(names(map_list[i]), start = -4)) %>%
        dplyr::select(V1, V2, identifier, method, season, sera_count)
      
    }
    
    antigen <- data.table::rbindlist(maps)
    df <- antigen %>%
      dplyr::mutate(season = as.numeric(season)) %>%
      dplyr::left_join(h3n2_vaccines,
                       by = "season")
    
      saveRDS(object = df,
              file = here::here("data/processed/h3_seq_distance/h3_seasons_distance.rds"))
    
    
    
    
    