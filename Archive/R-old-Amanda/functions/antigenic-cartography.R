########################################################################
#
#     epi()  ############
#
########################################################################

#Description: Selects the characters of a specific antigenic site from an amino acid string. Convert the string into a vector of characters, and select the characters that correspond to the given antigenic site. Re-condense the vector of characters into a string, and return that string.

epi <- function(sequence, 
                antigenic_site)
{
  # Arguments 
  # "sequence" = a character string of the amino acid sequence NOT starting from the Methionine start, but instead after the leader sequence has been removed. Please refer to Antigenic Evergreen residues for assistance.  
  
  # "antigenic_site" = a vector of integers indicating the amino acid residues of the antigenic site of interest
  
  vec <- s2c(sequence)
  site <- vec[antigenic_site]
  epit <- c2s(site)
  
  return(epit)
  
  #Value: Returns the antigenic site as a string
}

########################################################################
#
#     remove_underconstrained() #####
#
########################################################################

# Function for removing underconstrained values for antigenic cartography:



remove_underconstrained <- function (dataframe, 
                                     subtype, 
                                     time, 
                                     dimensions) 
{
  d <- dataframe
  
  #Arguments:
  
  # subtype = The influenza subtype for removing. Options include: "all", "H3N2", "H1N1", "typeb", "bvic", and "byam". The pre type B's are included in all of the groupings. B-vic = B-pre + B-vic; typeb = B-pre+B-vic+B-yam.
  
  if (subtype == "H3N2") {
    
    if (time == "pretiter") {
      pivoted <- d %>%
        filter(strain_type == "H3N2") %>%
        select(uniq_id, h3n2_vaccine_fullname, pretiter, strains_fullname) %>%
        mutate(sera = paste(h3n2_vaccine_fullname, uniq_id))  %>%
        select(-uniq_id, -h3n2_vaccine_fullname) %>%
        pivot_wider(names_from = sera, values_from = pretiter)
    }
    
    if (time == "postiter") {
      pivoted <- d %>%
        filter(strain_type == "H3N2") %>%
        select(uniq_id, h3n2_vaccine_fullname, postiter, strains_fullname) %>%
        mutate(sera = paste(h3n2_vaccine_fullname, uniq_id))  %>%
        select(-uniq_id, -h3n2_vaccine_fullname) %>%
        pivot_wider(names_from = sera, values_from = postiter)
    }
  }
  
  if (subtype == "H1N1") {
    
    if (time == "pretiter") {
      pivoted <- d %>%
        filter(strain_type == "H1N1") %>%
        select(uniq_id, h1n1_vaccine_fullname, pretiter, strains_fullname) %>%
        mutate(sera = paste(h1n1_vaccine_fullname, uniq_id))  %>%
        select(-uniq_id, -h1n1_vaccine_fullname) %>%
        pivot_wider(names_from = sera, values_from = pretiter)
    }
    
    if (time == "postiter") {
      pivoted <- d %>%
        filter(strain_type == "H1N1") %>%
        select(uniq_id, h1n1_vaccine_fullname, postiter, strains_fullname) %>%
        mutate(sera = paste(h1n1_vaccine_fullname, uniq_id))  %>%
        select(-uniq_id, -h1n1_vaccine_fullname) %>%
        pivot_wider(names_from = sera, values_from = postiter)
    }
  }
    
  if (subtype == "typeb") {
      
    if (time == "pretiter") {
      pivoted <- d %>%
        filter(strain_type %in% c("B-Pre", "B-Vic", "B-Yam")) %>%
        select(uniq_id, bvictoria_vaccine_fullname, yamagata_vaccine_fullname, pretiter, strains_fullname) %>%
        mutate(sera = paste(bvictoria_vaccine_fullname, yamagata_vaccine_fullname, uniq_id))  %>%
        select(-uniq_id, -bvictoria_vaccine_fullname, -yamagata_vaccine_fullname) %>%
        pivot_wider(names_from = sera, values_from = pretiter)
    }
    
    if (time == "postiter") {
      pivoted <- d %>%
        filter(strain_type %in% c("B-Pre", "B-Vic", "B-Yam")) %>%
        select(uniq_id, bvictoria_vaccine_fullname, yamagata_vaccine_fullname, postiter, strains_fullname) %>%
        mutate(sera = paste(bvictoria_vaccine_fullname, yamagata_vaccine_fullname, uniq_id))  %>%
        select(-uniq_id, -bvictoria_vaccine_fullname, -yamagata_vaccine_fullname) %>%
        pivot_wider(names_from = sera, values_from = postiter)
    }
  }
  
  if (subtype == "bvic") {
    
    if (time == "pretiter") {
      pivoted <- d %>%
        filter(strain_type %in% c("B-Pre", "B-Vic")) %>%
        select(uniq_id, bvictoria_vaccine_fullname, pretiter, strains_fullname) %>%
        mutate(sera = paste(bvictoria_vaccine_fullname, uniq_id))  %>%
        select(-c(uniq_id, bvictoria_vaccine_fullname)) %>%
        pivot_wider(names_from = sera, values_from = pretiter)
    }
    
    if (time == "postiter") {
      pivoted <- d %>%
        filter(strain_type %in% c("B-Pre", "B-Vic")) %>%
        select(uniq_id, bvictoria_vaccine_fullname, postiter, strains_fullname) %>%
        mutate(sera = paste(bvictoria_vaccine_fullname, uniq_id))  %>%
        select(-c(uniq_id, bvictoria_vaccine_fullname)) %>%
        pivot_wider(names_from = sera, values_from = postiter)
    }
  }
  
  if (subtype == "byam") {
    
    if (time == "pretiter") {
      pivoted <- d %>%
        filter(strain_type %in% c("B-Pre", "B-Yam")) %>%
        select(uniq_id, yamagata_vaccine_fullname, pretiter, strains_fullname) %>%
        mutate(sera = paste(yamagata_vaccine_fullname, uniq_id))  %>%
        select(-c(uniq_id, yamagata_vaccine_fullname)) %>%
        pivot_wider(names_from = sera, values_from = pretiter)
    } 
    
    if (time == "postiter") {
      pivoted <- d %>%
        filter(strain_type %in% c("B-Pre", "B-Yam")) %>%
        select(uniq_id, yamagata_vaccine_fullname, postiter, strains_fullname) %>%
        mutate(sera = paste(yamagata_vaccine_fullname, uniq_id))  %>%
        select(-c(uniq_id, yamagata_vaccine_fullname)) %>%
        pivot_wider(names_from = sera, values_from = postiter)
    }
  }
  
  if (subtype == "all") {
    
    if (time == "pretiter") {
      
      pivoted <- d %>%
        select(uniq_id, pretiter, strains_fullname)%>%
        mutate(sera = uniq_id)  %>%
        select(-c(uniq_id)) %>%
        pivot_wider(names_from = sera, values_from = pretiter)
    } 
    
    if (time == "postiter") {
      
      pivoted <- d %>%
        select(uniq_id, postiter, strains_fullname) %>%
        mutate(sera = uniq_id)  %>%
        select(-c(uniq_id)) %>%
        pivot_wider(names_from = sera, values_from = postiter)
    }
  }
    
    pivoted[is.na(pivoted)] <- 1
    pivoted_2 <- pivoted %>%
      mutate(across(.cols = where(is.numeric), ~ ifelse(. == "NULL", 1, .)))
    a <- pivoted_2
    
    
    #The Racmacs mapping program it's currently set up so that strains with < n + 1 measurable titers get this warning, where n is the number of dimensions. 
    
    while (sum(rowSums(a[,-1] >= 10) < (dimensions + 1)) > 0 || sum(colSums(a[,-1] >= 10) < dimensions + 1) > 0) {
      
      row_select <- rowSums(a[,-1] >= 10) > (dimensions + 1)
      col_select <- c(TRUE, colSums(a[,-1] >= 10) > (dimensions + 1))
      
      a <- a[row_select, col_select]
    }
    final <- a %>%
      mutate(across(.cols = everything(), as.character),
             across(.cols = everything(), ~ifelse(. == "1", "*", .)),
             across(.cols = everything(), ~ifelse(. == "5", "<10", .)))
    
    return(final)
  
}

# #######################################################################
# AcmapToDF Function ########

# Function for changing the Acmap objects into a dataframe for plotting. There is an option for merging with a key

## Variables for troubleshooting:

#map <- read.acmap(here::here(paste0("data/processed/", 
#                             "h1", 
#                             "_maps/", 
#                             "h1", 
#                             "_pre_sd_1d.ace")))

AcmapToDF <- function(map, MergeWithKey = FALSE, UniqIDKey = FALSE) {
  antigen <- as.data.frame(agCoords(map))
  antigen$identifier <- agNames(map)
  antigen$type <- "antigen"
  sera <- as.data.frame(srCoords(map))
  sera$identifier <- srNames(map)
  sera$type <- "sera"
  if ("V3" %in% colnames(antigen)) {
  a <- full_join(antigen, 
                 sera,
                 by = c("V1", "V2", "V3", "identifier", "type"))
  } else if ("V2" %in% colnames(antigen)) {
    a <- full_join(antigen, 
                   sera,
                   by = c("V1", "V2", "identifier", "type"))
  } else if ("V1" %in% colnames(antigen)) {
    a <- full_join(antigen, 
                   sera,
                   by = c("V1", "identifier", "type"))
  } else stop("Don't know what to join by")
  
  if (MergeWithKey == FALSE) {return(a)} else {
      b <- a %>%
      mutate(uniq_id = ifelse(type == "sera", word(identifier, -1), NA)) %>%
      full_join(UniqIDKey)
    return(b)
  }
}

#
#
######################### TI_desc_distance() #########################
#
#

# Function to order the strain panel by increasing antigenic distance. Layout is similar to Yang Ge's graphs. 

TI_desc_distance <- function(distance_df, vaccine_strain, antigenic_method, dose_of_interest, min_age = 0, max_age = 150, fill_by = "none") {
  
  #### ARGUMENTS ####
  
  # distance_df = A dataframe.  Contains the distances and information about the titers. Needs to contain a column of ages (`age`). A column of vaccine strain names (`vaccine_fullname`), distance measure methods (`method`), distances of those measures (`distance`) and the dose of the individual (`dose`). 
  
  # vaccine_strain = A string. The reference vaccine strain for comparison too. There must be a column in the distance_df that is labeled `vaccine_fullname` because this is where it filters.
  
  # antigenic_method = A string. Distance measure that you're interested in. 
  
  # dose = A string. Select between, "SD", "HD", or c("SD", "HD")
  
  # min_age = An integer. Minimum age that should be included in the plots. Default is set to 0. 
  
  # max_age = An integer. Maximum age that should be included in the plots. Default is set to 150.
  
  # fill_by = A column name that it should be separated by. Options include: "none", "gender", "dose"
  
  #This section determines the factor levels of the different viruses in the panel by ordering them by decreasing distance (chosen by `antigenic_method`), relative to the chosen `vaccine_strain`. These are not necessarily in order of year!
  
  fac_levels <- distance_df %>%
    filter(vaccine_fullname == vaccine_strain, 
           method == antigenic_method) %>%
    select(strain_short, distance) %>%
    distinct() %>%
    arrange(desc(distance)) %>%
    pull(strain_short)
  
  if (fill_by == "none") {
    
    plt <- distance_df %>%
      filter(age >= min_age & age <= max_age, 
             vaccine_fullname == vaccine_strain, 
             method == antigenic_method,
             dose %in% dose_of_interest,
             !is.na(distance)) %>%
      ggplot(aes(x = titerincrease, 
                 y = factor(strain_short, 
                            levels = fac_levels))) +
      geom_boxplot() +
      labs(title = paste(vaccine_strain, 
                         "Vaccine Titer Increase and", 
                         antigenic_method),
           subtitle = paste0("Age Range: ",
                             min_age,
                             " - ",
                             max_age,
                             "\nIncludes ",
                             dose_of_interest,
                             " Dose Vaccinations")) +
      xlab("Log 2 Titer Increase") +
      ylab(paste("Ordered by Increasing", 
                 antigenic_method))
  }
  
  if (fill_by == "gender") {
    
    plt <- distance_df %>%
      filter(age >= min_age & age <= max_age, 
             vaccine_fullname == vaccine_strain, 
             method == antigenic_method,
             dose %in% dose_of_interest,
             !is.na(distance)) %>%
      ggplot(aes(x = titerincrease, 
                 y = factor(strain_short, 
                            levels = fac_levels),
                 fill = gender)) +
      geom_boxplot() +
      labs(title = paste(vaccine_strain, 
                         "Vaccine Titer Increase and", 
                         antigenic_method),
           subtitle = paste0("Age Range: ",
                             min_age,
                             " - ",
                             max_age,
                             "\nIncludes ",
                             dose_of_interest,
                             " Dose Vaccinations")) +
      xlab("Log 2 Titer Increase") +
      ylab(paste("Ordered by Increasing", 
                 antigenic_method))
  }
  
  if (fill_by == "dose") {
    
    plt <- distance_df %>%
      filter(age >= min_age & age <= max_age, 
             vaccine_fullname == vaccine_strain, 
             method == antigenic_method,
             dose %in% dose_of_interest,
             !is.na(distance)) %>%
      ggplot(aes(x = titerincrease, 
                 y = factor(strain_short, 
                            levels = fac_levels),
                 fill = dose)) +
      geom_boxplot() +
      labs(title = paste(vaccine_strain, 
                         "Vaccine Titer Increase and", 
                         antigenic_method),
           subtitle = paste0("Age Range: ",
                             min_age,
                             " - ",
                             max_age,
                             "\nIncludes ",
                             dose_of_interest,
                             " Dose Vaccinations")) +
      xlab("Log 2 Titer Increase") +
      ylab(paste("Ordered by Increasing", 
                 antigenic_method))
    
  }
  
  
  return(plt)
}



########################################################################
#
#     CreateSaveAcmap() ##########
#
########################################################################

#Description: Creates Acmaps over `x` optimizations and saves the lowest stress optimization Acmap. 

CreateSaveAcmap <- function(dataframe, 
                            dimensions, 
                            optimizations, 
                            subtype, 
                            timing, 
                            specific_filtering, 
                            save = TRUE, 
                            subset = FALSE) {
  suppressWarnings({
    #Suppress the warnings generated by Racmacs. It raises a flag about unstable, but they mention that many points may make raise this warning.

map <- make.acmap(sr_names = colnames(dataframe)[-1], 
           ag_names = as.character(dataframe$strains_fullname), 
           titer_table = dataframe[,-1], 
           number_of_dimensions = dimensions, 
           number_of_optimizations = optimizations) %>%
  keepBestOptimization()
})

if (save == TRUE) {
  if (subset == FALSE) {
  save.acmap(map, filename = here::here(paste0("data/processed/", 
                                             subtype, 
                                             "_maps/", 
                                             subtype, 
                                             "_", 
                                             timing, 
                                             "_", 
                                             specific_filtering, 
                                             "_", 
                                             dimensions, 
                                             "d.ace"))) 
  } else {
    save.acmap(map, filename = here::here(paste0("data/processed/subset_analysis/cartography/", 
                                                 subtype, 
                                                 "_", 
                                                 specific_filtering, 
                                                 ".ace"))) 
  }} else { return(map)}

}


########################################################################
#
#     subset_Acmaps()  #######
#
########################################################################

## Create antigenic cartographies of a smaller subset of data. From the subset datasets, it filters the subsets to only SD recipients, removes the underconstrained points for the post-titer maps in 2 dimensions, and then creates and save the AcMaps after 100 optimizations stratified by subtype (H1N1 or H3N2).

# Create the files that should be pulled. Including the subset analysis and the full datasets
subset_Acmaps <- function(location = "data/processed/subset_analysis/datasets") {
  
  file_list <- paste0(location, "/", list.files(path = here::here(location))) #Pull all of the subset file names from the datasets folder
  
  subsets <- lapply(X = file_list, FUN = function(.) readRDS(file = here::here(.))) #Pull the actual dataframes
  names(subsets) <- str_match(file_list, "datasets/\\s*(.*?)\\s*.rds")[,2] #Name the subset based on the subset number
  
  subsets <- lapply(subsets, filter, dose == "SD")
  # For the Subset analysis only include SD and post vaccination titers:
  # Only looking at 2 dimensions
  
  for (i in 1:length(subsets)) {
    for (j in c("H1N1", "H3N2")) {
      
      reduced_dataset <- remove_underconstrained(subsets[[i]], #Remove the underconstrained values
                                                 subtype = j, 
                                                 time = "postiter", 
                                                 dimensions = 2)
      CreateSaveAcmap(dataframe = reduced_dataset, #Save the Acmaps to the data folder
                      dimensions = 2, 
                      optimizations = 100, 
                      subtype = j, 
                      timing = "post", 
                      specific_filtering = names(subsets)[i],
                      subset = TRUE)
    }
  }
  
}
