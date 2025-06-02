## Prepare the data ################

prepare_data <- function(dataframe,
                         distance_method = c("cart_2d_post", "p_epi", "year"), # Distance measures to include
                         strains_of_interest = c("H1N1", "H3N2"),
                         subset = FALSE) #Subtypes to include 
  {
  
  ## Variables for troubleshooting
  
  #dataframe <- readRDS(file = here::here("data/processed/subset_analysis/distance/distance_data.rds"))
  ##distance_method <- c("cart_2d_post", "p_epi", "year")
  #strains_of_interest <- c("H1N1", "H3N2")
  #subset <- TRUE
 df <-  dataframe %>%
   {
     if (subset == TRUE) {
       dplyr::select(.data = .,
                     c(subset, uniq_id, strain_type, season, age, dose, vac_short, prevactiter, postvactiter, titerincrease, method, distance, strain_short))
     } else {
     dplyr::select(.data = .,
                   c(uniq_id, strain_type, season, age, dose, vac_short, prevactiter, postvactiter, titerincrease, method, distance, strain_short))  }
   } %>%
    droplevels() %>%
    dplyr::filter(method %in% distance_method, 
                !is.na(distance),
                strain_type %in% strains_of_interest,
                season != 2020,
                !is.na(titerincrease)) %>%
    dplyr::mutate(method = factor(method, 
                                  levels = c("year", "p_epi", "cart_2d_post"), 
                                  labels = c("Year", "P-epitope", "Cartography")),
         dose = factor(dose, 
                       levels = c("SD", "HD")),
         strain_type = factor(strain_type, 
                              levels = c("H1N1", "H3N2")), #Will need to be changed to include type B
         vac_short = forcats::fct_cross(strain_type, vac_short)) %>% #Include the subtype in the vaccine name
   {
     if (subset == TRUE) {
       dplyr::group_by(.data = .,
                       subset, method, season, vac_short) 
     } else {
       dplyr::group_by(.data = .,
                       method, season, vac_short)
         
     }
   } %>%
  dplyr::mutate(distance_season = distance/max(distance)) %>% #Distance for comparing within subtype and season
   {
     if (subset == TRUE) {
       dplyr::group_by(.data = .,
                       subset, method, vac_short)
     } else {
       dplyr::group_by(.data = ., 
                       method, vac_short)
     }
   } %>%
  dplyr::mutate(distance_strain = distance/max(distance)) %>% #Distance for comparing within subtype and strain
   {
     if (subset == TRUE) {
       dplyr::group_by(.data = .,
                       subset, method, season)
     } else {
       dplyr::group_by(.data = .,
                       method, season)
     }
   } %>%
  dplyr::mutate(distance_vaccine = distance/max(distance)) %>% #Distance for comparing within a season/vaccine
  dplyr::ungroup() %>%
  tidyr::pivot_longer(cols = c(distance_season:distance_vaccine), 
               names_to = "distance_type",
               values_to = "distance_norm") %>%
  dplyr::mutate(distance_type = factor(distance_type,
                                levels = c("distance_season", "distance_strain", "distance_vaccine"),
                                labels = c("Season", "Strain", "Vaccine"))) %>%
    droplevels() %>%
    dplyr::distinct()
 
 return(df)

}

## Individual Plot ####################
individual_plot <- function(dataframe,
                            individual_to_plot = 1) {

## Variables for troubleshooting

#dataframe <- data
#individual_to_plot <- 101
 p <- dataframe %>%
  #Choose one individual to look at who received SD:
  #Looking at only the distance_type of Season, so only looking at one specific season
  dplyr::filter(distance_type == "Season",
         dose == "SD",
         uniq_id == unique(dataframe$uniq_id)[individual_to_plot]) %>%
  tidyr::pivot_longer(cols = c(prevactiter, postvactiter),
               names_to = 'timing',
               values_to = 'titer') %>%
  #Look Pre-/Post-vaccination titers:
   dplyr::mutate(timing = factor(timing, 
                                 levels = c('prevactiter', 'postvactiter'), 
                                 labels = c('Pre-vaccination', 'Post-vaccination'))) %>%
   ggplot2::ggplot(aes(x = distance_norm,
             y = titer,
             color = timing)) +
   ggplot2::geom_line() + 
   ggplot2::geom_point() +
   ggplot2::facet_grid(rows = vars(vac_short),
             cols = vars(method)) +
   ggplot2::theme_bw() +
   ggplot2::labs(color = "Timing") +
   ggplot2::xlab("Normalized Distance") +
   ggplot2::ylab(latex2exp::TeX("$log_2$(HAI Titer)")) +
   ggplot2::theme(legend.position = "bottom") +
  #Add labels to indicate the HAI virus:
   ggrepel::geom_text_repel(data = subset(dataframe, 
                                          uniq_id == unique(dataframe$uniq_id)[individual_to_plot] & distance_type == "Season"), 
                  aes(x = distance_norm, 
                      y = -1, 
                      label = strain_short), 
                  size = 2, 
                  angle = 90, 
                  inherit.aes = FALSE, 
                  direction = "x") 
 
 return(p)
}

## Linear Models ######################

linear_model <- function(dataframe,
                         normalization = c("Season", "Strain", "Vaccine", "Original"),
                         outcome = c("prepost", "titerincrease"), 
                         dose = c("SD", "SDHD"),
                         subset = FALSE) {
  
  #Set warnings for arguments
  if (normalization %in% c("Season", "Strain", "Vaccine", "Original") == FALSE) 
    stop(paste(normalization, ":normalization method does not exist. Please use Season, Strain, Vaccine or Original."))
  if (outcome %in% c("prepost", "titerincrease") == FALSE) 
    stop("Outcome does not exist. Please use: prepost or titerincrease")
  if (dose %in% c("SD", "SDHD") == FALSE)
    stop("Set dose to either: SD for all individuals with standard dose vaccination or SDHD for all indiviudals equal to or greater than 65 years of age, for both the SD and HD vaccines")
  
  # Downselect the data frame to the data of interest
    dataframe <- dataframe %>%
      {
        if (dose == "SD")
        dplyr::filter(.data = .,
               dose == "SD") #Keep only SD regardless of age
        else 
          dplyr::filter(.data = .,
                 age >=65) #Keep only elderly individuals
        } %>%
      dplyr::filter(distance_type == normalization) %>% #investigate one normalization scheme at a time
      droplevels()

  
  if (outcome == "prepost") {
    #If the outcome is pre- and post-vaccination, the dataframe needs to be pivoted. 
    
    mdl <- dataframe %>%
      tidyr::pivot_longer(cols = c(prevactiter, postvactiter),
                   names_to = "timing",
                   values_to = "titer") %>%
      dplyr::mutate(timing = factor(timing, 
                             levels = c("prevactiter", "postvactiter"), 
                             labels = c("Pre-vaccination", "Post-vaccination"))) %>%
      {
        if(normalization %in% c("Season", "Vaccine")) {
          if (subset == TRUE) 
            dplyr::group_by(.data = ., 
                     subset, strain_type, distance_type, method, vac_short, timing, dose, season) else
              dplyr::group_by(.data = ., 
                       strain_type, distance_type, method, vac_short, timing, dose, season)
         #Season and vaccine normalization requires the inclusion of the season variable
           } else {
             if (subset == TRUE) 
               dplyr::group_by(.data = ., subset, strain_type, distance_type, method, vac_short, timing, dose) else
                 dplyr::group_by(.data = ., strain_type, distance_type, method, vac_short, timing, dose)
           }
            #Strain normalization does not require the season variable since some strains are administered more than one season
           
      } %>%
      dplyr::summarize(count = dplyr::n_distinct(uniq_id), #Summarize the linear model results
                intercept = stats::lm(titer ~ distance_norm)$coefficients[["(Intercept)"]], #Starting
                slope = stats::lm(titer ~ distance_norm)$coefficients[["distance_norm"]],
                rse = summary(stats::lm(titer ~ distance_norm))$sigma,
                df1 = summary(stats::lm(titer ~ distance_norm))$df[[1]],
                dfresdiual = summary(stats::lm(titer ~ distance_norm))$df[[2]],
                rsqr = summary(stats::lm(titer ~ distance_norm))$r.squared,
                f_stat = summary(stats::lm(titer ~ distance_norm))$fstatistic[["value"]],
                f_ndf = summary(stats::lm(titer ~ distance_norm))$fstatistic[["numdf"]],
                f_ddf = summary(stats::lm(titer ~ distance_norm))$fstatistic[["dendf"]],
                max_x = max(distance_norm), #Maximum x-value. usually 1, but for vaccine-normalization it does change
                pred = seq(0, max_x, length.out = 11), #X-values for prediction fits. Use the maximum x value for the upper range
                conf_upr95 = as.numeric(stats::predict(stats::lm(titer ~ distance_norm), 
                                                newdata = data.frame(distance_norm = seq(0, max_x, length.out = 11)),
                                                interval = "confidence")[,3]),
                conf_lwr95 = as.numeric(stats::predict(stats::lm(titer ~ distance_norm), 
                                                newdata = data.frame(distance_norm = seq(0, max_x, length.out = 11)),
                                                interval = "confidence")[,2])) %>%
      dplyr::rowwise() %>%
      dplyr::mutate(dose_label = paste0("(n=", count, ")"), #Create a legend label
             dose_label = forcats::fct_cross(dose, dose_label, sep = " ")) %>% #Retain the factor levels of dose
      dplyr::ungroup()
    
  } else { #If outcome == titerincrease
    
    mdl <- dataframe %>%
      {
        if(normalization %in% c("Season", "Vaccine")) {
          if (subset == TRUE) {
            dplyr::group_by(.data = ., 
                     subset, strain_type, distance_type, method, vac_short, dose, season)
          }  else {
              dplyr::group_by(.data = .,
                       strain_type, distance_type, method, vac_short, dose, season)
          }} else {
            if (subset == TRUE) {
              dplyr::group_by(.data = ., 
                       subset, strain_type, distance_type, method, vac_short, dose)
            }  else {
              dplyr::group_by(.data = ., 
                       strain_type, distance_type, method, vac_short, dose)
            }
          }
        } %>%
          #include season variable for only season and vaccine normalization
      dplyr::summarize(count = dplyr::n_distinct(uniq_id),#Summarize the linear regression results
                intercept = stats::lm(titerincrease ~ distance_norm)$coefficients[["(Intercept)"]], #Starting
                slope = stats::lm(titerincrease ~ distance_norm)$coefficients[["distance_norm"]],
                rse = summary(stats::lm(titerincrease ~ distance_norm))$sigma,
                df1 = summary(stats::lm(titerincrease ~ distance_norm))$df[[1]],
                dfresdiual = summary(stats::lm(titerincrease ~ distance_norm))$df[[2]],
                rsqr = summary(stats::lm(titerincrease ~ distance_norm))$r.squared,
                f_stat = summary(stats::lm(titerincrease ~ distance_norm))$fstatistic[["value"]],
                f_ndf = summary(stats::lm(titerincrease ~ distance_norm))$fstatistic[["numdf"]],
                f_ddf = summary(stats::lm(titerincrease ~ distance_norm))$fstatistic[["dendf"]],
                max_x = max(distance_norm),
                pred = seq(0, max_x, length.out = 11),
                conf_upr95 = as.numeric(stats::predict(stats::lm(titerincrease ~ distance_norm), 
                                                newdata = data.frame(distance_norm = seq(0, max_x, length.out = 11)),
                                                interval = "confidence")[,3]),
                conf_lwr95 = as.numeric(stats::predict(stats::lm(titerincrease ~ distance_norm), 
                                                newdata = data.frame(distance_norm = seq(0, max_x, length.out = 11)),
                                                interval = "confidence")[,2])) %>%
      dplyr::ungroup() %>%
      dplyr::rowwise() %>%
      dplyr::mutate(dose_label = paste0("(n=", count, ")"), #Create a legend label
             dose_label = forcats::fct_cross(dose, dose_label, sep = " ")) %>% #Retain the factor levels of dose
      dplyr::ungroup()
  }
  
  return(mdl)
  
}



## Sample Size ###################################
# Determine the sample size of individuals in each grouping

sample_size_fun <- function(dataframe,
                            subset = FALSE) {
  
  if (unique(dataframe$distance_type) == "Strain") { #if By strain don't include season
    if (subset == TRUE) {
      sample_size <- dataframe %>%
        dplyr::select(subset, strain_type, vac_short, dose, count) %>%
        dplyr::distinct() %>%
        { 
          if("HD" %in% levels(.$dose)) {
            dplyr::mutate(.data = .,
                  count = paste0(dose, ": ", count)) %>%
              dplyr::group_by(subset, strain_type, vac_short) %>%
              dplyr::summarize(count = paste(count, collapse = "; ")) %>%
              dplyr::ungroup() 
         } else . 
          } %>%
        dplyr::mutate(count = paste0(vac_short, ": ", count)) %>% 
        dplyr::group_by(subset, strain_type) %>%
        dplyr::summarize(count = paste(count, collapse = "\n"))
    } else {
    sample_size <- dataframe %>%
      dplyr::select(strain_type, vac_short, dose, count) %>%
      dplyr::distinct() %>%
      { 
        if("HD" %in% levels(.$dose)) {
          dplyr::mutate(.data = .,
                count = paste0(dose, ": ", count)) %>%
            dplyr::group_by(strain_type, vac_short) %>%
            dplyr::summarize(count = paste(count, collapse = "; ")) %>%
            dplyr::ungroup() 
          } else . 
        } %>%
      dplyr::mutate(count = paste0(vac_short, ": ", count)) %>% 
      dplyr::group_by(strain_type) %>%
      dplyr::summarize(count = paste(count, collapse = "\n"))
 } } else {  #By season and vaccine are stratified by season not strain
    
   if (subset == TRUE) {
    sample_size <- dataframe %>%
      dplyr::select(subset, season, dose, count) %>%
      dplyr::distinct() %>%
      { 
        if("HD" %in% levels(.$dose)) {
          dplyr::mutate(.data = .,
                count = paste0(dose, ": ", count)) %>%
            dplyr::group_by(subset, season) %>%
            dplyr::summarize(count = paste(count, collapse = "; ")) %>%
            dplyr::ungroup() 
          } else . 
        } %>%
      dplyr::group_by(subset, season) %>%
      dplyr::summarize(count = paste(count, collapse = "\n"))
   } else {
     sample_size <- dataframe %>%
       dplyr::select(season, dose, count) %>%
       dplyr::distinct() %>%
       { 
         if("HD" %in% levels(.$dose)) {
           dplyr::mutate(.data = .,
                 count = paste0(dose, ": ", count)) %>%
             dplyr::group_by(season) %>%
             dplyr::summarize(count = paste(count, collapse = "; ")) %>%
             dplyr::ungroup() 
           } else . 
         } %>%
       dplyr::group_by(season) %>%
       dplyr::summarize(count = paste(count, collapse = "\n"))
   }
    
  }
  return(sample_size)
}


## Equations #########################
equation_labels <- function(dataframe,
                            subset = FALSE) {
  
  dataframe <- dataframe %>%
    dplyr::mutate(across(where(is.numeric), 
                         round, 
                         3))
  
  if ("timing" %in% colnames(dataframe)) { # This means it's pre/post
    dataframe <- dataframe %>%
      { 
        if (unique(dataframe$distance_type) == "Strain") {
          if (subset == TRUE) {
            dplyr::select(.data = ., 
                   subset, strain_type, dose, vac_short, method, timing, slope, intercept, rsqr) 
         } else {
            dplyr::select(.data = ., 
                   strain_type, dose, vac_short, method, timing, slope, intercept, rsqr) 
       }
          }  else {
          if (subset == TRUE) {
            dplyr::select(.data = ., 
                   subset, strain_type, season, dose, vac_short, method, timing, slope, intercept, rsqr) 
           } else {
              dplyr::select(.data = ., 
                     strain_type, season, dose, vac_short, method, timing, slope, intercept, rsqr) 
          }}
      } %>%
      dplyr::distinct() %>%
      dplyr::mutate(rsqr = ifelse(rsqr < 0.001,
                           "<0.001",
                           paste(rsqr))) %>%
      { 
        if("HD" %in% levels(.$dose)) {
          dplyr::mutate(.data = .,
               working = paste0("Titer=",
                                slope,
                                "*Distance+",
                                intercept,
                                " (Rsqr:",
                                rsqr,
                                ")")) %>%
          dplyr::select(-c(slope, intercept, rsqr)) %>%
          tidyr::pivot_wider(names_from = timing,
                      values_from = working) %>%
          dplyr::mutate(equation = paste0("Pre:",
                                   `Pre-vaccination`,
                                   "\nPost:",
                                   `Post-vaccination`),
                 .keep = "unused") %>%
          tidyr::pivot_wider(names_from = dose, 
                      values_from = equation) %>%
          dplyr::mutate(equation = paste0("SD:",
                                   SD,
                                   "\nHD:",
                                   HD),
                 .keep = "unused") 
          }
        else {
          tidyr::pivot_wider(data = .,
                      names_from = timing, 
                      values_from = c(slope:rsqr)) %>%
          dplyr::mutate(equation = paste0("Pre: Titer=", 
                                     `slope_Pre-vaccination`,
                                     "*Distance+",
                                     `intercept_Pre-vaccination`,
                                     "\nRsqr: ",
                                     `rsqr_Pre-vaccination`, 
                                     "\nPost: Titer=", 
                                     `slope_Post-vaccination`,
                                     "*Distance+",
                                     `intercept_Post-vaccination`,
                                     "\nRsqr: ",
                                     `rsqr_Post-vaccination`),
                   .keep = "unused")
        }
        }
  } else {
    dataframe <-  dataframe %>%
      { 
        if (unique(dataframe$distance_type) == "Strain") {
          if (subset == TRUE) {
            dplyr::select(.data = ., 
                   subset, strain_type, dose, vac_short, method, slope, intercept, rsqr)
          } else {
            dplyr::select(.data = ., 
                   strain_type, dose, vac_short, method, slope, intercept, rsqr) 
          }
        } else {
          if (subset == TRUE) {
            dplyr::select(.data = ., 
                   subset, strain_type, season, dose, vac_short, method, slope, intercept, rsqr) 
          } else {
            dplyr::select(.data = ., 
                   strain_type, season, dose, vac_short, method, slope, intercept, rsqr) 
          }
        }
      } %>%
      dplyr::distinct() %>%
      {
        if ("HD" %in% levels(.$dose)) {
            dplyr::mutate(.data = .,
                   equation = paste0("TI=", 
                                     slope,
                                     "*Distance+",
                                     intercept),
                   .keep = "unused") %>%
            dplyr::select(-rsqr) %>%
          tidyr::pivot_wider(names_from = dose, 
                             values_from = equation) %>%
            dplyr::mutate(equation = paste0("SD: ",
                                     SD,
                                     "\nHD: ",
                                     HD),
                   .keep = "unused")
        } else {
          dplyr::mutate(.data = .,
                       rsqr = ifelse(rsqr < 0.001,
                                     "<0.001",
                                     paste(rsqr)),
                       equation = paste0("TI=", 
                                         slope,
                                         "*Distance+",
                                         intercept,
                                         "\nRsqr: ",
                                         rsqr),
                       .keep = "unused") %>%
            dplyr::select(-rsqr) }
      }
  }
  return(dataframe)
}



## Mass Linear Regression Model Save Function ##########

save_linear_models <- function(dataframe = data,
                                norms = c("Strain", "Season", "Vaccine"),
                                doses = c("SD", "SDHD"),
                                outcomes = c("prepost", "titerincrease"),
                                strains = c("H1N1", "H3N2"),
                                seasons = 2014:2019,
                                save_plots = TRUE,
                                save_files = TRUE,
                                save_with_labels = FALSE) {
  
  #Variables for troubleshooting:
  #dataframe <- data
  #norms <- c("Strain", "Season", "Vaccine")
  #doses <- c("SD", "SDHD")
  #outcomes <- c("prepost", "titerincrease")
  #strains <- c("H1N1", "H3N2")
  #seasons <- 2014:2019
  #save_plots <- TRUE
  #save_files <- FALSE
  #save_with_labels <- FALSE
  
  
  ## Create the linear regression models: 
  mdl_list <- vector(mode = "list")
  for(i in norms) {
    for (j in outcomes) {
      for (k in doses) {
        x <- paste(i, j, k)
        mdl_list[[x]] <- linear_model(dataframe = data,
                                      normalization = i,
                                      outcome = j,
                                      dose = k)
      }
    }
  }
  
  print("Models are completed")
  
  ## Determine samples sizes for the plots:
  sample_size_list <- vector(mode = "list")
  
  for (i in 1:length(mdl_list)) {
    x <- names(mdl_list[i])
    sample_size_list[[x]] <- sample_size_fun(dataframe = mdl_list[[i]])
  }
  
  print("Sample sizes are completed")
  
  ##Determine the linear regression equations for the plots.
  
  slope_list <- vector(mode = "list")
  
  for (i in 1:length(mdl_list)) {
    x <- names(mdl_list[i])
    slope_list[[x]] <- equation_labels(dataframe = mdl_list[[i]])
  }
  
  print("Slopes are completed")
  
  if (save_files == TRUE) {
    saveRDS(mdl_list,
            file = here::here("data/processed/linear_regression_models/mdl_list.rds"))
    
    saveRDS(sample_size_list,
            file = here::here("data/processed/linear_regression_models/size_list.rds"))
    
    saveRDS(slope_list,
            file = here::here("data/processed/linear_regression_models/slope_list.rds"))
    print("Saved Lists")
  }  
  
}
  



## Season Level Plots##########

# This function creates the season-based antibody landscape plots for the manuscript. 

season_plots <- function(dataframe, 
                            # The data frame that contains the individual titers and distances
                         season_to_plot, 
                            # The season of interest, this is a character between 2014:2019
                         mdl_to_use, 
                            # String that contains: "{distance_type} {outcome of interest} {dose of interest}" in that order and form. Distance types available: Season, Vaccine, Strain. Outcome of interest available: prepost and titerincrease. Dose of interest available: SD or SDHD.
                         mdl = mdl_list,
                            # The list of linear models from the dataset. The mdl_to_use needs to match one of the named linear models component
                         slopes = slope_list,
                            # The list of slopes from the dataset. This data will be adding labels that contain the linear regression results to the upper right hand corner. The mdl_to_use needs to match one of the named slope_list component
                         sample_size = sample_size_list,
                            # The list of sample sizes from the dataset. The mdl_to_use needs to match one of the named slope_list component
                         include_equation_label = FALSE)
                            # Binary argument to turn on or off the addition of the equations in the upper corner of the plots.
{
## Variables for troubleshooting
 # dataframe <- data
#  season_to_plot <- 2014
#  mdl_to_use <- "Vaccine titerincrease SD"
#  mdl = mdl_list
#  slopes <- slope_list
#  sample_size = sample_size_list
#  include_equation_label = FALSE
  
  # Create the variables
  
  variables <- base::strsplit(mdl_to_use,
                        split = " ")[[1]] # Using the mdl_to_use separate into the individual settings
  
  # Use these to filter the data set to the desired data:
  normalization <- variables[1] # Used to filter the datasets
  
  outcome <- variables[2] #Used for the y-axis. If prepost, plots both pre and post, if titerincrease only titerincrease
  
  dose <- variables[3] #Used for grouping. If SD, only plots individuals with SD vaccination. If SDHD plots individuals greater than 64 and distinguished between vaccine dose.
  
  # Reduce the linear model, slopes, and samples size lists to the important ones.
  # Filter to include only the season of interest
  mdl <- mdl[[mdl_to_use]] %>%
    dplyr::filter(season == season_to_plot) %>%
    droplevels()
  
  slopes <- slopes[[mdl_to_use]] %>%
    dplyr::filter(season == season_to_plot) %>%
    droplevels()
  
  samples <- sample_size[[mdl_to_use]] %>%
    dplyr::filter(season == season_to_plot) %>%
    droplevels()
  
  # Reduce the dataframe to the SD only or the >= 65 age only.
  # Filter to include the distance of interest and season
  
    dataframe <- dataframe %>%
      dplyr::filter(distance_type == normalization,
             season == season_to_plot) %>%
      {
        if (dose == "SD") { #If SD need to remove all individuals who did not receive SD
          dplyr::filter(.data = ., dose == "SD") 
        } else if (dose == "SDHD") { #If SDHD need to remove all non elderly individuals
          dplyr::filter(.data = ., age >= 65)
        } else { #Produce a warning if the dose does not match
          stop(paste(dose, ": Not found. Only the SD only individuals all ages, and the SD HD individuals of 65+ years are to be analyzed using this code.")) 
        }
      } %>%
      droplevels() %>%
    dplyr::group_by(dose) %>%
    dplyr::mutate(dose_label = paste0("(n=", dplyr::n_distinct(uniq_id), ")"),
           dose_label = forcats::fct_cross(dose, dose_label, sep = " ")) %>%
    dplyr::ungroup()
  
  # Plot the prepost after pivoting and relabeling. Usually this plot is only for the standard dose only comparison. Not the HD/SD. Code isn't optimized for HD/SD comparison. 
  
  if (outcome == "prepost") {
    
    p <-  dataframe %>%
      tidyr::pivot_longer(cols = c(prevactiter, postvactiter),
                   names_to = "timing",
                   values_to = "titer") %>%
      dplyr::mutate(timing = factor(timing, 
                                    levels = c("prevactiter", "postvactiter"), 
                                    labels = c("Pre-vaccination", "Post-vaccination"))) %>%
      dplyr::arrange(distance_norm) %>%
      ggplot2::ggplot(aes(x = distance_norm, 
                 y = titer,
                 color = timing,
                 fill = timing, #Distinguish timing by color/fill
                 linetype = timing)) + # Distinguish different doses by linetype
      ggplot2::geom_jitter(size = 1, shape = ".", width = 0) + #Spreads out the points for readability. Only apply jitter in the Y-direction since the HAI titers are discrete and we don't want to apply jitter in the x-direction so as to not confuse the reader with regards to what the distance is
      ggplot2::geom_line(data = mdl,
                   aes(x = pred,
                       y = slope*pred+intercept,
                       color = timing)) + # Include the Linear regression best-fit line. Does not extend past the points used for linear regression
      ggplot2::geom_ribbon(data = mdl,
                  aes(x = pred,
                      y = pred*slope+intercept,
                      ymax = conf_upr95,
                      ymin = conf_lwr95),
                  alpha = 0.5,
                  color = NA) + # Include the 95% se confidence intervals of the models. do not include bounding lines
      ggplot2::facet_grid(rows = vars(vac_short, dose_label), 
                 cols = vars(method)) + #Stratify by vaccine received (H1N1 and H3N2 component) and the distance calculation method
      ggplot2::ylab(latex2exp::TeX("$log_2$(HAI Titer)")) +
      ggplot2::xlab("Distance from Vaccine Strain") +
      ggplot2::labs(title = paste(season_to_plot, "Season"),
           color = "Timing",
           linetype = "Timing") + 
      ggplot2::scale_y_continuous(breaks = seq(0, 10, 2),
                         limits = c(0, 10)) +
      ggplot2::theme_bw()+
      theme(legend.position = "bottom") +
      ggplot2::guides(color = guide_legend(override.aes=list(fill=NA,
                                                    alpha = 1)),
             linetype = guide_legend(override.aes=list(fill = NA,
                                                       alpha = 1)),
             fill = "none") + # Only allow the line color/ type to be shown in the legend.
      ggplot2::scale_color_manual(values = c("#56B4E9", "#E69F00")) + # Blue and Orange
      ggplot2::scale_fill_manual(values = c("#56B4E9", "#E69F00")) # Blue and Orange
    
    if (dose == "SDHD") 
      warning("season_plot() function is not optimized for plotting pre-/post-vaccination titers and SD and HD vaccinees.")
    
  }
  
  if (outcome == "titerincrease") {
    if (dose == "SDHD") {
    
    p <- dataframe %>%
      ggplot2::ggplot(aes(x = distance_norm, 
                 y = titerincrease,
                 color = dose_label,
                 fill = dose_label,
                 linetype = dose_label)) +
      ggplot2::geom_jitter(size = 1, shape = ".", width = 0) +
      ggplot2::geom_line(data = mdl,
                aes(x = pred,
                    y = slope*pred+intercept,
                    color = dose_label)) +
      ggplot2::geom_ribbon(data = mdl,
                  aes(x = pred,
                      y = pred*slope+intercept,
                      ymax = conf_upr95,
                      ymin = conf_lwr95),
                  alpha = 0.5,
                  color = NA) +
      ggplot2::facet_grid(rows = vars(vac_short), 
                 cols = vars(method)) +
      ggplot2::ylab(latex2exp::TeX("$log_2$(HAI Titer) Increase")) +
      ggplot2::xlab("Distance from Vaccine Strain") +
      ggplot2::labs(title = paste(season_to_plot, "Season"),
           color = "Dose",
           linetype = "Dose") +
      ggplot2::scale_y_continuous(breaks = seq(-1, 4, 1),
                         limits = c(-1, 4)) +
      ggplot2::theme_bw()+
      ggplot2::geom_hline(yintercept = 0) +
      ggplot2::theme(legend.position = "bottom") +
      ggplot2::guides(color = guide_legend(override.aes=list(fill=NA,
                                                    alpha = 1)),
             fill = "none") + # Only allow the line color/ type to be shown in the legend.
      ggplot2::scale_color_manual(values = c("#56B4E9", "#E69F00")) + # Blue and Orange
      ggplot2::scale_fill_manual(values = c("#56B4E9", "#E69F00")) # Blue and Orange
    
    } else if (dose == "SD") {
      
      p <- dataframe %>%
        ggplot2::ggplot(aes(x = distance_norm, 
                   y = titerincrease)) +
        ggplot2::geom_jitter(size = 1, shape = ".", width = 0) +
        ggplot2::geom_line(data = mdl,
                  aes(x = pred,
                      y = slope*pred+intercept)) +
        ggplot2::geom_ribbon(data = mdl,
                    aes(x = pred,
                        y = pred*slope+intercept,
                        ymax = conf_upr95,
                        ymin = conf_lwr95),
                    alpha = 0.5,
                    color = NA) +
        ggplot2::facet_grid(rows = vars(vac_short, dose_label), 
                   cols = vars(method)) +
        ggplot2::ylab(latex2exp::TeX("$log_2$(HAI Titer) Increase")) +
        ggplot2::xlab("Distance from Vaccine Strain") +
        ggplot2::labs(title = paste(season_to_plot, "Season")) +
        ggplot2::scale_y_continuous(breaks = seq(-1, 4, 1),
                           limits = c(-1, 4)) +
        ggplot2::theme_bw()+
        ggplot2::geom_hline(yintercept = 0,
                   linetype = "dashed") +
        ggplot2::guides(color = "none",
               fill = "none")# + # Only allow the line color/ type to be shown in the legend.
        #scale_color_manual(values = c("#56B4E9", "#E69F00")) + # Blue and Orange
       # scale_fill_manual(values = c("#56B4E9", "#E69F00")) # Blue and Orange
    }
    
    }
  if (include_equation_label == TRUE) {
    p <- p + 
      ggplot2::geom_label(data = slopes, 
                 aes(x = Inf, y = Inf, label = equation, vjust = 1, hjust = 1), 
                 size = 2,
                 inherit.aes = FALSE)
  }
  return(p)
}

## Season Level Subset Plots ##########

# This function creates the season-based antibody landscape plots for the manuscript but looking at the subsetted data

season_plots_subset <- function(dataframe, 
                         # The data frame that contains the individual titers and distances
                         season_to_plot, 
                         # The season of interest, this is a character between 2014:2019
                         mdl_to_use, 
                         # String that contains: "{distance_type} {outcome of interest} {dose of interest}" in that order and form. Distance types available: Season, Vaccine, Strain. Outcome of interest available: prepost and titerincrease. Dose of interest available: SD or SDHD.
                         mdl = mdl_list,
                         # The list of linear models from the dataset. The mdl_to_use needs to match one of the named linear models component
                         slopes = slope_list,
                         # The list of slopes from the dataset. This data will be adding labels that contain the linear regression results to the upper right hand corner. The mdl_to_use needs to match one of the named slope_list component
                         sample_size = sample_size_list)
                         # The list of sample sizes from the dataset. The mdl_to_use needs to match one of the named slope_list component
{
  # Data for troubleshooting
  ##dataframe <- data
  ##season_to_plot <- 2014
  ##mdl_to_use <- "Season titerincrease SDHD"
  ##mdl <- mdl_list
  ##slopes <- slope_list
  ##sample_size <- sample_size_list
  
  # Create the variables
  
  variables <- strsplit(mdl_to_use,
                        split = " ")[[1]] # Using the mdl_to_use separate into the individual settings
  
  # Use these to filter the data set to the desired data:
  normalization <- variables[1] # Used to filter the datasets
  
  outcome <- variables[2] #Used for the y-axis. If prepost, plots both pre and post, if titerincrease only titerincrease
  
  dose <- variables[3] #Used for grouping. If SD, only plots individuals with SD vaccination. If SDHD plots individuals greater than 64 and distinguished between vaccine dose.
  
  #Set warnings for arguments
  if (is.numeric(season_to_plot) == FALSE) 
    stop(paste(season_to_plot, ":season_to_plot needs to be a number."))
  if (outcome != "titerincrease") 
    stop(paste(outcome, ":Outcome does not exist. Please use only titerincrease"))
  if (dose %in% c("SD", "SDHD") == FALSE)
    stop("Set dose to either: SD for all individuals with standard dose vaccination or SDHD for all indiviudals equal to or greater than 65 years of age, for both the SD and HD vaccines")
  
  # Reduce the linear model, slopes, and samples size lists to the important ones.
  # Filter to include only the season of interest
  mdl <- mdl[[mdl_to_use]] %>%
    dplyr::filter(season == season_to_plot) %>%
    droplevels()
  
  slopes <- slopes[[mdl_to_use]] %>%
    dplyr::filter(season == season_to_plot) %>%
    droplevels()
  
  samples <- sample_size[[mdl_to_use]] %>%
    dplyr::filter(season == season_to_plot) %>%
    droplevels()
  
  # Reduce the dataframe to the SD only or the >= 65 age only.
  # Filter to include the distance of interest and season
  
  dataframe <- dataframe %>%
    dplyr::filter(distance_type == normalization,
                  season == season_to_plot) %>%
    {
      if (dose == "SD") { #If SD need to remove all individuals who did not receive SD
        dplyr::filter(.data = ., dose == "SD") 
      } else if (dose == "SDHD") { #If SDHD need to remove all non elderly individuals
        dplyr::filter(.data = ., age >= 65)
      } else { #Produce a warning if the dose does not match
        stop(paste(dose, ": Not found. Only the SD only individuals all ages, and the SD HD individuals of 65+ years are to be analyzed using this code.")) 
      }
    } %>%
    droplevels() %>%
    dplyr::group_by(dose) %>%
    dplyr::mutate(dose_label = paste0("(n=", dplyr::n_distinct(uniq_id), ")"),
                  dose_label = forcats::fct_cross(dose, dose_label, sep = " ")) %>%
    dplyr::ungroup()
  
  # Plot the prepost after pivoting and relabeling. Usually this plot is only for the standard dose only comparison. Not the HD/SD. Code isn't optimized for HD/SD comparison. 
    
    if (dose == "SD") {
    
    p <- dataframe %>%
      ggplot2::ggplot(aes(x = distance_norm, 
                 y = titerincrease,
                 color = subset,
                 fill = subset,
                 linetype = subset)) +
      ggplot2::geom_jitter(size = 1, shape = ".", width = 0) +
      ggplot2::geom_line(data = mdl,
                aes(x = pred,
                    y = slope*pred+intercept,
                    color = subset)) +
      ggplot2::geom_ribbon(data = mdl,
                  aes(x = pred,
                      y = pred*slope+intercept,
                      ymax = conf_upr95,
                      ymin = conf_lwr95),
                  alpha = 0.5,
                  color = NA) +
      ggplot2::facet_grid(rows = vars(vac_short), 
                 cols = vars(method)) +
      ggplot2::ylab(latex2exp::TeX("$log_2$(HAI Titer) Increase")) +
      ggplot2::xlab("Distance from Vaccine Strain") +
      ggplot2::labs(title = paste(season_to_plot, "Season", subset(dataframe,
                                                          season == season_to_plot)$dose_label[1]),
           color = "Subset",
           linetype = "Subset") +
      ggplot2::scale_y_continuous(breaks = seq(-1, 4, 1),
                         limits = c(-1, 4)) +
      ggplot2::theme_bw()+
      ggplot2::geom_hline(yintercept = 0) +
      ggplot2::theme(legend.position = "bottom") +
      ggplot2::guides(color = guide_legend(override.aes=list(fill=NA,
                                                    alpha = 1)),
             fill = "none") #+ # Only allow the line color/ type to be shown in the legend.
      #scale_color_manual(values = c("#56B4E9", "#E69F00")) + # Blue and Orange
      #scale_fillmanual(values = c("#56B4E9", "#E69F00")) # Blue and Orange
    
  return(p) 
    } 
    else 
      if (dose == "SDHD") {

      
      # Create the vector to save the strains to
      vaccine_virus <- levels(dataframe$vac_short) #Pull the vaccine strains to create individual plots for
      plt_list <- vector(mode = "list", 
                         length = length(vaccine_virus)) #Create the list for saving the plots to; make it the same length as the vaccine viruses
      names(plt_list) <- vaccine_virus #Name the plt_list after the vaccine strains
      
      
      for (i in 1:length(vaccine_virus)) {
        plt_list[[i]] <- dataframe %>%
          dplyr::filter(vac_short == vaccine_virus[i]) %>%
          ggplot2::ggplot(aes(x = distance_norm, 
                              y = titerincrease,
                              color = dose_label,
                              linetype = dose_label,
                              fill = dose_label)) +
          ggplot2::geom_jitter(size = 1, shape = ".", width = 0) +
          ggplot2::geom_line(data = subset(mdl, vac_short == vaccine_virus[i]),
                    aes(x = pred,
                        y = slope*pred+intercept)) +
          ggplot2::geom_ribbon(data = subset(mdl, vac_short == vaccine_virus[i]),
                      aes(x = pred,
                          y = pred*slope+intercept,
                          ymax = conf_upr95,
                          ymin = conf_lwr95),
                      alpha = 0.5,
                      color = NA) + 
          ggplot2::facet_grid(rows = vars(subset), 
                     cols = vars(method)) +
          ggplot2::ylab(latex2exp::TeX("$log_2$(HAI Titer) Increase")) +
          ggplot2::xlab("Distance from Vaccine Strain") +
          ggplot2::labs(title = paste("Season", season_to_plot, "---", vaccine_virus[i]),
               color = "Dose",
               fill = "Dose",
               linetype = "Dose") +
          ggplot2::scale_y_continuous(breaks = seq(-1, 4, 1),
                             limits = c(-1, 4)) +
          ggplot2::theme_bw()+
          ggplot2::geom_hline(yintercept = 0) +
          ggplot2::theme(legend.position = "bottom") +
          ggplot2::guides(color = guide_legend(override.aes=list(fill=NA,
                                                        alpha = 1)),
                 fill = "none") + # Only allow the line color/ type to be shown in the legend.
          ggplot2::scale_color_manual(values = c("#56B4E9", "#E69F00")) + # Blue and Orange
          ggplot2::scale_fill_manual(values = c("#56B4E9", "#E69F00")) # Blue and Orange
      } # Loop through vaccine strains to create individual plots
      return(plt_list)
    } #Creates Method ~ Subset, with 2 doses on each plot comparing SDHD results. Creates a list with each component being a vaccine strain
  }
  


## Subtype Level Plots##########

subtype_plots <- function(dataframe,
                         subtype_to_plot,
                         mdl_to_use,
                         mdl = mdl_list,
                         slopes = slope_list,
                         sample_size = sample_size_list,
                         include_equation_label = FALSE) {
  #Variables for troubleshooting:
  #dataframe <- data
  #subtype_to_plot <- "H3N2"
  #mdl_to_use <- "Strain titerincrease SD"
  #mdl <- mdl_list
  #slopes <- slope_list
  #sample_size <- sample_size_list

  
  # Create the variables 
  variables <- base::strsplit(mdl_to_use,
                        split = " ")[[1]]
  
  normalization <- variables[1]
  outcome <- variables[2]
  dose <- variables[3]
  
  # Reduce the lists to the important ones
  mdl <- mdl[[mdl_to_use]] %>%
    dplyr::filter(strain_type == subtype_to_plot) %>%
    droplevels()
  
  slopes <- slopes[[mdl_to_use]] %>%
    dplyr::filter(strain_type == subtype_to_plot) %>%
    droplevels()
  
  samples <- sample_size[[mdl_to_use]] %>%
    dplyr::filter(strain_type == subtype_to_plot) %>%
    droplevels()
  
  # Reduce the dataframe to the SD only or the >= 65 age only. 
  
    dataframe <- dataframe %>%
      dplyr::filter(distance_type == normalization, #Select the distance normalization
                    strain_type == subtype_to_plot) %>% #Select the subtype to be plotted
      {
        if (dose == "SD") { #If SD need to remove all individuals who did not receive SD
          dplyr::filter(.data = ., dose == "SD") 
        } else if (dose == "SDHD") { #If SDHD need to remove all non elderly individuals
          dplyr::filter(.data = ., age >= 65)
        } else { #Produce a warning if the dose does not match
            stop(paste(dose, ": Not found. Only the SD only individuals all ages, and the SD HD individuals of 65+ years are to be analyzed using this code.")) 
                 }
        } %>%
      droplevels() %>%
      dplyr::group_by(dose, vac_short) %>% # Add sample size labels to the corresponding dose for the legends
            dplyr::mutate(dose_label = paste0("(n=", dplyr::n_distinct(uniq_id), ")")) %>% # Total number of people in the dose 
      {
        if (dose == "SD") {
          dplyr::mutate(.data = .,
                        dose_label = forcats::fct_cross(dose, dose_label, sep = " ")) #Combine it with the dose
        } else if (dose == "SDHD") {
          dplyr::mutate(.data = .,
                        dose_label = forcats::fct_cross(dose, dose_label, sep = " ")) %>%
          dplyr::group_by(vac_short) %>%
          dplyr::arrange(dose) %>%
          dplyr::mutate(dose_label = paste(unique(dose_label), collapse = " | "))
        }
      } %>%
    ungroup()
  
  # Plot the prepost after pivoting and relabeling. Usually this plot is only for the standard dose only comparison. Not the HD/SD. 
  
  if (outcome == "prepost") {
    
    p <-  dataframe %>%
      tidyr::pivot_longer(cols = c(prevactiter, postvactiter),
                   names_to = "timing",
                   values_to = "titer") %>%
      dplyr::mutate(timing = factor(timing, levels = c("prevactiter", "postvactiter"), labels = c("Pre-vaccination", "Post-vaccination"))) %>%
      ggplot2::ggplot(aes(x = distance_norm, 
                 y = titer,
                 color = timing,#Distinguish timing by color/fill
                 fill = timing,
                 linetype = timing)) + # Distinguish different doses by linetype
      ggplot2::geom_jitter(size = 1, shape = ".", width = 0)  + #Spreads out the points for readability. Only apply jitter in the Y-direction since the HAI titers are discrete and we don't want to apply jitter in the x-direction so as to not confuse the reader with regards to what the distance is
      ggplot2::geom_line(data = mdl,
                aes(x = pred,
                    y = slope*pred+intercept,
                    color = timing,
                    linetype = timing)) +
      ggplot2::geom_ribbon(data = mdl,
                  aes(x = pred,
                      y = pred*slope+intercept,
                      ymax = conf_upr95,
                      ymin = conf_lwr95,
                      linetype = timing),
                  alpha = 0.5,
                  color = NA) + # Include the 95% se confidence intervals of the models. do not include bounding lines
      ggplot2::facet_grid(rows = vars(vac_short, dose_label), 
                 cols = vars(method)) + #Stratify by vaccine received (H1N1 and H3N2 component) and the distance calculation method
      ggplot2::ylab(latex2exp::TeX("$log_2$(HAI Titer)")) +
      ggplot2::xlab("Distance from Vaccine Strain") +
      ggplot2::labs(title = paste(subtype_to_plot),
           color = "Timing",
           fill = "Timing",
           linetype = "Timing") +
      ggplot2::scale_y_continuous(breaks = seq(0, 10, 2),
                         limits = c(0, 10)) +
      ggplot2::theme_bw()+
      ggplot2::theme(legend.position = "bottom") +
        guides(color = guide_legend(override.aes=list(fill=NA,
                                                      alpha = 1)),
               fill = "none") + # Only allow the line color/ type to be shown in the legend.
        scale_color_manual(values = c("#56B4E9", "#E69F00")) + # Blue and Orange
        scale_fill_manual(values = c("#56B4E9", "#E69F00")) # Blue and Orange
      
      if (dose == "SDHD") 
        warning("season_plot() function is not optimized for plotting pre-/post-vaccination titers and SD and HD vaccinees.")
      
    }
  
  if (outcome == "titerincrease") {
    
    if (dose == "SDHD") {
    mdl <- mdl %>%
      dplyr::group_by(vac_short) %>%
      dplyr::arrange(dose) %>%
      dplyr::mutate(dose_label = paste(unique(dose_label), collapse = " | ")) %>%
      dplyr::ungroup()
    
    p <- dataframe %>%
      ggplot2::ggplot(aes(x = distance_norm, 
                 y = titerincrease,
                 color = dose,
                 fill = dose,
                 linetype = dose)) +
      ggplot2::geom_jitter(size = 1, shape = ".", width = 0) +
      ggplot2::geom_line(data = mdl,
                aes(x = pred,
                    y = slope*pred+intercept,
                    color = dose)) +
      ggplot2::geom_ribbon(data = mdl,
                  aes(x = pred,
                      y = pred*slope+intercept,
                      ymax = conf_upr95,
                      ymin = conf_lwr95),
                  alpha = 0.5,
                  color = NA) + 
      ggplot2::facet_grid(rows = vars(vac_short, dose_label), 
                       cols = vars(method)) +
      ggplot2::ylab(latex2exp::TeX("$log_2$(HAI Titer) Increase")) +
      ggplot2::xlab("Distance from Vaccine Strain") +
      ggplot2::labs(title = paste(subtype_to_plot),
           color = "Dose",
           fill = "Dose",
           linetype = "Dose") +
      ggplot2::scale_y_continuous(breaks = seq(-1, 4, 1),
                         limits = c(-1, 4)) +
      ggplot2::theme_bw()+
      ggplot2::geom_hline(yintercept = 0) +
      ggplot2::theme(legend.position = "bottom") +
      ggplot2::guides(color = guide_legend(override.aes=list(fill=NA,
                                                    alpha = 1)),
             fill = "none") + # Only allow the line color/ type to be shown in the legend.
      ggplot2::scale_color_manual(values = c("#000080", "#8b0000")) + # Navy and Dark Red
      ggplot2::scale_fill_manual(values = c("#000080", "#8b0000")) # Navy and Dark Red
} else if (dose == "SD") {
  p <- dataframe %>%
    ggplot2::ggplot(aes(x = distance_norm, 
               y = titerincrease)) +
    ggplot2::geom_jitter(size = 1, 
                shape = ".", 
                width = 0) +
    ggplot2::geom_line(data = mdl,
              aes(x = pred,
                  y = slope*pred+intercept)) +
    ggplot2::geom_ribbon(data = mdl,
                aes(x = pred,
                    y = pred*slope+intercept,
                    ymax = conf_upr95,
                    ymin = conf_lwr95),
                alpha = 0.5,
                color = NA) + 
    ggplot2::facet_grid(rows = vars(vac_short, dose_label), 
               cols = vars(method)) +
    ggplot2::ylab(latex2exp::TeX("$log_2$(HAI Titer) Increase")) +
    ggplot2::xlab("Distance from Vaccine Strain") +
    ggplot2::scale_y_continuous(breaks = seq(-1, 4, 1),
                       limits = c(-1, 4)) +
    ggplot2::theme_bw()+
    ggplot2::geom_hline(yintercept = 0,
               linetype = "dashed")
}
      
  }
  
  if (include_equation_label == TRUE) {
    p <- p + ggplot2::geom_label(data = slopes, 
               aes(x = Inf, y = Inf, label = equation, vjust = 1, hjust = 1), 
               size = 1,
               inherit.aes = FALSE) 
  }
  return(p)
}





## Strain Level Subset Plots##########

strain_plots_subset <- function(dataframe,
                         subtype_to_plot,
                         mdl_to_use,
                         mdl = mdl_list,
                         slopes = slope_list,
                         sample_size = sample_size_list) {
  ## Variables for troubleshooting
  #dataframe <- data
  #subtype_to_plot <- "H1N1"
  #mdl_to_use <- "Strain titerincrease SDHD"
  #mdl <- mdl_list
  #slopes <- slope_list
  #sample_size <- sample_size_list
  
  
  # Create the variables 
  variables <- base::strsplit(mdl_to_use,
                        split = " ")[[1]]
  
  normalization <- variables[1]
  outcome <- variables[2]
  dose <- variables[3]
  
  # Reduce the lists to the important ones
  mdl <- mdl[[mdl_to_use]] %>%
    dplyr::filter(strain_type == subtype_to_plot) %>%
    droplevels()
  
  slopes <- slopes[[mdl_to_use]] %>%
    dplyr::filter(strain_type == subtype_to_plot) %>%
    droplevels()
  
  samples <- sample_size[[mdl_to_use]] %>%
    dplyr::filter(strain_type == subtype_to_plot) %>%
    droplevels()
  
  # Reduce the dataframe to the SD only or the >= 65 age only. 
  
  dataframe2 <- dataframe %>%
    dplyr::filter(distance_type == normalization, #Select the distance normalization
                  strain_type == subtype_to_plot) %>% #Select the subtype to be plotted
    {
      if (dose == "SD") { #If SD need to remove all individuals who did not receive SD
        dplyr::filter(.data = ., dose == "SD") 
      } else if (dose == "SDHD") { #If SDHD need to remove all non elderly individuals
        dplyr::filter(.data = ., age >= 65)
      } else { #Produce a warning if the dose does not match
        stop(paste(dose, ": Not found. Only the SD only individuals all ages, and the SD HD individuals of 65+ years are to be analyzed using this code.")) 
      }
    } %>%
    droplevels() %>%
    dplyr::group_by(subset, dose, vac_short) %>% # Add sample size labels to the corresponding dose for the legends
    dplyr::mutate(dose_label = paste0("(n=", n_distinct(uniq_id), ")")) %>% # Total number of people in the dose 
        dplyr::mutate(.data = .,
                      dose_label = forcats::fct_cross(dose, dose_label, sep = " ")) %>% #Combine it with the dose
    droplevels() %>%
    dplyr::ungroup()
  
  # Create the vector to save the strains to
vaccine_virus <- levels(dataframe2$vac_short) #Pull the vaccine strains to create individual plots for
plt_list <- vector(mode = "list", 
                   length = length(vaccine_virus)) #Create the list for saving the plots to; make it the same length as the vaccine viruses
names(plt_list) <- vaccine_virus #Name the plt_list after the vaccine strains
  
  # Plot the prepost after pivoting and relabeling. Usually this plot is only for the standard dose only comparison. Not the HD/SD. 
  
  if (outcome == "prepost") {
    warning("Prepost outcome is not fully optimized for plotting. Don't be surprised if they are not pretty graphs")
    
    p <-  dataframe %>%
      pivot_longer(cols = c(prevactiter, postvactiter),
                   names_to = "timing",
                   values_to = "titer") %>%
      mutate(timing = factor(timing, levels = c("prevactiter", "postvactiter"), labels = c("Pre-vaccination", "Post-vaccination"))) %>%
      ggplot(aes(x = distance_norm, 
                 y = titer,
                 color = timing,#Distinguish timing by color/fill
                 fill = timing,
                 linetype = dose_label)) + # Distinguish different doses by linetype
      geom_jitter(size = 1, shape = ".", width = 0) + #Spreads out the points for readability. Only apply jitter in the Y-direction since the HAI titers are discrete and we don't want to apply jitter in the x-direction so as to not confuse the reader with regards to what the distance is
      geom_line(data = mdl,
                aes(x = pred,
                    y = slope*pred+intercept,
                    color = timing)) +
      geom_ribbon(data = mdl,
                  aes(x = pred,
                      y = pred*slope+intercept,
                      ymax = conf_upr95,
                      ymin = conf_lwr95,
                      linetype = dose_label),
                  alpha = 0.5,
                  color = NA) + # Include the 95% se confidence intervals of the models. do not include bounding lines
      facet_grid(rows = vars(vac_short), 
                 cols = vars(method)) + #Stratify by vaccine received (H1N1 and H3N2 component) and the distance calculation method
      ylab(latex2exp::TeX("$log_2$(HAI Titer)")) +
      xlab("Distance from Vaccine Strain") +
      labs(title = paste(subtype_to_plot),
           color = "Timing",
           fill = "Timing")
    scale_y_continuous(breaks = seq(0, 10, 2),
                       limits = c(0, 10)) +
      theme_bw()+
      theme(legend.position = "bottom") +
      guides(color = guide_legend(override.aes=list(fill=NA,
                                                    alpha = 1)),
             fill = "none") + # Only allow the line color/ type to be shown in the legend.
      scale_color_manual(values = c("#56B4E9", "#E69F00")) + # Blue and Orange
      scale_fill_manual(values = c("#56B4E9", "#E69F00")) # Blue and Orange
    
    if (dose == "SDHD") 
      warning("season_plot() function is not optimized for plotting pre-/post-vaccination titers and SD and HD vaccinees.")
    
  } # Not optimized for plotting
  
# Plot the titer increase.
  if (outcome == "titerincrease") {
    if (dose == "SDHD") {
    
#For the SDHD it is very cluttered having the vaccine strains and the 3 different measures so create a list that contains each individual strain stratified by subset and method. This is to compare if the different subsets gave similar comparisons within the same strain
    
    for (i in 1:length(vaccine_virus)) {
    plt_list[[i]] <- dataframe2 %>%
      dplyr::filter(vac_short == vaccine_virus[i]) %>%
      ggplot2::ggplot(aes(x = distance_norm, 
                 y = titerincrease,
                 color = dose_label,
                 linetype = dose_label,
                 fill = dose_label)) +
      ggplot2::geom_jitter(size = 1, shape = ".", width = 0) +
      ggplot2::geom_line(data = subset(mdl, vac_short == vaccine_virus[i]),
                aes(x = pred,
                  y = slope*pred+intercept,
                    linetype = dose_label)) +
      ggplot2::geom_ribbon(data = subset(mdl, vac_short == vaccine_virus[i]),
                  aes(x = pred,
                      y = pred*slope+intercept,
                      ymax = conf_upr95,
                      ymin = conf_lwr95),
                  alpha = 0.5,
                  color = NA) + 
      ggplot2::facet_grid(rows = vars(subset), 
                 cols = vars(method)) +
      ggplot2::ylab(latex2exp::TeX("$log_2$(HAI Titer) Increase")) +
      ggplot2::xlab("Distance from Vaccine Strain") +
      ggplot2::labs(title = paste(vaccine_virus[i]),
           color = "Dose",
           fill = "Dose",
           linetype = "Dose") +
      ggplot2::scale_y_continuous(breaks = seq(-1, 4, 1),
                         limits = c(-1, 4)) +
      ggplot2::theme_bw()+
      ggplot2::geom_hline(yintercept = 0) +
      ggplot2::theme(legend.position = "bottom") +
      ggplot2::guides(color = guide_legend(override.aes=list(fill=NA,
                                                    alpha = 1)),
             fill = "none") + # Only allow the line color/ type to be shown in the legend.
      ggplot2::scale_color_manual(values = c("#000080", "#8b0000")) + # Navy and Dark Red
      ggplot2::scale_fill_manual(values = c("#000080", "#8b0000")) # Navy and Dark Red
    } # Loop through vaccine strains to create individual plots
      return(plt_list)
     #Creates Method ~ Subset, with 2 doses on each plot comparing SDHD results. Creates a list with each component being a vaccine strain
    } else 
      if (dose == "SD") {
      p <- dataframe2 %>%
        ggplot2::ggplot(aes(x = distance_norm, 
                            y = titerincrease,
                            color = subset,
                            fill = subset)) +
        ggplot2::geom_jitter(size = 1, shape = ".", width = 0) +
        ggplot2::geom_line(data = mdl,
                  aes(x = pred,
                      y = slope*pred+intercept)) +
        ggplot2::geom_ribbon(data = mdl,
                    aes(x = pred,
                        y = pred*slope+intercept,
                        ymax = conf_upr95,
                        ymin = conf_lwr95),
                    alpha = 0.5,
                    color = NA) + 
        ggplot2::facet_grid(rows = vars(vac_short), 
                   cols = vars(method)) +
        ggplot2::ylab(latex2exp::TeX("$log_2$(HAI Titer) Increase")) +
        ggplot2::xlab("Distance from Vaccine Strain") +
        ggplot2::labs(title = paste(subtype_to_plot),
             color = "Subset",
             fill = "Subset",
             linetype = "Subset") +
        ggplot2::scale_y_continuous(breaks = seq(-1, 4, 1),
                           limits = c(-1, 4)) +
        ggplot2::theme_bw()+
        ggplot2::geom_hline(yintercept = 0) +
        ggplot2::theme(legend.position = "bottom") +
        ggplot2::guides(color = guide_legend(override.aes=list(fill=NA,
                                                      alpha = 1)),
               fill = "none")# + # Only allow the line color/ type to be shown in the legend.
        #scale_color_manual(values = c("#56B4E9", "#E69F00")) + # Blue and Orange
        #scale_fill_manual(values = c("#56B4E9", "#E69F00"))
      
      return(p)
    } #Creates Method ~ Vac_short, with 5 subsets on each plot comparing SD results
  }
}

## Antibody Landscape Plots ######################

antibody_landscapes <- function(dataframe = data,
                                norms = c("Strain", "Season", "Vaccine"),
                                doses = c("SD", "SDHD"),
                                outcomes = c("prepost", "titerincrease"),
                                strains = c("H1N1", "H3N2"),
                                seasons = 2014:2019,
                                save_plots = TRUE,
                                save_files = TRUE,
                                save_with_labels = FALSE) {
  
  #Variables for troubleshooting:
  #dataframe <- data
  #norms <- c("Strain", "Season", "Vaccine")
  #doses <- c("SD", "SDHD")
  #outcomes <- c("prepost", "titerincrease")
  #strains <- c("H1N1", "H3N2")
  #seasons <- 2014:2019
  #save_plots <- TRUE
  #save_files <- FALSE
  #save_with_labels <- FALSE
  
  
  ## Create the linear regression models: 
  mdl_list <- vector(mode = "list")
  for(i in norms) {
    for (j in outcomes) {
      for (k in doses) {
        x <- paste(i, j, k)
        mdl_list[[x]] <- linear_model(dataframe = data,
                                      normalization = i,
                                      outcome = j,
                                      dose = k)
      }
    }
  }
  
  print("Models are completed")
  
  ## Determine samples sizes for the plots:
  sample_size_list <- vector(mode = "list")
  
  for (i in 1:length(mdl_list)) {
    x <- names(mdl_list[i])
    sample_size_list[[x]] <- sample_size_fun(dataframe = mdl_list[[i]])
  }
  
  print("Sample sizes are completed")
  
  ##Determine the linear regression equations for the plots.
  
  slope_list <- vector(mode = "list")
  
  for (i in 1:length(mdl_list)) {
    x <- names(mdl_list[i])
    slope_list[[x]] <- equation_labels(dataframe = mdl_list[[i]])
  }
  
  print("Slopes are completed")
  
  if (save_files == TRUE) {
    saveRDS(mdl_list,
            file = here::here("data/processed/linear_regression_models/mdl_list.rds"))
    
    saveRDS(sample_size_list,
            file = here::here("data/processed/linear_regression_models/size_list.rds"))
    
    saveRDS(slope_list,
            file = here::here("data/processed/linear_regression_models/slope_list.rds"))
    print("Saved Lists")
  }  
  
  
  ## Strain plots
  if ("Strain" %in% norms) {
    for (j in doses) {
      for(k in outcomes) {
        
        plt_list <- vector(mode = "list")
        
        for (i in c("H1N1", "H3N2")) {
          
          plt_list[[i]] <- subtype_plots(dataframe = data,
                                         subtype_to_plot = i,
                                         mdl_to_use = paste("Strain", k, j),
                                         mdl = mdl_list,
                                         slopes = slope_list,
                                         sample_size = sample_size_list,
                                         include_equation_label = FALSE)
          
        }
        
        if (save_plots == TRUE) {
          
          cowplot::plot_grid(plotlist = plt_list,
                             ncol = 1,
                             rel_heights = c(1,1.5),
                             labels = "AUTO"
          ) %>%
            cowplot::ggsave2(filename = here::here(paste0("figures/supplemental/typeA_", tolower(j), "_", k, "_strain.png")), 
                             width = 4,
                             height = 10,
                             scale = 1.5)
          print("Saved plots")
          
        }
        
        
        print(paste("Figures and files for Strain", k, j, "are completed"))
      }
      print("am i stopping here?")
    }
    print(paste(j, "dose is completed"))
  }
  print("strain plots are completed")
  
  ## Season plots
  if (sum(norms %in% c("Vaccine", "Season")) != 0) {
    for (k in norms[which(norms %in% c("Vaccine", "Season"))]) {
      for (j in doses) {
        for(l in outcomes) {
          
          overall_model <- paste(k, l, j)
          plt_list <- vector(mode = "list")
          
          for (i in seasons) {
            x <- paste0("season_", i)
            
            plt_list[[x]] <- season_plots(dataframe = data,
                                          season_to_plot = i,
                                          mdl_to_use = overall_model,
                                          mdl = mdl_list,
                                          slopes = slope_list,
                                          sample_size = sample_size_list,
                                          include_equation_label = FALSE)
          }
          
          if (save_plots == TRUE) {
            
            cowplot::plot_grid(plotlist = plt_list, 
                               ncol = 2,
                               labels = "AUTO"
            ) %>%
              cowplot::ggsave2(filename = here::here(paste0("figures/supplemental/typeA_", tolower(j), "_", l, "_", tolower(k), ".png")), 
                               width = 12,
                               height = 15)
            
            if (save_with_labels == TRUE & j == "SD" & k == "Season" & l == "titerincrease") {
              
              labels <- data %>%
                dplyr::select(season, vac_short, strain_short, method, distance_norm, distance_type) %>%
                dplyr::filter(season == 2017,
                       distance_type == k,
                       strain_short %in% c("NC/99", "SC/18", "Shan/93", "TX/77")) %>%
                distinct()
              
              plt_list[[4]] <- plt_list[[4]] +
                ggplot2::geom_label_repel(data = labels,
                                 aes(x = distance_norm,
                                     y = 0,
                                     label = strain_short),
                                 color = "black",
                                 inherit.aes = FALSE,
                                 size = 2,
                                 direction = "y",
                                 box.padding = 0.1) 
            }
            
            cowplot::ggsave2(filename = here::here(paste0("figures/manuscript/typeA_", tolower(j), "_inc_season17_", tolower(k), ".png")),
                             plot = plt_list[[4]],
                             width = 6,
                             height = 5)
            
          } 
          
          
        }
      }
    }
  }
}
