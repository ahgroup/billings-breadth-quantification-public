## Two Antigenic Unit Cutoff calculations with normalized distances #########

cutoff_calc <- function(dataframe,
                        normalization, # String. Normalization of distances to be used
                        dose_of_interest) {
  df <- dataframe %>%
    dplyr::filter(method == "Cartography", #Need
                  distance_type == normalization,
                  dose == dose_of_interest) %>%
    dplyr::select(season, vac_short, strain_type, distance_type, distance) %>%
    {
      if(normalization == "Strain") dplyr::group_by(.data = ., distance_type, strain_type, vac_short) else
        dplyr::group_by(.data = ., distance_type, season, strain_type, vac_short)
    } %>%
    dplyr::summarise(cutoff = 2/max(distance)) %>%
    ungroup()

  return(df)
}


## Compare the two Dose tables together ###########

dose_compare <- function(dataframe,
                         normalization)   {

sd <- cutoff_calc(data,
                    normalization = normalization,
                    dose_of_interest = "SD")

hd <- cutoff_calc(data,
                  normalization = normalization,
                  dose_of_interest = "HD")

if (identical(sd, hd) == TRUE) {
  print("The SD and HD cutoff dataframes are identical.")
  return(sd)
} else {
  stop("The SD and HD cutoff dataframes are NOT identical. ")
}

}

## Figure Labels ######################

fig_label <- function(text, region="figure", pos="topleft", cex=NULL, ...) {

  region <- match.arg(region, c("figure", "plot", "device"))
  pos <- match.arg(pos, c("topleft", "top", "topright",
                          "left", "center", "right",
                          "bottomleft", "bottom", "bottomright"))

  if(region %in% c("figure", "device")) {
    ds <- dev.size("in")
    # xy coordinates of device corners in user coordinates
    x <- graphics::grconvertX(c(0, ds[1]), from="in", to="user")
    y <- graphics::grconvertY(c(0, ds[2]), from="in", to="user")

    # fragment of the device we use to plot
    if(region == "figure") {
      # account for the fragment of the device that
      # the figure is using
      fig <- par("fig")
      dx <- (x[2] - x[1])
      dy <- (y[2] - y[1])
      x <- x[1] + dx * fig[1:2]
      y <- y[1] + dy * fig[3:4]
    }
  }

  # much simpler if in plotting region
  if(region == "plot") {
    u <- par("usr")
    x <- u[1:2]
    y <- u[3:4]
  }

  sw <- graphics::strwidth(text, cex=cex) * 60/100
  sh <- graphics::strheight(text, cex=cex) * 60/100

  x1 <- switch(pos,
               topleft     =x[1] + sw,
               left        =x[1] + sw,
               bottomleft  =x[1] + sw,
               top         =(x[1] + x[2])/2,
               center      =(x[1] + x[2])/2,
               bottom      =(x[1] + x[2])/2,
               topright    =x[2] - sw,
               right       =x[2] - sw,
               bottomright =x[2] - sw)

  y1 <- switch(pos,
               topleft     =y[2] - sh,
               top         =y[2] - sh,
               topright    =y[2] - sh,
               left        =(y[1] + y[2])/2,
               center      =(y[1] + y[2])/2,
               right       =(y[1] + y[2])/2,
               bottomleft  =y[1] + sh,
               bottom      =y[1] + sh,
               bottomright =y[1] + sh)

  old.par <- par(xpd=NA)
  on.exit(par(old.par))

  text(x1, y1, text, cex=cex, ...)
  return(invisible(c(x,y)))
}



## Weighting Example Plot #####################

weighting_example <- function(file = "figures/manuscript/example_weighting.png") {

png(here::here(file), width = 1280, height = 480, res= 300, pointsize =6)

par(mfrow=c(1,3))

plot(x = c(0, 1),
     y = c(1, 1),
     type = "l",
     main = "Unweighted",
     xlab = "Normalized Distance",
     ylab = "Weighting applied to fitted values",
     ylim = c(0, 1))
polygon(x = c(0, 1, 1, 0),
        y = c(1, 1, 0, 0),
        col = "grey")
fig_label("A", cex=2)

plot(x = c(0, 1),
     y = c(1, 0),
     type = "l",
     main = "Linear",
     xlab = "Normalized Distance",
     ylab = "")
polygon(x = c(0, 1, 0),
        y = c(1, 0, 0),
        col = "grey")
fig_label("B", cex=2)
plot(x = c(0, 0.5, 0.5, 1),
     y = c(1, 1, 0, 0),
     type = "l",
     main = "Antigenic Distance",
     xlab = "Normalized Distance",
     ylab = "")
polygon(x = c(0, 0.5, 0.5, 0),
        y = c(1, 1, 0, 0),
        col = "grey")
fig_label("C", cex=2)

dev.off()
}

## Area under the curve Calculations ############

#Split the data into appropriate groupings for linear regressions
#Create Normalized distance by season, vaccine, and method

auc_calc <- function(dataframe,
                     cutoff_dataframe = readRDS(file = here::here("results/tables/AU_cutoff.rds")),
                     normalization = c("Season", "Strain", "Vaccine"),
                     dose_of_interest = c("SD", "SDHD"),
                     example = FALSE,
                     season_of_interest = NULL,
                     strain_of_interest = NULL,
                     gradation = TRUE) {

# For Troubleshooting:
#  dataframe <- data
#  cutoff_dataframe = readRDS(file = here::here("tables/supplemental/AU_cutoff.rds"))
#  normalization = "Strain"
#  dose_of_interest <- "SD"
#  example = TRUE
#  season_of_interest = NULL
#  strain_of_interest <- "H1N1:CA/09"
#  gradation = TRUE
#Set warnings for arguments
  if (normalization %in% c("Season", "Strain", "Vaccine") == FALSE)
    stop(paste(normalization, ":normalization method does not exist. Please use Season, Strain, or Vaccine."))
  if (dose_of_interest %in% c("SD", "SDHD") == FALSE)
    stop(paste(dose, ": Dose does not exist. Set dose to either: SD for all individuals with standard dose vaccination or SDHD for all indiviudals equal to or greater than 65 years of age, for both the SD and HD vaccines"))
  #if(example == TRUE & normalization == "Strain")
  #  stop("Plots are by Season. Change normalization parameter to either Season or Vaccine.")

  df <- dataframe %>%
    {
      if (dose_of_interest == "SD") #Filter data for SD: all individuals who received SD only all ages
        dplyr::filter(.data = .,
               dose == "SD") else
                 dplyr::filter(.data = ., #Filter data for SDHD: 65+ years individuals who received SD or HD
                        age >= 65)
    } %>%
    dplyr::filter(distance_type == normalization) %>% # Select the normalization of interest
    dplyr::select(season, strain_type, vac_short, method, dose, titerincrease, distance_norm) #Select the necessary columns


    split <- df %>%
      {
        if (normalization == "Strain") dplyr::group_by(.data = .,
                                                       strain_type, vac_short, method, dose) # Don't group_by season
        else dplyr::group_by(.data = .,
                             season, strain_type, vac_short, method, dose) # Do include season
       } %>%
      dplyr::group_split() #Split it into different dataframes for linear regression/predictions


  final <- data.frame(matrix(ncol = 11)) # Create the final dataframe to be filled

  colnames(final) <- c("distance_norm", "pred", "w_anti", "w_linear", "season", "vac_short", "method", "strain_type", "cutoff", "distance_type", "dose") #Name the dataframe

  # Combine datasets with AU Cut offs
  # Run initial linear regression
  # Create prediction data from fits
  # Weight the fitted data by the different schemes
  # Result: Data frame that contains all season, vaccine, normalized distances, and weighted fitted titers

  for (i in 1:length(split)) {
    if (normalization == "Strain") {
      limit <- split[[i]] %>%
        dplyr::select(vac_short, method, strain_type, dose) %>% #Retain all columns that match the AU cutoff df.
        unique() %>%
        dplyr::left_join(subset(cutoff_dataframe, distance_type == normalization)) #Include only the AU cutoff rows that match the normalization of interest
    } else { # Includes the season variable
      limit <- split[[i]] %>%
        dplyr::select(season, vac_short, method, strain_type, dose) %>%
        unique() %>%
        dplyr::left_join(subset(cutoff_dataframe, distance_type == normalization))
    }

    fit <-  stats::lm(data = split[[i]], formula = titerincrease ~ distance_norm) #Fit the linear regressions for predictions

    predict_data <- data.frame(distance_norm = c(seq(from = 0,
                                                     to = max(split[[i]]$distance_norm),
                                                     length.out = 100))) # Create the new data frame from 0 to the maximum x value available. For season and strain normalizations this should be 1, but for the vaccine the maximum will vary.

    predict_data$pred <- stats::predict.lm(fit,
                                    newdata = predict_data) #Predict the data using the fit and dataframe
    final <- predict_data %>%
      dplyr::mutate(w_anti = ifelse(distance_norm < limit$cutoff,
                                    pred*1,
                                    pred*0), # Weighting scheme of AU cutoff
             w_linear = pred*(-1*distance_norm+1)) %>% # Weighting scheme of linear decline with increasing distance
      cbind(limit) %>% #Combine with the AU cutoff value
      rbind(final) %>% #Append it to the previous final
      dplyr::filter(!is.na(distance_norm)) #Remove the NA from the first row
  }

  help <- final %>%
    {
      if (normalization == "Strain") #Set the group_by settings for finding the AUC
       {dplyr::group_by(.data = .,
                 strain_type, vac_short, method, dose) } else {
                   dplyr::group_by(.data = .,
                            season, strain_type, vac_short, method, dose) }
    } %>%
    dplyr::arrange(distance_norm) %>% #These are trapezoidal AUC, need to arrange in ascending order the distance
    dplyr::summarize(auc_pred = pracma::trapz(x = distance_norm, #Unweighted AUC
                               y = pred),
              auc_anti = pracma::trapz(x = distance_norm, #Antigenic Unit AUC
                               y = w_anti),
              auc_linear = pracma::trapz(x = distance_norm, #Linear AUC
                                 y = w_linear)) %>%
    dplyr::mutate(dplyr::across(tidyselect:::where(is.numeric), round, 3)) %>% #Round all the numbers to 3 decimal places
    dplyr::ungroup()


  if (example == FALSE) { #Will produce the dataframe with AUCs
    if (dose_of_interest == "SD") {
      help <- help %>%
        dplyr::select(-dose)
    }
    return(help)

  } else if (example == TRUE) { #Produces plots of the AUC

    if (normalization == "Strain") {

      labels <- help %>%
        dplyr::filter(vac_short == strain_of_interest) %>% #Select strain of interest
        tidyr::pivot_longer(cols = starts_with("auc"),
                            names_to = "weighting",
                            values_to = "auc") %>%
        dplyr::mutate(weighting = factor(weighting,
                                         levels = c("auc_pred", "auc_linear", "auc_anti"),
                                         labels = c("Unweighted", "Linear", "2 Antigenic Units")),
                      label = paste("Area: ", auc)) %>%
        dplyr::select(vac_short, method, dose, weighting, label)

      if (gradation == TRUE) {
      #Prepare the dataframe
      final_2 <- final %>%
        dplyr::filter(vac_short == strain_of_interest) %>%
        dplyr::mutate(unweighted = "Unweighted",
                      linear = "Linear",
                      au2 = "2 Antigenic Units") %>%
        tidyr::pivot_longer(cols = c(unweighted:au2),
                            names_to = "scheme",
                            values_to = "weighting") %>%
        dplyr::select(-scheme) %>%
        dplyr::rename(titerincrease = pred) %>%
        dplyr::mutate(gradation = ifelse(weighting == "Unweighted",
                                         1,
                                         ifelse(weighting == "2 Antigenic Units" & distance_norm < cutoff,
                                                1,
                                                ifelse(weighting == "2 Antigenic Units" & distance_norm >= cutoff,
                                                       0,
                                                       1-distance_norm))),
                      weighting = factor(weighting, levels = c("Unweighted", "Linear", "2 Antigenic Units")))

      #Make the plot

      p <- final_2 %>%
        ggplot2::ggplot(aes(x = distance_norm,
                   y = titerincrease)) +
        ggplot2::geom_segment(aes(xend = distance_norm,
                         yend = 0,
                         color = gradation),
                     size =2)+
        ggplot2::geom_line() +
        ggplot2::geom_jitter(data = subset(df, vac_short == strain_of_interest),
                    aes(x = distance_norm,
                        y = titerincrease),
                    size = 1,
                    shape = ".",
                    width = 0,
                    color = "black") + #Spreads out the points for readability. Only apply jitter in the Y-direction since the HAI titers are discrete and we don't want to apply jitter in the x-direction so as to not confuse the reader with regards to what the distance is
        ggplot2::geom_hline(yintercept = 0,
                   linetype = "dashed") +
        ggplot2::theme_bw() +
        ggplot2::facet_grid(rows = vars(weighting),
                   cols = vars(method)) +
        ggplot2::geom_label(data = labels,
                   aes(x = Inf,
                       y = Inf,
                       label = label,
                       hjust = 1,
                       vjust = 1),
                   size = 3) +
        ggplot2::xlab("Normalized Distance from CA/09 H1N1") +
        ggplot2::ylab(latex2exp::TeX("$log_2$(HAI Titer) Increase")) +
        ggplot2::scale_color_gradient(low = "white",
                             high = "blue") +
        ggplot2::labs(color = "Weight") +
        ggplot2::theme(legend.position = "bottom")  +
        ggplot2::ylim(c(-1, 4))
    } else {
      #Prepare the dataframe
      final_2 <- final %>%
        dplyr::filter(vac_short == strain_of_interest) %>%
        tidyr::pivot_longer(cols = pred:w_linear,
                            names_to = "scheme",
                            values_to = "titerincrease") %>%
        dplyr::select(-scheme) %>%
        dplyr::rename(titerincrease = pred) %>%
        dplyr::mutate(gradation = ifelse(weighting == "Unweighted",
                                         1,
                                         ifelse(weighting == "2 Antigenic Units" & distance_norm < cutoff,
                                                1,
                                                ifelse(weighting == "2 Antigenic Units" & distance_norm >= cutoff,
                                                       0,
                                                       1-distance_norm))),
                      weighting = factor(weighting, levels = c("Unweighted", "Linear", "2 Antigenic Units")))

      #Make the plot

      p <- final_2 %>%
        ggplot2::ggplot(aes(x = distance_norm,
                   y = titerincrease)) +
        ggplot2::geom_segment(aes(xend = distance_norm,
                         yend = 0,
                         color = gradation),
                     size =2)+
        ggplot2::geom_line() +
        ggplot2::geom_jitter(data = subset(df, vac_short == strain_of_interest),
                    aes(x = distance_norm,
                        y = titerincrease),
                    size = 1,
                    shape = ".",
                    width = 0,
                    color = "black") + #Spreads out the points for readability. Only apply jitter in the Y-direction since the HAI titers are discrete and we don't want to apply jitter in the x-direction so as to not confuse the reader with regards to what the distance is
        ggplot2::geom_hline(yintercept = 0) +
        ggplot2::theme_bw() +
        ggplot2::facet_grid(rows = vars(weighting),
                   cols = vars(method)) +
        ggplot2::geom_label(data = labels,
                   aes(x = Inf,
                       y = Inf,
                       label = label,
                       hjust = 1,
                       vjust = 1),
                   size = 3) +
        ggplot2::xlab("Normalized Distance from CA/09 H1N1") +
        ggplot2::ylab(latex2exp::TeX("$log_2$(HAI Titer) Increase")) +
        ggplot2::scale_color_gradient(low = "white",
                             high = "black") +
        ggplot2::labs(color = "Weight") +
        ggplot2::theme(legend.position = "bottom")  +
        ggplot2::ylim(c(-1, 4))
    }


    } else if (normalization %in% c("Season", "Vaccine")) {

    #Create label for the plots. These include the values of the AUC shown in the corner
    labels <- help %>%
      dplyr::filter(season == season_of_interest) %>% #Select season of interest
      dplyr::mutate(label = paste0("Unweighted: ",
                            auc_pred,
                            "\nLinear: ",
                            auc_linear,
                            "\nAU: ",
                            auc_anti)) %>%
      dplyr::select(season, dose, strain_type, vac_short, method, label)

    #Create the plots
    p <- final %>%
      dplyr::filter(season == season_of_interest) %>%
      tidyr::pivot_longer(cols = c(w_anti, w_linear, pred),
                   names_to = "weighting",
                   values_to = "titerincrease") %>%
      dplyr::mutate(weighting = factor(weighting, levels = c("pred", "w_linear", "w_anti"), labels = c("Unweighted", "Linear", "2 Antigenic Units"))) %>%
      ggplot2::ggplot(aes(x = distance_norm,
                 y = titerincrease,
                 color = dose)) +
      ggplot2::geom_line(aes(linetype = weighting)) +
      ggplot2::geom_hline(yintercept = 0) +
      ggplot2::theme_bw() +
      ggplot2::facet_grid(rows = vars(vac_short),
                 cols = vars(method)) +
      ggplot2::geom_label(data = subset(labels, dose == "SD"),
                 aes(x = Inf,
                     y = Inf,
                     label = label,
                     hjust = 1,
                     vjust = 1),
                 size = 2) +
      ggplot2::labs(title = paste0("Application of Weighting to the ",
                          season_of_interest,
                          " Season"),
           subtitle = paste0(dose_of_interest, " individuals\n", normalization, "-based Normalized Distance"),
           linetype = "Weighting",
           color = "Dose")+
      ggplot2::xlab("Normalized Distance") +
      ggplot2::ylab(latex2exp::TeX("$log_2$(HAI Titer) Increase")) +
      ggplot2::theme(legend.position = "bottom")

    if (dose_of_interest == "SD") {
      p <- p +
        ggplot2::scale_color_manual(values = "black") +
        ggplot2::guides(color = "none")
    } else {
      p <-  p +
        ggplot2::scale_color_manual(values = c("black", "red")) +
        ggplot2::geom_label(data = subset(labels, dose == "HD"),
                   aes(x = -Inf,
                       y = Inf,
                       label = label,
                       hjust = 0,
                       vjust = 1),
                   size = 2)
    }
    }
    return(p)

  }
}

## Titer Increase Calculator ####
titer_increase_cal <- function(dataframe,
                               dose_of_interest,
                               normalization) {
  df <- dataframe %>%
    {
      if (dose_of_interest == "SD") #Filter data for SD: all individuals who received SD only all ages
        dplyr::filter(.data = .,
                      dose == "SD") else
                        dplyr::filter(.data = ., #Filter data for SDHD: 65+ years individuals who received SD or HD
                                      age >= 65)
    } %>%
    dplyr::filter(distance_type == normalization) %>% # Select the normalization of interest
    dplyr::select(subset, season, strain_type, vac_short, method, dose, titerincrease) #Select the necessary columns


  split <- df %>%
    {
      if (normalization == "Strain") dplyr::group_by(.data = .,
                                                     subset, strain_type, vac_short, method, dose) # Don't group_by season
      else dplyr::group_by(.data = .,
                           subset, season, strain_type, vac_short, method, dose) # Do include season
    } %>%
    summarize(mean_ti = mean(titerincrease)) %>%
    {
      if (dose_of_interest == "SD") dplyr::select(.data = .,
                                                  -dose)
      else .
    }

  return(split)
}
