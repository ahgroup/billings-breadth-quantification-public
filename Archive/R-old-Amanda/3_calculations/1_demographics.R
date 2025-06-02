# Purpose of this code:

# This code is to visualize the basic demographics of the cohort

# Prerequisties:

# Need to run:

#* ../code/1_cleaning/1_cleaning_code.R


# Products of this code:

## Tables

#There are 4 tables produced here. They are all saved as RDS files to the figures/demographics folder. The first is information of the vaccine components, and the last two are of the study cohort demographics.


#Load the packages:

library(tidyverse) #data wrangling
library(Racmacs)
library(here) #locating files
source(here::here("R/functions/antigenic-cartography.R"))


#Load the data:
d <- readRDS(file = here::here("data/processed/clean_data.rds")) %>%
  dplyr::select(uniq_id, season, age, dose) %>%
  dplyr::distinct()

#Individual that received HD but is actually filtered out. 
d %>%
  dplyr::filter(dose == "HD", age == 61) %>%
  dplyr::pull(uniq_id)

### Independent cohort factors:

# >= 65 ####
#For all individuals greater than or equal to 65 years of age

elderly <- d %>%
  dplyr::filter(age >= 65)

age <- elderly %>%
  dplyr::group_by(season, dose) %>%
  dplyr::summarize(age_med = median(age),
            age_low = min(age),
            age_high = max(age),
            n = n()) %>%
  dplyr::mutate(column = paste(season, dose, sep = "_")) %>%
  dplyr::ungroup() %>%
  dplyr::select(-c(season, dose)) %>%
  dplyr::mutate(age_range = paste(age_low, "-", age_high, sep = ""))

demo_hd <- age %>%
  dplyr::select(column, n) %>%
  tidyr::pivot_wider(names_from = column, values_from = n) %>%
  base::rbind(tidyr::pivot_wider(data = base::subset(age, 
                                         select = c(column, age_med)), 
                           names_from = column, 
                           values_from = age_med)) %>%
  base::rbind(tidyr::pivot_wider(data = base::subset(age, 
                                         select = c(column, age_range)), 
                           names_from = column, 
                           values_from = age_range))

rownames(demo_hd) <- c("Total Individuals", "Median Age (Years)", "Age Range")
#Rownames on tibbles is deprecated but it works, for now

saveRDS(demo_hd, file = here::here("results/tables/cohort_demographics_elderly.rds"))
  
## SD ####
## For all individuals who received the standard dose vaccination, regardless of age

d <- d %>% 
  dplyr::filter(dose == "SD") %>%
  dplyr::select(-dose)

age <- d %>%
  dplyr::group_by(season) %>%
  dplyr::summarize(age_med = median(age),
            age_low = min(age),
            age_high = max(age),
            n = n()) %>%
  dplyr::mutate(column = paste(season, sep = "_"),
         .keep = "unused") %>%
  dplyr::ungroup() %>%
  dplyr::mutate(age_range = paste(age_low, "-", age_high, sep = ""))

demo <- age %>%
  dplyr::select(column, n) %>%
  tidyr::pivot_wider(names_from = column, values_from = n) %>%
  base::rbind(tidyr::pivot_wider(data = base::subset(age, 
                                               select = c(column, age_med)), 
                                 names_from = column, 
                                 values_from = age_med)) %>%
  base::rbind(tidyr::pivot_wider(data = base::subset(age, 
                                               select = c(column, age_range)), 
                                 names_from = column, 
                                 values_from = age_range))

rownames(demo) <- c("Total Individuals", "Median Age (Years)", "Age Range")
#Rownames on tibbles is deprecated but it works, for now

  working_save <- demo %>%
    tibble::rownames_to_column("rowname") %>%
    tibble::column_to_rownames("rowname")
  
  saveRDS(demo, file = here::here("results/tables/cohort_demographics_standard_dose.rds"))

# Analysis ####

## Data ####

overall <- readRDS(file = here::here("data/processed/h1_seq_distance/h1_distance_measure_key.rds")) %>%
    dplyr::filter(method %in% c("cart_2d_pre", "cart_2d_post")) %>%
    dplyr::mutate(season = "overall",
         sera_count = ifelse(method == "cart_2d_pre", #Values from sample size section below
                             1374,
                             1526)) %>%
    dplyr::rename(identifier = strains_fullname) %>%
    dplyr::select(identifier, method, season, sera_count, h1n1_vaccine_fullname, distance)

h1_distance <- readRDS(file = here::here("data/processed/h1_seq_distance/h1_seasons_distance.rds")) %>%
  dplyr::group_by(season, 
           method) %>%
  dplyr::mutate(distance = sqrt((V1-V1[which(identifier==paste(h1n1_vaccine_fullname))])^2 +
                           (V2-V2[which(identifier==paste(h1n1_vaccine_fullname))])^2))  %>%
  dplyr::ungroup() %>%
  dplyr::select(-c(V1:V2)) %>%
  base::rbind(overall)

overall <- readRDS(file = here::here("data/processed/h3_seq_distance/h3_distance_measure_key.rds")) %>%
  dplyr::filter(method %in% c("cart_2d_pre", "cart_2d_post")) %>%
  dplyr::mutate(season = "overall",
         sera_count = ifelse(method == "cart_2d_pre", #Values from sample size section below
                             1831,
                             1911)) %>%
  dplyr::rename(identifier = strains_fullname) %>%
  dplyr::select(identifier, method, season, sera_count, h3n2_vaccine_fullname, distance)

h3_distance <- readRDS(file = here::here("data/processed/h3_seq_distance/h3_seasons_distance.rds")) %>%
  dplyr::group_by(season, 
           method) %>%
  dplyr::mutate(distance = sqrt((V1-V1[which(identifier==paste(h3n2_vaccine_fullname))])^2 +
                           (V2-V2[which(identifier==paste(h3n2_vaccine_fullname))])^2))  %>%
  dplyr::ungroup() %>%
  dplyr::select(-c(V1:V2)) %>%
  base::rbind(overall)

### H1N1 ####

plt <- h1_distance %>%
  tidyr::pivot_wider(names_from = method, 
              values_from = c(distance, sera_count)) %>%
  dplyr::mutate(season = paste0(season,
                         "; Pre: ",
                         sera_count_cart_2d_pre,
                         "; Post: ",
                         sera_count_cart_2d_post),
         .keep = "unused") %>%
  ggplot2::ggplot(aes(x = distance_cart_2d_pre,
             y = distance_cart_2d_post)) +
  ggplot2::geom_point(aes(color = h1n1_vaccine_fullname)) +
  ggplot2::geom_smooth(method = "lm") +
  ggplot2::facet_wrap(facets = vars(season),
             scales = "free") +
  ggplot2::labs(title = "H1N1 relative cartographic distances of pre- and post-vaccination by season",
       subtitle = "Facets: Season; Pre-vaccination sera sample size; Post-vaccination sera sample size",
       caption = "Blue line = linear regression with 95% confidence intervals\nSera used was standard dose individuals only\nR = Pearson's correlation coefficient\nBlack line = Line of identity (y=x)",
       color = "Vaccine") +
  ggplot2::xlab("Pre-vaccination Distances") +
  ggplot2::ylab("Post-vaccination Distances") +
  ggpubr::stat_cor(size = 2) +
  ggplot2::theme_bw() +
  ggplot2::theme(legend.position = "bottom") +
  ggplot2::geom_abline(intercept = 0,
              slope = 1,
              color = "black") +
  ggplot2::ylim(c(0,6)) +
  ggplot2::xlim(c(0,6))

ggplot2::ggsave(plt, 
         filename = here::here("results/figures/h1_prepost_mapDistance.png"),
         width = 7,
         height = 9.5)


### H3N2 ####


plt <- h3_distance %>%
  tidyr::pivot_wider(names_from = method, 
              values_from = c(distance, sera_count)) %>%
  dplyr::mutate(season = paste0(season,
                         "; Pre: ",
                         sera_count_cart_2d_pre,
                         "; Post: ",
                         sera_count_cart_2d_post),
         .keep = "unused") %>%
  ggplot2::ggplot(aes(x = distance_cart_2d_pre,
             y = distance_cart_2d_post)) +
  ggplot2::geom_point(aes(color = h3n2_vaccine_fullname)) +
  ggplot2::geom_smooth(method = "lm") +
  ggplot2::facet_wrap(facets = vars(season)) +
  ggplot2::labs(title = "H3N2 relative cartographic distances of pre- and post-vaccination by season",
       subtitle = "Facets: Season; Pre-vaccination sera sample size; Post-vaccination sera sample size",
       caption = "Blue line = linear regression with 95% confidence intervals\nSera used was standard dose individuals only\nR = Pearson's correlation coefficient\nBlack line = Line of identity (y=x)",
       color = "Vaccine") +
  ggplot2::xlab("Pre-vaccination Distances") +
  ggplot2::ylab("Post-vaccination Distances") +
  ggpubr::stat_cor(size = 2) +
  ggplot2::theme_bw()+
  ggplot2::theme(legend.position = "bottom") +
  ggplot2::guides(color = guide_legend(nrow = 2)) +
  ggplot2::geom_abline(intercept = 0,
              slope = 1,
              color = "black") +
  ggplot2::ylim(c(0,6.5)) +
  ggplot2::xlim(c(0,6.5))

  ggplot2::ggsave(plt, 
         filename = here::here("results/figures/h3_prepost_mapDistance.png"),
         width = 7,
         height = 9.5)

## Sample Size tables ####
  
overall_sample <- AcmapToDF(map = Racmacs::read.acmap(filename = here::here("data/processed/h3_maps/h3_pre_sd_2d.ace"))) %>%
  dplyr::mutate(sera_count_cart_2d_pre = sum(type == "sera"),
         virus_count_cart_2d_pre = sum(type == "antigen"),
         season = "Overall",
         subtype = "H3N2")  %>%
    dplyr::select(season, subtype, virus_count_cart_2d_pre, sera_count_cart_2d_pre) %>%
    dplyr::distinct()

overall_sample <- AcmapToDF(map = Racmacs::read.acmap(filename = here::here("data/processed/h3_maps/h3_post_sd_2d.ace"))) %>%
  dplyr::mutate(sera_count_cart_2d_post = sum(type == "sera"),
         virus_count_cart_2d_post = sum(type == "antigen"))  %>%
  dplyr::select(virus_count_cart_2d_post, sera_count_cart_2d_post) %>%
  dplyr::distinct() %>%
  base::cbind(overall_sample)

overall_sample.2 <- AcmapToDF(map = Racmacs::read.acmap(filename = here::here("data/processed/h1_maps/h1_pre_sd_2d.ace"))) %>%
  dplyr::mutate(sera_count_cart_2d_pre = sum(type == "sera"),
         virus_count_cart_2d_pre = sum(type == "antigen"),
         season = "Overall",
         subtype = "H1N1")  %>%
  dplyr::select(season, subtype, virus_count_cart_2d_pre, sera_count_cart_2d_pre) %>%
  dplyr::distinct()

overall_sample <- AcmapToDF(map = Racmacs::read.acmap(filename = here::here("data/processed/h1_maps/h1_post_sd_2d.ace"))) %>%
  dplyr::mutate(sera_count_cart_2d_post = sum(type == "sera"),
         virus_count_cart_2d_post = sum(type == "antigen"))  %>%
  dplyr::select(virus_count_cart_2d_post, sera_count_cart_2d_post) %>%
  dplyr::distinct() %>%
  base::cbind(overall_sample.2) %>%
  dplyr::full_join(overall_sample,
                   by = c("virus_count_cart_2d_post", "sera_count_cart_2d_post", "season", "subtype", "virus_count_cart_2d_pre",
                          "sera_count_cart_2d_pre")) %>%
  dplyr::select(season, subtype, virus_count_cart_2d_pre, virus_count_cart_2d_post, sera_count_cart_2d_pre, sera_count_cart_2d_post) 

sample_size <- h1_distance %>%
  dplyr::select(-distance) %>%
  dplyr::filter(season != "overall") %>%
  dplyr::group_by(season, method) %>%
  dplyr::mutate(virus_count = dplyr::n_distinct(identifier),
         subtype = "H1N1",
         .keep = "unused") %>%
  dplyr::ungroup() %>%
  dplyr::distinct() %>%
  tidyr::pivot_wider(names_from = method, 
                     values_from = c(sera_count, virus_count)) %>%
  dplyr::select(season, subtype, virus_count_cart_2d_pre, virus_count_cart_2d_post, sera_count_cart_2d_pre, sera_count_cart_2d_post) %>% 
  base::rbind(overall_sample)

h3_distance %>%
  dplyr::select(-distance) %>%
  dplyr::filter(season != "overall") %>%
  dplyr::group_by(season, method) %>%
  dplyr::mutate(virus_count = dplyr::n_distinct(identifier),
         subtype = "H3N2",
         .keep = "unused") %>%
  dplyr::ungroup() %>%
  dplyr::distinct() %>%
  tidyr::pivot_wider(names_from = method, 
                     values_from = c(sera_count, virus_count)) %>%
  dplyr::select(season, subtype, virus_count_cart_2d_pre, virus_count_cart_2d_post, sera_count_cart_2d_pre, sera_count_cart_2d_post)  %>%
  base::rbind(sample_size) %>%
  dplyr::arrange(season, subtype) %>%
  saveRDS(.,
          file = here::here("results/tables/SeasonBasedSampleSize.rds"))

