# 12/29/2022

# Purpose of this code:

#The purpose of this code is to compare the antigenic cartographic distances of the @fonville2015, @smith2004 and the current work investigated here.@smith2004 maps were created with naive ferret sera. @fonville2015 maps were created with naive ferret sera and naive human sera. Current work maps were created with preimmune human sera from either before and after vaccination with the Fluzone vaccine.This will also compare the with and without scaling of the maps.

# Load Packages:

library(Racmacs) #Manipulating Acmaps
library(ggrepel) #plotting
library(ggpubr) #Conducting correlations and adding to plots
library(tidyverse) #Data manipulation

#Load functions:

source(here::here("R/functions/antigenic-cartography.R")) #Antigenic cartography functions

fon_fer <- read.acmap(filename = here::here("data/processed/h3_maps/fonville_ferret_ds5.ace"))
fon_hum <- read.acmap(filename = here::here('data/processed/h3_maps/fonville_human_ds4.ace'))
skar <- read.acmap(filename = here::here('data/processed/h3_maps/h3_post_sd_2d.ace'))
skar_pre <-read.acmap(filename = here::here('data/processed/h3_maps/h3_pre_sd_2d.ace'))
smith <- read.acmap(filename = here::here('data/processed/h3_maps/smith_2004.ace'))
smith_ds <- read.acmap(here::here("data/processed/h3_maps/smith_ferret_check.ace"))
skar_orig <- skar


#Confirmation that the antigen/sera names are identical for positioning:
#agNames(skar_pre)
#agNames(skar)
#agNames(fon_fer)
#agNames(fon_hum)
#agNames(smith)

#Change the skar names to match those of smith and fonville
agNames(skar)[c(1:3, 5:10, 12, 15)] <- c("HK/1/68", "PC/1/73", 'TE/1/77', 'SI/2/87', "SD/9/93", "NA/933/95", "SY/5/97", "PM/2007/99", "FU/411/02", "WN/67/05", "VI/361/11")

agNames(skar_orig)[c(1:3, 5:10, 12, 15)] <- c("HK/1/68", "PC/1/73", 'TE/1/77', 'SI/2/87', "SD/9/93", "NA/933/95", "SY/5/97", "PM/2007/99", "FU/411/02", "WN/67/05", "VI/361/11")
agNames(skar_pre)[c(1:3, 5:10, 12, 15)] <- c("HK/1/68", "PC/1/73", 'TE/1/77', 'SI/2/87', "SD/9/93", "NA/933/95", "SY/5/97", "PM/2007/99", "FU/411/02", "WN/67/05", "VI/361/11")

#Align the maps to smith, 2004 map
#Include scaling or the skar maps. 
fon_fer <- Racmacs::realignMap(map = fon_fer,
                      target_map = smith)
fon_hum <- Racmacs::realignMap(map = fon_hum,
                      target_map = smith)
skar <- Racmacs::realignMap(map = skar,
                   target_map = smith)
skar_pre <- Racmacs::realignMap(map= skar_pre,
                       target_map = smith)
smith_ds <- Racmacs::realignMap(map = smith_ds,
                       target_map = smith)

# Scale the human study maps
fon_humS <- Racmacs::realignMap(map = fon_hum,
                       target_map = smith,
                       scaling = TRUE,
                       translation = TRUE)
skarS <- Racmacs::realignMap(map = skar,
                    target_map = smith,
                    scaling = TRUE,
                    translation = TRUE)
skar_preS <- Racmacs::realignMap(map= skar_pre,
                        target_map = smith,
                        scaling = TRUE,
                        translation = TRUE)

rmsd_tab <- matrix(data = NA,
                   nrow = 6,
                   ncol = 6)
colnames(rmsd_tab) <- c("Map 1", "Map 2", "# Ag", "Translation = T", "Scaling = T", "Both = T")
rmsd_tab[,1] <- c(rep("Smith", 3), rep("Fonville Ferret", 2), "Fonville Human")
rmsd_tab[,2] <- c("Fonville Ferret", "Fonville Human", "Skarlupka", "Fonville Human", "Skarlupka", "Skarlupka")

studies <- list(fon_fer, fon_hum, skar)

for (i in 1:length(studies)) {
  rmsd_tab[i,"# Ag"] <- sum(!is.na(Racmacs::procrustesData(map = studies[[i]],
                                                  comparison_map = smith,
                                                  translation = TRUE,
                                                  scaling = FALSE)$ag_dist))
  
  rmsd_tab[i,"Translation = T"] <- round(Racmacs::procrustesData(map = studies[[i]],
                                                        comparison_map = smith,
                                                        translation = TRUE,
                                                        scaling = FALSE)$ag_rmsd,
                                         digits = 3)
  
  rmsd_tab[i,"Scaling = T"] <- round(Racmacs::procrustesData(map = studies[[i]],
                                                    comparison_map = smith,
                                                    translation = FALSE,
                                                    scaling = TRUE)$ag_rmsd,
                                     digits = 3)
  
  rmsd_tab[i,"Both = T"] <- round(Racmacs::procrustesData(map = studies[[i]],
                                                 comparison_map = smith,
                                                 translation = TRUE,
                                                 scaling = TRUE)$ag_rmsd,
                                  digits = 3)
}

for (i in 1:2) {
  rmsd_tab[i+3,"# Ag"] <- sum(!is.na(Racmacs::procrustesData(map = studies[[i+1]],
                                                    comparison_map = studies[[1]],
                                                    translation = TRUE,
                                                    scaling = FALSE)$ag_dist))
  
  rmsd_tab[i+3,"Translation = T"] <- round(Racmacs::procrustesData(map = studies[[i+1]],
                                                          comparison_map = studies[[1]],
                                                          translation = TRUE,
                                                          scaling = FALSE)$ag_rmsd,
                                           digits = 3)
  
  rmsd_tab[i+3,"Scaling = T"] <- round(Racmacs::procrustesData(map = studies[[i+1]],
                                                      comparison_map = studies[[1]],
                                                      translation = FALSE,
                                                      scaling = TRUE)$ag_rmsd,
                                       digits = 3)
  
  rmsd_tab[i+3,"Both = T"] <- round(Racmacs::procrustesData(map = studies[[i+1]],
                                                   comparison_map = studies[[1]],
                                                   translation = TRUE,
                                                   scaling = TRUE)$ag_rmsd,
                                    digits = 3)
}


rmsd_tab[6,"# Ag"] <- sum(!is.na(Racmacs::procrustesData(map = studies[[3]],
                                                comparison_map = studies[[2]],
                                                translation = TRUE,
                                                scaling = FALSE)$ag_dist))

rmsd_tab[6,"Translation = T"] <- round(Racmacs::procrustesData(map = studies[[3]],
                                                      comparison_map = studies[[2]],
                                                      translation = TRUE,
                                                      scaling = FALSE)$ag_rmsd,
                                       digits = 3)

rmsd_tab[6,"Scaling = T"] <- round(Racmacs::procrustesData(map = studies[[3]],
                                                  comparison_map = studies[[2]],
                                                  translation = FALSE,
                                                  scaling = TRUE)$ag_rmsd,
                                   digits = 3)

rmsd_tab[6,"Both = T"] <- round(Racmacs::procrustesData(map = studies[[3]],
                                               comparison_map = studies[[2]],
                                               translation = TRUE,
                                               scaling = TRUE)$ag_rmsd,
                                digits = 3)


  saveRDS(rmsd_tab,
          file = here::here("results/tables/rmsd_procrustes.rds"))



fon_fer <- AcmapToDF(fon_fer) %>%
  dplyr::filter(type == "antigen") %>%
  dplyr::mutate(study = "fonville_ferret")
fon_hum <- AcmapToDF(fon_hum) %>%
  dplyr::filter(type == "antigen") %>%
  dplyr::mutate(study = "fonville_human")
smith <- AcmapToDF(smith) %>%
  dplyr::filter(type == "antigen") %>%
  dplyr::mutate(study = "smith")
skar <- AcmapToDF(skar) %>%
  dplyr::filter(type == "antigen") %>%
  dplyr::mutate(study = "current")
fon_humS <- AcmapToDF(fon_humS) %>%
  dplyr::filter(type == "antigen") %>%
  dplyr::mutate(study = "fonville_human_scaling")
skarS <- AcmapToDF(skarS) %>%
  dplyr::filter(type == "antigen") %>%
  dplyr::mutate(study = "current_scaling")
skar_preS <- AcmapToDF(skar_preS) %>%
  dplyr::filter(type == "antigen") %>%
  dplyr::mutate(study = "pre_scaling")
skar_pre <- AcmapToDF(skar_pre) %>%
  dplyr::filter(type == "antigen") %>%
  dplyr::mutate(study = "pre")
smith_ds <- AcmapToDF(smith_ds) %>%
  dplyr::filter(type == "antigen") %>%
  dplyr::mutate(study = "smith_ds")


tab <- base::rbind(fon_fer, fon_hum, smith, skar, fon_humS, skarS, skar_preS, skar_pre, smith_ds) %>%
  dplyr::select(-type)

table <-  tab %>%
  dplyr::filter(identifier %in% unique(smith$identifier)) %>%
  dplyr::filter(identifier %in% unique(base::subset(tab, study != "smith")$identifier))

relative <- unique(tab$identifier)
index <- c(5, 7, 8, 10, 12, 11, 13, 15, 14, 1, 2, 3, 4, 6, 9)
relative <- relative[order(index)]

procrus <- smith %>%
  dplyr::select(-type) %>%
  dplyr::rename(V1_smith = V1,
         V2_smith = V2,
         smith = study) %>%
  dplyr::full_join(base::subset(tab, study != "smith"),
            by = "identifier")

plt <- procrus %>%
  dplyr::filter(study %in% c("current", "fonville_ferret", "fonville_human")) %>%
  dplyr::mutate(study = factor(study, 
                               levels = c("fonville_ferret", "fonville_human", "current"), 
                               labels = c("Fonville Ferret", "Fonville Human", "Skarlupka"))) %>%
  ggplot2::ggplot() +
  ggplot2::geom_point(data = smith[,1:2],
             aes(x = V1,
                 y = V2),
             color = "red",
             pch = 1,
             alpha = 0.5) +
  ggplot2::geom_point(aes(x = V1_smith,
                 y = V2_smith))+
  ggplot2::geom_segment(aes(x = V1_smith,
                   y = V2_smith,
                   xend = V1,
                   yend = V2),
               arrow = arrow(length = unit(0.25, "cm"))) +
  ggplot2::facet_wrap(vars(study)) +
  ggrepel::geom_text_repel(aes(label = identifier,
                      x = V1_smith,
                      y = V2_smith),
                  size = 3) +
  ggplot2::ylab("Antigenic Unit") +
  ggplot2::xlab("Antigenic Unit") +
  ggplot2::theme_bw() 

suppressWarnings(cowplot::ggsave2(plt,
        filename = here::here("results/figures/smith_procrustes.png"),
        width = 9,
        height = 4,
        scale = 0.75))
  
  tab <- table %>%
    dplyr::group_by(study) %>%
    dplyr::mutate(index = sum(identifier == "PM/2007/99")) %>%
    dplyr::filter(index > 0) %>%
    dplyr::select(-index) %>%
    dplyr::group_by(study) %>%
    dplyr::mutate(distance = sqrt((V1-V1[which(identifier == "PM/2007/99")])^2 + (V2-V2[which(identifier == "PM/2007/99")])^2)) %>%
    dplyr::select(-c(V1, V2)) %>%
    tidyr::pivot_wider(names_from = study,
                values_from = distance)  %>%
    tidyr::pivot_longer(cols = -c(identifier, smith),
                 names_to = "study",
                 values_to = "distance") %>%
    dplyr::filter(study %in% c("fonville_ferret", "fonville_human", "current")) %>%
    dplyr::mutate(study = factor(study, 
                                 levels = c("fonville_ferret", "fonville_human", "current"), 
                                 labels = c("Fonville Ferret", "Fonville Human", "Skarlupka")))
  
  plt_list <- vector(mode = "list")
  for (i in c("Fonville Ferret", "Fonville Human", "Skarlupka")) {
    x <- ifelse(i == "Skarlupka", 25, 11)
    plt_list[[i]] <- tab %>%
      dplyr::filter(study == i) %>%
      ggplot2::ggplot(aes(x = smith, 
                 y = distance, 
                 label = identifier)) +
      ggplot2::geom_point() +
      ggplot2::geom_abline(intercept = 0, slope = 1, color = "black") +
      ggplot2::geom_smooth(method = "lm") +
      ggrepel::geom_label_repel(alpha = 0.5) +
      ggplot2::theme_bw() +
      ggplot2::xlab("Distance (Smith)") +
      ggplot2::ylab(paste0("Distance (",
                  i,
                  ")")) +
      ggplot2::ylim(c(0, x)) +
      ggplot2::xlim(c(0, x)) +
      ggpubr::stat_cor()
  }
  
  cowplot::plot_grid(plt_list[[1]], plt_list[[2]], plt_list[[3]],
                     nrow = 1,
                     labels = "AUTO") %>%
    cowplot::ggsave2(filename = here::here("results/figures/distance_corr.png"),
            width = 6,
            height = 2,
            scale = 2)



rsqr_tab <- matrix(nrow = length(relative),
                   ncol = 7)
rsqr_tab[,1] <- relative


split_tab <- table %>%
  dplyr::filter(study %in% c("smith", "fonville_ferret", "fonville_human", "current")) %>%
  dplyr::group_split(study)

distance_list <- vector(mode = "list")
for (i in 1:length(split_tab)) {
  a <- split_tab[[i]]
  b <- data.frame(relative_virus = rep(a$identifier, each = length(a$identifier)),
                  V1_rel = rep(a$V1, each = length(a$identifier)),
                  V2_rel = rep(a$V2, each = length(a$identifier)))
  
  distance_list[[i]] <- a %>%
    base::cbind(b) %>%
    dplyr::group_by(relative_virus) %>%
    dplyr::mutate(distance = sqrt((V1-V1_rel)^2 + (V2-V2_rel)^2)) %>%
    dplyr::select(study, identifier, relative_virus, distance)
}

for (i in 1:3) {
  distance_list[[i]] <- dplyr::left_join(distance_list[[i]],
                                  distance_list[[4]],
                                  by = c("relative_virus", "identifier")) %>%
    dplyr::group_by(study.x, relative_virus) %>%
    dplyr::summarize(r_squared = stats::cor(distance.x, 
                              distance.y, 
                              method = "pearson"),
              p_value = stats::cor.test(distance.x,
                                 distance.y,
                                 alternative = "two.sided",
                                 method = "pearson")$p.value) %>%
    dplyr::ungroup()
}

distance_list <- head(distance_list, -1)
c <- data.table::rbindlist(distance_list) %>%
  tidyr::pivot_wider(names_from = study.x,
              values_from = c(r_squared, p_value)) %>%
  dplyr::select(relative_virus, r_squared_fonville_ferret, p_value_fonville_ferret, r_squared_fonville_human, p_value_fonville_human, r_squared_current, p_value_current)


  c[match(relative, c$relative_virus),] %>%
    saveRDS(file = here::here("results/tables/rsqr.rds"))





