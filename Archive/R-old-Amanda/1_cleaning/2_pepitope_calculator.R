# Sequence-Based Distance Calculations

## PURPOSE: This file takes the Uniprot Accession IDs and pulls the full HA
## sequence. Then epitope residues for each antigenic site are extracted from
## the HA amino acid FASTA. RDS files are generated for the separate epitopes
## and then one final key that contains all of the variables. (xxx_key.rds).

## Subtypes: H1N1 and H3N2

# These parameters are to help with quickly rerunning code:

# save_plots turns on the saving of figures and tables. TRUE will save the
# images, FALSE will create the plots within the file still but will not save
# the final image. The images are viewable when knit under both settings.

save_plots <- TRUE

save_files <- TRUE

#The difference matrices of the different epitopes
save_extra <- TRUE


# Packages:
library(seqinr) #Fasta manipulation
library(Biostrings) #Fasta manipulation
library(stringr) #String manipulation
library(readxl) #Working with excel files
library(here) #Referencing files
library(tidyr) #Data manipulation
library(Rcpi) #pulling sequences from Uniprot
library(dplyr) #Data manipulation

# Global Environment Functions:
source(file = here::here("R/functions/2_pepitope.R"))

## Obtaining full HA sequences: ####

# This site:
# https://www.antigenic-cartography.org/surveillance/evergreen/HAnumbering/ was
# used to determine the cut offs of the HA1 and the leader sequence. UniProt
# used used for the HA sequences.

# Data
accession_no <- readxl::read_excel(here::here("data/raw/accession_no.xlsx"))

# Retrieve sequences from UniProt. Use accession numbers to pull amino acid
# sequence and the formal uniprot accession number used.

full <- accession_no %>%
  dplyr::filter(!is.na(uniprot_ac)) %>%
  dplyr::mutate(
    full_uniprot_ac = names(unlist(Rcpi::getSeqFromUniProt(uniprot_ac))),
    ha_sequence = unlist(Rcpi::getSeqFromUniProt(uniprot_ac))) %>%
  dplyr::right_join(accession_no, by = c("subtype", "full_strain_name", "analysis_name", "uniprot_ac")) %>%
  dplyr::mutate(ha_sequence = dplyr::coalesce(ha_sequence.x,
                                              ha_sequence.y),
                gisaid_ac = dplyr::coalesce(gisaid_ac.x,
                                            gisaid_ac.y),
                full_length = dplyr::coalesce(full_length.x,
                                              full_length.y),
                short_name = dplyr::coalesce(short_name.x,
                                             short_name.y),
                factor_order = dplyr::coalesce(factor_order.x,
                                               factor_order.y)) %>%
  dplyr::select(!ends_with(c(".x", ".y")))


sub <- full %>%
  dplyr::group_by(subtype) %>%
  dplyr::group_split()


#H1N1 Influenza
h1 <- base::list2DF(sub[[1]])
#H3N2 Influenza
h3 <- base::list2DF(sub[[2]])

#The full accession numbers and amino acid sequences are saved here to be recalled later, in case there are issues replicating the download process from Uniprot
base::saveRDS(object = full,
              file = here::here("data/processed/virus_protein_information.rds"))

## H1 Influenza ####
### Defining the epitopes ####
#H1 HA epitope residues with A/California/04/2009 numbering.

subtype <- "H1N1"

site_ha1 <- c(1:326)

site_a <- c(118, 120:122, 126:129, 132:135, 137, 139:143, 146, 147, 149, 165, 252, 253)

site_b <- c(124, 125, 152:157, 160, 162, 183:187, 189:191, 193:196)

site_c <- c(34:38, 40, 41, 43:45, 269:274, 276:278, 283, 288, 292, 295, 297, 298, 302, 303, 305:310)

site_d <- c(89, 94:96, 113, 117, 163, 164, 166:174, 176, 179, 198, 200, 202, 204:216, 222:227, 235, 237, 241, 243:245)

site_e <- c(47, 48, 50, 51, 53, 54, 56:58, 66, 68:75, 78:80, 82:86, 102, 257:261, 263, 267)

site_all_epitopes <- c(site_a, site_b, site_c, site_d, site_e)


epitope_tbl <- tidyr::tibble(`Subtype` = subtype,
                             `Antigenic Site` = c("A", 'B', 'C', 'D', 'E'),
                             `Epitope Residues` = c(paste(site_a, collapse = ", "),
                                                    paste(site_b, collapse = ", "),
                                                    paste(site_c, collapse = ", "),
                                                    paste(site_d, collapse = ", "),
                                                    paste(site_e, collapse = ", ")))

# Take the string, expand it into individual lists of characters then select the epitope sites then recollapse.
### Alignment ####
#The raw_ha_0 does not account for gaps that need to be added to align the sequences.
seqs_df <- h1 %>%
  dplyr::mutate(raw_ha_0 = stringr::str_sub(ha_sequence, 18, -1),
                length = stringr::str_length(raw_ha_0),
                ha_0 = ifelse(length == 548,
                              gsub('^(.{129})(.*)$', '\\1-\\2', raw_ha_0),
                              paste(raw_ha_0)),
                validate_length = stringr::str_length(ha_0))


epitope_list <- list(site_ha1, site_all_epitopes, site_a, site_b, site_c, site_d, site_e)
names(epitope_list) <- c("site_ha1", "site_all_epitopes", "site_a", "site_b", "site_c", "site_d", "site_e")

starting <- ncol(seqs_df)

for (i in 1:length(epitope_list)) {
  seqs_df <- combin(epitope_site = epitope_list[[i]],
                    dataframe = seqs_df,
                    reference_column = "ha_0")
  names(seqs_df)[starting+i] <- names(epitope_list)[i]
}

### Difference and P-value matrices ####

#### Entire HA sequence: HA0 ####
ha_0_difference_matrix <- difference_matrix(data = subset(seqs_df,
                                                          select = c(analysis_name, ha_0)))
pha0 <- cbind(ha_0_difference_matrix[,1:2],
              ha_0_difference_matrix[,-c(1:2)]/549)

#### HA head sequence: HA1 ####
ha_1 <- difference_matrix(data = subset(seqs_df,
                                        select = c(analysis_name, site_ha1)))
pha_1 <- cbind(ha_1[,1:2],
               ha_1[,-c(1:2)]/length(site_ha1))

#### All epitopes ####
ha_all_epitopes <- difference_matrix(data = subset(seqs_df,
                                                   select = c(analysis_name, site_all_epitopes)))
pall <- cbind(ha_all_epitopes[,1:2],
              ha_all_epitopes[,-c(1:2)]/length(site_all_epitopes))

#### Site A ####
ha_site_a <- difference_matrix(data = subset(seqs_df,
                                             select = c(analysis_name, site_a)))
pa <- cbind(ha_site_a[,1:2],
            ha_site_a[,-c(1:2)]/length(site_a))

#### Site B ####
ha_site_b <- difference_matrix(data = subset(seqs_df,
                                             select = c(analysis_name, site_b)))
pb <- cbind(ha_site_b[,1:2],
            ha_site_b[,-c(1:2)]/length(site_b))

#### Site C ####
ha_site_c <- difference_matrix(data = subset(seqs_df,
                                             select = c(analysis_name, site_c)))
pc <- cbind(ha_site_c[,1:2],
            ha_site_c[,-c(1:2)]/length(site_c))

#### Site D ####
ha_site_d <- difference_matrix(data = subset(seqs_df,
                                             select = c(analysis_name, site_d)))
pd <- cbind(ha_site_d[,1:2],
            ha_site_d[,-c(1:2)]/length(site_d))

#### Site E ####
ha_site_e <- difference_matrix(data = subset(seqs_df,
                                             select = c(analysis_name, site_e)))
pe <- cbind(ha_site_e[,1:2],
            ha_site_e[,-c(1:2)]/length(site_e))

#### P-epitope ####

sites <- list(pha0, pha_1, pall, pa, pb, pc, pd, pe)

names(sites) <- c("pha0", "pha_1", "pall", "pa", "pb", "pc", "pd", "pe")

sites <- lapply(sites,`[`,-2)

for (i in 1:length(sites)) {
  sites[[i]] <- tidyr::pivot_longer(sites[[i]],
                                    cols = -1,
                                    names_to = "strain_2",
                                    values_to = names(sites[i]))
}

p_epitope <- Reduce(
  function(x, y, ...) merge(x, y, all = TRUE, ...),
  sites
)

p_epitope$p_epi_site <- colnames(subset(p_epitope, select = pa:pe))[apply(subset(p_epitope, select = pa:pe),1,which.max)]

p_epitope$p_epi <- apply(subset(p_epitope, select = pa:pe), 1, max)

col_order <- c("analysis_name", "strain_2", "p_epi_site", "p_epi", "pha0", "pha_1", "pall", "pa", "pb", "pc", "pd", "pe")

p_epitope <- p_epitope[, col_order]

#### Anderson, 2018 Distance ####

p_anderson <- p_epitope %>%
  dplyr::select(c(analysis_name, strain_2, pa:pe)) %>%
  dplyr::mutate_if(is.numeric,
                   function(x) {x*20}) %>%
  dplyr::rowwise() %>%
  dplyr::mutate(anderson = sum(c(pa, pb, pc, pd, pe))/5)

working <- p_epitope %>%
  dplyr::select(-p_epi_site) %>%
  dplyr::left_join(p_anderson[,c(1:2,8)],
                   by = c("analysis_name", "strain_2")) %>%
  tidyr::pivot_longer(cols = p_epi:anderson,
                      names_to = "method",
                      values_to = "distance") %>%
  tidyr::pivot_wider(names_from = strain_2,
                     values_from = distance)

### Save Files for H1 ####

if (save_files == TRUE) {

  saveRDS(working, here::here("data/processed/h1_seq_distance/h1_key.rds"))

  if(save_extra == TRUE) {

    saveRDS(ha_0_difference_matrix,
            file = here("data/processed/h1_seq_distance/h1_ha0_difference.rds"))
    saveRDS(pha0,
            file = here("data/processed/h1_seq_distance/h1_ha0_pvalue.rds"))
    saveRDS(ha_1,
            file = here("data/processed/h1_seq_distance/h1_ha1_difference.rds"))
    saveRDS(pha_1,
            file = here("data/processed/h1_seq_distance/h1_ha1_pvalue.rds"))
    saveRDS(ha_all_epitopes,
            file = here("data/processed/h1_seq_distance/h1_all_epitopes_difference.rds"))
    saveRDS(pall,
            file = here("data/processed/h1_seq_distance/h1_all_epitopes_pvalue.rds"))
    saveRDS(ha_site_a,
            file = here("data/processed/h1_seq_distance/h1_siteA_difference.rds"))
    saveRDS(pa,
            file = here("data/processed/h1_seq_distance/h1_siteA_pvalue.rds"))
    saveRDS(ha_site_b,
            file = here("data/processed/h1_seq_distance/h1_siteB_difference.rds"))
    saveRDS(pb,
            file = here("data/processed/h1_seq_distance/h1_siteB_pvalue.rds"))
    saveRDS(ha_site_c,
            file = here("data/processed/h1_seq_distance/h1_siteC_difference.rds"))
    saveRDS(pc,
            file = here("data/processed/h1_seq_distance/h1_siteC_pvalue.rds"))
    saveRDS(ha_site_d,
            file = here("data/processed/h1_seq_distance/h1_siteD_difference.rds"))
    saveRDS(pd,
            file = here("data/processed/h1_seq_distance/h1_siteD_pvalue.rds"))
    saveRDS(ha_site_e,
            file = here("data/processed/h1_seq_distance/h1_siteE_difference.rds"))
    saveRDS(pe,
            file = here("data/processed/h1_seq_distance/h1_siteE_pvalue.rds"))
    saveRDS(p_epitope,
            file = here("data/processed/h1_seq_distance/h1_pepitope_pvalue.rds"))
    saveRDS(p_anderson,
            file = here("data/processed/h1_seq_distance/h1_anderson_distance.rds"))

  }

}

## H3 Influenza ####

### Defining the Epitopes
# H3 HA epitope residues from Deem p-epitope calculator

site_ha1 <- c(1:328)

site_a <- c(122, 124, 126, 130:133, 135, 137:138, 140, 142:146, 150, 152, 168)

site_b <- c(128, 129, 155:160, 163, 165, 186:190, 192:194, 196:198)

site_c <- c(44:48, 50, 51, 53, 54, 273, 275:276, 278:280, 294, 297, 299, 300, 304:305, 307:312)

site_d <- c(96, 102, 103, 117, 121, 167, 170:177, 179, 182, 201, 203, 207:209, 212:219, 226:230, 238, 240, 242, 244, 246:248)

site_e <- c(57, 59, 62, 63, 67, 75, 78, 80:83, 86:88, 91, 92, 94, 109, 260:262, 265)

site_all_epitopes <- c(site_a, site_b, site_c, site_d, site_e)

### Collapse & Alignment ####

epitope_tbl <- tidyr::tibble(`Subtype` = "H3N2",
                             `Antigenic Site` = c("A", 'B', 'C', 'D', 'E'),
                             `Epitope Residues` = c(paste(site_a, collapse = ", "),
                                                    paste(site_b, collapse = ", "),
                                                    paste(site_c, collapse = ", "),
                                                    paste(site_d, collapse = ", "),
                                                    paste(site_e, collapse = ", "))) %>%
  rbind(epitope_tbl) %>%
  dplyr::arrange(Subtype, `Antigenic Site`)

if (save_files == TRUE) {
  base::saveRDS(object = epitope_tbl,
                file = here::here("Results", "Tables", "antigenic_sites.rds"))

}


# Alignment
#Select only full length HA sequences:
seqs_df <- h3 %>%
  dplyr::mutate(ha_0 = stringr::str_sub(ha_sequence, 17, -1),
                length = stringr::str_length(ha_0))


epitope_list <- list(site_ha1, site_all_epitopes, site_a, site_b, site_c, site_d, site_e)
names(epitope_list) <- c("site_ha1", "site_all_epitopes", "site_a", "site_b", "site_c", "site_d", "site_e")
starting <- ncol(seqs_df)

for (i in 1:length(epitope_list)) {
  seqs_df <- combin(epitope_site = epitope_list[[i]],
                    dataframe = seqs_df,
                    reference_column = "ha_0")
  names(seqs_df)[starting+i] <- names(epitope_list)[i]
}

seqs_df <- seqs_df %>%
  select(subtype, analysis_name, length, ha_0, site_ha1:site_e)


### Difference and P-value matrices ####

#### Entire HA sequence: HA0 ####
ha_0_difference_matrix <- difference_matrix(data = subset(seqs_df,
                                                          length == 550,
                                                          select = c(analysis_name, ha_0)))
pha0 <- cbind(ha_0_difference_matrix[,1:2],
              ha_0_difference_matrix[,-c(1:2)]/549)

####HA Head: HA1 ####
ha_1 <- difference_matrix(data = subset(seqs_df,
                                        select = c(analysis_name, site_ha1)))
pha_1 <- base::cbind(ha_1[,1:2],
                     ha_1[,-c(1:2)]/length(site_ha1))

#### All Sites ####
ha_all_epitopes <- difference_matrix(data = subset(seqs_df,
                                                   select = c(analysis_name, site_all_epitopes)))
pall <- base::cbind(ha_all_epitopes[,1:2],
                    ha_all_epitopes[,-c(1:2)]/length(site_all_epitopes))

#### Site A ####
ha_site_a <- difference_matrix(data = subset(seqs_df,
                                             select = c(analysis_name, site_a)))
pa <- base::cbind(ha_site_a[,1:2],
                  ha_site_a[,-c(1:2)]/length(site_a))

#### Site B ####
ha_site_b <- difference_matrix(data = subset(seqs_df,
                                             select = c(analysis_name, site_b)))
pb <- base::cbind(ha_site_b[,1:2],
                  ha_site_b[,-c(1:2)]/length(site_b))

#### Site C ####
ha_site_c <- difference_matrix(data = subset(seqs_df,
                                             select = c(analysis_name, site_c)))
pc <- base::cbind(ha_site_c[,1:2],
                  ha_site_c[,-c(1:2)]/length(site_c))

#### Site D ####
ha_site_d <- difference_matrix(data = subset(seqs_df,
                                             select = c(analysis_name, site_d)))
pd <- base::cbind(ha_site_d[,1:2],
                  ha_site_d[,-c(1:2)]/length(site_d))

#### Site E ####
ha_site_e <- difference_matrix(data = subset(seqs_df,
                                             select = c(analysis_name, site_e)))
pe <- base::cbind(ha_site_e[,1:2],
                  ha_site_e[,-c(1:2)]/length(site_e))

#### P-epitope ####
sites <- list(pha0, pha_1, pall, pa, pb, pc, pd, pe)
names(sites) <- c("pha0", "pha_1", "pall", "pa", "pb", "pc", "pd", "pe")
sites <- base::lapply(sites,`[`,-2)

for (i in 1:length(sites)) {
  sites[[i]] <- tidyr::pivot_longer(sites[[i]],
                                    cols = -1,
                                    names_to = "strain_2",
                                    values_to = names(sites[i]))
}

p_epitope <- Reduce(
  function(x, y, ...) merge(x, y, all = TRUE, ...),
  sites
)

p_epitope$p_epi_site <- colnames(subset(p_epitope, select = pa:pe))[apply(subset(p_epitope, select = pa:pe),1,which.max)]
p_epitope$p_epi <- apply(subset(p_epitope, select = pa:pe), 1, max)

col_order <- c("analysis_name", "strain_2", "p_epi_site", "p_epi", "pha0", "pha_1", "pall", "pa", "pb", "pc", "pd", "pe")
p_epitope <- p_epitope[, col_order]

#### Anderson, 2018 Distance ####

p_anderson <- p_epitope %>%
  dplyr::select(c(analysis_name, strain_2, pa:pe)) %>%
  dplyr::mutate_if(is.numeric,
                   function(x) {x*20}) %>%
  dplyr::rowwise() %>%
  dplyr::mutate(anderson = sum(c(pa, pb, pc, pd, pe))/5)

working <- p_epitope %>%
  dplyr::select(-p_epi_site) %>%
  dplyr::left_join(p_anderson[,c(1:2,8)],
                   by = c("analysis_name", "strain_2")) %>%
  tidyr::pivot_longer(cols = p_epi:anderson,
                      names_to = "method",
                      values_to = "distance") %>%
  tidyr::pivot_wider(names_from = strain_2,
                     values_from = distance)

### Save Files for H3 ####

if (save_files == TRUE) {
  saveRDS(working, here::here("data/processed/h3_seq_distance/h3_key.rds"))

  if (save_extra == TRUE) {

    saveRDS(ha_0_difference_matrix, file = here("data/processed/h3_seq_distance/h3_ha0_difference.rds"))

    saveRDS(pha0, file = here("data/processed/h3_seq_distance/h3_ha0_pvalue.rds"))

    saveRDS(ha_1, file = here("data/processed/h3_seq_distance/h3_ha1_difference.rds"))

    saveRDS(pha_1, file = here("data/processed/h3_seq_distance/h3_ha1_pvalue.rds"))

    saveRDS(ha_all_epitopes, file = here("data/processed/h3_seq_distance/h3_all_epitopes_difference.rds"))

    saveRDS(pall, file = here("data/processed/h3_seq_distance/h3_all_epitopes_pvalue.rds"))

    saveRDS(ha_site_a, file = here("data/processed/h3_seq_distance/h3_siteA_difference.rds"))

    saveRDS(pa, file = here("data/processed/h3_seq_distance/h3_siteA_pvalue.rds"))

    saveRDS(ha_site_b, file = here("data/processed/h3_seq_distance/h3_siteB_difference.rds"))

    saveRDS(pb, file = here("data/processed/h3_seq_distance/h3_siteB_pvalue.rds"))

    saveRDS(ha_site_c, file = here("data/processed/h3_seq_distance/h3_siteC_difference.rds"))

    saveRDS(pc, file = here("data/processed/h3_seq_distance/h3_siteC_pvalue.rds"))

    saveRDS(ha_site_d, file = here("data/processed/h3_seq_distance/h3_siteD_difference.rds"))

    saveRDS(pd, file = here("data/processed/h3_seq_distance/h3_siteD_pvalue.rds"))

    saveRDS(ha_site_e, file = here("data/processed/h3_seq_distance/h3_siteE_difference.rds"))

    saveRDS(pe, file = here("data/processed/h3_seq_distance/h3_siteE_pvalue.rds"))

    saveRDS(p_epitope,
            file = here("data/processed/h3_seq_distance/h3_pepitope_pvalue.rds"))

    saveRDS(p_anderson,
            file = here("data/processed/h3_seq_distance/h3_anderson_distance.rds"))
  }

}

