## Functions for 1_cleaning_code.R

## Cleaning the overall data:

cleaning_all <- function(file = "data/raw/data.rds",
                         save_file = TRUE) {
  
  df <- base::readRDS(file = here::here(file)) %>%
    dplyr::mutate(uniq_id = ifelse(uniq_id == "uga2018_id_265 (c022)", "uga2018_id_265", uniq_id)) %>% #This individual is coded differently. Change to match
    dplyr::distinct()
  
  if(save_file == TRUE) {
    base::saveRDS(object = df,
            file = here::here("data/processed/clean_data.rds"))
  }
  
  return(df)
}

virus_clean <- function(file = "data/raw/accession_no.xlsx") {
  df <- readxl::read_excel(here::here(file)) %>%
    dplyr::arrange(factor_order) %>%
    dplyr::select(c(full_strain_name, analysis_name, short_name, uniprot_ac, ha_sequence, gisaid_ac, full_length)) %>%
    dplyr::mutate(uniprot_ac = ifelse(is.na(uniprot_ac), NA, paste0("UniProt: ", uniprot_ac)),
           gisaid_ac = ifelse(is.na(gisaid_ac), NA, paste0("Gisaid: ", gisaid_ac)),
           ha_sequence = ifelse(is.na(ha_sequence), NA, paste0("HA sequence: ", ha_sequence)),
           ha_source = dplyr::coalesce(uniprot_ac, gisaid_ac, ha_sequence)) %>%
    dplyr::select(-c(uniprot_ac, gisaid_ac, ha_sequence)) %>%
    dplyr::mutate(full_strain_name = stringr::str_to_title(full_strain_name),
           full_strain_name = stringr::str_replace_all(string = full_strain_name,
                                                       pattern = "_",
                                                       replacement = "/"),
           full_strain_name = stringr::str_replace(string = full_strain_name,
                                                   pattern = "(?<=[:digit:])n(?=[:digit:])",
                                                   replacement = "N"))
  return(df)
  
}

## Pull the vaccine strains used for Selection: ###########

### Vaccine data list:

pull_vaccine_strains <- function(dataframe = individuals,
                                 season_of_interest = c(2014:2019),
                                 dose_of_interest = c("SD", "HD"),
                                 subtypes_of_interest = c("H1N1", "H3N2", "B-Pre", "B-Yam", "B-Vic")) {
  # Variables for troubleshooting
  ##dataframe <- individuals
  ##season_of_interest <- 2014:2019
  ##dose_of_interest <- c("SD", "HD")
  ##subtypes_of_interest <- c("H1N1", "H3N2", "B-Pre", "B-Yam", "B-Vic")
  
  vaccine_strains <- dataframe %>%
  dplyr::filter(season %in% season_of_interest,
                dose %in% dose_of_interest,
                strain_type %in% subtypes_of_interest) %>%
  dplyr::select(strain_type, dplyr::ends_with("vaccine_fullname")) %>%
  dplyr::distinct() %>%
  tidyr::pivot_longer(cols = dplyr::ends_with("vaccine_fullname"),
                      names_to = "vaccine_subtype",
                      values_to = "vaccine_fullname") %>%
  dplyr::filter(vaccine_fullname != "None") %>% 
  dplyr::select(vaccine_fullname) %>%
  dplyr::distinct() %>%
  dplyr::pull() %>%
  base::droplevels()
  
  return(vaccine_strains)

}


## Subset the data for sensitivity analysis ######################

subset_the_data <- function(dataframe = virus, #Virus dataframe that contains all of the viruses
                            number_of_subsets = 5, # Number of different dataframes
                            number_of_viruses = 10, #Number of different viruses to pull per straintype
                            starting_seed = 123,
                            individuals_df = individuals) #Starting seed to use for randomized sampling
{ 
  #Variables for Troubleshooting: 
  ##dataframe <- virus
  ##number_of_viruses <- 10
  ##number_of_subsets <- 5
  ##starting_seed <- 123
  #individuals_df <- individuals
  
  #Pull the vaccine strains from the individuals dataframe, these need to be included in the analysis:
  vaccine_list <- pull_vaccine_strains(dataframe = individuals_df)
  
  base::saveRDS(object = vaccine_list,
          file = here::here("data/processed/vaccine_viruses.rds"))
  
    
  
  # Create the dataframe that is to be filled
  subset_vector <- base::data.frame(subset = NA, 
                              subtype = NA, 
                              analysis_name = NA, 
                              factor_order = NA) #There are three 
  
  #Fill the dataframe  
  for (i in 1:number_of_subsets) { #Run this for x amount of times
    
    set.seed((starting_seed+(1-i))) # Set the seed for the randomized sampling
    
    subset_vector <- dataframe %>%
      dplyr::filter(!(analysis_name %in% vaccine_list)) %>%
      dplyr::group_by(subtype) %>%
      dplyr::slice_sample(n = number_of_viruses) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(subset = i) %>%
      dplyr::select(subset, subtype, analysis_name, factor_order) %>%
      base::rbind(subset_vector) %>%
      dplyr::filter(!is.na(subset)) %>%
      dplyr::distinct()
  }
  
  
  
  #Save the subsetted data for all of the individuals for future analysis                            
  for (i in 1:number_of_subsets) {
    #Pull the vector to use for filtering
    vector <- subset_vector %>%
      dplyr::filter(subset == i) %>%
      dplyr::pull(analysis_name)
    
    #Filter the individuals
    individuals_df %>%
      dplyr::filter(strains_fullname %in% c(vector, paste(vaccine_list))) %>%
      base::saveRDS(file = here::here(paste0("data/processed/subset_analysis/datasets/subset_", i, ".rds")))
  }
  
  return(subset_vector)
  
}

  
