# remove contamination
# function to remove ASVs identified in the control respective to a sample
remove_contamination <- function(data, sample, map_sample = sample_control_map_long, ...){
  # make helper function to extract specific controls
  extract_controls <- function(x = sample_control_map_long,
                               asvs_original = data, 
                               sample = sample, 
                               type = "Extraction_control", ...){
    # get IDs
    ids <- x %>%
      filter(Sample_name == sample) %>% 
      filter(Control_type == type) %>% 
      pull(Control_ID)
    
    # get ASVs
    control_ASVs <- asvs_original %>% 
      ungroup() %>% 
      filter(Sample %in% ids) %>%
      filter(Abundance > 0) %>% 
      pull(ASV) %>% 
      unique()
    return(control_ASVs)
  }  
  
  # Obtain control sample IDs
  
  extraction_control <- extract_controls(sample = sample, type = "Extraction_control")
  filtration_control <- extract_controls(sample = sample, type = "Filtration_control")
  PCR_control <- extract_controls(sample = sample, type = "PCR_control")
  
  # ASVs in all controls
  asvs_in_control <- c(extraction_control, filtration_control, PCR_control) %>% 
    unique()
  
  #
  data %>%
    filter(Sample == sample) %>% 
    filter(!ASV %in% asvs_in_control)
}