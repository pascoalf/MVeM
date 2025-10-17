# remove contamination
# function to remove ASVs identified in the control respective to a sample
remove_contamination <- function(data, sample, 
                                 map_sample = sample_control_map_long, 
                                 threshold = NULL, option = NULL, 
                                 output = "standard", ...){
  # parameters warning
  if(!is.null(option)){
    if(!option %in% c("automatic", "threshold")){
      stop("The option argument can be either 'automatic' or 'threshold'")
    }
  }
  if(!output %in% c("standard", "contaminants")){
    stop("The output argument can be either 'standard' or 'contaminants'")
  }
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

  # Safe ASVs -- too abundant in original sample to be removed  
  # Threshold option
  if(is.null(option)){
    data %>%
      filter(Sample == sample) %>% 
      filter(!ASV %in% asvs_in_control)
    asvs_in_control.df <- data.frame(ASV = asvs_in_control)
  } else if(option == "threshold"){
    safe_ASVs <- data %>% 
      filter(Sample == sample) %>%
      filter(Abundance >= threshold) %>% 
      pull(ASV) %>% 
      unique()
    safe_ASVs.df <- data.frame(ASV = safe_ASVs)    
    # Remove safe ASVs from contaminant list
    asvs_in_control.df <- data.frame(ASV = asvs_in_control) %>% 
      anti_join(safe_ASVs.df, by = "ASV")
  } else if(option == "automatic"){
    # unsupervised clustering
    ucluster <- suppressWarnings(define_rb(data))
    safe_ASVs <- ucluster %>% 
      filter(Sample == sample) %>%
      filter(Classification != "Rare") %>% 
      pull(ASV) %>% 
      unique()
    safe_ASVs.df <- data.frame(ASV = safe_ASVs)
    # Remove safe ASVs from contaminant list
    asvs_in_control.df <- data.frame(ASV = asvs_in_control) %>% 
      anti_join(safe_ASVs.df, by = "ASV")    
  }

  # Select ASVs with zero counts on data
  emptyASVs <- data %>% 
    filter(Sample == sample) %>%
    filter(Abundance == 0) %>% 
    pull(ASV)
  #
  if(output == "standard"){
    no_cont_table <- data %>%
      filter(Sample == sample) %>% 
      filter(!ASV %in% asvs_in_control.df$ASV)
    return(no_cont_table)
  } else if(output == "contaminants"){
    asvs_in_control.df <- 
      asvs_in_control.df %>% 
      filter(!ASV %in% emptyASVs)
    names(asvs_in_control.df) <- sample
    return(asvs_in_control.df)
  }
}
