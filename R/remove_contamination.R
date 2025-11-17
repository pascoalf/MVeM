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
  if(!output %in% c("standard", "contaminants", "saved")){
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
  
  # Obtain control sample ASVss
  extraction_control <- extract_controls(sample = sample, type = "Extraction_control")
  filtration_control <- extract_controls(sample = sample, type = "Filtration_control")
  PCR_control <- extract_controls(sample = sample, type = "PCR_control")
  
  # ASVs in all controls
  asvs_in_control <- c(extraction_control, filtration_control, PCR_control) %>% 
    unique()

  # Safe ASVs -- too abundant in original sample to be removed  
  # Threshold option
  if(is.null(option)){
    asvs_in_control.df <- data.frame(ASV = asvs_in_control)
  } else if(option == "threshold"){
    safe_ASVs <- data %>% 
      filter(Sample == sample) %>%
      mutate(relativeAbundance = Abundance*100/sum(Abundance)) %>% 
      filter(relativeAbundance >= threshold) %>% 
      pull(ASV) %>% 
      unique()
    safe_ASVs.df <- data.frame(ASV = safe_ASVs)    
    # Remove safe ASVs from contaminant list
    asvs_in_control.df <- data.frame(ASV = asvs_in_control) %>% 
      anti_join(safe_ASVs.df, by = "ASV")
  } else if(option == "automatic"){
    # Select ASVs to remove based on ulrb clusters
    # by temporary combination of env and control sample
    remove_from_control <- function(type = "Extraction_control", ...){
      # Get sample and respective control
      control_sample <- map_sample %>%
        filter(Control_type == type,
               Sample_name == sample) %>% 
        select(Sample_name, Control_ID) %>% 
        distinct()

      # Combine sample and control temporarily
      temp1 <- data %>% ungroup() %>% 
        filter(Sample %in% control_sample[1,]) %>% 
        mutate(Source = ifelse(Sample %in% control_sample[1,2], "Control", "Sample")) %>% 
        mutate(Sample = "temp1") %>% 
        filter(Abundance > 0)

      # Cluster using ulrb
        if(dim(temp1)[1] == 3){
              temp1_ucluster <- suppressWarnings(
                define_rb(
                  temp1, 
                  classification_vector = c("Rare", "Abundant"),
                  simplified = TRUE))
        } else if(dim(temp1)[1] < 3){
          # Manual input
          temp1_ucluster <- temp1 %>% 
            mutate(Classification = "Undetermined")
          } else {
            temp1_ucluster <- suppressWarnings(
              define_rb(temp1, 
                        simplified = TRUE))
            }

      # Reformat after clustering
      temp1_ucluster <- temp1_ucluster %>% 
        select(Source, ASV, Classification) %>% 
        distinct()
      ## Extract all possibilities
      # rare in control
      c_asv_rare <- temp1_ucluster %>% filter(Source == "Control", Classification == "Rare") %>% pull(ASV) %>% unique()
      # undetermined in control
      c_asv_und <- temp1_ucluster %>% filter(Source == "Control", Classification == "Undetermined") %>% pull(ASV) %>% unique()
      # abundant in control
      c_asv_abu <- temp1_ucluster %>% filter(Source == "Control", Classification == "Abundant") %>% pull(ASV) %>% unique()
      # rare in env
      e_asv_rare <- temp1_ucluster %>% filter(Source == "Sample", Classification == "Rare") %>% pull(ASV) %>% unique()
      # undetermined in control
      e_asv_und <- temp1_ucluster %>% filter(Source == "Sample", Classification == "Undetermined") %>% pull(ASV) %>% unique()
      # abundant in control
      e_asv_abu <- temp1_ucluster %>% filter(Source == "Sample", Classification == "Abundant") %>% pull(ASV) %>% unique()
      
      ## Decide to remove or save each ASV
      asvs_output <- temp1_ucluster %>% 
        mutate(output = case_when(ASV %in% c_asv_rare & ASV %in% e_asv_rare ~ "Remove",
                                  ASV %in% c_asv_und & ASV %in% e_asv_und ~ "Remove",
                                  ASV %in% c_asv_und ~ "Remove",
                                  ASV %in% c_asv_abu ~ "Remove",
                                  ASV %in% c_asv_rare & ASV %in% e_asv_und ~ "Save",
                                  ASV %in% c_asv_rare & ASV %in% e_asv_abu ~ "Save",
                                  TRUE ~ NA))
      
      # output is a data frame with ASVs to save
      remove_asvs <- asvs_output %>% 
        filter(ASV %in% extract_controls(sample = sample, type = type)) %>% 
        filter(output == "Remove")
      # alternative output
      save_asvs <- asvs_output %>% 
        filter(ASV %in% extract_controls(sample = sample, type = type)) %>% 
        filter(output != "Remove")
      # in case there is nothing
      if(dim(remove_asvs)[1] == 0){
        remove_asvs <- data.frame(Source = NA, ASV = NA, Classification = NA, output = NA)
        save_asvs <- data.frame(Source = NA, ASV = NA, Classification = NA, output = NA)
      }
      if(output == "saved"){
        return(save_asvs)
      } else {
        return(remove_asvs)        
      }
    }
    # for each control type
    asvs_in_control.df <- map(.x = c("Extraction_control", "Filtration_control","PCR_control"),
        .f = ~remove_from_control(type = .x)) %>% 
      bind_rows() %>% ungroup() %>% 
      select(ASV) %>% 
      distinct() %>% 
      filter(!is.na(ASV))
  }
  
  # Select ASVs with zero counts in data
  emptyASVs <- data %>% 
    filter(Sample == sample) %>%
    filter(Abundance == 0) %>% 
    pull(ASV)
  #
  if(output == "standard"){
    no_cont_table <- data %>%
      filter(Sample == sample) %>% 
      filter(!ASV %in% asvs_in_control.df$ASV) %>% 
      filter(Abundance > 0)
    return(no_cont_table)
  } else if(output == "contaminants"){
    asvs_in_control.df <- 
      asvs_in_control.df %>% 
      filter(!ASV %in% emptyASVs)
    names(asvs_in_control.df) <- sample
    return(asvs_in_control.df)
  } else if(output == "saved"){
    return(asvs_in_control.df)
    warning("This option is under test!")
  }
}
