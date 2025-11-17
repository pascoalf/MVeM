# Load function to remove contamination, based on control map
source("./R/remove_contamination.R")

# Store sample names in a vector
sample_names <- sample_control_map_df$Sample_name %>% unique() 

# Remove contamination for all samples and re-merge in a single data frame
abundance_table_no_cont <- map(.x = sample_names, 
                               .f = ~remove_contamination(data = mi_table,
                                                          sample = .x)) %>% 
  bind_rows()

# To obtain a list of the ASVs that were considered contaminants in each sample
list_of_contaminants <- map(.x = sample_names, 
                            .f = ~remove_contamination(data = mi_table,
                                                       sample = .x, 
                                                       output = "contaminants")) %>% 
  bind_rows()

#
write.csv(abundance_table_no_cont, "./rc2_all_samples.csv")
write.csv(list_of_contaminants, "./rc2_all_contaminants.csv")
