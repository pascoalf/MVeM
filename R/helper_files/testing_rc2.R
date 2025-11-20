# Load function to remove contamination, based on control map
source("./R/remove_contamination.R")

# Store sample names in a vector
sample_names <- sample_control_map_df$Sample_name %>% unique() 

# Remove contamination for all samples and re-merge in a single data frame
abundance_table_no_cont <- map(.x = sample_names[1:5], 
                               .f = ~remove_contamination(data = mi_table,
                                                          sample = .x, 
                                                          output = "standard",
                                                          option = "automatic")) %>% 
  bind_rows()

# To obtain a list of the ASVs that were considered contaminants in each sample
list_of_contaminants <- map(.x = sample_names, 
                            .f = ~remove_contamination(data = mi_table,
                                                       sample = .x, 
                                                       output = "contaminants", 
                                                       option = "automatic")) %>% 
  bind_rows()

#
write.csv(abundance_table_no_cont, "./rc2_all_samples.csv")
write.csv(list_of_contaminants, "./rc2_all_contaminants.csv")


source("./R/remove_contamination.R"); map(.x = sample_names, 
    .f = ~remove_contamination(data = mi_table,
                               sample = .x, 
                               output = "saved", 
                               option = "automatic"))


remove_contamination(data = mi_table, 
                     sample = "M2-3-16S_S1_L001_R1_001", 
                     map_sample = sample_control_map_long,
                     output = "saved",
                     option = "automatic")
#


