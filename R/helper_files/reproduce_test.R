#
source("./R/remove_contamination.R")

table_1 <- read.csv("./results/Table_1.csv", sep = ";")

#
table_1_long <- table_1 %>% 
  pivot_longer(cols = names(table_1)[3:31],
    names_to = "Sample",
               values_to = "Abundance")


# Store sample names in a vector
sample_names <- sample_control_map_df$Sample_name %>% unique()

# Example without threshold
example_automatic <- map(.x = sample_names, 
                         .f = ~remove_contamination(data = abundance_table_long_filtered, 
                                                    sample = .x, 
                                                    output = "standard",
                                                    option = "automatic")) %>% 
  bind_rows()

# If you wanto to verify which ASVs were considered contaminants without thresholds
contaminants_automatic <- map(.x = sample_names,
                              .f = ~remove_contamination(data = abundance_table_long_filtered,
                                                         sample = .x,
                                                         option = "automatic",
                                                         output = "contaminants")) %>% 
  bind_rows()

#
saved_automatic <- map(.x = sample_names,
                       .f = ~remove_contamination(data = abundance_table_long_filtered,
                                                  sample = .x,
                                                  option = "automatic",
                                                  output = "saved")) %>% 
  bind_rows()
