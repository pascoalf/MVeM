# contamination bug
library(tidyr)
library(dplyr)
library(purrr)

mi_table <- read.csv("./results/abundance_table_long_filtered.csv")

View(mi_table)

#
sample_names <- sample_control_map_df$Sample_name %>% unique() 

#
#mi_abundance_table_no_cont <- 
mi_check <- map(.x = sample_names, .f = ~remove_contamination(data = mi_table, sample = .x, 
                                                  map_sample = sample_control_map_long,
                                                  threshold = 10000,
                                                  option = "automatic")) %>% 
  bind_rows()


map(.x = sample_names, .f = ~remove_contamination(data = mi_table, sample = .x, 
                                                  map_sample = sample_control_map_long,
                                                  threshold = 200, option = "automatic",
                                                  output = "contaminants"))

##  
