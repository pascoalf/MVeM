# contamination bug
library(tidyr)
library(dplyr)
library(purrr)

mi_table <- read.csv("./results/abundance_table_long_filtered.csv")


sample_control_map_df <- readxl::read_xlsx("refs/sample_control_map_example_complete.xlsx")

# Convert to long format 
sample_control_map_long <- sample_control_map_df %>% 
  pivot_longer(cols = c("Extraction_control", 
                        "Filtration_control", 
                        "PCR_control"),
               values_to = "Control_ID",
               names_to = "Control_type")

View(mi_table)

#
sample_names <- sample_control_map_df$Sample_name %>% unique() 

#
#mi_abundance_table_no_cont <- 
mi_check <- map(.x = sample_names, .f = ~remove_contamination(data = mi_table, sample = .x, 
                                                  map_sample = sample_control_map_long,
                                                  threshold = 1000,
                                                  option = "automatic")) %>% 
  bind_rows()


## check all options
# no threshold 
map(.x = sample_names, .f = ~remove_contamination(data = mi_table, sample = .x, 
                                                  map_sample = sample_control_map_long,
                                                  threshold = NULL,
                                                  output = "contaminants",
                                                  option = NULL))
map(.x = sample_names, .f = ~remove_contamination(data = mi_table, sample = .x, 
                                                  map_sample = sample_control_map_long,
                                                  threshold = NULL,
                                                  output = "standard",
                                                  option = NULL))
map(.x = sample_names, .f = ~remove_contamination(data = mi_table, sample = .x, 
                                                  map_sample = sample_control_map_long,
                                                  threshold = 1000,
                                                  output = "standard",
                                                  option = "threshold"))
map(.x = sample_names, .f = ~remove_contamination(data = mi_table, sample = .x, 
                                                  map_sample = sample_control_map_long,
                                                  threshold = 1000,
                                                  output = "contaminants",
                                                  option = "threshold"))
map(.x = sample_names, .f = ~remove_contamination(data = mi_table, sample = .x, 
                                                  map_sample = sample_control_map_long,
                                                  threshold = NULL,
                                                  output = "contaminants",
                                                  option = NULL))

##  
