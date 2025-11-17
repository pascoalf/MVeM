# contamination bug
library(tidyr)
library(dplyr)
library(purrr)
library(ulrb)
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
                                                  threshold = 0.01,
                                                  output = "contaminants",
                                                  option = "threshold"))
map(.x = sample_names, .f = ~remove_contamination(data = mi_table, sample = .x, 
                                                  map_sample = sample_control_map_long,
                                                  threshold = 4539,
                                                  output = "contaminants",
                                                  option = "threshold"))
map(.x = sample_names, .f = ~remove_contamination(data = mi_table, sample = .x, 
                                                  map_sample = sample_control_map_long,
                                                  output = "contaminants",
                                                  option = "automatic"))

map(.x = sample_names, .f = ~remove_contamination(data = mi_table, sample = .x, 
                                                  map_sample = sample_control_map_long,
                                                  output = "contaminants"))

## Manually verify
mi_table %>% 
  filter(Sample == sample_names[11]) %>% 
  View()
  

### check what happens if we use remove contaminants from control samples
mi_data2 <- read.csv("~/Downloads/abundance_table_long_filtered.csv")

head(mi_data2$Sample)

mi_data2_nozero <- mi_data2 %>% filter(Abundance > 0)

map(.x = unique(mi_data2$Sample), .f = ~remove_contamination(data = mi_data2, 
                                                  sample = .x, 
                                                  map_sample = sample_control_map_long,
                                                  output = "contaminants",
                                                  option = "automatic"))

remove_contamination(data = mi_data2, 
                     sample = "M3-2-16S_S1_L001_R1_001", 
                     map_sample = sample_control_map_long,
                     output = "contaminants",
                     option = "automatic")


## testing new remove contaminants
temp1_test <- remove_contamination(data = mi_data2, 
                     sample = "M3-2-16S_S1_L001_R1_001", 
                     map_sample = sample_control_map_long,
                     output = "contaminants",
                     option = "automatic") %>% as.data.frame()

write.csv(temp1_test, "temp1_test.csv")



source("R/remove_contamination.R");
remove_contamination(data = mi_data2, 
                     sample = "M1-1-16S_S1_L001_R1_001", 
                     map_sample = sample_control_map_long,
                     output = "contaminants",
                     option = "automatic")


