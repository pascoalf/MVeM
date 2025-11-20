# load packages
library(dplyr)
library(tidyr)
library(ulrb)
library(stringr)
library(purrr)

# Load blast results
all_hits <- read.table("./results/blast_results_taxonomy", header = FALSE, sep = "\t", # change file path as needed
                     col.names = c("Query accession", "Query sequence length",
                                   "Subject seq-id",    "Subject accession",
                                   "Subject sequence length",   "evalue", "Bit Score",
                                   "Raw Score", "Alignment length", "Percentage of identical matches",
                                   "Number of identical matches",   "Number of mismatches",
                                   "Number of positive scoring matches", "Total number of gaps",
                                   "Taxonomy ID", "Scientific name",    "Subject blast name", "Subject common name",
                                   "Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"))

# Ban list
ban_list <- read.table("./refs/ban_list.txt", header = FALSE) %>% 
  rename(Genus = V1,
         Species = V2) %>% 
  mutate(binomial_name = paste(Genus, Species)) %>% 
  pull(binomial_name)

# Target genes
target_genes <- read.table("refs/accession_mitochondrial_list.txt", header = TRUE) ## last accessed 16 Oct 2025

# Filter valid hits
filtered_hits <- all_hits %>%
  filter(Alignment.length >= 190,
         !Scientific.name %in% ban_list,
         Subject.accession %in% target_genes$Subject.accession,
         # Remove environmental sample rows
         !grepl("environmental sample", Species, ignore.case = TRUE)) %>%
         # Normalize to first two words for species-level matching
         mutate(Scientific.name = sub("^([A-Za-z]+\\s+[A-Za-z]+).*", "\\1", Species))

# Obtain top hits and remove environmental samples hits before summarizing
top_hits <- filtered_hits %>%
  group_by(Query.accession) %>%
  filter(Bit.Score == max(Bit.Score)) %>%
  filter(Percentage.of.identical.matches == max(Percentage.of.identical.matches)) %>% 
  ungroup()

source("check_ties.R")
source("assign_LCA.R")

# Best hits, with LCA
taxonomic_assignments <- top_hits %>%
  group_by(Query.accession) %>% 
  nest() %>% 
  mutate(LCA = map(.x = data, 
                   .f = ~assign_LCA(.x)[1])) %>% 
  mutate(taxa = map(.x = data, .f = ~unique(.x$Species))) %>% 
  mutate(isTie = map(.x = taxa, .f = ~ifelse(length(unique(.x)) == 1, FALSE, TRUE))) %>%
  mutate(Level = map(.x = data, 
                     .f = ~assign_LCA(.x)[2])) %>% 
  unnest(c(LCA, Level, data, isTie)) %>% 
  group_by(Query.accession) %>% 
  slice_head(n = 1) %>% 
  mutate(FinalAssignment = ifelse(isTRUE(isTie), LCA, Species)) %>% 
  mutate(Species = ifelse(!is.na(Level), NA, Species)) %>% 
  mutate(Genus = case_when(Level == "Genus" ~ NA, TRUE ~ Genus)) %>% 
  mutate(Family = case_when(Level == "Family" ~ NA, TRUE ~ Family)) %>% 
  select(Query.accession,
         FinalAssignment, 
         Bit.Score, evalue, Alignment.length,
         Percentage.of.identical.matches,
         Number.of.identical.matches,
         Number.of.mismatches,
         Domain, Phylum, Class,
         Order, Family, Genus,
         Species)

tax_assign_merged <- taxonomic_assignments %>% 
  select(ASV = Query.accession,
         FinalAssignment, 
         Bit.Score, evalue, Alignment.length,
         Percentage.of.identical.matches,
         Number.of.identical.matches,
         Number.of.mismatches,
         Domain, Phylum, Class,
         Order, Family, Genus,
         Species) %>% 
  filter(!is.na(FinalAssignment)) %>%  # Remove unassigned ASVs
  left_join(ASVs.df, by = "ASV") # ASVs.df was made in the DADA2 section

# View results in your R session
View(tax_assign_merged)

# Save taxonomic assignments into memory
write.csv(tax_assign_merged, "results/taxonomic_assignments.csv", row.names = FALSE)

filt_tax_assignments <- tax_assign_merged %>%
  filter(Class %in% c("Mammalia", "Actinopteri", "Chondrichthyes"))

# Save final taxonomic assignments into memory
write.csv(filt_tax_assignments, "results/taxonomic_assignments_filtered.csv", row.names = FALSE)

# Create abundance table in long format
abundance_table_long <- ASV_table %>% # ASV_table was made in DADA2 section
  prepare_tidy_data(sample_names = row.names(ASV_table), samples_in = "rows") %>% 
  rename(Sequence = Taxa_id) %>% 
  left_join(ASVs.df, by = "Sequence") %>% 
  left_join(filt_tax_assignments, by = "ASV")

# Creates abundance table in wide format
abundance_table_wide <- abundance_table_long %>% 
  filter(!is.na(FinalAssignment)) %>%
  pivot_wider(names_from = Sample, values_from = Abundance)

table_1 <- abundance_table_wide %>%
  select(ASV, FinalAssignment,
         18:last_col(),
         Domain, Phylum, Class, Order, Family, Genus, Species) %>%
  arrange(ASV)
colnames(table_1) <- gsub("-16S_S1_L001_R1_001", "", colnames(table_1))

# Save wide format abundance table
write.csv(table_1, "results/Table_1.csv", row.names = FALSE)

# Filter local low abundance (< 0.01%)
abundance_table_long_filtered <- abundance_table_long %>% 
  group_by(Sample) %>% 
  mutate(relativeAbundance = Abundance*100/sum(Abundance)) %>% 
  mutate(Abundance = ifelse(Abundance == 1, 0, Abundance),
         Abundance = ifelse(relativeAbundance > 0.01, Abundance, 0),
         Abundance = ifelse(is.na(Abundance), 0, Abundance)) %>% 
  select(-Sequence.y, -relativeAbundance) %>% 
  rename(Sequence = Sequence.x)

# Load sample_control_map
sample_control_map_df <- readxl::read_xlsx("refs/sample_control_map_example_complete.xlsx")

# Convert to long format 
sample_control_map_long <- sample_control_map_df %>% 
  pivot_longer(cols = c("Extraction_control", 
                        "Filtration_control", 
                        "PCR_control"),
               values_to = "Control_ID",
               names_to = "Control_type")

# Load function to remove contamination, based on control map
source("./R/remove_contamination.R")

# Store sample names in a vector
sample_names <- sample_control_map_df$Sample_name %>% unique() 

# Remove contamination for all samples and re-merge in a single data frame
abundance_table_no_cont <- map(.x = sample_names, 
                               .f = ~remove_contamination(data = abundance_table_long_filtered,
                                                          sample = .x)) %>% 
  bind_rows()

# To obtain a list of the ASVs that were considered contaminants in each sample
list_of_contaminants <- map(.x = sample_names, 
                            .f = ~remove_contamination(data = abundance_table_long_filtered,
                                                       sample = .x, 
                                                       output = "contaminants")) %>% 
  bind_rows()
  
# Convert to wide format
abundance_table_no_cont_wide <- abundance_table_no_cont %>% 
  pivot_wider(names_from = Sample, values_from = Abundance)

table_2 <- abundance_table_no_cont_wide %>%
  select(ASV, FinalAssignment,
         18:last_col(),
         Domain, Phylum, Class, Order, Family, Genus, Species) %>%
  filter(!is.na(FinalAssignment)) %>% 
  mutate(across(where(is.numeric), ~replace_na(.x, 0))) %>%
  arrange(ASV)
colnames(table_2) <- gsub("-16S_S1_L001_R1_001", "", colnames(table_2))

# Save the output table
write.csv(table_2, file = "results/Table_2.csv", row.names = FALSE)
