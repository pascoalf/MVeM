# load packages
library(dplyr)
library(tidyr)
library(ulrb)
library(stringr)
library(purrr)

# load blast results
all_hits <- read.csv("../eDNA/blast_results_atlantida", header = FALSE, # change file path as needed
                     col.names = c("Query accession", "Query sequence length",
                                   "Subject seq-id",    "Subject accession",
                                   "Subject sequence length",   "evalue", "Bit Score",
                                   "Raw Score", "Alignment length", "Percentage of identical matches",
                                   "Number of identical matches",   "Number of mismatches",
                                   "Number of positive scoring matches", "Total number of gaps",
                                   "Taxonomy ID", "Scientific name",    "Subject blast name", 
                                   "Subject common name"))

# ban list
ban_list <- read.table("./refs/ban_list.txt", header = FALSE) %>% 
  rename(Genus = V1,
         Species = V2) %>% 
  mutate(binomial_name = paste(Genus, Species)) %>% 
  pull(binomial_name)

# target genes
target_genes <- read.table("../eDNA/lista_accessions_16S.seq", header = FALSE)
# some data cleaning
target_genes <- target_genes %>% 
  rename(Subject.accession = V1) %>% 
  mutate(Subject.accession = str_remove(Subject.accession, "\\.\\d+"))

filtered_hits <- all_hits %>%
  filter(Alignment.length >= 190,
         !Scientific.name %in% ban_list,
         Subject.accession %in% target_genes$Subject.accession,
         Subject.blast.name %in% c("bony fishes", "whales & dolphins", "sharks & rays")) %>% 

# Obtain top hits and Remove environmental samples hits before summarizing
top_hits <- filtered_hits %>%
  # Remove environmental sample rows before anything else
  filter(!grepl("environmental sample", Scientific.name, ignore.case = TRUE)) %>% 
  
  # Normalize to first two words for species-level matching
  mutate(Scientific.name = sub("^([A-Za-z]+\\s+[A-Za-z]+).*", "\\1", Scientific.name)) %>%
  
  group_by(Query.accession) %>%
  filter(Bit.Score == max(Bit.Score)) %>%
  filter(Percentage.of.identical.matches == max(Percentage.of.identical.matches)) %>% 
  ungroup()


source("check_ties.R")
source("assign_LCA.R")

# Add reference for families
delphinidae_family <- read.table("delphinidae_family.txt"); names(delphinidae_family) <- "Genus"
pleuronectidae_family <- read.table("pleuronectidae_family.txt"); names(pleuronectidae_family) <- "Genus"
ziphiidae_family <- read.table("ziphiidae_family.txt"); names(ziphiidae_family) <- "Genus"
salmonidae_family <- read.table("Salmonidae_family.txt"); names(salmonidae_family) <- "Genus"
mugilidae_family <- read.table("Mugilidae_family.txt"); names(mugilidae_family) <- "Genus"

# best hits, with LCA
taxonomic_assignments <- top_hits %>%
  group_by(Query.accession) %>% 
  nest() %>% 
  mutate(LCA = map(.x = data, 
                   .f = ~assign_LCA(.x))) %>% 
  mutate(taxa = map(.x = data, .f = ~unique(.x$Scientific.name))) %>% 
  mutate(isTie = map(.x = taxa, .f = ~ifelse(length(unique(.x)) == 1, FALSE, TRUE))) %>% 
  unnest(c(LCA, data, isTie)) %>% 
  group_by(Query.accession) %>% 
  slice_head(n = 1) %>% 
  mutate(FinalAssignment = ifelse(isTRUE(isTie), LCA, Scientific.name)) %>% 
  select(Query.accession,
         Scientific.name, 
         LCA, FinalAssignment, 
         Bit.Score, evalue, Alignment.length,
         Percentage.of.identical.matches,
         Number.of.identical.matches,
         Number.of.mismatches)

# view results in your R session
View(taxonomic_assignments)

# Save final assignments into memory
#write.csv(taxonomic_assignments, "taxonomic_assignments_atlantida.csv")

# transform blast results to compatible format
ASV_ncbi <- taxonomic_assignments %>% 
  select(ASV = Query.accession, FinalAssignment) %>% 
  filter(!is.na(FinalAssignment)) %>%  # Remove unassigned ASVs
  left_join(ASVs.df, by = "ASV") # ASVs.df was made in the DADA2 section

# Create abundance table (long format)
abundance_table_long <- ASV_table %>% # ASV_table was made in DADA2 section
  prepare_tidy_data(sample_names = row.names(ASV_table), samples_in = "rows") %>% 
  rename(Sequence = Taxa_id) %>% 
  left_join(ASVs.df, by = "Sequence") %>% 
  left_join(ASV_ncbi, by = "ASV")

# Filter local low abundance (< 0.01%)
total_reads <- abundance_table_long %>%
  group_by(Sample) %>%
  summarise(total = sum(Abundance, na.rm = TRUE))

abundance_table_long_filtered <- abundance_table_long %>%
  left_join(total_reads, by = "Sample") %>%
  mutate(freq = Abundance / total * 100,
         Abundance = ifelse(freq < 0.01, 0, Abundance)) %>%
  select(-total, -freq)

## Create control map, connecting samples to their controls
sample_control_map <- list(
  "M1-1-16S_S1_L001_R1_001" = c("CE1-16S_S1_L001_R1_001", "CF1-1-16S_S1_L001_R1_001"),
  "M1-2-16S_S1_L001_R1_001" = c("CE1-16S_S1_L001_R1_001", "CF1-2-16S_S1_L001_R1_001"),
  "M1-3-16S_S1_L001_R1_001" = c("CE1-16S_S1_L001_R1_001", "CF1-3-16S_S1_L001_R1_001"),
  "M2-1-16S_S1_L001_R1_001" = c("CE2-16S_S1_L001_R1_001", "CF2-1-16S_S1_L001_R1_001"),
  "M2-1-NZY-16S_S1_L001_R1_001" = c("CE2-16S_S1_L001_R1_001", "CF2-1-16S_S1_L001_R1_001"),
  "M2-2-16S_S1_L001_R1_001" = c("CE2-16S_S1_L001_R1_001", "CF2-2-16S_S1_L001_R1_001"),
  "M2-2-NZY-16S_S1_L001_R1_001" = c("CE2-16S_S1_L001_R1_001", "CF2-2-16S_S1_L001_R1_001"),
  "M2-3-16S_S1_L001_R1_001" = c("CE2-16S_S1_L001_R1_001", "CF2-3-16S_S1_L001_R1_001"),
  "M2-3-NZY-16S_S1_L001_R1_001" = c("CE2-16S_S1_L001_R1_001", "CF2-3-16S_S1_L001_R1_001"),
  "M3-1-16S_S1_L001_R1_001" = c("CE3-16S_S1_L001_R1_001", "CF2-1-16S_S1_L001_R1_001"),
  "M3-2-16S_S1_L001_R1_001" = c("CE3-16S_S1_L001_R1_001", "CF2-2-16S_S1_L001_R1_001"),
  "M3-3-16S_S1_L001_R1_001" = c("CE3-16S_S1_L001_R1_001", "CF2-3-16S_S1_L001_R1_001")
)

## Identify ASVs present in control samples
asvs_in_controls <- abundance_table_long_filtered %>%
  filter(Sample %in% unlist(sample_control_map),
         Abundance > 0) %>%
  distinct(Sample, ASV)

## Remove those ASVs from their corresponding environmental samples
for (sample_name in names(sample_control_map)) {
  controls <- sample_control_map[[sample_name]]
  
  contaminant_asvs <- asvs_in_controls %>%
    filter(Sample %in% controls) %>%
    pull(ASV) %>%
    unique()
  
  abundance_table_long_filtered <- abundance_table_long_filtered %>%
    mutate(Abundance = ifelse(Sample == sample_name & ASV %in% contaminant_asvs, 0, Abundance))
}

# Save final long-format abundance table
write.csv(abundance_table_long_filtered, file = "abundance_table_long_atlantida.csv", row.names = FALSE)

# Create wide-format abundance table
abundance_table_wide <- abundance_table_long_filtered %>% 
  pivot_wider(names_from = Sample, values_from = Abundance)

# Save wide-format table
write.csv(abundance_table_wide, file = "abundance_table_wide_atlantida.csv", row.names = FALSE)
