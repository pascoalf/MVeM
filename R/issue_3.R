# Check issue 3
# Ban list
all_hits <- read.csv("tmp_files/all_hits.csv", sep = ";")

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


# assign_LCA function
assign_LCA <- function(x){
  # make possible LCAs
  dom_LCA <- x %>% pull(Domain) %>% unique()
  phyl_LCA <- x %>%  pull(Phylum) %>% unique()
  class_LCA <- x %>% pull(Class) %>% unique()
  ord_LCA <- x %>%  pull(Order) %>% unique()
  fam_LCA <- x %>% pull(Family) %>% unique()
  genus_LCA <- x %>%  pull(Genus) %>% unique()
  species_LCA <- x %>%pull(Species) %>% unique()
  #
  taxonomic_levels <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")
  #
  if(length(dom_LCA) > 1){
    LCA <- "Uncertain"
    Level <- taxonomic_levels
  } else if(length(phyl_LCA) > 1){
    LCA <- dom_LCA
    Level <- taxonomic_levels[2:7]
  } else if(length(class_LCA) > 1){
    LCA <- phyl_LCA
    Level <- taxonomic_levels[3:7]
  } else if(length(ord_LCA) > 1){
    LCA <- class_LCA
    Level <- taxonomic_levels[4:7]
  } else if(length(fam_LCA) > 1){
    LCA <- ord_LCA
    Level <- taxonomic_levels[5:7]
  } else if(length(genus_LCA) > 1){
    LCA <- fam_LCA
    Level <- taxonomic_levels[6:7]
  } else if(length(species_LCA) > 1){
    LCA <- paste(genus_LCA, "sp.")
    Level <- taxonomic_levels[7]
  } else {
    LCA <- species_LCA
    Level <- NA
  }
  return(c(LCA, Level))
}

## check ASV_0027
top_hits %>% 
  filter(Query.accession == "ASV_0027") %>% 
  assign_LCA()
  

taxonomic_levels <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")

# Best hits, with LCA
#taxonomic_assignments <- 
test1 <- top_hits %>%
  group_by(Query.accession) %>% 
  nest() %>% 
  mutate(LCA = map(.x = data, 
                   .f = ~assign_LCA(.x)[1])) %>% 
  mutate(taxa = map(.x = data, .f = ~unique(.x$Species))) %>% 
  mutate(isTie = map(.x = taxa, .f = ~ifelse(length(unique(.x)) == 1, FALSE, TRUE))) %>%
  mutate(Level = map(.x = data, 
                     .f = ~assign_LCA(.x)[-1])) %>%  
  unnest(c(LCA, data, isTie)) %>% 
  group_by(Query.accession) %>% 
  slice_head(n = 1) %>% 
  mutate(FinalAssignment = ifelse(isTRUE(isTie), LCA, Species)) %>% 
  mutate(Species = ifelse(!is.na(Level), NA, Species)) %>% 
  mutate(Genus = case_when("Genus" %in% Level[[1]] ~ NA, TRUE ~ Genus)) %>% 
  mutate(Family = case_when("Family" %in% Level[[1]] ~ NA, TRUE ~ Family)) %>% 
  mutate(Order = case_when("Order" %in% Level[[1]] ~ NA, TRUE ~ Order)) %>% 
  mutate(Class = case_when("Class" %in% Level[[1]] ~ NA, TRUE ~ Class)) %>% 
  mutate(Phylum = case_when("Phylum" %in% Level[[1]] ~ NA, TRUE ~ Phylum)) %>%
  select(Query.accession,
         FinalAssignment, 
         Bit.Score, evalue, Alignment.length,
         Percentage.of.identical.matches,
         Number.of.identical.matches,
         Number.of.mismatches,
         Domain, Phylum, Class,
         Order, Family, Genus,
         Species)
###########

# Make data frame with unique ASVs ID and Sequence
ASVs.df <- ASV_table %>% 
  colnames() %>% 
  as.data.frame() %>% 
  rename(Sequence = ".") %>% 
  distinct() %>% 
  {n <- nrow(.)
  mutate(., ASV = paste0("ASV_", sprintf(paste0("%0", nchar(n), "d"), row_number()))) 
  } %>%
  arrange(ASV)
##
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


filt_tax_assignments <- tax_assign_merged %>%
  filter(Class %in% c("Mammalia", "Actinopteri", "Chondrichthyes"))


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


#############
table1_tmp <- readxl::read_xlsx("tmp_files/Table_1.xlsx")

table1_tmp_long <- table1_tmp %>% 
  pivot_longer(cols = names(table1_tmp)[3:322],
               names_to = "Sample",
               values_to = "Abundance")

table1_tmp_long_filtered <- table1_tmp_long %>% 
  group_by(Sample) %>% 
  mutate(relativeAbundance = Abundance*100/sum(Abundance)) %>% 
  mutate(Abundance = ifelse(Abundance == 1, 0, Abundance),
         Abundance = ifelse(relativeAbundance > 0.01, Abundance, 0),
         Abundance = ifelse(is.na(Abundance), 0, Abundance)) 

# Load function to remove contamination, based on control map
source("./R/remove_contamination.R")

# Store sample names in a vector
sample_names <- sample_control_map_df$Sample_name %>% unique() 


# Example without threshold
example_automatic <- map(.x = sample_names, 
                         .f = ~remove_contamination(data = table1_tmp_long_filtered, 
                                                    sample = .x, 
                                                    output = "standard",
                                                    option = "automatic")) %>% 
  bind_rows()

# If you want to verify which ASVs were considered contaminants without thresholds
contaminants_automatic <- map(.x = sample_names,
                              .f = ~remove_contamination(data = abundance_table_long_filtered,
                                                         sample = .x,
                                                         option = "automatic",
                                                         output = "contaminants")) %>% 
  bind_rows()

# If you want to verify which ASVs were saved
saved_automatic <- map(.x = sample_names,
                       .f = ~remove_contamination(data = abundance_table_long_filtered,
                                                  sample = .x,
                                                  option = "automatic",
                                                  output = "saved")) %>% 
  bind_rows()

