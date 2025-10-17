library(dplyr)
library(stringr)
mito_ref <- read.table("./refs/Accessions_Chordata_mithocondrion.txt", header = FALSE)
g16S_ref <- read.table("./refs/gene_16_list.txt", header = FALSE)

#
accession_mito <- mito_ref %>% 
  rbind(g16S_ref)

write.table(accession_mito, "./refs/accession_mitochondrial_list.txt", row.names = FALSE)


target_genes <- read.table("refs/accession_mitochondrial_list.txt", header = TRUE) ## last accessed 16 Oct 2025

target_genes <- target_genes %>% 
  mutate(Subject.accession = str_remove(Subject.accession, "\\.\\d+"))
