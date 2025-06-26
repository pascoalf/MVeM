# make FASTA file
library(seqinr)
library(dplyr)

# load ASV table
ASV_table <- read.table("./input/ASV_table1.tsv")

# make data frame with unique ASVs ID and Sequence
ASVs.df <- ASV_table %>% 
    colnames() %>% 
    as.data.frame() %>% 
    rename(Sequence = ".") %>% 
    distinct() %>% 
  mutate(ASV = paste0("ASV_", row_number(.)))

# Make FASTA file 
write.fasta(sequences = as.list(ASVs.df$Sequence), 
            names = ASVs.df$ASV, 
            "./input/ASV.fasta",
            as.string = TRUE)

