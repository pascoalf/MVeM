# make FASTA file
library(seqinr)
library(dplyr)

# load ASV table
ASV_table <- read.table("./results/ASV_table.tsv")

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

# Make FASTA file 
write.fasta(sequences = as.list(ASVs.df$Sequence), 
            names = ASVs.df$ASV, 
            "./results/ASV.fasta",
            as.string = TRUE)

