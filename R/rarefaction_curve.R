# rarefaction curve
library(vegan)
library(dplyr)

# remove unnecessary columns
ASV_matrix.1 <- abundance_table_wide %>% 
  select(-Sequence, -Scientific.name) 

#
asc_col <- ASV_matrix.1$ASV
ASV_matrix.1$ASV <- NULL
rownames(ASV_matrix.1) <- asc_col

#  
ASV_matrix <- ASV_matrix.1 %>% t()
# rarefaction curve
rarecurve(ASV_matrix, 
          step = 500, 
          xlab = "Sequencing depth",
          ylab = "Number of ASVs")
