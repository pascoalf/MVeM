# rarefaction curve
library(vegan)
library(dplyr)

# Load ASV table (from DADA2 output)
# The ASV.table is a TSV file where samples are rows and ASV sequences are columns
ASV_rarefaction <- read.delim("ASV_table1_eDNA.tsv", header = TRUE, row.names = 1, sep = "\t", check.names = FALSE)

# Replace NA's with 0
ASV_rarefaction[is.na(ASV_rarefaction)] <- 0

# Ensure all entries are numeric (in case they were read as characters)
ASV_rarefaction <- apply(ASV_rarefaction, 2, as.numeric)
rownames(ASV_rarefaction) <- rownames(read.delim("ASV_table1_eDNA.tsv", header = TRUE, sep = "\t", check.names = FALSE, row.names = 1))

# Replace sample names to shorter version
rownames(ASV_rarefaction) <- str_remove(rownames(ASV_rarefaction), "-16S_S1_L001_R1_001")

# Rarefaction curve
rarecurve(
  ASV_rarefaction,
  step = 500,
  xlab = "Sequencing depth",
  ylab = "Number of ASVs"
)
