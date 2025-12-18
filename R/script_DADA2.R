# Packages used
library(dada2); packageVersion("dada2") ## we used 1.22
library(ShortRead)
library(seqinr) # to make FASTA file
library(dplyr)
library(ggplot2)

path <- "./path_to_directory" # CHANGE ME to the directory containing the fastq files after unzipping (and after Cutadapt trimming if applied).
# Verify files in path
list.files(path)

# Forward and reverse fastq file names have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
# CHANGE according to your file names
# note: fastq.gz files usually don't need to be decompressed for this step
fnFs <- sort(list.files(path, pattern = "_R1_001.fastq", full.names = TRUE))
fnFs
fnRs <- sort(list.files(path, pattern = "_R2_001.fastq", full.names = TRUE))
fnRs

# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.namesF <- sapply(strsplit(basename(fnFs), "[.]"), `[`, 1)
sample.namesF
sample.namesR <- sapply(strsplit(basename(fnRs), "[.]"), `[`, 1)
sample.namesR

# Place filtered files in filtered/ subdirectory
filtFs <- file.path(path, "filtered", paste0(sample.namesF, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.namesR, "_R_filt.fastq.gz"))
names(filtFs) <- sample.namesF
names(filtRs) <- sample.namesR

# Example for 5 samples
# Quality of forward reads
plotQualityProfile(fnFs[1:5])
# Quality of reverse reads
plotQualityProfile(fnRs[1:5])

# To inspect many samples at once and save the plots
QProfile_Fw <- plotQualityProfile(fnFs, aggregate = TRUE)
QProfile_Fw
ggsave("./results/QProfile_Fw.tiff", 
       plot = QProfile_Fw, 
       width = 15, 
       height = 12, 
       dpi = 600)

QProfile_Rv <- plotQualityProfile(fnRs, aggregate = TRUE)
QProfile_Rv
ggsave("./results/QProfile_Rv.tiff", 
       plot = QProfile_Rv, 
       width = 15, 
       height = 12, 
       dpi = 600)

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, 
                     truncLen = c(170,150), ## change according to quality profiles 
                     maxN = 0, maxEE = c(2,2), truncQ = 2, rm.phix = TRUE, 
                     compress = TRUE, multithread = FALSE, # On Windows set multithread=FALSE
                     ## OPTIONAL: if you need to remove primers at this stage, you can use trimLeft
                     #trimLeft = c(nchar("AGACGAGAAGACCCTATG"),                      
                     #            nchar("GGATTGCGCTGTTATCCC"))
                     )

# Learn the Error Rates
errF <- learnErrors(filtFs, multithread = FALSE)
errR <- learnErrors(filtRs, multithread = FALSE)
plotErrors(errF, nominalQ=TRUE)
plotErrors(errR, nominalQ=TRUE)

derepFs <- derepFastq(filtFs, verbose = TRUE)
names(derepFs) <- sample.namesF
derepRs <- derepFastq(filtRs, verbose = TRUE)
names(derepRs) <- sample.namesR

# Identify unique sequences
dadaFs <- dada(derepFs, err = errF, multithread = FALSE)
dadaRs <- dada(derepRs, err = errR, multithread = FALSE)

dadaFs[[1]]
dadaRs[[1]]

# Merge paired reads
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose = TRUE)
# Inspect the merger data.frame from the first sample
head(mergers[[1]])

# Construct sequence table
seqtab <- makeSequenceTable(mergers)
dim(seqtab)

# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))
hist(nchar(getSequences(seqtab)), main = "Distribution of Sequence lengths")

# Remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method = "consensus", multithread = FALSE, verbose = TRUE)

# Check percentage of non-chimeric sequences
sum(seqtab.nochim)/sum(seqtab)

# Track reads through the pipeline
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)

colnames(track) <- c("Input reads", "Filtered reads", "Denoised Fw reads", "Denoised Rv reads", "Merged reads", "Non-chimeric reads")
rownames(track) <- sample.namesF ## sample.namesF is just to indicate the sample ID
rownames(track) <- gsub("-16S_S1_L001_R1_001", "", rownames(track)) # Change according to your sample names
head(track)

# Save summary table
write.csv(track, file = './results/Summary_table.csv')

# Change object name
ASV_table <- seqtab.nochim

# Create .csv file
write.table(ASV_table, file='./results/ASV_table.tsv', quote = FALSE, sep = '\t', col.names = NA)

