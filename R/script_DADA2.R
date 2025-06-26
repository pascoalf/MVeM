library(dada2); packageVersion("dada2")
library(ShortRead)

setwd("C:/Users/cetus/Desktop/ATLANTIDA_eDNA")
getwd()

path <- "C:/Users/cetus/Desktop/ATLANTIDA_eDNA" # CHANGE ME to the directory containing the fastq files after unzipping.
list.files(path)

# Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
fnFs <- sort(list.files(path, pattern="_R1_001.fastq", full.names = TRUE))
fnFs
fnRs <- sort(list.files(path, pattern="_R2_001.fastq", full.names = TRUE))
fnRs

# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.namesF <- sapply(strsplit(basename(fnFs), "[.]"), `[`, 1)
sample.namesF
sample.namesR <- sapply(strsplit(basename(fnRs), "[.]"), `[`, 1)
sample.namesR

#Inspect read quality profiles
plotQualityProfile(fnFs[1:5])
plotQualityProfile(fnRs[1:5])

#F limit=240 bps, Rlimit=210 bps

# Place filtered files in filtered/ subdirectory
filtFs <- file.path(path, "filtered", paste0(sample.namesF, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.namesR, "_R_filt.fastq.gz"))
names(filtFs) <- sample.namesF
names(filtRs) <- sample.namesR


out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(240,210),
              maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
              compress=TRUE, multithread=FALSE) # On Windows set multithread=FALSE
head(out)

#Learn the Error Rates
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)
plotErrors(errF, nominalQ=TRUE)
plotErrors(errR, nominalQ=TRUE)

#Sample Filtering
dadaFs <- dada(filtFs, err=errF, multithread=TRUE)
dadaRs <- dada(filtRs, err=errR, multithread=TRUE)

dadaFs[[1]]
dadaRs[[1]]

#Merge paired reads
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)
# Inspect the merger data.frame from the first sample
head(mergers[[1]])

#Construct sequence table
seqtab <- makeSequenceTable(mergers)
dim(seqtab)

# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))
# Remove non-target-length sequences from your sequence table
seqtab2 <- seqtab[,nchar(colnames(seqtab)) %in% 200:300]
table(nchar(getSequences(seqtab2)))

#Remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab2, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)
table(nchar(getSequences(seqtab.nochim)))

sum(seqtab.nochim)/sum(seqtab2)

#Track reads through the pipeline
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))

# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)

#Create .csv file
write.table(seqtab.nochim, file='ASV_table1.tsv', quote=FALSE, sep='\t', col.names = NA)

#Download FASTA ASV
sq <- getSequences(seqtab.nochim)
writeFasta(sq, "ASV.fa")