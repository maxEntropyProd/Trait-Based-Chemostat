# Re-analyzing SES data using DADA2
# Updated by ANB on 8 January 2021
# From Fall 2018
# Retrieved data from VAMPS2: https://vamps2.mbl.edu/

# Libraries
library('dada2');packageVersion('dada2')

# DADA2 PROCESSING----
##REMEMBER TO CHANGE THIS DEPENDING ON THE ANALYSIS YOU ARE CONDUCTING:
path = "/Users/abulseco/Dropbox/MBL_Postdoc/Experiments/2018_Pilot_SESProject/Microbial_Data/Downloads_from_VAMPS/JJV_MEG_Bv4v5" 
setwd("~/Dropbox/MBL_Postdoc/Experiments/2018_Pilot_SESProject/Microbial_Data/Downloads_from_VAMPS/JJV_MEG_Bv4v5")

#sort files into R1 (forward) and R2 (reverse)
fnFs <- sort(list.files(path, pattern="R1.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="R2.fastq", full.names = TRUE))
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
head(sample.names)
View(sample.names)

# Inspect the quality plots
plotQualityProfile(fnFs[1:2]) 
plotQualityProfile(fnRs[1:2])

# Create a directory for the filtered data
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

# In this next step we will quality filter and remove primers. 

# It is common for chimeric sequences to be a majority (even large majority) of the total number of 
# unique sequence variants inferred by the dada(...) algorithm, but they should not be a majority of 
# the sequencing reads. That is, there may be many chimeras, but they are relatively low abundance.
# The most common reason that far too many reads are flagged as chimeric is that primer sequences were 
# not removed prior to starting the dada2 workflow. The ambiguous nucleotides in universal primer sequences 
# will be interpreted as real variation by the dada2 pipeline, and this interferes with the chimera algorithm. 
# In most cases, if you see over 25% of your sequencing reads being flagged chimeric, remove the primers from your 
# reads and restart the workflow with the clean data.

# Documentation on filterAndTrim: https://rdrr.io/bioc/dada2/man/filterAndTrim.html

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(250, 225),
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE, trimLeft = c(19, ),
                     compress=TRUE, multithread=TRUE)

head(out)
out
saveRDS(out,"/Users/abulseco/Dropbox/MBL_Postdoc/Experiments/2018_Pilot_SESProject/Microbial_Data/Downloads_from_VAMPS/JJV_MEG_Bv4v5/phyloseq_files_TRIMMED/out.rds")
# Use this when you need to read it back in
# out <- readRDS("/Users/abulseco/Dropbox/MBL_Postdoc/Experiments/2018_Pilot_SESProject/Microbial_Data/Downloads_from_VAMPS/JJV_MEG_Bv4v5/phyloseq_files/out.rds")

# Inspect the filtered profiles
plotQualityProfile(filtFs[1:2]) 
plotQualityProfile(filtRs[1:2])

# Learn errors and dereplicate
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)
saveRDS(errF,"/Users/abulseco/Dropbox/MBL_Postdoc/Experiments/2018_Pilot_SESProject/Microbial_Data/Downloads_from_VAMPS/JJV_MEG_Bv4v5/phyloseq_files_TRIMMED/errF.rds")
saveRDS(errR,"/Users/abulseco/Dropbox/MBL_Postdoc/Experiments/2018_Pilot_SESProject/Microbial_Data/Downloads_from_VAMPS/JJV_MEG_Bv4v5/phyloseq_files_TRIMMED/errR.rds")
# errF <- readRDS("/Users/abulseco/Dropbox/MBL_Postdoc/Experiments/2018_Pilot_SESProject/Microbial_Data/Downloads_from_VAMPS/JJV_MEG_Bv4v5/errF.rds")
# errR <- readRDS("/Users/abulseco/Dropbox/MBL_Postdoc/Experiments/2018_Pilot_SESProject/Microbial_Data/Downloads_from_VAMPS/JJV_MEG_Bv4v5/errR.rds")

# Plot the errors
plotErrors(errF, nominalQ=TRUE)
plotErrors(errR, nominalQ=TRUE)

# Dereplicate and save objects
derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)
saveRDS(derepFs,"/Users/abulseco/Dropbox/MBL_Postdoc/Experiments/2018_Pilot_SESProject/Microbial_Data/Downloads_from_VAMPS/JJV_MEG_Bv4v5/phyloseq_files_TRIMMED/derepFs.rds")
saveRDS(derepRs,"/Users/abulseco/Dropbox/MBL_Postdoc/Experiments/2018_Pilot_SESProject/Microbial_Data/Downloads_from_VAMPS/JJV_MEG_Bv4v5/phyloseq_files_TRIMMED/derepRs.rds")
# derepFs <- readRDS("/Users/abulseco/Dropbox/MBL_Postdoc/Experiments/2018_Pilot_SESProject/Microbial_Data/Downloads_from_VAMPS/JJV_MEG_Bv4v5/phyloseq_files_TRIMMED/derepFs.rds")
# derepRs <- readRDS("/Users/abulseco/Dropbox/MBL_Postdoc/Experiments/2018_Pilot_SESProject/Microbial_Data/Downloads_from_VAMPS/JJV_MEG_Bv4v5/phyloseq_files_TRIMMED/derepRs.rds")

# RUN DADA
dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
dadaRs <- dada(derepRs, err=errR, multithread=TRUE)
saveRDS(dadaFs, "/Users/abulseco/Dropbox/MBL_Postdoc/Experiments/2018_Pilot_SESProject/Microbial_Data/Downloads_from_VAMPS/JJV_MEG_Bv4v5/phyloseq_files_TRIMMED/dadaFs.rds")
saveRDS(dadaRs, "/Users/abulseco/Dropbox/MBL_Postdoc/Experiments/2018_Pilot_SESProject/Microbial_Data/Downloads_from_VAMPS/JJV_MEG_Bv4v5/phyloseq_files_TRIMMED/dadaRs.rds")
# dadaFs <- readRDS("/Users/abulseco/Dropbox/MBL_Postdoc/Experiments/2018_Pilot_SESProject/Microbial_Data/Downloads_from_VAMPS/JJV_MEG_Bv4v5/phyloseq_files_TRIMMED/dadaFs.rds")
# dadaRs <- readRDS("/Users/abulseco/Dropbox/MBL_Postdoc/Experiments/2018_Pilot_SESProject/Microbial_Data/Downloads_from_VAMPS/JJV_MEG_Bv4v5/phyloseq_files_TRIMMED/dadaRs.rds")

# MERGE PAIRS
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)
head(mergers[[1]])
saveRDS(mergers, "/Users/abulseco/Dropbox/MBL_Postdoc/Experiments/2018_Pilot_SESProject/Microbial_Data/Downloads_from_VAMPS/JJV_MEG_Bv4v5/phyloseq_files_TRIMMED/mergers.rds")
# mergers <- readRDS("/Users/abulseco/Dropbox/MBL_Postdoc/Experiments/2018_Pilot_SESProject/Microbial_Data/Downloads_from_VAMPS/JJV_MEG_Bv4v5/phyloseq_files_TRIMMED/mergers.rds")

# Make a sequence table
seqtab <- makeSequenceTable(mergers)
dim(seqtab) # 26 32422
table(nchar(getSequences(seqtab)))
saveRDS(seqtab, "/Users/abulseco/Dropbox/MBL_Postdoc/Experiments/2018_Pilot_SESProject/Microbial_Data/Downloads_from_VAMPS/JJV_MEG_Bv4v5/phyloseq_files_TRIMMED/seqtab.rds")
# seqtab <- readRDS("/Users/abulseco/Dropbox/MBL_Postdoc/Experiments/2018_Pilot_SESProject/Microbial_Data/Downloads_from_VAMPS/JJV_MEG_Bv4v5/phyloseq_files_TRIMMED/seqtab.rds")

# Then remove chimreas
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim) # 26 1702
sum(seqtab.nochim)/sum(seqtab) # 0.8707548
saveRDS(seqtab.nochim, "/Users/abulseco/Dropbox/MBL_Postdoc/Experiments/2018_Pilot_SESProject/Microbial_Data/Downloads_from_VAMPS/JJV_MEG_Bv4v5/phyloseq_files_TRIMMED/seqtab.nochim.rds")
# seqtab <- readRDS("/Users/abulseco/Dropbox/MBL_Postdoc/Experiments/2018_Pilot_SESProject/Microbial_Data/Downloads_from_VAMPS/JJV_MEG_Bv4v5/phyloseq_files/seqtab.nochim.rds")

# Track your sequences through the DADA2 pipeline
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)
track

# Assign taxonomy
taxa <- assignTaxonomy(seqtab.nochim, "silva_nr_v138_train_set.fa.gz", multithread=TRUE)
head(taxa)
saveRDS(taxa, "/Users/abulseco/Dropbox/MBL_Postdoc/Experiments/2018_Pilot_SESProject/Microbial_Data/Downloads_from_VAMPS/JJV_MEG_Bv4v5/phyloseq_files_TRIMMED/taxa.rds")

## PREPARE FOR PHYLOSEQ IMPORT
# PREPARING FOR PHYLOSEQ IMPORT
# ==================================================================

# Constructing a sample table that we can use for Phyloseq import
seqs <- getSequences(seqtab.nochim)
asv_seqs <- colnames(seqtab.nochim)
asv_headers <- vector(dim(seqtab.nochim)[2], mode="character")

for (i in 1:dim(seqtab.nochim)[2]) {
  asv_headers[i] <- paste(">ASV", i, sep="_")
}

# making and writing out a fasta of our final ASV seqs:
asv_fasta <- c(rbind(asv_headers, asv_seqs))
write(asv_fasta, "ASVs.fa")

# count table:
asv_tab <- t(seqtab.nochim)
row.names(asv_tab) <- sub(">", "", asv_headers)
write.table(asv_tab, "ASVs_counts.tsv", sep="\t", quote=F, col.names=NA)

# tax table:
asv_tax <- taxa
row.names(asv_tax) <- sub(">", "", asv_headers)
write.table(asv_tax, "ASVs_taxonomy.tsv", sep="\t", quote=F, col.names=NA)

# Then to make the tree, use ASV.fa file in FastTree with python. 