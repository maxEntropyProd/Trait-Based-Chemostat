# Final Phyloseq file
# Use this file to import data into Phyloseq 
# Imported data analyzed "2021_01_SES_DADA2_Reprocess.R"

# Import Libraries
library(phyloseq)
library(data.table)
library(ggplot2)
library(dplyr)

# PHYLOSEQ ANALYSIS
# ===========================================================================
# Loading into Phyloseq----
mat = read.table("ASVs_counts.tsv", header = TRUE, sep = "\t", row.names = 1)
tax = read.table("ASVs_taxonomy.tsv", header = TRUE, row.names = 1)
meta = read.table("sample_info_final.txt", header = TRUE, row.names = 1)
tree = read_tree("asv-famsa-trimal.tre")

mat = as.matrix(mat)
tax = as.matrix(tax)

OTU = otu_table(mat, taxa_are_rows = TRUE)
TAX = tax_table(tax)
META = sample_data(meta)
tail(META)
head(TAX)
head(OTU)

phy = phyloseq(OTU,TAX, META, tree)
phy

# Get to know your phyloseq object
ntaxa(phy) # 1700, this is less than before but makes more sense I think
nsamples(phy) # 26
sample_variables(phy)
taxa_names(phy)
sample_names(phy)
sample_sums(phy)

# Get some statistics on your dataset
num.reads = sum(sample_sums(phy))
lowest.sam = sort(sample_sums(phy)) 
mean.seq = mean(sample_sums(phy))  
std.seq = sd(sample_sums(phy))/sqrt(26) 
med.seq = median(sample_sums(phy))
phy.summary <- data.frame(num.reads, mean.seq, std.seq, med.seq)
phy.summary      

seq.dt = data.table(as(sample_data(phy), "data.frame"),
                    TotalReads = sample_sums(phy), keep.rownames = TRUE)
seq.dt # information per sample after filtering out taxa 

# Visualize dataset properties
readsumsdf = data.frame(nreads = sort(taxa_sums(phy), TRUE), sorted = 1:ntaxa(phy), type = "Nodes")
readsumsdf = rbind(readsumsdf, data.frame(nreads = sort(sample_sums(phy),TRUE),sorted = 1:nsamples(phy), type = "Samples"))
sample_sums(phy)
title = "Total number of reads"
p = ggplot(readsumsdf, aes(x=sorted, y = nreads)) + geom_bar(stat="identity")
p + ggtitle(title) + 
  facet_wrap(~type, 1, scales = "free") + 
  scale_y_log10() +
  pretty.theme()

## FILTERING & TRANSFORMATIONS ===========================================================
# Filtering out taxa 
# Getting rid of unwanted taxonomic groups
# We need to make sure we are actually using the correct file for filtered out abundances
# See "stacked barplot" file to confirm

phy %>%
  subset_taxa(Family != "Mitochondria" &
                Class != "Chloroplast") -> phy.f
phy #1702 taxa prior to filtering out above taxa
phy.f # 1420 taxa remain

# Making a phyloseq object, filtered, without pond sample
phy.f.nopond = subset_samples(phy.f, chemostat_ID !="POND")
phy.f.nopond
sample_names(phy.f.nopond)

# Conduct any transformations you are interested in
per = transform_sample_counts(phy.f, function (x) x/sum(x)*100)
lowabundnames = filter_taxa(per, function(x) sum(x) > 0.1)
lowabundnames = filter_taxa(per, function(x) mean(x) > 0.1) ## This filters out anything that has less than a mean of 0.2% relative abundance acrocss all samples
per.above1 = prune_taxa(lowabundnames, per) 
ntaxa(per.above1) # 90 remain? 

# Remove pond samples from relative abundance
per = transform_sample_counts(phy.f, function (x) x/sum(x)*100)
per_nopond = subset_samples(per, chemostat_ID !="POND")

# Filter out pond and then calculate a new relative abundance? 
per = transform_sample_counts(phy.f.nopond, function (x) x/sum(x)*100)
lowabundnames = filter_taxa(per, function(x) mean(x) > 0.1) ## OTUs with a mean greater than 0.1 are kept.
per.nopond.above1= prune_taxa(lowabundnames, per) 
ntaxa(per.nopond.above1) # 88 remain - seems low

# And then try this other way of filtering out samples
per = transform_sample_counts(phy.f, function (x) x/sum(x)*100) #1445
per_abundf = filter_taxa(per, function(x) sum(x) > 10, TRUE)
ntaxa(per_abundf) 
per_abundf # 34 taxa remain

## REDOING FILTERING TO 5% BASED ON DISCUSSION WITH JOE
per = transform_sample_counts(phy.f, function (x) x/sum(x)*100) #1445
per_abundf = filter_taxa(per, function(x) sum(x) >= 5, TRUE)
ntaxa(per_abundf) 
per_abundf # 59 taxa remain
