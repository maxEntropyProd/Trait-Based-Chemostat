# phyloseq object generation
# Bulseco et al., 'microbial growth strategies' 

# This code assumes you have run "DADA2_Pipeline.R"

# SETUP ENVIRONMENT----
## Import Libraries----
library(phyloseq); library(data.table); library(ggplot2); library(dplyr)
library(data.table)

# PHYLOSEQ ANALYSIS----
# Loading into Phyloseq----
mat = read.table("input_files/ASVs_counts.tsv", header = TRUE, sep = "\t", row.names = 1)
tax = read.table("input_files/ASVs_taxonomy.tsv", header = TRUE, row.names = 1)
meta = read.table("input_files/sample_info_final.txt", header = TRUE, row.names = 1)
tree = read_tree("input_files/asv-famsa-trimal.tre") # comment out if not using tree 

# Set count matrix and tax files as matrices
mat = as.matrix(mat)
tax = as.matrix(tax)

OTU = otu_table(mat, taxa_are_rows = TRUE)
TAX = tax_table(tax)
META = sample_data(meta)
tail(META)
head(TAX)
head(OTU)

# Create phyloseq object by combining count matrix, tax file, and metadata
phy = phyloseq(OTU,TAX, META, tree)
phy

# Explore----
# Get to know your phyloseq object
ntaxa(phy) # 1702, this is less than before but makes more sense I think
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
seq.dt 

## Filter----
# Getting rid of unwanted taxonomic groups

# Remove those that match Mitochondria and chloroplasts 
phy %>%
  subset_taxa(Family != "Mitochondria" &
                Class != "Chloroplast") -> phy.f
phy #1702 taxa prior to filtering out above taxa
phy.f # 1420 taxa remain

# rarefy for when necessary
depth <- min(sample_sums(phy.f))
depth # 18006
phy.f.rare <- rarefy_even_depth(phy.f, sample.size = depth, rngseed = 123, replace = FALSE)
phy.f # 1420 taxa
phy.f.rare # 1217 taxa (lose 203 taxa)
summary(sample_sums(phy.f.rare))# confirming rarefying worked

# Making a phyloseq object, filtered, without pond sample
# this phyloseq object is not rarefied
phy.f.nopond = subset_samples(phy.f, chemostat_ID !="POND")
phy.f.nopond
sample_names(phy.f.nopond) # Should remove pond 

# Filter out pond on the rarefied dataset
phy.f.rare.nopond = subset_samples(phy.f.rare, chemostat_ID !="POND")
phy.f.rare.nopond
sample_names(phy.f.rare.nopond) # Should remove pond 
summary(sample_sums(phy.f.rare.nopond))# confirming rarefying worked

# Conduct any transformations you are interested in
per = transform_sample_counts(phy.f, function (x) x/sum(x)*100)
lowabundnames = filter_taxa(per, function(x) sum(x) > 0.1)
lowabundnames = filter_taxa(per, function(x) mean(x) > 0.1) ## This filters out anything that has less than a mean of 0.2% relative abundance acrocss all samples
per.above1 = prune_taxa(lowabundnames, per) 
ntaxa(per.above1) # 90 remain? 

# Remove pond samples from relative abundance
per = transform_sample_counts(phy.f, function (x) x/sum(x)*100)
per_nopond = subset_samples(per, chemostat_ID !="POND")
per_nopond # should be 25 samples (no pond) 

# Interested in subsetting by chemostat ID
# Reminder: per_nopond is where relative abundance has been calculated 
# And pond samples have been removed
per_MC1 <- subset_samples(per_nopond, chemostat_ID=="MC1")
per_MC2 <- subset_samples(per_nopond, chemostat_ID=="MC2")
per_MC1 # 12 samples
per_MC2 # 13 samples

# Save phyloseq objects----
saveRDS(phy, "phyloseq_objects/phy_noFiltering.rds")
saveRDS(phy.f, "phyloseq_objects/phy_taxaFilter.rds")
saveRDS(per, "phyloseq_objects/phy_taxaFilter_relAbund.rds")
saveRDS(per_nopond, "phyloseq_objects/phy_taxaFilter_relAbund_noPond.rds")
saveRDS(per_MC1, "phyloseq_objects/phy_taxaFilter_relAbund_noPond_MC1.rds")
saveRDS(per_MC2, "phyloseq_objects/phy_taxaFilter_relAbund_noPond_MC2.rds")
saveRDS(phy.f.rare, "phyloseq_objects/phy_taxaFilter_rarefied.rds")
saveRDS(phy.f.rare.nopond, "phyloseq_objects/phy_taxaFilter_rarefied_noPond.rds")


