#!/usr/bin/env Rscript
library(phyloseq)
library(ape) # importing trees
library(tidyverse)
library(vegan)

#### Load data ####
# Change file paths as necessary
metafp <- "MICB_421_Soil_Metadata.tsv"
meta <- read_delim(metafp, delim="\t")

otufp <- "feature-table.txt"
otu <- read_delim(file = otufp, delim="\t", skip=1)

taxfp <- "taxonomy.tsv"
tax <- read_delim(taxfp, delim="\t")

phylotreefp <- "tree.nwk"
phylotree <- read.tree(phylotreefp)

#### Format OTU table ####
# OTU tables should be a matrix
# with rownames and colnames as OTUs and sampleIDs, respectively
# Note: tibbles do not allow rownames so if you imported with read_delim, change back

# save everything except first column (OTU ID) into a matrix
otu_mat <- as.matrix(otu[,-1])
# Make first column (#OTU ID) the rownames of the new matrix
rownames(otu_mat) <- otu$`#OTU ID`
# Use the "otu_table" function to make an OTU table
OTU <- otu_table(otu_mat, taxa_are_rows = TRUE) 
class(OTU)

#### Format sample metadata ####
# Save everything except sampleid as new data frame
samp_df <- as.data.frame(meta[,-1])
# Make sampleids the rownames
View(meta)
rownames(samp_df)<- meta$'#SampleID'
# Make phyloseq sample data with sample_data() function
SAMP <- sample_data(samp_df)
class(SAMP)

#### Formatting taxonomy ####
# Convert taxon strings to a table with separate taxa rank columns
tax_mat <- tax %>% select(-Confidence)%>%
  separate(col=Taxon, sep="; "
           , into = c("Domain","Phylum","Class","Order","Family","Genus","Species")) %>%
  as.matrix() # Saving as a matrix
# Save everything except feature IDs 
tax_mat <- tax_mat[,-1]
# Make sampleids the rownames
rownames(tax_mat) <- tax$`Feature ID`
# Make taxa table
TAX <- tax_table(tax_mat)
class(TAX)

mpt <- phyloseq(OTU, SAMP, TAX, phylotree)

#### Looking at phyloseq object #####
# View components of phyloseq object with the following commands
otu_table(mpt)
sample_data(mpt)
tax_table(mpt)
phy_tree(mpt)


#### Accessor functions ####
# These functions allow you to see or summarise data

# If we look at sample variables and decide we only want some variables, we can view them like so:
sample_variables(mpt)
# colnames(sample_data(atamaca))
get_variable(mpt, c("Region","Environmental.Source","Soil.Classification","Compaction.Treatment")) # equivalent to "select" in tidyverse

## Let's say we want to filter OTU table by sample. 
# We can first view sample names:
sample_names(mpt)
# How many samples do we have?
nsamples(mpt)
# What is the sum of reads in each sample?
sample_sums(mpt)
# Save the sample names of the 3 samples with the most reads
getsamps <- names(sort(sample_sums(mpt), decreasing = TRUE)[1:3])
# filter to see taxa abundances for each sample
get_taxa(mpt, getsamps) 

## Conversely, let's say we want to compare OTU abundance of most abundant OTU across samples
# Look at taxa names
taxa_names(mpt)
# How many taxa do we have?
ntaxa(mpt)
# What is the total read count for each taxa?
taxa_sums(mpt)
# Let's find the top 3 most abundant taxa
gettaxa <- names(sort(taxa_sums(mpt), decreasing = TRUE)[1:3] )
get_sample(mpt, gettaxa)


######### ANALYZE ##########
# Remove non-bacterial sequences, if any
mpt_filt <- subset_taxa(mpt,  Domain == "d__Bacteria" & Class!="c__Chloroplast" & Family !="f__Mitochondria")
# Remove ASVs that have less than 2 counts total
mpt_filt_nolow <- filter_taxa(mpt_filt, function(x) sum(x)>2, prune = TRUE)
# Remove samples with less than 2 reads
mpt_final <- prune_samples(sample_sums(mpt_filt_nolow)>2, mpt_filt_nolow)
# Remove samples where month is na
#mpt_final <- subset_samples(mpt_filt_nolow_samps, !is.na(month) )

# Rarefy samples
# rngseed sets a random number. If you want to reproduce this exact analysis, you need
# to set rngseed the same number each time
# t transposes the table to use rarecurve function
# cex decreases font size
rarecurve(t(as.data.frame(otu_table(mpt_final))), cex=0.1)
