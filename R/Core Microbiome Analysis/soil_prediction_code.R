#!/usr/bin/env Rscript
library(phyloseq)
library(ape) # importing trees
library(tidyverse)
library(vegan)
library(microbiome)
library("ggVennDiagram")
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

##### Saving #####
save(mpt_final, file="mpt_final.RData")

## Core microbiome analysis
# Convert to relative abundance
mpt_RA <- transform_sample_counts(mpt_final, fun=function(x) x/sum(x))


# Filter dataset by soil classification
mpt_brun <- subset_samples(mpt_RA, str_detect(Soil.Classification, "Brunisolic Gray Luvisol"))
mpt_orth_hfp <- subset_samples(mpt_RA, str_detect(Soil.Classification, "Orthic Humo-Ferric Podzol"))
mpt_gley_edb <- subset_samples(mpt_RA, str_detect(Soil.Classification, "Gleyed Eluviated Dystric Brunisol"))
mpt_orth_gl <- subset_samples(mpt_RA, str_detect(Soil.Classification, "Orthic Gray Luvisol"))
mpt_gley_gl <- subset_samples(mpt_RA, str_detect(Soil.Classification, "Gleyed Gray Luvisol"))
mpt_orth_db <- subset_samples(mpt_RA, str_detect(Soil.Classification, "Orthic Dystric Brunisol"))
mpt_gley_db <- subset_samples(mpt_RA, str_detect(Soil.Classification, "Gleyed Dystric Brunisol"))
mpt_mes_uh <- subset_samples(mpt_RA, str_detect(Soil.Classification, "Mesic Ultic Haploxeralfs"))
mpt_aq_gl <- subset_samples(mpt_RA, str_detect(Soil.Classification, "Aquic Glossudalfs"))
soil_na <- subset_samples(mpt_RA, is.na(Soil.Classification))

# What ASVs are found in more than 70% of samples in each antibiotic usage category?
# trying changing the prevalence to see what happens
prev_threshold <- 0.05
brun_ASVs <- core_members(mpt_brun, detection=0, prevalence = prev_threshold)
orth_hfp_ASVs <- core_members(mpt_orth_hfp, detection=0, prevalence = prev_threshold)
gley_edb_ASVs <- core_members(mpt_gley_edb, detection=0, prevalence = prev_threshold)
orth_gl_ASVs <- core_members(mpt_orth_gl, detection=0, prevalence = prev_threshold)
gley_gl_ASVs <- core_members(mpt_gley_gl, detection=0, prevalence = prev_threshold)
orth_db_ASVs <- core_members(mpt_orth_db, detection=0, prevalence = prev_threshold)
gley_db_ASVs <- core_members(mpt_gley_db, detection=0, prevalence = prev_threshold)
mes_uh_ASVs <- core_members(mpt_mes_uh, detection=0, prevalence = prev_threshold)
aq_gl_ASVs <- core_members(mpt_aq_gl, detection=0, prevalence = prev_threshold)
soil_na_ASVs <- core_members(soil_na, detection=0, prevalence = prev_threshold)

# can plot those ASVs' relative abundance
prune_taxa(brun_ASVs,mpt_RA) %>% 
  plot_bar(fill="Genus") + 
  facet_wrap(.~`Soil.Classification`, scales ="free")
prune_taxa(orth_hfp_ASVs,mpt_RA) %>% 
  plot_bar(fill="Genus") + 
  facet_wrap(.~`Soil.Classification`, scales ="free")
prune_taxa(gley_edb_ASVs,mpt_RA) %>% 
  plot_bar(fill="Genus") + 
  facet_wrap(.~`Soil.Classification`, scales ="free")
prune_taxa(orth_gl_ASVs,mpt_RA) %>% 
  plot_bar(fill="Genus") + 
  facet_wrap(.~`Soil.Classification`, scales ="free")
prune_taxa(gley_gl_ASVs,mpt_RA) %>% 
  plot_bar(fill="Genus") + 
  facet_wrap(.~`Soil.Classification`, scales ="free")
prune_taxa(orth_db_ASVs,mpt_RA) %>% 
  plot_bar(fill="Genus") + 
  facet_wrap(.~`Soil.Classification`, scales ="free")
prune_taxa(gley_db_ASVs,mpt_RA) %>% 
  plot_bar(fill="Genus") + 
  facet_wrap(.~`Soil.Classification`, scales ="free")
prune_taxa(mes_uh_ASVs,mpt_RA) %>% 
  plot_bar(fill="Genus") + 
  facet_wrap(.~`Soil.Classification`, scales ="free")
prune_taxa(aq_gl_ASVs,mpt_RA) %>% 
  plot_bar(fill="Genus") + 
  facet_wrap(.~`Soil.Classification`, scales ="free")
prune_taxa(soil_na_ASVs,mpt_RA) %>% 
  plot_bar(fill="Genus") + 
  facet_wrap(.~`Soil.Classification`, scales ="free")
