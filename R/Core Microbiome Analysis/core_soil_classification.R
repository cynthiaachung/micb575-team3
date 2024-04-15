#!/usr/bin/env Rscript
library(phyloseq)
library(ape) # importing trees
library(tidyverse)
library(vegan)
library(microbiome)
library("ggVennDiagram")
library(dplyr)
library(patchwork)
library(UpSetR)



# load dephyloseq function
dephyloseq = function(phylo_obj){ 
  meta = as.data.frame(as.matrix(phylo_obj@sam_data)) 
  metacols = ncol(meta)+1
  otu = as.data.frame(t(as.matrix(phylo_obj@otu_table)))
  mo = merge(meta, otu, by=0) 
  tax = as.data.frame(phylo_obj@tax_table) %>% select(-c(Class, asv_id))
  tax = tax %>% rownames_to_column(var="asv_id")
  mo = mo %>% pivot_longer(cols = -c(1:metacols), names_to = "asv_id", values_to="asv_abundance") 
  mot = full_join(mo, tax)
  output = mot 
}


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
otu_mat <- as.matrix(otu[,-1]) # remove index row (remove first row)
# Make first column (#OTU ID) the rownames of the new matrix, ASVs as rows
rownames(otu_mat) <- otu$`#OTU ID` #replace index and replace from column called OTU ID
# Use the "otu_table" function to make an OTU table
OTU <- otu_table(otu_mat, taxa_are_rows = TRUE) #create phyloseq object from otu
class(OTU) #it is a phyloseq object


#### Format sample metadata ####
# Save everything except sampleid as new data frame
samp_df <- as.data.frame(meta[,-1])
# Make sampleids the rownames
rownames(samp_df)<- meta$'#SampleID'
# Make phyloseq sample data with sample_data() function
SAMP <- sample_data(samp_df)
class(SAMP)


#### Formatting taxonomy ####
# Convert taxon strings to a table with separate taxa rank columns
tax_mat <- tax %>% select(-Confidence)%>% #remove confidence column
  separate(col=Taxon, sep="; "
           , into = c("Domain","Phylum","Class","Order","Family","Genus","Species")) %>%
  as.matrix() # convert all this into a matrix
# Save everything except feature IDs 
tax_mat <- tax_mat[,-1]
# Make sampleids the rownames
rownames(tax_mat) <- tax$`Feature ID`
# Make taxa table
TAX <- tax_table(tax_mat)
class(TAX)


#### Create phyloseq object ####
# Merge all into a phyloseq object!
phylobj <- phyloseq(OTU, SAMP, TAX, phylotree)



#### Looking at phyloseq object #####
# View components of phyloseq object with the following commands
# Can see the individual files created
otu_table(phylobj)
sample_data(phylobj)
tax_table(phylobj)
phy_tree(phylobj)

#### Accessor functions ####
# These functions allow you to see or summarise data

# If we look at sample variables and decide we only want some variables, we can view them like so:
sample_variables(phylobj) 
# colnames(sample_data(atamaca))
# knows to look at metadata part of phylobj
get_variable(phylobj, c("Soil.Classification","Compaction.Treatment","Tree.Cover")) # equivalent to "select" in tidyverse

## Let's say we want to filter OTU table by sample. 
# We can first view sample names: tells us all the different samples listed in OTU table
sample_names(phylobj) 
# How many samples do we have?
nsamples(phylobj) 
# What is the sum of reads in each sample?
sample_sums(phylobj) 
# Save the sample names of the 3 samples with the most reads
getsamps <- names(sort(sample_sums(phylobj), decreasing = TRUE)[1:3]) 
# filter to see taxa abundances for each sample
get_taxa(phylobj, getsamps) #recreate taxa table with the ones with the most reads


## Conversely, let's say we want to compare OTU abundance of most abundant OTU across samples
# Look at taxa names
taxa_names(phylobj)
# How many taxa do we have?
ntaxa(phylobj)
# What is the total read count for each taxa?
taxa_sums(phylobj)
# Let's find the top 3 most abundant taxa
gettaxa <- names(sort(taxa_sums(phylobj), decreasing = TRUE)[1:3] )
get_sample(phylobj, gettaxa)


######### ANALYZE ##########
# filtering step, removes from all components
# Remove non-bacterial sequences, if any and exclude chloroplast class and mitochondria family
phylobj_filt <- subset_taxa(phylobj,  Domain == "d__Bacteria" & Class!="c__Chloroplast" & Family !="f__Mitochondria")
#------------------Ran this (optional) check results-------------------------------
# Remove ASVs that have less than 5 counts total #optional
# applies a new function, go through sum of ASVs if greater than 5 keep, otherwise remove
phylobj_filt_nolow <- filter_taxa(phylobj_filt, function(x) sum(x)>5, prune = TRUE)
#---------------------------__________---------------------------------------------
# Remove samples with less than 100 reads
phylobj_filt_nolow_samps <- prune_samples(sample_sums(phylobj_filt_nolow)>100, phylobj_filt_nolow)
#------------------Ran this (optional) check results-------------------------------
# Remove samples where Soil Classification and Compaction level is na
phylobj_final <- subset_samples(phylobj_filt_nolow_samps, !is.na(Soil.Classification) )
phylobj_final <- subset_samples(phylobj_filt_nolow_samps, !is.na(Compaction.Treatment) )
#---------------------------__________---------------------------------------------

# saved as RData files, hard to save as human readable files
save(phylobj_final, file="phylobj_final.RData") 

#### Taxonomy bar plots ####
#### Soil Classification ####
taxa_sums <- taxa_sums(phylobj_final)
taxa_names <- taxa_names(phylobj_final)
tax_table <- as.data.frame(tax_table(phylobj_final))
tax_table_rare <- as.data.frame(tax_table(phylobj_final))
meta <- sample_data(phylobj_final)


# Calculate prevalence of each taxa
taxa_presence <- taxa_sums(phylobj_final)
## summarize at Class level
phylobj_class <- tax_glom(phylobj_final, taxrank = "Class")
# calculate the number of reads in each sample.
phylobj_class@sam_data$read_depth = sample_sums(phylobj_class)

# Filter dataset by soil classification
phylobj_brun <- subset_samples(phylobj_class, str_detect(Soil.Classification, "Brunisolic Gray Luvisol"))
phylobj_orth_hfp <- subset_samples(phylobj_class, str_detect(Soil.Classification, "Orthic Humo-Ferric Podzol"))
phylobj_gley_edb <- subset_samples(phylobj_class, str_detect(Soil.Classification, "Gleyed Eluviated Dystric Brunisol"))
phylobj_orth_gl <- subset_samples(phylobj_class, str_detect(Soil.Classification, "Orthic Gray Luvisol"))
phylobj_gley_gl <- subset_samples(phylobj_class, str_detect(Soil.Classification, "Gleyed Gray Luvisol"))
phylobj_orth_db <- subset_samples(phylobj_class, str_detect(Soil.Classification, "Orthic Dystric Brunisol"))
phylobj_gley_db <- subset_samples(phylobj_class, str_detect(Soil.Classification, "Gleyed Dystric Brunisol"))
phylobj_mes_uh <- subset_samples(phylobj_class, str_detect(Soil.Classification, "Mesic Ultic Haploxeralfs"))
phylobj_aq_gl <- subset_samples(phylobj_class, str_detect(Soil.Classification, "Aquic Glossudalfs"))
phylobj_soil_na <- subset_samples(phylobj_class, is.na(Soil.Classification))

# convert to relative abundance
brun_RA <- transform_sample_counts(phylobj_brun, function(x) x/sum(x)) #relative
orth_hfp_RA <- transform_sample_counts(phylobj_orth_hfp, function(x) x/sum(x)) #relative
gley_edb_RA <- transform_sample_counts(phylobj_gley_edb, function(x) x/sum(x)) #relative
orth_gl_RA <- transform_sample_counts(phylobj_orth_gl, function(x) x/sum(x)) #relative
gley_gl_RA <- transform_sample_counts(phylobj_gley_gl, function(x) x/sum(x)) #relative
orth_db_RA <- transform_sample_counts(phylobj_orth_db, function(x) x/sum(x)) #relative
gley_db_RA <- transform_sample_counts(phylobj_gley_db, function(x) x/sum(x)) #relative
mes_uh_RA <- transform_sample_counts(phylobj_mes_uh, function(x) x/sum(x)) #relative
aq_gl_RA <- transform_sample_counts(phylobj_aq_gl, function(x) x/sum(x)) #relative
soil_na_RA <- transform_sample_counts(phylobj_soil_na, function(x) x/sum(x)) #relative

# put phyloseq object to a matrix and calculate relative abundance
# Convert to core members
det = 0
prev = 0.7
brun_df <- psmelt(brun_RA)
brun_core <- core_members(phylobj_brun, detection = det, prevalence = prev)
brun_df <- brun_df %>% filter(OTU %in% brun_core)
#brun_df$relative_abundance <- brun_df$read_depth/sum(brun_df$read_depth)
orth_hfp_core <- core_members(phylobj_orth_hfp, detection = det, prevalence = prev)
orth_hfp_df <- psmelt(orth_hfp_RA)
orth_hfp_df <- orth_hfp_df %>% filter(OTU %in% orth_hfp_core)
#orth_hfp_df$relative_abundance <- orth_hfp_df$read_depth/sum(orth_hfp_df$read_depth)
gley_edb_core <- core_members(phylobj_gley_edb, detection = det, prevalence = prev)
gley_edb_df <- psmelt(gley_edb_RA)
gley_edb_df <- gley_edb_df %>% filter(OTU %in% gley_edb_core)
#gley_edb_df$relative_abundance <- gley_edb_df$read_depth/sum(gley_edb_df$read_depth)
orth_gl_core <- core_members(phylobj_orth_gl, detection = det, prevalence = prev)
orth_gl_df <- psmelt(orth_gl_RA)
orth_gl_df <- orth_gl_df %>% filter(OTU %in% orth_gl_core)
#orth_gl_df$relative_abundance <- orth_gl_df$read_depth/sum(orth_gl_df$read_depth)
gley_gl_core <- core_members(phylobj_gley_gl, detection = det, prevalence = prev)
gley_gl_df <- psmelt(gley_gl_RA)
gley_gl_df <- gley_gl_df %>% filter(OTU %in% gley_gl_core)
#gley_gl_df$relative_abundance <- gley_gl_df$read_depth/sum(gley_gl_df$read_depth)
orth_db_core <- core_members(phylobj_orth_db, detection = det, prevalence = prev)
orth_db_df <- psmelt(orth_db_RA)
orth_db_df <- orth_db_df %>% filter(OTU %in% orth_db_core)
#orth_db_df$relative_abundance <- orth_db_df$read_depth/sum(orth_db_df$read_depth)
gley_db_core <- core_members(phylobj_gley_db, detection = det, prevalence = prev)
gley_db_df <- psmelt(gley_db_RA)
gley_db_df <- gley_db_df %>% filter(OTU %in% gley_db_core)
#gley_db_df$relative_abundance <- gley_db_df$read_depth/sum(gley_db_df$read_depth)
mes_uh_core <- core_members(phylobj_mes_uh, detection = det, prevalence = prev)
mes_uh_df <- psmelt(mes_uh_RA)
mes_uh_df <- mes_uh_df %>% filter(OTU %in% mes_uh_core)
#mes_uh_df$relative_abundance <- mes_uh_df$read_depth/sum(mes_uh_df$read_depth)
aq_gl_core <- core_members(phylobj_aq_gl, detection = det, prevalence = prev)
aq_gl_df <- psmelt(aq_gl_RA)
aq_gl_df <- aq_gl_df %>% filter(OTU %in% aq_gl_core)
#aq_gl_df$relative_abundance <- aq_gl_df$read_depth/sum(aq_gl_df$read_depth)
soil_na_core <- core_members(phylobj_soil_na, detection = det, prevalence = prev)
soil_na_df <- psmelt(soil_na_RA)
soil_na_df <- soil_na_df %>% filter(OTU %in% soil_na_core)
#soil_na_df$relative_abundance <- soil_na_df$read_depth/sum(soil_na_df$read_depth)

#core microbiome matrix
max_length <- max(length(brun_core), length(orth_hfp_core), 
                  length(gley_edb_core), length(orth_gl_core),
                  length(gley_gl_core), length(orth_db_core),
                  length(gley_db_core), length(mes_uh_core),
                  length(aq_gl_core))
core_mic_list <- data.frame(brun_core, orth_hfp_core, gley_edb_core,
                            orth_gl_core, gley_gl_core, orth_db_core,
                            mes_uh_core, aq_gl_core)

#plot
first_venn <- ggVennDiagram(x = core_mic_list)
first_venn
ggsave("first_venn.png", units = "px", 
       width = 3500, height = 3500)
