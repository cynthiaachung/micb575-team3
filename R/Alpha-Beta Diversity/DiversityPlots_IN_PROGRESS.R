# Alpha-Beta Diversity Plots


#!/usr/bin/env Rscript
# NOTE: to install phyloseq, please use the following code instead of the usual "install.packages" function:
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("phyloseq")


library(phyloseq)
library(ape) # importing trees
library(tidyverse)
library(vegan)
library(picante)
library(patchwork)



#### Load data ####
# Change file paths as necessary
metafp <- "qiime_export/metadata.tsv"
meta <- read_delim(metafp, delim="\t") #imports as tibbles

otufp <- "qiime_export/feature-table.txt"
otu <- read_delim(file = otufp, delim="\t", skip=1)

taxfp <- "qiime_export/taxonomy.tsv"
tax <- read_delim(taxfp, delim="\t")

phylotreefp <- "qiime_export/tree.nwk"
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
save(phylobj_final, file="phylobj_final.RData") #before rarefyed

# Rarefy samples
# rngseed sets a random number. If you want to reproduce this exact analysis, you need
# to set rngseed the same number each time
# t transposes the table to use rarecurve function
# cex decreases font size
# rarefaction, 1) rarefaction curve, using final filtered file, cex = font size
rarecurve(t(as.data.frame(otu_table(phylobj_final))), cex=0.1)
# random subsampling so data can look different each time, so we set a random number rngseed = 1 (or any number) will sample every single time
phylobj_rare <- rarefy_even_depth(phylobj_final, rngseed = 1, sample.size = 2526) 
phylobj_rare2 <- rarefy_even_depth(phylobj_final, rngseed = 2, sample.size = 1500)
###################################Keep 2526 if you have enough samples otherwise go lower
# using R's plotting, not ggplot
# after this the console says some samples where removed

##### Saving #####
# saved as RData files, hard to save as human readable files
save(phylobj_rare, file="phylobj_rare.RData")
save(phylobj_rare2, file="phylobj_rare2.RData")  




#### Load in RData ####
load("phylobj_rare.RData")
load("phylobj_final.RData")

#### Alpha diversity ######
plot_richness(phylobj_rare) 
plot_richness(phylobj, measures = c("Shannon","Chao1")) # only interested in Shannon and Chao1

# save into object called gg_richness
# define x axis
# see only shannon and Fisher
# relabeled x axis
# have them be boxplots not dots
gg_richness <- plot_richness(phylobj_rare, x = "Soil.Classification", measures = c("Shannon","Fisher")) +
  xlab("Soil Classification") +
  geom_boxplot()
gg_richness

gg_richness_compac <- plot_richness(phylobj_rare, x = "Compaction.Treatment", measures = c("Shannon","Fisher")) +
  xlab("Compaction.Treatment") +
  geom_boxplot()
gg_richness_compac

# save plot, 
# jpeg or png, png has higher resolution
ggsave(filename = "plot_richness.png"
       , gg_richness
       , height=4, width=6)

ggsave(filename = "plot_richness_compac.png"
       , gg_richness_compac
       , height=4, width=6)

estimate_richness(phylobj_rare)
#####################explain in Github Readme why we control for different variables###!!!!!!!######## paranova analysis, and control for other parameters like site, herbecide...
#####################start here##################################################################################################################

x <- get_variable(phylobj, c("Site","Soil.Classification","LTSP.Treatment","Compaction.Treatment","Tree.Cover")) # equivalent to "select" in tidyverse

obj1 <- subset_samples(phylobj_rare, Ecozone == "IDFBC")
print(sample_data(obj1))
print(otu_table(obj1))
gg_richness_filtered <- plot_richness(obj1, x = "Soil.Classification", measures = c("Shannon","Fisher")) +
  xlab("Soil Classification") +
  geom_boxplot()
gg_richness_filtered

#-----------------------------------------------------------------not getting what I want--------------------------------------------------------------------
# Extract the unique sites from the Site column
all_sites <- unique(sample_data(phylobj_rare)$Site)

# Print the unique sites
print(all_sites)

one_site = subset_samples(phylobj_rare, Herbicide.Use == "0")

gg_richness <- plot_richness(one_site, x = "Soil.Classification", measures = c("Shannon")) +
  xlab("Soil Classification") +
  geom_boxplot() +
  ggtitle(paste("Richness plot for No herbacide"))
gg_richness


# Create an empty list to store the richness plots for each site
richness_plots <- list()

# Loop through each specific site
for (site in all_sites) {
  # Subset the sample data to include only the specific site
  specific_site_data <- subset_samples(phylobj_rare, Site == site)
  # Create the richness plot for the specific site
  gg_richness <- plot_richness(specific_site_data, x = "Soil.Classification", measures = c("Shannon")) +
    xlab("Soil Classification") +
    geom_boxplot() +
    ggtitle(paste("Richness plot for Site:", site))
  # Store the richness plot in the list
  richness_plots[[site]] <- gg_richness
}

# Set the filename to the site name with the desired file extension
filename <- paste("Richness_", gsub(" ", "_", site), ".png", sep = "")

# Save the plot with the filename
#ggsave(filename = filename, plot = gg_richness, height = 4, width = 6)

# Print or save the richness plots for each site
#for (site in all_sites) {
#  print(richness_plots[[site]])
#  ggsave(filename = paste("Richness_", gsub(" ", "_", site), ".png", sep = "")
#         , gg_richness
#         , height=4, width=6)
#}

#-------------------------------------------------------------------------------------------------------------------------------------


# phylogenetic diversity

# calculate Faith's phylogenetic diversity as PD
# measured phylogenetic distances
# t transposes it, that is the way it is read...
phylo_dist <- pd(t(otu_table(phylobj_rare)), phy_tree(phylobj_rare),
                 include.root=F) #F, don't want it to be rooted
# calculate phylogenetic distance using OTU table
?pd

# add PD to metadata table
sample_data(phylobj_rare)$PD <- phylo_dist$PD

# plot any metadata category against the PD
plot_pd <- ggplot(sample_data(phylobj_rare), aes(Soil.Classification, PD, fill = Soil.Classification)) + 
  geom_boxplot() +
  xlab("Soil Classification") +
  ylab("Phylogenetic Diversity") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# view plot
plot_pd

ggsave("plot_pd.png"
       , plot_pd
       , height=8, width=10)

############ PD for compaction
plot_pd_compac <- ggplot(sample_data(phylobj_rare), aes(Compaction.Treatment, PD, fill = Compaction.Treatment)) + 
  geom_boxplot() +
  xlab("Compaction Treatment") +
  ylab("Phylogenetic Diversity") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# view plot
plot_pd_compac

ggsave("plot_pd_compac.png"
       , plot_pd_compac
       , height=4, width=5)


#### Beta diversity #####
bc_dm <- distance(phylobj_rare, method="bray") #look at bc_dm it is like a pyramid, pairwise
# check which methods you can specify
#?distance

#create pcoa coordinate system
pcoa_bc <- ordinate(phylobj_rare, method="PCoA", distance=bc_dm)

#plot the coordination
#color based on body site
plot_ordination(phylobj_rare, pcoa_bc, color = "Soil.Classification", shape="Soil.Classification")

#save to a variable
gg_pcoa_compac <- plot_ordination(phylobj_rare, pcoa_bc, color = "Compaction.Treatment") + #, shape="Compaction.Treatment") +
  labs(pch="Compaction Treatment", col = "Compaction Treatment") + #rename channels
  stat_ellipse()
gg_pcoa_compac

ggsave("plot_pcoa_compac.png"
       , gg_pcoa_compac
       , height=5, width=8)

gg_pcoa <- plot_ordination(phylobj_rare, pcoa_bc, color = "Soil.Classification") +
  labs(pch="Soil Classification", col = "Soil Classification") + #rename channels
  stat_ellipse(level = 0.95)
gg_pcoa

ggsave("plot_pcoa.png"
       , gg_pcoa
       , height=5, width=8)

#-------------------------------------------------------------------------------------------------------------------------------------

#### Taxonomy bar plots ####
#### Soil Classification ####
taxa_sums <- taxa_sums(phylobj_final)
taxa_names <- taxa_names(phylobj_final)
tax_table <- as.data.frame(tax_table(phylobj_final))
tax_table_rare <- as.data.frame(tax_table(phylobj_rare))
meta <- sample_data(phylobj_final)

# Better bar plot of taxonomy Class ###########
#plot_bar(phylobj_rare, fill="Class") 
#?plot_bar

# Convert to relative abundance instead of absolute, to make it more comparable between samples, create new phyloseq object
phylobj_RA <- transform_sample_counts(phylobj_rare, function(x) x/sum(x)) #relative

# To remove black bars, "glom" by phylum first, grouped based on phylum
phylobj_class <- tax_glom(phylobj_RA, taxrank = "Class", NArm=FALSE)

# Make dataframe from phyloseq object using psmelt
df <- psmelt(phylobj_class)

# Calculate average abundances
average_abundance <- df %>%
  group_by(Soil.Classification, Class) %>%
  summarise(AverageAbundance = mean(Abundance)) %>%
  ungroup()

# Plot it
gg_taxa_class <- ggplot(average_abundance, aes(x = Soil.Classification, y = AverageAbundance, fill = Class)) +
  geom_bar(stat = "identity") + #, position = "stack") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "Soil Classification", y = "Average Class Abundance", title = "Average Microbial Class Abundance by Soil Classification")

ggsave("plot_taxaclass.png"
       , gg_taxa_class
       , height=8, width =12)

#-------------------------------------------------------------------------------------------------------------------------------------

# Better bar plot of taxonomy Phylum ###########
#plot_bar(phylobj_rare, fill="Phylum") 
#?plot_bar

# Convert to relative abundance instead of absolute, to make it more comparable between samples, create new phyloseq object
phylobj_RA <- transform_sample_counts(phylobj_rare, function(x) x/sum(x)) #relative

# To remove black bars, "glom" by phylum first, grouped based on phylum
phylobj_phylum <- tax_glom(phylobj_RA, taxrank = "Phylum", NArm=FALSE)

# Make dataframe from phyloseq object using psmelt
df_phylum <- psmelt(phylobj_phylum)

# Calculate average abundances
average_abundance <- df_phylum %>%
  group_by(Soil.Classification, Phylum) %>%
  summarise(AverageAbundance = mean(Abundance)) %>%
  ungroup()

# Plot it
gg_taxa_phylum <- ggplot(average_abundance, aes(x = Soil.Classification, y = AverageAbundance, fill = Phylum)) +
  geom_bar(stat = "identity") + #, position = "stack") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "Soil Classification", y = "Average Phylum Abundance", title = "Average Microbial Phylum Abundance by Soil Classification")
gg_taxa_phylum

ggsave("plot_taxaphylum.png"
       , gg_taxa_phylum
       , height=8, width =12)



#-------------------------------------------------------------------------------------------------------------------------------------

# Better bar plot of taxonomy Domain ###########
#plot_bar(phylobj_rare, fill="Domain") 
#?plot_bar

# Convert to relative abundance instead of absolute, to make it more comparable between samples, create new phyloseq object
phylobj_RA <- transform_sample_counts(phylobj_rare, function(x) x/sum(x)) #relative

# To remove black bars, "glom" by phylum first, grouped based on phylum
phylobj_domain <- tax_glom(phylobj_RA, taxrank = "Domain", NArm=FALSE)

# Make dataframe from phyloseq object using psmelt
df_domain <- psmelt(phylobj_domain)

# Calculate average abundances
average_abundance <- df_domain %>%
  group_by(Soil.Classification, Domain) %>%
  summarise(AverageAbundance = mean(Abundance)) %>%
  ungroup()

# Plot it
gg_taxa_domain <- ggplot(average_abundance, aes(x = Soil.Classification, y = AverageAbundance, fill = Domain)) +
  geom_bar(stat = "identity") + #, position = "stack") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "Soil Classification", y = "Average Domain Abundance", title = "Average Microbial Domain Abundance by Soil Classification")
gg_taxa_domain

ggsave("plot_taxadomain.png"
       , gg_taxa_domain
       , height=8, width =12)









#-------------------------------------------------------------------------------------------------------------------------------------

#### Taxonomy bar plots ####
#### Compaction Treatment ####

# Better bar plot of taxonomy Class ###########
#plot_bar(phylobj_rare, fill="Class") 
#?plot_bar

# Convert to relative abundance instead of absolute, to make it more comparable between samples, create new phyloseq object
phylobj_RA <- transform_sample_counts(phylobj_rare, function(x) x/sum(x)) #relative

# To remove black bars, "glom" by phylum first, grouped based on phylum
phylobj_class <- tax_glom(phylobj_RA, taxrank = "Class", NArm=FALSE)

# Make dataframe from phyloseq object using psmelt
df <- psmelt(phylobj_class)

# Calculate average abundances
average_abundance <- df %>%
  group_by(Compaction.Treatment, Class) %>%
  summarise(AverageAbundance = mean(Abundance)) %>%
  ungroup()

# Plot it
gg_taxa_class_compac <- ggplot(average_abundance, aes(x = Compaction.Treatment, y = AverageAbundance, fill = Class)) +
  geom_bar(stat = "identity") + #, position = "stack") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "Compaction Treatment", y = "Average Class Abundance", title = "Average Microbial Class Abundance by Compaction Treatment")
gg_taxa_class_compac

ggsave("plot_taxaclass_compac.png"
       , gg_taxa_class_compac
       , height=8, width =12)

#-------------------------------------------------------------------------------------------------------------------------------------

# Better bar plot of taxonomy Phylum ###########
#plot_bar(phylobj_rare, fill="Phylum") 
#?plot_bar

# Convert to relative abundance instead of absolute, to make it more comparable between samples, create new phyloseq object
phylobj_RA <- transform_sample_counts(phylobj_rare, function(x) x/sum(x)) #relative

# To remove black bars, "glom" by phylum first, grouped based on phylum
phylobj_phylum <- tax_glom(phylobj_RA, taxrank = "Phylum", NArm=FALSE)

# Make dataframe from phyloseq object using psmelt
df_phylum <- psmelt(phylobj_phylum)

# Calculate average abundances
average_abundance <- df_phylum %>%
  group_by(Compaction.Treatment, Phylum) %>%
  summarise(AverageAbundance = mean(Abundance)) %>%
  ungroup()

# Plot it
gg_taxa_phylum_compac <- ggplot(average_abundance, aes(x = Compaction.Treatment, y = AverageAbundance, fill = Phylum)) +
  geom_bar(stat = "identity") + #, position = "stack") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "Compaction Treatment", y = "Average Phylum Abundance", title = "Average Microbial Phylum Abundance by Compaction Treatment")
gg_taxa_phylum_compac

ggsave("plot_taxaphylum_compac.png"
       , gg_taxa_phylum_compac
       , height=8, width =12)


