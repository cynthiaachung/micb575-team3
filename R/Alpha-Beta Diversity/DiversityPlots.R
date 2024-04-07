# Alpha-Beta Diversity Plots


#!/usr/bin/env Rscript
# NOTE: to install phyloseq, please use the following code instead of the usual "install.packages" function:
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("phyloseq")




#------------------------------------------Run section every time------------------------------------------------
library(phyloseq)
library(ape) # importing trees
library(tidyverse)
library(vegan)
library(picante)
library(patchwork)
library(viridis) # Load the viridis package for color pallets
library(gridExtra)



#-----------------------------Run below to create phyloseq object, or ignore if you already have it in environment-----------------------------
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






#------------------------------------------Ignore section below, not relevant to final plots------------------------------------------------
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







#------------------------------------------Run below for phylseq object creation------------------------------------------------
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
#Kept rarefaction parameter of 2526 if you have enough samples otherwise go lower, but gave enough samples so good enough
# using R's plotting, not ggplot
# after this the console says some samples where removed

##### Saving #####
# saved as RData files, hard to save as human readable files
save(phylobj_rare, file="phylobj_rare.RData")
save(phylobj_rare2, file="phylobj_rare2.RData")  





#----------------Run below if you already created phyloseq object, and now want to load it, ignore if already in environment---------
#### Load in RData ####
load("phylobj_rare.RData")
load("phylobj_final.RData")
#-----------------------------------------------------------------------------------------------------------------------------------------



###################################################  Alpha diversity Plots ##############################################




#-----------------------------------------------------------------Figure 1--------------------------------------------------------------------
# Figure 1: Alpha Diversity Richness - Controlled for LTSP Treatment
# Extract the unique LTSP.Treatment from the LTSP.Treatment column
all_ltsp <- unique(sample_data(phylobj_rare)$LTSP.Treatment)

# Print the unique sites
print(all_ltsp)
# Result: "OM2" "OM1" "REF" "OM3"

# filter out NAs from soil classification (Worked!)
phylobj_rare <- subset_samples(phylobj_rare, !is.na(Soil.Classification) & Soil.Classification != "NA")

# Reorder LTSP.Treatment factor levels
sample_data(phylobj_rare)$LTSP.Treatment <- factor(sample_data(phylobj_rare)$LTSP.Treatment,
                                                   levels = c("REF", "OM1", "OM2", "OM3"))

# *******Change LTSP_Treatment == to REF, OM1, OM2, OM3 to generate desired plots********
#LTSP_Treatment = "REF"
#filtered = subset_samples(phylobj_rare, LTSP.Treatment == LTSP_Treatment)

# Run it! THis is the good one gg_richness_horiz
gg_richness_horiz <- plot_richness(phylobj_rare, x = "Soil.Classification", measures = c("Shannon")) +
  geom_boxplot(aes(fill = Soil.Classification)) +
  #uncomment below to do color blind palette
  scale_fill_viridis(discrete = TRUE, option = "D") + # Use a color-blind friendly palette
  facet_grid(. ~ LTSP.Treatment, scales = "free_x", space = "free_x") +
  xlab("Soil Classification") +
  ggtitle("Alpha Diversity Richness for Soil Classification Controlled by LTSP Treatment") +
  theme(
    strip.background = element_blank(),
    strip.placement = "outside",
    strip.text = element_text(size = 14), # Make facet labels larger
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.title.x.bottom = element_text(size = 14),
    axis.title.y = element_text(size = 14), # Adjust y-axis title size
    axis.text.y = element_text(size = 14, angle = 0, hjust = 0.5), # Make y-axis labels horizontal and adjust alignment
    legend.text = element_text(size = 14), # Make legend text larger
    legend.title = element_text(size = 14), # Make legend title larger
    panel.spacing = unit(2, "lines"), # Increase panel spacing for more distinction
    panel.border = element_rect(colour = "black", fill = NA, size = 1), # Add/enhance panel borders
    plot.title = element_text(size = 20, hjust = 0.5) # Make plot title larger and centered
  )
gg_richness_horiz

ggsave("Figure_1_AlphaDiversityRichness_horizontal_panels_colorblind.png"
       , gg_richness_horiz
       , height=9, width=14)




#----------------------------------------------Figure 1 Kruskal Wallis analysis-------------------------------------------------------
# First Need to extract information
alphadiv <- estimate_richness(phylobj_rare)
samp_dat <- sample_data(phylobj_rare)
samp_dat_wdiv <- data.frame(samp_dat, alphadiv)
#View(samp_dat_wdiv)

# Treatment levels to iterate over
treatments <- c("REF", "OM1", "OM2", "OM3")

# Loop through each treatment
for(LTSP_Treatment in treatments) {
  # Filter the data for the current treatment
  filtered_div_data <- subset(samp_dat_wdiv, LTSP.Treatment == LTSP_Treatment)
  
  # Perform the Kruskal-Wallis test for Shannon diversity
  kruskal_obs <- kruskal.test(Shannon ~ `Soil.Classification`, data = filtered_div_data)
  
  # Extract the p-value
  p_value <- kruskal_obs$p.value
  
  # Check if the result is significant and print the result
  if(p_value < 0.05) {
    cat("For", LTSP_Treatment, ": The difference in Shannon diversity across soil classifications is statistically significant. P-value:", p_value, "\n")
  } else {
    cat("For", LTSP_Treatment, ": The difference in Shannon diversity across soil classifications is not statistically significant. P-value:", p_value, "\n")
  }
}





#-----------------------------------Supplemental 1 below, run entire section--------------------------------------------
# phylogenetic diversity, PD
# Figure 1 Supplemental, Faith's phylogenetic diversity
gg_faith_pd <- ggplot(sample_data(phylobj_rare), aes(x = Soil.Classification, y = PD, fill = Soil.Classification)) +
  geom_boxplot() + # Create boxplot
  facet_grid(. ~ LTSP.Treatment, scales = "free_x", space = "free_x") + # Faceting by treatment
  xlab("Soil Classification") + # Customize X-axis label
  ylab("Faith's Phylogenetic Diversity") + # Customize Y-axis label
  ggtitle("Faith's Phylogenetic Diversity across Soil Classification and LTSP Treatment") + # Add title
  theme(
    strip.background = element_blank(),
    strip.placement = "outside",
    strip.text = element_text(size = 14),
    axis.title.x.bottom = element_text(size = 14),
    axis.title.x = element_blank(), # Hide x-axis title
    axis.text.x = element_blank(), # Hide x-axis labels
    axis.title.y = element_text(size = 14),
    axis.text.y = element_text(size = 14, angle = 0, hjust = 0.5),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 14),
    panel.spacing = unit(2, "lines"),
    panel.border = element_rect(colour = "black", fill = NA, size = 1),
    plot.title = element_text(size = 20, hjust = 0.5)
  )
gg_faith_pd 

ggsave("Supplemental_1.png"
       , gg_faith_pd 
       , height=10, width=15)

gg_faith_pd_colorblind <- gg_faith_pd + scale_fill_viridis(discrete = TRUE, option = "D")

ggsave("Supplemental_1_colorblind.png"
       , gg_faith_pd_colorblind 
       , height=10, width=15)

# KW stats for Faith's PD
alphadiv <- estimate_richness(phylobj_rare)
samp_dat <- sample_data(phylobj_rare)
samp_dat_wdiv <- data.frame(samp_dat, alphadiv)
#View(samp_dat_wdiv)

# Treatment levels to iterate over
treatments <- c("REF", "OM1", "OM2", "OM3")

# Loop through each treatment
for(LTSP_Treatment in treatments) {
  # Filter the data for the current treatment
  filtered_div_data <- subset(samp_dat_wdiv, LTSP.Treatment == LTSP_Treatment)
  
  # Perform the Kruskal-Wallis test
  kruskal_obs <- kruskal.test(PD ~ `Soil.Classification`, data = filtered_div_data)
  
  # Extract the p-value
  p_value <- kruskal_obs$p.value
  
  # Check if the result is significant and print the result
  if(p_value < 0.05) {
    cat("For ", LTSP_Treatment, ": The difference in Faith's Phylogenetic Diversity across soil classifications is statistically significant. P-value:", p_value, "\n")
  } else {
    cat("For ", LTSP_Treatment, ": The difference in Faith's Phylogenetic Diversity across soil classifications is not statistically significant. P-value:", p_value, "\n")
  }
}




######################################################  Beta diversity Plots ####################################################

#-----------------------------------------------------Figure 2 Beta diversity-----------------------------------------------------
#### Beta diversity #####
#bc_dm <- distance(phylobj_rare, method="bray") #look at bc_dm it is like a pyramid, pairwise
# check which methods you can specify
#?distance

# filter out NAs from soil classification (not working actually, just removed them from plot)
phylobj_rare <- subset_samples(phylobj_rare, !is.na(Soil.Classification) & Soil.Classification != 'NA')

#create pcoa coordinate system
pcoa_bc <- ordinate(phylobj_rare, method="PCoA", distance=bc_dm)

# soil plots
gg_pcoa <- plot_ordination(phylobj_rare, pcoa_bc, color = "Soil.Classification") +
  labs(pch="Soil Classification", col = "Soil Classification") + 
  stat_ellipse(level = 0.95) +
  scale_color_discrete(na.translate = FALSE) + # This will exclude NA values from the legend and the plot
  theme(text = element_text(size = 18), # Changes global text size
        axis.text = element_text(size = 18), # Changes axis text size
        legend.text = element_text(size = 14)) + # Changes legend text size
  ggtitle("Beta Diversity of Soil Classifications - Bray Curtis PCOA Plot") 
gg_pcoa

ggsave("Figure_2_PCOA_soil.classification.png"
       , gg_pcoa
       , height=5, width=10)



#--------------------------------Figure 2 Beta diversity group soil types then stats-----------------------------------------------------

##combine soil categories

phylobj_rare_soilgroups <- phylobj_rare

# Access sample data
sample_data_df <- data.frame(sample_data(phylobj_rare_soilgroups))

sample_data_df <- sample_data_df %>%
  filter(!is.na(Soil.Classification)) # This removes rows with NA in Soil.Classification

# Change soil classifications names to group them
sample_data_df$Soil.Classification[sample_data_df$Soil.Classification %in% c("Aquic Glossudalfs")] <- "Soil Group 1"
sample_data_df$Soil.Classification[sample_data_df$Soil.Classification %in% c("Brunisolic Gray Luvisol", "Orthic Gray Luvisol, Gleyed Gray Luvisol")] <- "Soil Group 2"
sample_data_df$Soil.Classification[sample_data_df$Soil.Classification %in% c("Gleyed Dystric Brunisol","Mesic Ultic Haploxeralfs", "Orthic Dystric Brunisol", "Orthic Humo-Ferric Podzol")] <- "Soil Group 3"

# Updating phyloseq object with new sample data
sample_data(phylobj_rare_soilgroups) <- sample_data(sample_data_df)

pcoa_bc <- ordinate(phylobj_rare_soilgroups, method="PCoA", distance=bc_dm)

# plot and save
gg_pcoa <- plot_ordination(phylobj_rare_soilgroups, pcoa_bc, color = "Soil.Classification") +
  labs(pch="Soil Classification", col = "Soil Classification") + 
  stat_ellipse(level = 0.95) +
  scale_color_discrete(na.translate = FALSE) + # This will exclude NA values from the legend and the plot
  theme(text = element_text(size = 18), # Changes global text size
        axis.text = element_text(size = 18), # Changes axis text size
        legend.text = element_text(size = 18)) + # Changes legend text size
  ggtitle("Soil Classifications Grouped") +
  theme(plot.title = element_text(hjust = 0.5)) # Centers the title
gg_pcoa

ggsave("Figure_2.2_SoilGrouped.png"
       , gg_pcoa
       , height=5, width=8)

# colorblind plot and save
gg_pcoa_colorblind <- plot_ordination(phylobj_rare_soilgroups, pcoa_bc, color = "Soil.Classification") +
  labs(pch="Soil Classification", col = "Soil Classification") + 
  stat_ellipse(level = 0.95) +
  #scale_color_discrete(na.translate = FALSE) + # This will exclude NA values from the legend and the plot
  scale_color_viridis(discrete = TRUE, option = "D", na.value = NA, begin = 0, end = 1, direction = 1, na.translate = FALSE) + # Using viridis color scale
  theme(text = element_text(size = 18), # Changes global text size
        axis.text = element_text(size = 18), # Changes axis text size
        legend.text = element_text(size = 18)) + # Changes legend text size
  ggtitle("Soil Classifications Grouped") +
  theme(plot.title = element_text(hjust = 0.5)) # Centers the title
gg_pcoa_colorblind

ggsave("Figure_2.2_SoilGrouped_colorblind.png"
       , gg_pcoa_colorblind
       , height=5, width=8)


### PERMANOVA (Permutational ANOVA) ####
# non-parametric version of ANOVA
# Takes a distance matrix, which can be calculated with any kind of metric you want
# e.g. Bray, Jaccard, Unifrac
# Need the package, "vegan"
samp_dat_wdiv <- data.frame(sample_data(phylobj_rare_soilgroups), estimate_richness(phylobj_rare_soilgroups))

dm_bray <- vegdist(t(otu_table(phylobj_rare_soilgroups)), method="bray")
permanova_result <- adonis2(dm_bray ~ Soil.Classification, data=samp_dat_wdiv)
permanova_result 





#-------------------------------------------------------------------Supplemental 2, run this-----------------------------------------------------------------

# Better bar plot of taxonomy Genus ###########
#plot_bar(phylobj_rare, fill="Genus") 
#?plot_bar

# filter out NAs from soil classification 
phylobj_rare <- subset_samples(phylobj_rare, !is.na(Soil.Classification) & Soil.Classification != "NA")

# Convert to relative abundance instead of absolute, to make it more comparable between samples, create new phyloseq object
phylobj_RA <- transform_sample_counts(phylobj_rare, function(x) x/sum(x)) #relative

# To remove black bars, "glom" by phylum first, grouped based on phylum
phylobj_genus <- tax_glom(phylobj_RA, taxrank = "Genus", NArm=FALSE)

# Make dataframe from phyloseq object using psmelt
df <- psmelt(phylobj_genus)

# Calculate average abundances
average_abundance <- df %>%
  group_by(Soil.Classification, Genus) %>%
  summarise(AverageAbundance = mean(Abundance)) %>%
  ungroup()

# Better plot
genus_levels <- na.omit(unique(average_abundance$Genus))
gg_taxa_genus <- ggplot(average_abundance, aes(x = Soil.Classification, y = AverageAbundance, fill = Genus)) +
  geom_bar(stat = "identity") +
  scale_fill_viridis(discrete = TRUE, na.value = "transparent", breaks = genus_levels) +  # Explicitly specify breaks to exclude NA
  labs(x = "Soil Classification", y = "Average Genus Abundance", title = "Taxonomy Bar Plot for Soil Classifications - Genus") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
    axis.text.y = element_text(size = 14),
    axis.title = element_text(size = 16),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 16),
    plot.title = element_text(size = 20, hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(size = 18),
    strip.text.x = element_text(size = 16)
  )

gg_taxa_genus

ggsave("Supplemental_2_TaxaBarPlot_colorblind.png"
       , gg_taxa_genus
       , height=8, width =12)


# Older plot
# Plot it
gg_taxa_genus <- ggplot(average_abundance, aes(x = Soil.Classification, y = AverageAbundance, fill = Genus)) +
  geom_bar(stat = "identity") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "Soil Classification", y = "Average Genus Abundance", title = "Average Microbial Genus Abundance by Soil Classification")
gg_taxa_genus






#--------------------------------------------Figure 3, run it!----------------------------------------
# try again, use figure 4 as reference, same type of plot
# figure 3 equivalent ---> figure 4
# Soil classification ---> LTSP treatment
# Plot compaction treatment ----> compaction treatment 

## Combine soil categories

phylobj_rare_soilgroups <- phylobj_rare

# Access sample data
sample_data_df <- data.frame(sample_data(phylobj_rare_soilgroups))

sample_data_df <- sample_data_df %>%
  filter(!is.na(Soil.Classification)) # This removes rows with NA in Soil.Classification

# Change soil classifications names to group them
sample_data_df$Soil.Classification[sample_data_df$Soil.Classification %in% c("Aquic Glossudalfs")] <- "Soil Group 1"
sample_data_df$Soil.Classification[sample_data_df$Soil.Classification %in% c("Brunisolic Gray Luvisol", "Orthic Gray Luvisol, Gleyed Gray Luvisol")] <- "Soil Group 2"
sample_data_df$Soil.Classification[sample_data_df$Soil.Classification %in% c("Gleyed Dystric Brunisol","Mesic Ultic Haploxeralfs", "Orthic Dystric Brunisol", "Orthic Humo-Ferric Podzol")] <- "Soil Group 3"

# Updating phyloseq object with new sample data
sample_data(phylobj_rare_soilgroups) <- sample_data(sample_data_df)

gg_richness_compac_bysoil <- plot_richness(phylobj_rare_soilgroups, x = "Compaction.Treatment", measures = c("Shannon")) +
  geom_boxplot(aes(fill = Compaction.Treatment)) +
  facet_grid(. ~ Soil.Classification, scales = "free_x", space = "free_x") +
  xlab("Compaction Treatment") +
  ggtitle("Alpha Diversity Richness for Compaction Treatment - Controlled by Soil Group") +
  theme(
    strip.background = element_blank(),
    strip.placement = "outside",
    strip.text = element_text(size = 18), # Make facet labels larger
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.title.x.bottom = element_text(size = 18),
    axis.title.y = element_text(size = 18), # Adjust y-axis title size
    axis.text.y = element_text(size = 18, angle = 0, hjust = 0.5), # Make y-axis labels horizontal and adjust alignment
    legend.text = element_text(size = 18), # Make legend text larger
    legend.title = element_text(size = 18), # Make legend title larger
    panel.spacing = unit(2, "lines"), # Increase panel spacing for more distinction
    panel.border = element_rect(colour = "black", fill = NA, size = 1), # Add/enhance panel borders
    plot.title = element_text(size = 20, hjust = 0.5) # Make plot title larger and centered
  )
gg_richness_compac_bysoil

ggsave("Figure_3_AlphaDiversityRichness_compac_vs_soilgroup.png"
       , gg_richness_compac_bysoil
       , height=8, width=15)

gg_richness_compac_bysoil_colorblind <- gg_richness_compac_bysoil + scale_fill_viridis(discrete = TRUE, option = "D")
gg_richness_compac_bysoil_colorblind

ggsave("Figure_3_AlphaDiversityRichness_compac_vs_soilgroup_colorblind.png"
       , gg_richness_compac_bysoil_colorblind
       , height=8, width=15)


#### Figure 3, stats

# First Need to extract information
alphadiv <- estimate_richness(phylobj_rare_soilgroups)
samp_dat <- sample_data(phylobj_rare_soilgroups)
samp_dat_wdiv <- data.frame(samp_dat, alphadiv)
#View(samp_dat_wdiv)

# Soil groups to iterate over, all three groups
soilgroups <- c("Soil Group 1", "Soil Group 2", "Soil Group 3")

# Loop through each treatment
for(soilgroup in soilgroups) {
  # Filter the data for the current treatment
  filtered_div_data <- subset(samp_dat_wdiv, Soil.Classification == soilgroup)
  
  # Perform the Kruskal-Wallis test for Shannon diversity
  kruskal_obs <- kruskal.test(Shannon ~ `Compaction.Treatment`, data = filtered_div_data)
  
  # Extract the p-value
  p_value <- kruskal_obs$p.value
  
  # Check if the result is significant and print the result
  if(p_value < 0.05) {
    cat("For", soilgroup, ": The difference in Shannon diversity across compaction treatment is statistically significant. P-value:", p_value, "\n")
  } else {
    cat("For", soilgroup, ": The difference in Shannon diversity across compaction treatment is not statistically significant. P-value:", p_value, "\n")
  }
}





#------------------------------------------------Supplemental 3 - PD for compaction--------------------------------------------
plot_pd_compac <- ggplot(sample_data(phylobj_rare), aes(Compaction.Treatment, PD, fill = Compaction.Treatment)) + 
  geom_boxplot() +
  xlab("Compaction Treatment") +
  ylab("Phylogenetic Diversity") +
  ggtitle("Faith's Phylogenetic Diversity across Soil Classification and LTSP Treatment") + # Add title
  theme(
    axis.title.x.bottom = element_text(size = 18),
    axis.title.x = element_blank(), # Hide x-axis title
    axis.text.x = element_blank(), # Hide x-axis labels
    axis.title.y = element_text(size = 18),
    axis.text.y = element_text(size = 18, angle = 0, hjust = 0.5),
    legend.text = element_text(size = 18),
    legend.title = element_text(size = 18),
    panel.spacing = unit(2, "lines"),
    panel.border = element_rect(colour = "black", fill = NA, size = 1),
    plot.title = element_text(size = 20, hjust = 0.5)
  )
plot_pd_compac

# save it
ggsave("Supplemental_3.png"
       , plot_pd_compac
       , height=8, width=13)

# Colorblind plot
plot_pd_compac_colorblind <- plot_pd_compac + scale_fill_viridis(discrete = TRUE, option = "D")
plot_pd_compac_colorblind

# save it
ggsave("Supplemental_3_colorblind.png"
       , plot_pd_compac_colorblind
       , height=8, width=13)


# KW stats for Faith's PD
alphadiv <- estimate_richness(phylobj_rare)
samp_dat <- sample_data(phylobj_rare)
samp_dat_wdiv <- data.frame(samp_dat, alphadiv)
#View(samp_dat_wdiv)

# Perform the Kruskal-Wallis test
kruskal_obs <- kruskal.test(PD ~ `Compaction.Treatment`, data = samp_dat_wdiv)
  
# Extract the p-value
p_value <- kruskal_obs$p.value
  
# Check if the result is significant and print the result
if(p_value < 0.05) {
  cat("The difference in Faith's Phylogenetic Diversity across Compaction Treatment is statistically significant. P-value:", p_value, "\n")
  } else {
  cat("The difference in Faith's Phylogenetic Diversity across Compaction Treatment is not statistically significant. P-value:", p_value, "\n")
}






#-----------------------------------------------Figure 4 Alpha diversity compaction----------------------------------------
gg_richness_compac_horiz <- plot_richness(phylobj_rare, x = "Compaction.Treatment", measures = c("Shannon")) +
  geom_boxplot(aes(fill = Compaction.Treatment)) +
  facet_grid(. ~ LTSP.Treatment, scales = "free_x", space = "free_x") +
  xlab("Compaction Treatment") +
  ggtitle("Alpha Diversity Richness for Compaction Treatment Controlled by Organic Matter Removal Level") +
  theme(
    strip.background = element_blank(),
    strip.placement = "outside",
    strip.text = element_text(size = 18), # Make facet labels larger
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.title.x.bottom = element_text(size = 18),
    axis.title.y = element_text(size = 18), # Adjust y-axis title size
    axis.text.y = element_text(size = 18, angle = 0, hjust = 0.5), # Make y-axis labels horizontal and adjust alignment
    legend.text = element_text(size = 18), # Make legend text larger
    legend.title = element_text(size = 18), # Make legend title larger
    panel.spacing = unit(2, "lines"), # Increase panel spacing for more distinction
    panel.border = element_rect(colour = "black", fill = NA, size = 1), # Add/enhance panel borders
    plot.title = element_text(size = 20, hjust = 0.5) # Make plot title larger and centered
  )
gg_richness_compac_horiz

ggsave("Figure_4_AlphaDiversityRichness_compac_vs_LTSP.png"
       , gg_richness_compac_horiz
       , height=8, width=15)

gg_richness_compac_horiz_colorblind <- gg_richness_compac_horiz + scale_fill_viridis(discrete = TRUE, option = "D")
gg_richness_compac_horiz_colorblind

ggsave("Figure_4_AlphaDiversityRichness_compac_vs_LTSP_colorblind.png"
       , gg_richness_compac_horiz_colorblind
       , height=8, width=15)

#-----------------------------------------------Figure 4 stats----------------------------------------

# First Need to extract information
alphadiv <- estimate_richness(phylobj_rare)
samp_dat <- sample_data(phylobj_rare)
samp_dat_wdiv <- data.frame(samp_dat, alphadiv)
#View(samp_dat_wdiv)

# Treatment levels to iterate over, only OM1 and OM2
treatments <- c("OM1", "OM2")

# Loop through each treatment
for(LTSP_Treatment in treatments) {
  # Filter the data for the current treatment
  filtered_div_data <- subset(samp_dat_wdiv, LTSP.Treatment == LTSP_Treatment)
  
  # Perform the Kruskal-Wallis test for Shannon diversity
  kruskal_obs <- kruskal.test(Shannon ~ `Compaction.Treatment`, data = filtered_div_data)
  
  # Extract the p-value
  p_value <- kruskal_obs$p.value
  
  # Check if the result is significant and print the result
  if(p_value < 0.05) {
    cat("For", LTSP_Treatment, ": The difference in Shannon diversity across compaction treatment is statistically significant. P-value:", p_value, "\n")
  } else {
    cat("For", LTSP_Treatment, ": The difference in Shannon diversity across compaction treatment is not statistically significant. P-value:", p_value, "\n")
  }
}











#-----------------------------------------------------------------Taxonomy plots----------------------------------------

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

#----------------------------------------------------------------------don't need---------------------------------------------------------------

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



#------------------------------------------------------------------don't need-------------------------------------------------------------------

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




#--------------------------------------------------------------------don't need-----------------------------------------------------------------

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

#----------------------------------------------------------------don't need---------------------------------------------------------------------

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

#-------------------------------------------------------don't need------------------------------------------------------------------------------


# Better bar plot of taxonomy Genus ###########
#plot_bar(phylobj_rare, fill="Genus") 
#?plot_bar

# Convert to relative abundance instead of absolute, to make it more comparable between samples, create new phyloseq object
phylobj_RA <- transform_sample_counts(phylobj_rare, function(x) x/sum(x)) #relative

# To remove black bars, "glom" by phylum first, grouped based on phylum
phylobj_genus <- tax_glom(phylobj_RA, taxrank = "Genus", NArm=FALSE)

# Make dataframe from phyloseq object using psmelt
df <- psmelt(phylobj_genus)

# Calculate average abundances
average_abundance <- df %>%
  group_by(Compaction.Treatment, Genus) %>%
  summarise(AverageAbundance = mean(Abundance)) %>%
  ungroup()

# Plot it
gg_taxa_genus_compac <- ggplot(average_abundance, aes(x = Compaction.Treatment, y = AverageAbundance, fill = Genus)) +
  geom_bar(stat = "identity") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "Compaction Treatment", y = "Average Genus Abundance", title = "Average Microbial Genus Abundance by Compaction Treatment")
gg_taxa_genus_compac

ggsave("plot_taxa_genus_compac.png"
       , gg_taxa_genus_compac
       , height=8, width =12)




#------------------------------------------Don't want plots below, data is Not controlled for anything------------------------------------------------
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
#########explain in Github Readme why we control for different variables###!!!!!!!######## paranova analysis, and control for other parameters like site, herbecide...

x <- get_variable(phylobj, c("Site","Soil.Classification","LTSP.Treatment","Compaction.Treatment","Tree.Cover")) # equivalent to "select" in tidyverse

obj1 <- subset_samples(phylobj_rare, Ecozone == "IDFBC")
print(sample_data(obj1))
print(otu_table(obj1))
gg_richness_filtered <- plot_richness(obj1, x = "Soil.Classification", measures = c("Shannon","Fisher")) +
  xlab("Soil Classification") +
  geom_boxplot()
gg_richness_filtered



#-------------------------------------------------------Ignore section below, not used------------------------------------------------------------------------------
# phylogenetic diversity, PD

# calculate Faith's phylogenetic diversity as PD
# measured phylogenetic distances
# t transposes it, that is the way it is read...
phylo_dist <- pd(t(otu_table(phylobj_rare)), phy_tree(phylobj_rare),
                 include.root=F) #F, don't want it to be rooted
# calculate phylogenetic distance using OTU table
#?pd

# add PD to metadata table
sample_data(phylobj_rare)$PD <- phylo_dist$PD

# plot any metadata category against the PD
plot_pd <- ggplot(sample_data(phylobj_rare), aes(Soil.Classification, PD, fill = Soil.Classification)) + 
  geom_boxplot() +
  xlab("Soil Classification") +
  ylab("Phylogenetic Diversity") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
plot_pd

ggsave("plot_pd.png"
       , plot_pd
       , height=8, width=10)






#-----------------------------------------------------------------don't use, old way to plot figure 3, looked bad----------------------------------------
# First subset data by soil groups
#sample_data_df$Soil.Classification[sample_data_df$Soil.Classification %in% c("Aquic Glossudalfs")] <- "Soil_Group1"
#sample_data_df$Soil.Classification[sample_data_df$Soil.Classification %in% c("Brunisolic Gray Luvisol", "Orthic Gray Luvisol, Gleyed Gray Luvisol")] <- "Soil_Group2"
#sample_data_df$Soil.Classification[sample_data_df$Soil.Classification %in% c("Gleyed Dystric Brunisol","Mesic Ultic Haploxeralfs", "Orthic Dystric Brunisol", "Orthic Humo-Ferric Podzol")] <- "Soil_Group3"

soilgroup_1 <- subset_samples(phylobj_rare, Soil.Classification %in% c("Aquic Glossudalfs"))
soilgroup_2 <- subset_samples(phylobj_rare, Soil.Classification %in% c("Brunisolic Gray Luvisol", "Orthic Gray Luvisol, Gleyed Gray Luvisol"))
soilgroup_3 <- subset_samples(phylobj_rare, Soil.Classification %in% c("Gleyed Dystric Brunisol","Mesic Ultic Haploxeralfs", "Orthic Dystric Brunisol", "Orthic Humo-Ferric Podzol"))

gg_richness_compac_soil1 <- plot_richness(soilgroup_1, x = "Compaction.Treatment", measures = c("Shannon")) +
  xlab("Compaction Treatment") +
  geom_boxplot(aes(fill = Compaction.Treatment)) +
  ggtitle(paste("Soil Group 1"))
gg_richness_compac_soil1

gg_richness_compac_soil2 <- plot_richness(soilgroup_2, x = "Compaction.Treatment", measures = c("Shannon")) +
  xlab("Compaction Treatment") +
  geom_boxplot(aes(fill = Compaction.Treatment)) +
  ggtitle(paste("Soil Group 2"))
gg_richness_compac_soil2

gg_richness_compac_soil3 <- plot_richness(soilgroup_3, x = "Compaction.Treatment", measures = c("Shannon")) +
  xlab("Compaction Treatment") +
  geom_boxplot(aes(fill = Compaction.Treatment)) +
  ggtitle(paste("Soil Group 3"))
gg_richness_compac_soil3

library(grid)
# Arrange the plots side by side
combined_plot <- grid.arrange(
  gg_richness_compac_soil1, 
  gg_richness_compac_soil2, 
  gg_richness_compac_soil3, 
  ncol = 3,
  top = textGrob("Richness Plot for Compaction Treatment by Soil Groups", gp = gpar(fontface = "bold", fontsize = 20))
)
combined_plot

ggsave("Figure3_richness_compac_soilgroups.png"
       , combined_plot
       , height=10, width=20)






#-----------------------------------------------------------don't need these--------------------------------------------------
# other types of plots for soil classification
# don't run, not used
gg_richness <- plot_richness(phylobj_rare, x = "Soil.Classification", measures = c("Shannon")) +
  xlab("Soil Classification") +
  geom_boxplot(aes(fill = Soil.Classification)) +
  facet_wrap(~ LTSP.Treatment, scales = "free_x") + # Facet by LTSP_Treatment
  ggtitle(paste("Alpha Diversity Richness for Soil Classification Contolled by LTSP Treatment"))
gg_richness

# don't run, not used
gg_richness_stacked <- plot_richness(phylobj_rare, x = "Soil.Classification", measures = c("Shannon")) +
  geom_boxplot(aes(fill = Soil.Classification)) +
  facet_grid(LTSP.Treatment ~ ., scales = "free_x", space = "free_x") + # Facet in a single column
  theme(axis.text.x = element_text(angle = 45, hjust = 1), # Adjust text angle if needed for readability
        strip.background = element_blank(), # Optional: removes the background of facet labels for a cleaner look
        strip.placement = "outside") + # Place facet labels above each panel
  xlab("Soil Classification") +
  ggtitle("Alpha Diversity Richness for Soil Classification Controlled by LTSP Treatment") +
  theme(
    # Hide x-axis titles and text for all but the bottom panel
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    # To ensure x-axis title and text appear only at the bottom
    axis.title.x.bottom = element_text(),
    axis.text.x.bottom = element_text()
  )
gg_richness_stacked

# don't run, not used
gg_richness_vertically <- plot_richness(phylobj_rare, x = "Soil.Classification", measures = c("Shannon")) +
  geom_boxplot(aes(fill = Soil.Classification)) +
  facet_grid(LTSP.Treatment ~ ., scales = "free_x", space = "free_x") +
  xlab("Soil Classification") +
  ggtitle("Alpha Diversity Richness for Soil Classification Controlled by LTSP Treatment") +
  theme(
    strip.background = element_blank(),
    strip.placement = "outside",
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.title.x.bottom = element_text(),
    #axis.text.x.bottom = element_text(angle = 90, hjust = 1, size = 14), # Adjust x-axis label size
    axis.title.y = element_text(size = 14), # Adjust y-axis title size
    axis.text.y = element_text(size = 14), # Adjust y-axis label size
    panel.spacing = unit(1, "lines")
  )
gg_richness_vertically




#----------------------------------------compaction pcoa plot not used---------------------------------------------------
# compaction plots
gg_pcoa_compac <- plot_ordination(phylobj_rare, pcoa_bc, color = "Compaction.Treatment") + 
  labs(pch="Compaction Treatment", col = "Compaction Treatment") + 
  stat_ellipse()
gg_pcoa_compac

ggsave("plot_pcoa_compac.png"
       , gg_pcoa_compac
       , height=5, width=8)




