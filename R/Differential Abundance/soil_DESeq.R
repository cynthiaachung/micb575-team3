# if you didn't install the DESeq2 package, run the following
BiocManager::install("DESeq2")

#!/usr/bin/env Rscript
library(tidyverse)
library(phyloseq)
library(DESeq2)


#### Load data ####
load("desktop/M1/MICB_575/project_2/R exports/soil_rare.RData")
load("R/soil_final.RData")


#### DESeq ####
soil_deseq <- phyloseq_to_deseq2(soil_final, ~`Compaction.Treatment`)
DESEQ_soil <- DESeq(soil_deseq)

otu_table(soil_final)
sample_data(soil_final)

## NOTE: If you get a zeros error, then you need to add '1' count to all reads
soil_plus1 <- transform_sample_counts(soil_final, function(x) x+1) 
soil_deseq <- phyloseq_to_deseq2(soil_plus1, design = ~`Compaction.Treatment`)
DESEQ_soil <- DESeq(soil_deseq)


#Creating varibale comparing compaction treatment 0 to REF
res0 <- results(DESEQ_soil, tidy=TRUE, 
               #this will ensure that No is your reference group
               contrast = c("Compaction.Treatment","C0","REF"))

#Creating varibale comparing compaction treatment 1 to REF
res1 <- results(DESEQ_soil, tidy=TRUE, 
                #this will ensure that No is your reference group
                contrast = c("Compaction.Treatment","C1","REF"))

#Creating varibale comparing compaction treatment 2 to REF
res2 <- results(DESEQ_soil, tidy=TRUE, 
                #this will ensure that No is your reference group
                contrast = c("Compaction.Treatment","C2","REF"))


View(res)



#C0
## Volcano plot: effect size VS significance

ggplot(res0) +
  geom_point(aes(x=log2FoldChange, y=-log10(padj)))
  

## Make variable to color by whether it is significant + large change
vol_plot0 <- res0 %>%
  mutate(significant = padj<0.01 & abs(log2FoldChange)>2) %>%
  ggplot() +
  geom_point(aes(x=log2FoldChange, y=-log10(padj), col=significant)) +
  geom_vline(xintercept = 0, linetype="dotted")


#C1
## Volcano plot: effect size VS significance
ggplot(res1) +
  geom_point(aes(x=log2FoldChange, y=-log10(padj)))

## Make variable to color by whether it is significant + large change
vol_plot1 <- res1 %>%
  mutate(significant = padj<0.01 & abs(log2FoldChange)>2) %>%
  ggplot() +
  geom_point(aes(x=log2FoldChange, y=-log10(padj), col=significant))+
  geom_vline(xintercept = 0, linetype="dotted")


#C2
## Volcano plot: effect size VS significance
ggplot(res2) +
  geom_point(aes(x=log2FoldChange, y=-log10(padj)))

## Make variable to color by whether it is significant + large change
vol_plot2 <- res2 %>%
  mutate(significant = padj<0.05 & abs(log2FoldChange)>2) %>%
  ggplot() +
  geom_point(aes(x=log2FoldChange, y=-log10(padj), col=significant))+
  geom_vline(xintercept = 0, linetype="dotted")


vol_plot2

#Save files
ggsave(filename="vol_plot2.png",vol_plot2)
# To get table of results
sigASVs0 <- res0 %>% 
  filter(padj<0.01 & abs(log2FoldChange)>2) %>%
  dplyr::rename(ASV=row)
View(sigASVs0)

sigASVs1 <- res1 %>% 
  filter(padj<0.01 & abs(log2FoldChange)>2) %>%
  dplyr::rename(ASV=row)
View(sigASVs1)

sigASVs2 <- res1 %>% 
  filter(padj<0.01 & abs(log2FoldChange)>2) %>%
  dplyr::rename(ASV=row)
View(sigASVs2)


# Get only asv names
sigASVs_vec0 <- sigASVs0 %>%
  pull(ASV)

# Prune phyloseq file
mpt_DESeq0 <- prune_taxa(sigASVs_vec0,soil_final)
sigASVs0_names <- tax_table(DESEQ_soil) %>% 
  as.data.frame() %>%
  rownames_to_column(var="ASV") %>%
  right_join(sigASVs0) %>%
  arrange(log2FoldChange) %>%
  mutate(Species = make.unique(Species)) %>%
  mutate(Species = factor(Species, levels=unique(Species)))

sigASVs0_names <- tax_table(mpt_DESeq0) 
sigASVs0_names <- as.data.frame(sigASVs0_names)
sigASVs0_names <- sigASVs0_names %>%
  rownames_to_column(var="ASV") %>%
  right_join(sigASVs0) %>%
  arrange(log2FoldChange) %>%
  mutate(Species = make.unique(Species)) %>%
  mutate(Species = factor(Species, levels=unique(Species)))

#Create bar plot of ASV names that are significantly changed 
ggplot(sigASVs0_names) +
  geom_bar(aes(x=Species, y=log2FoldChange), stat="identity")+
  geom_errorbar(aes(x=Species, ymin=log2FoldChange-lfcSE, ymax=log2FoldChange+lfcSE)) +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))


# Get only asv names for res1
sigASVs_vec1 <- sigASVs1 %>%
  pull(ASV)

# Prune phyloseq file
mpt_DESeq1 <- prune_taxa(sigASVs_vec1,soil_final)
sigASVs1_names <- tax_table(DESEQ_soil) %>% as.data.frame() %>%
  rownames_to_column(var="ASV") %>%
  right_join(sigASVs1) %>%
  arrange(log2FoldChange) %>%
  mutate(Genus = make.unique(Genus)) %>%
  mutate(Genus = factor(Genus, levels=unique(Genus)))

sigASVs1_names <- tax_table(mpt_DESeq1) 
sigASVs1_names <- as.data.frame(sigASVs1_names)
sigASVs1_names <- sigASVs1_names %>%
  rownames_to_column(var="ASV") %>%
  right_join(sigASVs1) %>%
  arrange(log2FoldChange) %>%
  mutate(Genus = make.unique(Genus)) %>%
  mutate(Genus = factor(Genus, levels=unique(Genus)))

#Create bar plot of ASV names that are significantly changed 
ggplot(sigASVs1_names) +
  geom_bar(aes(x=Genus, y=log2FoldChange), stat="identity")+
  geom_errorbar(aes(x=Genus, ymin=log2FoldChange-lfcSE, ymax=log2FoldChange+lfcSE)) +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))

#species
sigASVs_vec1 <- sigASVs1 %>%
  pull(ASV)

# Prune phyloseq file
mpt_DESeq1 <- prune_taxa(sigASVs_vec1,soil_final)
sigASVs1_names <- tax_table(DESEQ_soil) %>% as.data.frame() %>%
  rownames_to_column(var="ASV") %>%
  right_join(sigASVs1) %>%
  arrange(log2FoldChange) %>%
  mutate(Species = make.unique(Species)) %>%
  mutate(Species = factor(Species, levels=unique(Species)))

sigASVs1_names <- tax_table(mpt_DESeq1) 
sigASVs1_names <- as.data.frame(sigASVs1_names)
sigASVs1_names <- sigASVs1_names %>%
  rownames_to_column(var="ASV") %>%
  right_join(sigASVs1) %>%
  arrange(log2FoldChange) %>%
  mutate(Species = make.unique(Species)) %>%
  mutate(Species = factor(Species, levels=unique(Species)))

#Create bar plot of ASV names that are significantly changed 
ggplot(sigASVs1_names) +
  geom_bar(aes(x=Species, y=log2FoldChange), stat="identity")+
  geom_errorbar(aes(x=Species, ymin=log2FoldChange-lfcSE, ymax=log2FoldChange+lfcSE)) +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))

# Get only asv names for res2
sigASVs_vec2 <- sigASVs2 %>%
  pull(ASV)

# Prune phyloseq file
mpt_DESeq2 <- prune_taxa(sigASVs_vec1,soil_final)
sigASVs2_names <- tax_table(DESEQ_soil) %>% as.data.frame() %>%
  rownames_to_column(var="ASV") %>%
  right_join(sigASVs2) %>%
  arrange(log2FoldChange) %>%
  mutate(Genus = make.unique(Genus)) %>%
  mutate(Genus = factor(Genus, levels=unique(Genus)))

sigASVs2_names <- tax_table(mpt_DESeq2) 
sigASVs2_names <- as.data.frame(sigASVs2_names)
sigASVs2_names <- sigASVs2_names %>%
  rownames_to_column(var="ASV") %>%
  right_join(sigASVs2) %>%
  arrange(log2FoldChange) %>%
  mutate(Genus = make.unique(Genus)) %>%
  mutate(Genus = factor(Genus, levels=unique(Genus)))

#Create bar plot of ASV names that are significantly changed 
ggplot(sigASVs2_names) +
  geom_bar(aes(x=Genus, y=log2FoldChange), stat="identity")+
  geom_errorbar(aes(x=Genus, ymin=log2FoldChange-lfcSE, ymax=log2FoldChange+lfcSE)) +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))




