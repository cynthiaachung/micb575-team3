library(tidyverse)
library(phyloseq)
library(indicspecies)

#Loading phyloseq object 
load("C:/Users/nmhow/OneDrive/Desktop/soil_final.RData")
View(soil_final)
soil_final

#
soil_new <- subset_samples(soil_final,!is.na(Soil.Classification))  
mpt_genus <- tax_glom(soil_new, "Genus", NArm = FALSE)
mpt_genus_RA <- transform_sample_counts(mpt_genus, fun=function(x) x/sum(x))

isa_mpt <- multipatt(t(otu_table(mpt_genus_RA)), cluster = sample_data(mpt_genus_RA)$`Soil.Classification`)
summary(isa_mpt)

tax_table(soil_new) %>% as.data.frame() %>% rownames_to_column(var="ASV")

isa_mpt$sign %>%
  rownames_to_column(var="ASV") %>%
  left_join(taxtable) %>%
  filter(p.value<0.05) %>% View()



