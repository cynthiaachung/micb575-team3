library(tidyverse)
library(phyloseq)
library(indicspecies)

#Loading phyloseq object 
load("C:/Users/nmhow/OneDrive/Desktop/soil_final.RData")
View(soil_final)
soil_final

#Using taxglom 
soil_new <- subset_samples(soil_final,!is.na(Compaction.Treatment))  
mpt_genus <- tax_glom(soil_new, "Genus", NArm = FALSE)
mpt_genus_RA <- transform_sample_counts(mpt_genus, fun=function(x) x/sum(x))

isa_mpt <- multipatt(t(otu_table(mpt_genus_RA)), cluster = sample_data(mpt_genus_RA)$`Compaction.Treatment`)
summary(isa_mpt)

tax_table(soil_new) %>% as.data.frame() %>% rownames_to_column(var="ASV")

isa<-isa_mpt$sign %>%
  rownames_to_column(var="ASV") %>%
  left_join(taxtable) %>%
  filter(p.value<0.05) 
  

write.csv(isa, "C:/Users/nmhow/OneDrive/Documents/UBC - Masters/MICB 575/micb575-team3/isa.csv", row.names=FALSE)

 

soil_new2 <- subset_samples(soil_final,!is.na(Soil.Classification))  
mpt_genus <- tax_glom(soil_new2, "Genus", NArm = FALSE)
mpt_genus_RA <- transform_sample_counts(mpt_genus, fun=function(x) x/sum(x))

isa_mpt <- multipatt(t(otu_table(mpt_genus_RA)), cluster = sample_data(mpt_genus_RA)$`Soil.Classification`)
summary(isa_mpt)

tax_table(soil_new) %>% as.data.frame() %>% rownames_to_column(var="ASV")

isa2<-isa_mpt$sign %>%
  rownames_to_column(var="ASV") %>%
  left_join(taxtable) %>%
  filter(p.value<0.05) 


write.csv(isa2, "C:/Users/nmhow/OneDrive/Documents/UBC - Masters/MICB 575/micb575-team3/isa2.csv", row.names=FALSE)




