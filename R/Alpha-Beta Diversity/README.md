# Plots Ready for Slides 🏄🏽‍♀️
Figure 1:
- Figure 1 and stats DONE
- Supplemental 1 with stats DONE

Figure 2:
- Figure 2 plot and stats DONE
- Supplemental 2 DONE (don't need stats)

Figure 3:
- Figure 3 Plot and stats
- Supplemental 3 and stats done

Figure 4:
- Figure 4 plot and stats DONE

For figure 5, see other folder

**Note below are only the color blind friendly plots, the normal ones can also be found in the files section if needed**

___________________________________________________________________________________

# Figure 1: Soil Classification Alpha Diversity Richness - Controlled for LTSP Treatment

![Figure 1 color blind](https://github.com/cynthiaachung/micb575-team3/blob/main/R/Alpha-Beta%20Diversity/Figure_1_AlphaDiversityRichness_horizontal_panels_colorblind.png)

## K-W stats analysis

For REF : The difference in Shannon diversity across soil classifications is statistically significant. P-value: 0.03456481 

For OM1 : The difference in Shannon diversity across soil classifications is statistically significant. P-value: 0.0001358425 

For OM2 : The difference in Shannon diversity across soil classifications is statistically significant. P-value: 5.283175e-05 

For OM3 : The difference in Shannon diversity across soil classifications is statistically significant. P-value: 0.03308952 

Note from prof: Talk about trends


## Supplemental: Phylogenetic diversity soil classification **DONE**

![Supplemental 1 color blind](https://github.com/cynthiaachung/micb575-team3/blob/main/R/Alpha-Beta%20Diversity/Supplemental_1_colorblind.png)

## K-W analysis

For  REF : The difference in Faith's Phylogenetic Diversity across soil classifications is statistically significant. P-value: 0.004645887 

For  OM1 : The difference in Faith's Phylogenetic Diversity across soil classifications is statistically significant. P-value: 6.034932e-05 

For  OM2 : The difference in Faith's Phylogenetic Diversity across soil classifications is statistically significant. P-value: 0.000103259 

For  OM3 : The difference in Faith's Phylogenetic Diversity across soil classifications is statistically significant. P-value: 0.01202364 

___________________________________________________________________________________

# Figure 2 Beta Diversity Bray-Curtis PCOA plot - Soil Classification:

![Beta](https://github.com/cynthiaachung/micb575-team3/blob/main/R/Alpha-Beta%20Diversity/Figure_2_PCOA_soil.classification.png)

## Permanova Stats

**Result: Pr(>F) values for Soil.Classification is 0.001 meaning there is a significant difference**

![permanova](https://github.com/cynthiaachung/micb575-team3/blob/main/R/Alpha-Beta%20Diversity/Screen%20Shot%202024-04-04%20at%2010.07.11%20PM.png)

## Soils grouped if plot needed:
![Soils grouped](https://github.com/cynthiaachung/micb575-team3/blob/main/R/Alpha-Beta%20Diversity/Figure_2.2_SoilGrouped.png)

## Soils grouped color blind pallet, if plot needed:
![Soils grouped colorblind](https://github.com/cynthiaachung/micb575-team3/blob/main/R/Alpha-Beta%20Diversity/Figure_2.2_SoilGrouped_colorblind.png)

## Supplemental 2: Taxa Bar Plot Genus Soil Classification **DONE**

![Taxa-Genus Soil Classification](https://github.com/cynthiaachung/micb575-team3/blob/main/R/Alpha-Beta%20Diversity/Supplemental_2_TaxaBarPlot_colorblind.png)

___________________________________________________________________________________

# Figure 3: Alpha Diversity Richness Compaction, controlled by soil group

## Separated by soil groups 

![Figure 3](https://github.com/cynthiaachung/micb575-team3/blob/main/R/Alpha-Beta%20Diversity/Figure_3_AlphaDiversityRichness_compac_vs_soilgroup_colorblind.png)

## KW STATS

For Soil Group 1 : The difference in Shannon diversity across compaction treatment is not statistically significant. P-value: 0.6341708 

For Soil Group 2 : The difference in Shannon diversity across compaction treatment is not statistically significant. P-value: 0.2529662 

For Soil Group 3 : The difference in Shannon diversity across compaction treatment is not statistically significant. P-value: 0.2782093 

## Supplemental 3 Faith's PD, compaction treatment

![Supplemental 3](https://github.com/cynthiaachung/micb575-team3/blob/main/R/Alpha-Beta%20Diversity/Supplemental_3_colorblind.png)

## Supplemental 3 stats KW

The difference in Faith's Phylogenetic Diversity across Compaction Treatment is not statistically significant. P-value: 0.2585947 

___________________________________________________________________________________

# Figure 4: Alpha Diversity compaction controlled for LTSP treatment

## Controlled for LTSP

color blind pallet:

![Figure 4](https://github.com/cynthiaachung/micb575-team3/blob/main/R/Alpha-Beta%20Diversity/Figure_4_AlphaDiversityRichness_compac_vs_LTSP_colorblind.png)

### K-W stats analysis 

For OM1 : The difference in Shannon diversity across compaction treatment is not statistically significant. P-value: 0.324537 

For OM2 : The difference in Shannon diversity across compaction treatment is statistically significant. P-value: 0.006780397 

___________________________________________________________________________________

# Figure 5 Volcano Plots: 
see other folder

___________________________________________________________________________________
# OLD DATA BELOW, NOT USED

## Rarefaction Curve with R

Looks different from QIIME. So a Phyloseq object with the previous rarefaction parameter (2526) is created and another one using 1500 arbitrarly to lose less samples
Both phyloseq objects are used to create curves
![Rarefaction Curve Generated by R](https://github.com/cynthiaachung/micb575-team3/blob/main/R/Alpha-Beta%20Diversity/Rarefaction%20Curve%20from%20R)
Rarefaction parameters:
- Rare1 = 2526
- Rare2 = 1500

_I will be using Rare1 data for the plots below, unless I see very little data shown in the plots, then I will change to Rare2 and do the plots all again, but for now, only Rare1 is used_

## Alpha Diversity

### We Must First Control for other Variables!!!
**Here is how the RICHNESS plots looks without controlling for other variables, there is too much variability, so we should control for other variables**

*Compaction Treatment - Richness (Not Controlled) - Rare1*


![Richness Not Filtered](https://github.com/cynthiaachung/micb575-team3/blob/main/R/Alpha-Beta%20Diversity/plot_richness_compac_not_filtered.png)

*Soil Classification - Richness (Not Controlled) - Rare1*


![Richness Not Filtered](https://github.com/cynthiaachung/micb575-team3/blob/main/R/Alpha-Beta%20Diversity/plot_richness_not_filtered.png)

**Explanation:**
According to the Richness plots above, we can see a large number of points outside the boxplots, indicating that there is a large range of alpha diversity measures within each soil classification OR compaction treatment. This suggests that we need to explore controlling for other variables that make affect alpha diversity measures such as "Site", "Hebacide Use"...



### Richness and Evenness
The goal here is to first select for other variables before plotting alpha diversity richness plots

(note: I am still having issues controlling for variables... working on these plots so to be continued...)



### Phylogentic Diversity

*Soil Classification - Phylogenetic Diversity (Not Controlled) - Rare1*![Phylogenetic Diversity](https://github.com/cynthiaachung/micb575-team3/blob/main/R/Alpha-Beta%20Diversity/plot_pd.png)
**Explanation:**  Visually we can see that the box plots for soil type phylogenetic diversity have some variability in phylogenetic diversity, most notably the group "Orthic Gray Luvisol, Gleyed Gray Luvisol". However we need to do statistical analysis to see if this difference is significant.

(stats to be continued...)

*Compaction Treatment - Phylogenetic Diversity (Not Controlled) - Rare1*![Phylogenetic Diversity](https://github.com/cynthiaachung/micb575-team3/blob/main/R/Alpha-Beta%20Diversity/plot_pd_compac.png)
**Explanation:**  Visually we can see that the box plots for compaction treatment phylogenetic diversity have very little difference in their median values. However we need to do statistical analysis to check if there might be any significant difference.

(stats to be continued...)

## Beta Diversity
------>>>>>>>>>>>>>Do paranova
### PCOA Plot - Bray-Curtis (Richness and Abundance/Eveness)

*Soil Classification - Beta Diversity PCOA plots*![Beta Diversity - PCOA plots](https://github.com/cynthiaachung/micb575-team3/blob/main/R/Alpha-Beta%20Diversity/plot_pcoa.png)

**Explanation:**  This plot really highlights that: 
- Soil Classification "Aquic Glossudalfs" has microbial communities that are different from the other groups.
- "Brunisolic Gray Luvisol" and "Orthis Gray Luvisol, Gleyed Gray Luvisol" are very close to eachother and are both farther away than other groups
- The rest of the soil types each have a specific cluster but overlap somewhat
- It is worth to note that "NA" forms a tight cluster indicating that they may all be similar soil types
(stats to be continued...)


*Compaction Treatment - Beta Diversity PCOA plots*![Beta Diversity - PCOA plots](https://github.com/cynthiaachung/micb575-team3/blob/main/R/Alpha-Beta%20Diversity/plot_pcoa_compac.png)

**Explanation:**  This plot shows visually that: 
- compaction treatment C0 and REF are very similar and compaction treatments C1 and C2 are farther away in terms of beta diversity.
- Important to note that there are a few red (C0) and purple (REF) points that are positioned close the C1 and C2.
- There are also not that many data points for C1 and C2
(stats to be continued...)


## Taxonomy bar plots

---> look at Genus

### Soil Classification

*Soil Classification - Taxa Bar Plots - Phylum*![Taxa-Phylum Plots](https://github.com/cynthiaachung/micb575-team3/blob/main/R/Alpha-Beta%20Diversity/plot_taxaphylum.png)


*Soil Classification - Taxa Bar Plots - Class*![Taxa-Class Plots](https://github.com/cynthiaachung/micb575-team3/blob/main/R/Alpha-Beta%20Diversity/plot_taxaclass.png)

**For fun, showing that they are all bacteria!**

*Soil Classification - Taxa Bar Plots - Domain*![Taxa-Domain Plots](https://github.com/cynthiaachung/micb575-team3/blob/main/R/Alpha-Beta%20Diversity/plot_taxadomain.png)

(stats to come...)


### Compaction Treatment

*Compaction Treatment - Taxa Bar Plots - Phylum*![Taxa-Phylum Plots](https://github.com/cynthiaachung/micb575-team3/blob/main/R/Alpha-Beta%20Diversity/plot_taxaphylum_compac.png)

*Compaction Treatment - Taxa Bar Plots - Class*![Taxa-Phylum Plots](https://github.com/cynthiaachung/micb575-team3/blob/main/R/Alpha-Beta%20Diversity/plot_taxaclass_compac.png)

(stats to come...)
