**See "Report Figures" Folder for all the figures used in the paper including formatted versions**

# Figure 1: Soil classifications are grouped into three distinct clusters

![Beta](https://github.com/cynthiaachung/micb575-team3/blob/main/R/Alpha-Beta%20Diversity/Figure_2_PCOA_soil.classification.png)

Figure 1. Soil classifications are grouped into three distinct clusters. The Bray-Curtis distances between the seven soil classification was plotted on a PCoA plot. Bray-Curtis provides a measurement of beta-diversity with respect to abundance and evenness, and is thus a metric for community dissimilarity. The soil classifications formed three distinct clusters: Soil Group 1 of Aquic Glassudalfs; Soil Group 2 of Brunisolic Gray Luvisol, Orthic Gray Luvisol, Gleyed Gray Luvisol; and Soil Group 3 of the remaining soil types. The ellipses were generated around each  soil group at a level of 0.95. The PERMANOVA test determined the significance between groups to be 0.001 (Pr(>F)). The viridis package was used to generate an inclusive color palette viewable to all readers.



## Permanova Stats

**Result: Pr(>F) values for Soil.Classification is 0.001 meaning there is a significant difference**

![permanova](https://github.com/cynthiaachung/micb575-team3/blob/main/R/Alpha-Beta%20Diversity/Screen%20Shot%202024-04-04%20at%2010.07.11%20PM.png)

## Soils grouped if plot needed:
![Soils grouped](https://github.com/cynthiaachung/micb575-team3/blob/main/R/Alpha-Beta%20Diversity/Figure_2.2_SoilGrouped.png)

## Soils grouped color blind pallet, if plot needed:
![Soils grouped colorblind](https://github.com/cynthiaachung/micb575-team3/blob/main/R/Alpha-Beta%20Diversity/Figure_2.2_SoilGrouped_colorblind.png)


## Supplemental 1: Organic matter removal significantly impacts phylogenetic distance of the soil microbial community

![Supplemental 1 color blind](https://github.com/cynthiaachung/micb575-team3/blob/main/R/Alpha-Beta%20Diversity/Supplemental_1_colorblind.png)

Supplemental 1. Organic matter removal significantly impacts phylogenetic distance of the soil microbial community. Faith’s phylogenetic diversity was measured across seven  soil classifications among the four organic removal levels (REF, OM1, OM2, OM3). Faith’s Phylogenetic diversity provides a measurement of alpha-diversity with respect to richness, abundance, and phylogenetic distance; the quantified metric makes up the y-axis. A Kruskal-Wallis statistical analysis was performed and significance was determined as p < 0.05 as indicated by the asterisk (*).

## K-W statistics analysis

For  REF : The difference in Faith's Phylogenetic Diversity across soil classifications is statistically significant. P-value: 0.004645887 

For  OM1 : The difference in Faith's Phylogenetic Diversity across soil classifications is statistically significant. P-value: 6.034932e-05 

For  OM2 : The difference in Faith's Phylogenetic Diversity across soil classifications is statistically significant. P-value: 0.000103259 

For  OM3 : The difference in Faith's Phylogenetic Diversity across soil classifications is statistically significant. P-value: 0.01202364 

___________________________________________________________________________________

# Figure 2: Organic matter removal influences changes in soil microbial diversity dependent on soil type

![Figure 1 color blind](https://github.com/cynthiaachung/micb575-team3/blob/main/R/Alpha-Beta%20Diversity/Figure_1_AlphaDiversityRichness_horizontal_panels_colorblind.png)

Figure 2. Organic matter removal influences changes in soil microbial diversity dependent on soil type. Shannon’s diversity was measured across seven soil types (color coded) to the reference plots (REF) where no organic matter was removed and with increasing organic matter removal (OM1-OM3). Shannon’s diversity provides a measurement of alpha-diversity with respect to richness and abundance; the quantified metric makes up the y-axis. A Kruskal-Wallis statistical analysis was performed for each soil classification’s Shannon’s diversity within each LTSP Treatment type. Significance was determined at p <0.05 as indicated by the asterisk (*). The viridis package was used to generate an inclusive color palette viewable to all readers.

## K-W statistics analysis

For REF : The difference in Shannon diversity across soil classifications is statistically significant. P-value: 0.03456481 

For OM1 : The difference in Shannon diversity across soil classifications is statistically significant. P-value: 0.0001358425 

For OM2 : The difference in Shannon diversity across soil classifications is statistically significant. P-value: 5.283175e-05 

For OM3 : The difference in Shannon diversity across soil classifications is statistically significant. P-value: 0.03308952 

Note from prof: Talk about trends



## Supplemental 2: Average genus abundance in each soil classification

![Taxa-Genus Soil Classification](https://github.com/cynthiaachung/micb575-team3/blob/main/R/Alpha-Beta%20Diversity/Supplemental_2_TaxaBarPlot_colorblind.png)

Supplemental 2. Average genus abundance in each soil classification. Orthic Gray Luvisol and Gleyed Gray Luvisol have been grouped into one. 5 genuses are found in the soil types with all soil classification containing Nocardia (dark purple), Pandoraea (dark blue), and Zymomonas (yellow). Three soil types, Aquic Glossudalfs, Orthic Dystric Brunisol, and Orthic Humo-Ferric Podzol contain Treponema (light green). Gleyed Dystric Brunisol is the only soil classification which contains Tepidbacter (dark green). The average genus abundance makes up the y-axis.

___________________________________________________________________________________

# Figure 3: Compaction level does not influence soil microbial diversity within soil groups

## Separated by soil groups 

![Figure 3](https://github.com/cynthiaachung/micb575-team3/blob/main/R/Alpha-Beta%20Diversity/Figure_3_AlphaDiversityRichness_compac_vs_soilgroup_colorblind.png)

Figure 3. Compaction level does not influence soil microbial diversity within soil groups.    The impact of compaction level (color coded) on the microbial diversity of each soil group was measured by Shannon’s Diversity. Shannon’s diversity provides a measurement of alpha-diversity with respect to richness and abundance; the quantified metric makes up the y-axis. Soil Group 1 and Soil Group 3 only have data for compaction treatment levels C0 (dark purple) and REF (yellow), while the Soil Group 2 contains C0, C1, and C2 (green). A Kruskal-Wallis statistical analysis was done within each of the three soil groups, significance was determined to be p > 0.05. There is no significant microbial diversity difference due to  compaction treatment within any soil groups. The viridis package was used to generate an inclusive color palette viewable to all readers.

## K-W statistics analysis

For Soil Group 1 : The difference in Shannon diversity across compaction treatment is not statistically significant. P-value: 0.6341708 

For Soil Group 2 : The difference in Shannon diversity across compaction treatment is not statistically significant. P-value: 0.2529662 

For Soil Group 3 : The difference in Shannon diversity across compaction treatment is not statistically significant. P-value: 0.2782093 

## Supplemental 3: Compaction Treatment does not have an impact on biodiversity

![Supplemental 3](https://github.com/cynthiaachung/micb575-team3/blob/main/R/Alpha-Beta%20Diversity/Supplemental_3_colorblind.png)

Supplemental 3. Compaction Treatment does not have an impact on biodiversity. Phylogenetic diversity was measured across each compaction treatment (C0-C2) as well as the reference plot. Significance was determined at p <0.05 and there was no significance between compaction treatments found. 

## K-W statistics analysis

The difference in Faith's Phylogenetic Diversity across Compaction Treatment is not statistically significant. P-value: 0.2585947 

___________________________________________________________________________________

# Figure 5: Removal of all organic matter (OM2) exacerbates the negative impact of compaction on soil microbial diversity

## Controlled for LTSP

color blind pallet:

![Figure 4](https://github.com/cynthiaachung/micb575-team3/blob/main/R/Alpha-Beta%20Diversity/Figure_4_AlphaDiversityRichness_compac_vs_LTSP_colorblind.png)

Figure 5. Removal of all organic matter (OM2) exacerbates the negative impact of compaction on soil microbial diversity. Shannon’s Diversity of across each compaction treatment level (C0-C2) was measured at all organic matter (OM) removal levels (REF, OM1, OM2, OM3). Shannon’s diversity provides a measurement of alpha-diversity with respect to richness and abundance; the quantified metric makes up the y-axis. The REF (yellow) only has data for compaction treatment REF, organic matter removal OM3 only has data for compaction treatment C0 (dark purple) and OM1 and OM2 only has data for CO, C1 (blue) and C2 (green) compaction treatment. Kruskal-Wallis was performed within the OM1 and OM2 LTSP treatment types to determine if there is a significant different in the alpha-diversity between the compaction treatment levels. A Wilcoxon rank-sum tests was performed between compaction treatment levels C0-C1, C1-C2, and C0-C2 within the OM2 LTSP treatment group. Significance was determined at p < 0.05 indicated by asterisk (*). The viridis package was used to generate an inclusive color palette viewable to all readers.

## K-W statistics analysis

For OM1 : The difference in Shannon diversity across compaction treatment is not statistically significant. P-value: 0.324537 

For OM2 : The difference in Shannon diversity across compaction treatment is statistically significant. P-value: 0.006780397 

### Wilcoxon Rank Sum 

![WX stats](https://github.com/cynthiaachung/micb575-team3/blob/main/R/Alpha-Beta%20Diversity/Screen%20Shot%202024-04-14%20at%205.39.20%20PM.png)

___________________________________________________________________________________

# Figure 4 Volcano Plots: 
see "Differential Abundance" folder
