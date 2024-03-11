
# setup -------------------------------------------------------------------
library(dplyr)
library(tidyverse)
library(corrplot)
library(ggplot2)
library(tidyr)
library(nnet)

# import data -------------------------------------------------------------
data <- read.csv("data/metadata.csv", header = TRUE) %>%
  # select variables of interest
  select(X.SampleID, Sample.Accession, Sample.Alias, Ecozone, Region, Site, Country, Collection.Date,
         Soil.Classification, Tree.Cover, LTSP.Treatment, Compaction.Treatment, pH, CN.Ratio) %>%
  # rename variables
  rename(sample_id = X.SampleID,
         sample_accession = Sample.Accession,
         sample_alias = Sample.Alias,
         location_ecozone = Ecozone,
         location_region = Region,
         location_site = Site,
         location_country = Country,
         collection_date = Collection.Date,
         info_soil = Soil.Classification,
         info_tree = Tree.Cover,
         info_ph = pH,
         info_cnratio = CN.Ratio,
         treatment_ltsp = LTSP.Treatment,
         treatment_compaction = Compaction.Treatment)

# metadata wrangling ------------------------------------------------------
# method 1: separating data into columns with binary
data_separated <- data %>%
  # tree cover
  mutate(tree_fir_doug = ifelse(str_detect(info_tree, "Douglas fir"), "1", "0"), # douglas fir
         tree_pine_lodge = ifelse(str_detect(info_tree, "Lodgepole pine"), "1", "0"), # lodgepole pine
         tree_fir_subalp = ifelse(str_detect(info_tree, "Subalpine fir"), "1", "0"), # subalpine fir
         tree_spruce_int = ifelse(str_detect(info_tree, fixed("Interior Spruce", ignore_case=TRUE)), "1", "0"), # interior spruce
         tree_spruce_bl = ifelse(str_detect(info_tree, "Black Spruce"), "1", "0"), # black spruce
         tree_pine_pond = ifelse(str_detect(info_tree, "Ponderosa pine"), "1", "0"), # ponderosa pine
         tree_pine_sug = ifelse(str_detect(info_tree, "sugar pine"), "1", "0"), # sugar pine
         tree_fir_wh = ifelse(str_detect(info_tree, "white fir"), "1", "0"), # white fir
         tree_seq_gi = ifelse(str_detect(info_tree, "giant sequoia"), "1", "0"), # giant sequoia
         tree_pine_jack = ifelse(str_detect(info_tree, "Jack Pine"), "1", "0"), # jack pine
         tree_fir_bals = ifelse(str_detect(info_tree, "Balsam fir"), "1", "0"), # balsam fir
         tree_birch_wh = ifelse(str_detect(info_tree, "White birch"), "1", "0"), # white birch
         tree_pine_red = ifelse(str_detect(info_tree, "Red Pine"), "1", "0"), # red pine
         tree_pine_lob = ifelse(str_detect(info_tree, "Loblolly Pine"), "1", "0"), # loblolly pine
         tree_beauberr = ifelse(str_detect(info_tree, "Beautyberry"), "1", "0"), # beautyberry
         tree_yaup = ifelse(str_detect(info_tree, "Yaupon"), "1", "0"), # yaupon 
         tree_sweetgum = ifelse(str_detect(info_tree, "Sweetgum"), "1", "0"), #sweetgum
         tree_oak = ifelse(str_detect(info_tree, "Oaks"), "1", "0"), # oaks
         tree_waxmyr = ifelse(str_detect(info_tree, "Wax Myrtle"), "1", "0"), # wax myrtle
  ) %>%
  # soil classification
  mutate(soil_brun = ifelse(str_detect(info_soil, "Brunisolic Gray Luvisol"), "1", "0"), # brunisolic gray luvisol
         soil_orth_hfp = ifelse(str_detect(info_soil, "Orthic Humo-Ferric Podzol"), "1", "0"),
         soil_gley_edb = ifelse(str_detect(info_soil, "Gleyed Eluviated Dystric Brunisol"), "1", "0"),
         soil_orth_gl = ifelse(str_detect(info_soil, "Orthic Gray Luvisol"), "1", "0"),
         soil_gley_gl = ifelse(str_detect(info_soil, "Gleyed Gray Luvisol"), "1", "0"),
         soil_orth_db = ifelse(str_detect(info_soil, "Orthic Dystric Brunisol"), "1", "0"),
         soil_gley_db = ifelse(str_detect(info_soil, "Gleyed Dystric Brunisol"), "1", "0"),
         soil_mes_uh = ifelse(str_detect(info_soil, "Mesic Ultic Haploxeralfs"), "1", "0"),
         soil_aq_gl = ifelse(str_detect(info_soil, "Aquic Glossudalfs"), "1", "0"),
         soil_na = ifelse(is.na(info_soil), "1", "0")
  ) %>%
  # ltsp
  mutate(ltsp_om1 = ifelse(treatment_ltsp=="OM1", "1", "0"),
         ltsp_om2 = ifelse(treatment_ltsp=="OM2", "1", "0"),
         ltsp_om3 = ifelse(treatment_ltsp=="OM3", "1", "0"),
         ltsp_ref = ifelse(treatment_ltsp=="REF", "1", "0"),
  ) %>%
  # compaction
  mutate(compact_1 = ifelse(treatment_compaction=="C1", "1", "0"),
         compact_2 = ifelse(treatment_compaction=="C2", "1", "0"),
         compact_3 = ifelse(treatment_compaction=="C3", "1", "0"),
         compact_ref = ifelse(treatment_compaction=="REF", "1", "0"),
  ) %>%
  # remove unsorted columns
  subset(select = -c(info_tree, info_soil, treatment_ltsp, treatment_compaction)) %>%
  # changes columns to numeric
  mutate_at(c(11:47), as.numeric)

# method 2: repeating rows for unique observations  
data_repeated <- data %>%
  # to ensure consistent capitalization
  mutate(info_tree = str_to_title(info_tree),
         info_soil = str_to_title(info_soil)) %>%
  # tree cover
  mutate(info_tree = strsplit(as.character(info_tree), ", ")) %>%
  unnest(info_tree) %>%
  # soil classification
  mutate(info_soil = strsplit(as.character(info_soil), ", ")) %>%
  unnest(info_soil)

# heat maps ---------------------------------------------------------------
# soil and tree
## frequency table
heat_soiltree <- data_repeated %>%
  group_by(info_soil, info_tree) %>%
  summarize(n = n(), .groups = "drop")
## heat map
ggplot(heat_soiltree, aes(x = info_soil, y = info_tree, fill = n)) +
  geom_tile(colour = "white",
            size=0.5) +
  labs(title = "Heatmap of Soil Classification versus Tree Cover",
       x = "Soil Classification",
       y = "Tree Cover",
       fill = "Count") +
  scale_fill_gradient(low = "white",
                      high = "blue",
                      na.value = "transparent",
                      limits = c(0,155)) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 1, hjust=1),
        axis.text.y = element_text(hjust=1))
## saves plot
ggsave("output/heatmap/heat_soiltree.png", width = 25, height = 20, units = "cm")

# soil and ltsp
## frequency table
heat_soilltsp <- data_repeated %>%
  group_by(info_soil, treatment_ltsp) %>%
  summarize(n = n(), .groups = "drop")
## heat map
ggplot(heat_soilltsp, aes(x = info_soil, y = treatment_ltsp, fill = n)) +
  geom_tile(colour = "white",
            size=0.5) +
  labs(title = "Heatmap of Soil Classification versus LTSP Treatment",
       x = "Soil Classification",
       y = "LTSP Treatment",
       fill = "Count") +
  scale_fill_gradient(low = "white",
                      high = "blue",
                      na.value = "transparent",
                      limits = c(0,300)) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 1, hjust=1),
        axis.text.y = element_text(hjust=1))
## saves plot
ggsave("output/heatmap/heat_soilltsp.png", width = 25, height = 20, units = "cm")

# soil and compaction
## frequency table
heat_soilcompact <- data_repeated %>%
  group_by(info_soil, treatment_compaction) %>%
  summarize(n = n(), .groups = "drop")
## heat map
ggplot(heat_soilcompact, aes(x = info_soil, y = treatment_compaction, fill = n)) +
  geom_tile(colour = "white",
            size=0.5) +
  labs(title = "Heatmap of Soil Classification versus Compaction",
       x = "Soil Classification",
       y = "Compaction",
       fill = "Count") +
  scale_fill_gradient(low = "white",
                      high = "blue",
                      na.value = "transparent",
                      limits = c(0,1000)) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 1, hjust=1),
        axis.text.y = element_text(hjust=1))
## saves plot
ggsave("output/heatmap/heat_soilcompact.png", width = 25, height = 20, units = "cm")

# soil and ph
## frequency table
heat_soilph <- data_repeated %>%
  mutate(ph_cat = (as.factor(ifelse(info_ph <= 3, "≤3",
                                    ifelse(info_ph < 4, "4", 
                                           ifelse(info_ph < 5, "5", 
                                                  ifelse(info_ph <= 6, "6", ">6")))))),
         ph_cat = factor(ph_cat, levels = c("≤3", "4", "5", "6", ">6"))
  ) %>%
  group_by(info_soil, ph_cat) %>%
  summarize(n = n(), .groups = "drop")
## heat map
ggplot(heat_soilph, aes(x = info_soil, y = ph_cat, fill = n)) +
  geom_tile(colour = "white",
            size=0.5) +
  labs(title = "Heatmap of Soil Classification versus pH",
       x = "Soil Classification",
       y = "pH",
       fill = "Count") +
  scale_fill_gradient(low = "white",
                      high = "blue",
                      na.value = "transparent",
                      limits = c(0,575)) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 1, hjust=1),
        axis.text.y = element_text(hjust=1))
## saves plot
ggsave("output/heatmap/heat_soilph.png", width = 25, height = 20, units = "cm")

# soil and cn ratio
## frequency table
heat_soilcn <- data_repeated %>%
  mutate(cn_cat = (as.factor(ifelse(info_cnratio <= 10, "≤10",
                                    ifelse(info_cnratio <= 20, "10-20", 
                                           ifelse(info_cnratio <= 30, "20-30", 
                                                  ifelse(info_cnratio <= 40, "30-40", 
                                                         ifelse(info_cnratio <= 50, "40-50", ">50"))))))),
         cn_cat = factor(cn_cat, levels = c("≤10", "10-20", "20-30", "30-40", "40-50", ">50"))
  ) %>%
  group_by(info_soil, cn_cat) %>%
  summarize(n = n(), .groups = "drop")
## heat map
ggplot(heat_soilcn, aes(x = info_soil, y = cn_cat, fill = n)) +
  geom_tile(colour = "white",
            size=0.5) +
  labs(title = "Heatmap of Soil Classification versus Carbon-Nitrogen Ratio",
       x = "Soil Classification",
       y = "CN Ratio",
       fill = "Count") +
  scale_fill_gradient(low = "white",
                      high = "blue",
                      na.value = "transparent",
                      limits = c(0,325)) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 1, hjust=1),
        axis.text.y = element_text(hjust=1))
## saves plot
ggsave("output/heatmap/heat_soilcn.png", width = 25, height = 20, units = "cm")

# multinomial regression analysis -----------------------------------------
# specify reference group
data_multinom <- data_repeated %>%
  mutate(info_soil = as.factor(info_soil)) %>%
  mutate(info_soil = relevel(info_soil, ref = "Orthic Gray Luvisol")) %>%
  # create categories for ph
  mutate(ph_cat = (as.factor(ifelse(info_ph <= 3, "≤3",
                                    ifelse(info_ph < 4, "4", 
                                           ifelse(info_ph < 5, "5", 
                                                  ifelse(info_ph <= 6, "6", ">6")))))),
         ph_cat = factor(ph_cat, levels = c("≤3", "4", "5", "6", ">6"))
  ) %>%
  # create categories for cn ratio
  mutate(cn_cat = (as.factor(ifelse(info_cnratio <= 10, "≤10",
                                    ifelse(info_cnratio <= 20, "10-20", 
                                           ifelse(info_cnratio <= 30, "20-30", 
                                                  ifelse(info_cnratio <= 40, "30-40", 
                                                         ifelse(info_cnratio <= 50, "40-50", ">50"))))))),
         cn_cat = factor(cn_cat, levels = c("≤10", "10-20", "20-30", "30-40", "40-50", ">50"))
  )

# soil and tree
## run the multinomial model, tree cover is independent variable and soil classification is the dependent variable 
multinom_soiltree <- multinom(info_soil ~ info_tree, data = data_multinom)
## check summary
summary(multinom_soiltree)
## extract coefficients and convert to data frame
coefmatrix_soiltree <- coef(multinom_soiltree)
coefmatrix_soiltree <- as.data.frame(coefmatrix_soiltree) %>%
  # reset row names to make them a column in the dataframe
  mutate(info_tree = rownames(coefmatrix_soiltree))
## melts the data frame to make it suitable for ggplot2
coefmatrix_soiltree_melt <- melt(coefmatrix_soiltree, id.vars = "info_tree") %>%
  # removes intercept column
  subset(variable != "(Intercept)") %>%
  # removes prefix
  mutate(variable = sub("^info_tree", "", variable))
## create barplot
ggplot(coefmatrix_soiltree_melt, aes(x = info_tree, y = value, fill = variable)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  labs(title = "Coefficients from Multinomial Logistic Regression for Tree Cover and Soil Classification",
       subtitle = "Reference Soil Type: Orthic Gray Luvisol",
       x = "Soil Type", y = "Coefficient", fill = "Tree Type") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 90, vjust = 1, hjust=1),
        axis.text.y = element_text(hjust=1))
## saves plot
ggsave("output/multinomial/multi_soiltree.png", width = 25, height = 20, units = "cm")

# soil and ltsp
## run the multinomial model, ltsp treatment is independent variable and soil classification is the dependent variable 
multinom_soilltsp <- multinom(info_soil ~ treatment_ltsp, data = data_multinom)
## check summary
summary(multinom_soilltsp)
## extract coefficients and convert to data frame
coefmatrix_soilltsp <- coef(multinom_soilltsp)
coefmatrix_soilltsp <- as.data.frame(coefmatrix_soilltsp) %>%
  # reset row names to make them a column in the dataframe
  mutate(treatment_ltsp = rownames(coefmatrix_soilltsp))
## melts the data frame to make it suitable for ggplot2
coefmatrix_soilltsp_melt <- melt(coefmatrix_soilltsp, id.vars = "treatment_ltsp") %>%
  # removes intercept column
  subset(variable != "(Intercept)") %>%
  # removes prefix
  mutate(variable = sub("^treatment_ltsp", "", variable))
## create barplot
ggplot(coefmatrix_soilltsp_melt, aes(x = treatment_ltsp, y = value, fill = variable)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  labs(title = "Coefficients from Multinomial Logistic Regression for LTSP Treatment and Soil Classification",
       subtitle = "Reference Soil Type: Orthic Gray Luvisol",
       x = "Soil Type", y = "Coefficient", fill = "LTSP Treatment") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 90, vjust = 1, hjust=1),
        axis.text.y = element_text(hjust=1))
## saves plot
ggsave("output/multinomial/multi_soilltsp.png", width = 25, height = 20, units = "cm")

# soil and compaction
## run the multinomial model, compaction is independent variable and soil classification is the dependent variable 
multinom_soilcompact <- multinom(info_soil ~ treatment_compaction, data = data_multinom)
## check summary
summary(multinom_soilcompact)
## extract coefficients and convert to data frame
coefmatrix_soilcompact <- coef(multinom_soilcompact)
coefmatrix_soilcompact <- as.data.frame(coefmatrix_soilcompact) %>%
  # reset row names to make them a column in the dataframe
  mutate(treatment_compaction = rownames(coefmatrix_soilcompact))
## melts the data frame to make it suitable for ggplot2
coefmatrix_soilcompact_melt <- melt(coefmatrix_soilcompact, id.vars = "treatment_compaction") %>%
  # removes intercept column
  subset(variable != "(Intercept)") %>%
  # removes prefix
  mutate(variable = sub("^treatment_compaction", "", variable))
## create barplot
ggplot(coefmatrix_soilcompact_melt, aes(x = treatment_compaction, y = value, fill = variable)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  labs(title = "Coefficients from Multinomial Logistic Regression for Compaction and Soil Classification",
       subtitle = "Reference Soil Type: Orthic Gray Luvisol",
       x = "Soil Type", y = "Coefficient", fill = "Compaction") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 90, vjust = 1, hjust=1),
        axis.text.y = element_text(hjust=1))
## saves plot
ggsave("output/multinomial/multi_soilcompact.png", width = 25, height = 20, units = "cm")

# soil and ph
## run the multinomial model, ph is independent variable and soil classification is the dependent variable 
multinom_soilph <- multinom(info_soil ~ ph_cat, data = data_multinom)
## check summary
summary(multinom_soilph)
## extract coefficients and convert to data frame
coefmatrix_soilph <- coef(multinom_soilph)
coefmatrix_soilph <- as.data.frame(coefmatrix_soilph) %>%
  # reset row names to make them a column in the dataframe
  mutate(ph_cat = rownames(coefmatrix_soilph))
## melts the data frame to make it suitable for ggplot2
coefmatrix_soilph_melt <- melt(coefmatrix_soilph, id.vars = "ph_cat") %>%
  # removes intercept column
  subset(variable != "(Intercept)") %>%
  # removes prefix
  mutate(variable = sub("^ph_cat", "", variable))
## create barplot
ggplot(coefmatrix_soilph_melt, aes(x = ph_cat, y = value, fill = variable)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  labs(title = "Coefficients from Multinomial Logistic Regression for pH and Soil Classification",
       subtitle = "Reference Soil Type: Orthic Gray Luvisol",
       x = "Soil Type", y = "Coefficient", fill = "pH Category") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 90, vjust = 1, hjust=1),
        axis.text.y = element_text(hjust=1))
## saves plot
ggsave("output/multinomial/multi_soilph.png", width = 25, height = 20, units = "cm")

# soil and cn ratio
## run the multinomial model, ph is independent variable and soil classification is the dependent variable 
multinom_soilcn <- multinom(info_soil ~ cn_cat, data = data_multinom)
## check summary
summary(multinom_soilcn)
## extract coefficients and convert to data frame
coefmatrix_soilcn <- coef(multinom_soilcn)
coefmatrix_soilcn <- as.data.frame(coefmatrix_soilcn) %>%
  # reset row names to make them a column in the dataframe
  mutate(cn_cat = rownames(coefmatrix_soilcn))
## melts the data frame to make it suitable for ggplot2
coefmatrix_soilcn_melt <- melt(coefmatrix_soilcn, id.vars = "cn_cat") %>%
  # removes intercept column
  subset(variable != "(Intercept)") %>%
  # removes prefix
  mutate(variable = sub("^cn_cat", "", variable))
## create barplot
ggplot(coefmatrix_soilcn_melt, aes(x = cn_cat, y = value, fill = variable)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  labs(title = "Coefficients from Multinomial Logistic Regression for CN Ratio and Soil Classification",
       subtitle = "Reference Soil Type: Orthic Gray Luvisol",
       x = "Soil Type", y = "Coefficient", fill = "CN Ratio Category") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 90, vjust = 1, hjust=1),
        axis.text.y = element_text(hjust=1))
## saves plot
ggsave("output/multinomial/multi_soilcn.png", width = 25, height = 20, units = "cm")

# archive -----------------------------------------------------------------
# variable subsets ---
## soil
soil <- data_wrangled %>%
  select((starts_with("soil")))
## trees
tree <- data_wrangled %>%
  select((starts_with("tree")))
## ltsp
ltsp <- data_wrangled %>%
  select((starts_with("ltsp")))
## compaction
compaction <- data_wrangled %>%
  select((starts_with("compact")))
## ph
ph <- data_wrangled$info_ph
## cn ratio
cnratio <- data_wrangled$info_cnratio

# correlation matrices ---
# soil and tree
## correlation matrix
corr_soiltree <- cor(tree, soil, use = "complete.obs") %>%
  round(2)
## correlation plot
corrplot_soiltree <- corrplot(corr_soiltree, 
                              type = "upper",
                              method = 'color',
                              addCoef.col = 'black',
                              tl.col = "black")

# soil and ltsp
## correlation matrix
corr_soilltsp <- cor(soil, ltsp, use = "complete.obs") %>%
  round(2)
## correlation plot
corrplot_soilltsp <- corrplot(corr_soilltsp, 
                              type = "full",
                              method = 'color',
                              addCoef.col = 'black',
                              tl.col = "black")

# soil and compaction
## correlation matrix
corr_soilcompact <- cor(soil, compaction, use = "complete.obs") %>%
  round(2)
## correlation plot
corrplot_soilcompact <- corrplot(corr_soilcompact, 
                                 type = "full",
                                 method = 'color',
                                 addCoef.col = 'black',
                                 tl.col = "black")

# soil and ph
## correlation matrix
corr_soilph <- cor(soil, ph, use = "complete.obs") %>%
  round(2)
## correlation plot
corrplot_soilph <- corrplot(corr_soilph, 
                            type = "full",
                            method = 'color',
                            addCoef.col = 'black',
                            tl.col = "black")

# soil and cn ratio
## correlation matrix
corr_soilcn <- cor(soil, cnratio, use = "complete.obs") %>%
  round(2)
## correlation plot
corrplot_soilcn <- corrplot(corr_soilcn, 
                            type = "full",
                            method = 'color',
                            addCoef.col = 'black',
                            tl.col = "black")

