library(dplyr)
library(ggplot2)
library(ggpubr)
library(superheat)
library(readr)
library(readxl)
library(tidyverse)
library(writexl)
library(reshape2)
library(superheat)
library(tidyverse)
library(cluster)
library(factoextra)
library(gridExtra)
library(pvclust)
library(dendextend)
library(devtools)
library(ComplexHeatmap)

setwd("C:\\Users\\erinm\\OneDrive\\Desktop\\GLRI Project\\Milwauke\\Data for Papier\\Mixtures MS\\Effect Evaluation\\Effect Data")

effectdata <- read_excel("C:\\Users\\erinm\\OneDrive\\Desktop\\GLRI Project\\Milwauke\\Data for Papier\\Effect Data\\Effect_Data_for_all_Analyses.xlsx") 
View(effectdata)

effect_data1 <- effectdata %>% select("Assay Name", "Assay Type", "Year", "Site", "Replicate", "Measurement", "Effect", "Units", "Sample Type")
View(effect_data1)

#2017####
#only significant ####
list(unique(effect_data1$Site))
effect2017  <- effect_data1 %>% filter(Year == "2017", !Site %in% c("KKL", "BK", "MQ-BK", "Blank"))

sigdiff_effect2017 <- effect2017 %>% filter(Measurement %in% c("E2-EQ", "17b_estradiol", "Testosterone", "Aryl hydrocarbon Receptor (AhR)",
                                                               "Estrogen Receptor alpha (ERa)", "Estrogen Response Element (ERE)",
                                                               "Antioxidant Response Element (ARE)-binding Nuclear factor (erythroid-derived 2)-like 2 (NRF2) (NRF2/ARE)",
                                                               "Peroxisome proliferator-activated receptor-gamma (PPARg)", "Peroxisome proliferator activating receptor (PPRE)",
                                                               "Pregnane-X-Receptor (PXR)", "CYP1A1_intestine", "CYP1A1_liver", "CYP2AD6_intestine",
                                                               "CYP2N13_intestine", "CYP3A_intestine", "CYP3A_liver", "UGT1A1_intestine", "UGT1A1_liver"))

sigdiff_effect2017$Measurement <- ifelse(sigdiff_effect2017$Measurement == "17b_estradiol" & sigdiff_effect2017$`Sample Type` == "Male Fathead Minnow Serum", "E2_M_plasma", sigdiff_effect2017$Measurement)
sigdiff_effect2017$Measurement <- ifelse(sigdiff_effect2017$Measurement == "17b_estradiol" & sigdiff_effect2017$`Sample Type` == "Female Fathead Minnow Serum", "E2_F_plasma", sigdiff_effect2017$Measurement)
sigdiff_effect2017$Measurement <- ifelse(sigdiff_effect2017$Measurement == "Testosterone" & sigdiff_effect2017$`Sample Type` == "Male Fathead Minnow Serum", "T_M_plasma", sigdiff_effect2017$Measurement)
sigdiff_effect2017$Measurement <- ifelse(sigdiff_effect2017$Measurement == "Aryl hydrocarbon Receptor (AhR)", "AhR_bioactivity", sigdiff_effect2017$Measurement)
sigdiff_effect2017$Measurement <- ifelse(sigdiff_effect2017$Measurement =="Estrogen Receptor alpha (ERa)", "ERa_bioactivity", sigdiff_effect2017$Measurement)
sigdiff_effect2017$Measurement <- ifelse(sigdiff_effect2017$Measurement =="Estrogen Response Element (ERE)", "ERE_bioactivity", sigdiff_effect2017$Measurement)
sigdiff_effect2017$Measurement <- ifelse(sigdiff_effect2017$Measurement =="Antioxidant Response Element (ARE)-binding Nuclear factor (erythroid-derived 2)-like 2 (NRF2) (NRF2/ARE)", "ARE_bioactivity", sigdiff_effect2017$Measurement)
sigdiff_effect2017$Measurement <- ifelse(sigdiff_effect2017$Measurement =="Peroxisome proliferator-activated receptor-gamma (PPARg)", "PPARg_bioactivity", sigdiff_effect2017$Measurement)
sigdiff_effect2017$Measurement <- ifelse(sigdiff_effect2017$Measurement =="Peroxisome proliferator activating receptor (PPRE)", "PPRE_bioactivity", sigdiff_effect2017$Measurement)
sigdiff_effect2017$Measurement <- ifelse(sigdiff_effect2017$Measurement =="Pregnane-X-Receptor (PXR)" & sigdiff_effect2017$`Assay Name` == "Attagene cis-Factorial", "PXRcis_bioactivity", sigdiff_effect2017$Measurement)
sigdiff_effect2017$Measurement <- ifelse(sigdiff_effect2017$Measurement =="Pregnane-X-Receptor (PXR)" & sigdiff_effect2017$`Assay Name` == "Attagene trans-Factorial", "PXRtrans_bioactivity", sigdiff_effect2017$Measurement)

list(unique(sigdiff_effect2017$Measurement))

mean_sdeffect2017 <- sigdiff_effect2017 %>% group_by(Site, Measurement)%>% summarize(meanEffect = mean(Effect))

mean_sdeffect2017$Measurement <- factor(mean_sdeffect2017$Measurement, levels= c("AhR_bioactivity", "PXRcis_bioactivity", "PXRtrans_bioactivity",
                                                                                 "ARE_bioactivity", "PPARg_bioactivity", "PPRE_bioactivity",
                                                                                 "CYP1A1_intestine", "CYP1A1_liver", "CYP2AD6_intestine", "CYP2N13_intestine", 
                                                                                 "CYP3A_intestine","CYP3A_liver",
                                                                                 "UGT1A1_intestine", "UGT1A1_liver",
                                                                                 "E2-EQ", "ERa_bioactivity", "ERE_bioactivity",
                                                                                 "E2_F_plasma", "E2_M_plasma", "T_M_plasma"))


wide_sdeffect2017 <- mean_sdeffect2017 %>% select("Site", "Measurement", "meanEffect") %>% spread(key = "Measurement", value = "meanEffect") %>%
  replace(is.na(.), 0) %>% column_to_rownames("Site") 
View(wide_sdeffect2017)

cluster_sdeffect2017 <- scale(wide_sdeffect2017, center = TRUE, scale = TRUE)
View(cluster_sdeffect2017)

t_cluster_sdeffect2017 <- t(cluster_sdeffect2017)

View(t_cluster_sdeffect2017)
set.seed(3)
pvclust_sdeffect2017 <- pvclust(t_cluster_sdeffect2017[,1:8], method.hclust = "average", method.dist = "correlation", nboot = 10000, parallel = FALSE, iseed = NULL)
clusterplot_sdeffect2017 <- plot(pvclust_sdeffect2017, hang = -1, cex = 0.5,  print.pv = "au")

pvrect(pvclust_sdeffect2017, alpha = 0.80, pv = "au", max.only = TRUE)
finalcluster_sdeffect2017 <- pvpick(pvclust_sdeffect2017, alpha = 0.80, pv = "au", max.only = TRUE)
finalcluster_sdeffect2017

#Heat map generation
#2017
dend_effect_2017 <- as.dendrogram(pvclust_sdeffect2017)
colour_dend_effect_2017 <- color_branches(dend_effect_2017, k = 2)

heatmap_effect2017 <- as.matrix(t(cluster_sdeffect2017))
names(heatmap_effect2017)
View(heatmap_effect2017)
list(unique(sigdiff_effect2017$Measurement))

row_ha = rowAnnotation(`Effect Type` = c("Xenobiotic Response-Related", "Xenobiotic Response-Related","Xenobiotic Response-Related",
                                       "Xenobiotic Response-Related", "Xenobiotic Response-Related", "Xenobiotic Response-Related",
                                       "Xenobiotic Response-Related", "Xenobiotic Response-Related",
                                       "Xenobiotic Response-Related", "Xenobiotic Response-Related",
                                       "Xenobiotic Response-Related", "Xenobiotic Response-Related",
                                       "Xenobiotic Response-Related", "Xenobiotic Response-Related",
                                       "Endocrine-Related", "Endocrine-Related",
                                       "Endocrine-Related", "Endocrine-Related",
                                       "Endocrine-Related", "Endocrine-Related"), col = list(`Effect Type` = c("Xenobiotic Response-Related" = "black", "Endocrine-Related" = "darkred")))
View(t(cluster_sdeffect2017))
column_ha = HeatmapAnnotation(cluster = c("1", "2", "2", "1", "2", "2", "2", "1"), col = list(cluster = c("1" = "azure2", "2" = "azure3")))
htef2017 <- Heatmap(heatmap_effect2017, col = c("white", "gold", "darkblue"),
        cluster_columns = dend_effect_2017, cluster_rows = FALSE, 
        row_title = "Measured Effects", column_title = "Site",
        column_title_side = "bottom", row_title_side = "right", cluster_row_slices = FALSE,
        width = unit(15, "cm"), height = unit(7, "cm"), row_names_gp = grid::gpar(fontsize = 7), show_heatmap_legend = TRUE,
        name = " ", top_annotation = column_ha, left_annotation = row_ha)

draw(htef2017, heatmap_legend_side = "left")


#2018####
#only significant####
list(unique(effect_data1$Site))
effect2018  <- effect_data1 %>% filter(Year == "2018", !Site %in% c("BK", "MQ-BK", "Blank"))
View(effect2018)
sigdiff_effect2018 <- effect2018 %>% filter(Measurement %in% c("E2-EQ", "Adrenoceptor Beta 1 (ADRB1)", "Aryl hydrocarbon Receptor (AhR)",
                                                               "Estrogen Receptor alpha (ERa)", "Estrogen Response Element (ERE)", "Glucocorticoid Receptor (GR)",
                                                               "Melanocortin 1 Receptor (MC1R)", "Antioxidant Response Element (ARE)-binding Nuclear factor (erythroid-derived 2)-like 2 (NRF2) (NRF2/ARE)",
                                                               "Peroxisome proliferator-activated receptor-a  (PPARa)", "Peroxisome proliferator-activated receptor-gamma (PPARg)",
                                                               "Peroxisome proliferator activating receptor (PPRE)", "Prostaglandin D2 Receptor (PTGDR)",
                                                               "Prostaglandin E Receptor 2 (PTGER2)", "Prostaglandin I2 Receptor (PTGIR)", 
                                                               "Pregnane-X-Receptor (PXR)", 
                                                               "Retinoid X receptor-b (RXRb)", "CYP1A1_liver"))

sigdiff_effect2018$Measurement <- ifelse(sigdiff_effect2018$Measurement == "17b_estradiol" & sigdiff_effect2018$`Sample Type` == "Male Fathead Minnow Serum", "E2_M_plasma", sigdiff_effect2018$Measurement)
sigdiff_effect2018$Measurement <- ifelse(sigdiff_effect2018$Measurement == "17b_estradiol" & sigdiff_effect2018$`Sample Type` == "Female Fathead Minnow Serum", "E2_F_plasma", sigdiff_effect2018$Measurement)
sigdiff_effect2018$Measurement <- ifelse(sigdiff_effect2018$Measurement == "Testosterone" & sigdiff_effect2018$`Sample Type` == "Male Fathead Minnow Serum", "T_M_plasma", sigdiff_effect2018$Measurement)
sigdiff_effect2018$Measurement <- ifelse(sigdiff_effect2018$Measurement == "Aryl hydrocarbon Receptor (AhR)", "AhR_bioactivity", sigdiff_effect2018$Measurement)
sigdiff_effect2018$Measurement <- ifelse(sigdiff_effect2018$Measurement =="Estrogen Receptor alpha (ERa)", "ERa_bioactivity", sigdiff_effect2018$Measurement)
sigdiff_effect2018$Measurement <- ifelse(sigdiff_effect2018$Measurement =="Estrogen Response Element (ERE)", "ERE_bioactivity", sigdiff_effect2018$Measurement)
sigdiff_effect2018$Measurement <- ifelse(sigdiff_effect2018$Measurement =="Antioxidant Response Element (ARE)-binding Nuclear factor (erythroid-derived 2)-like 2 (NRF2) (NRF2/ARE)", "ARE_bioactivity", sigdiff_effect2018$Measurement)
sigdiff_effect2018$Measurement <- ifelse(sigdiff_effect2018$Measurement =="Peroxisome proliferator-activated receptor-gamma (PPARg)", "PPARg_bioactivity", sigdiff_effect2018$Measurement)
sigdiff_effect2018$Measurement <- ifelse(sigdiff_effect2018$Measurement =="Peroxisome proliferator-activated receptor-a  (PPARa)", "PPARa_bioactivity", sigdiff_effect2018$Measurement)
sigdiff_effect2018$Measurement <- ifelse(sigdiff_effect2018$Measurement =="Peroxisome proliferator activating receptor (PPRE)", "PPRE_bioactivity", sigdiff_effect2018$Measurement)
sigdiff_effect2018$Measurement <- ifelse(sigdiff_effect2018$Measurement =="Pregnane-X-Receptor (PXR)" & sigdiff_effect2018$`Assay Name` == "Attagene cis-Factorial", "PXRcis_bioactivity", sigdiff_effect2018$Measurement)
sigdiff_effect2018$Measurement <- ifelse(sigdiff_effect2018$Measurement =="Pregnane-X-Receptor (PXR)" & sigdiff_effect2018$`Assay Name` == "Attagene trans-Factorial", "PXRtrans_bioactivity", sigdiff_effect2018$Measurement)
sigdiff_effect2018$Measurement <- ifelse(sigdiff_effect2018$Measurement =="Adrenoceptor Beta 1 (ADRB1)", "ADRB1_bioactivity", sigdiff_effect2018$Measurement)
sigdiff_effect2018$Measurement <- ifelse(sigdiff_effect2018$Measurement =="Glucocorticoid Receptor (GR)", "GR_bioactivity", sigdiff_effect2018$Measurement)
sigdiff_effect2018$Measurement <- ifelse(sigdiff_effect2018$Measurement =="Melanocortin 1 Receptor (MC1R)", "MC1R_bioactivity", sigdiff_effect2018$Measurement)
sigdiff_effect2018$Measurement <- ifelse(sigdiff_effect2018$Measurement =="Prostaglandin D2 Receptor (PTGDR)", "PTGDR_bioactivity", sigdiff_effect2018$Measurement)
sigdiff_effect2018$Measurement <- ifelse(sigdiff_effect2018$Measurement =="Prostaglandin E Receptor 2 (PTGER2)", "PTGER2_bioactivity", sigdiff_effect2018$Measurement)
sigdiff_effect2018$Measurement <- ifelse(sigdiff_effect2018$Measurement =="Prostaglandin I2 Receptor (PTGIR)", "PTGIR_bioactivity", sigdiff_effect2018$Measurement)
sigdiff_effect2018$Measurement <- ifelse(sigdiff_effect2018$Measurement =="Retinoid X receptor-b (RXRb)", "RXRb_bioactivity", sigdiff_effect2018$Measurement)

list(unique(sigdiff_effect2018$Measurement))


mean_sdeffect2018 <- sigdiff_effect2018 %>% group_by(Site, Measurement) %>% summarize(meanEffect = mean(Effect)) %>% filter(Measurement != "PXRtrans_bioactivity")
mean_sdeffect2018$Measurement <- factor(mean_sdeffect2018$Measurement, levels = c("AhR_bioactivity", "PXRcis_bioactivity", "ADRB1_bioactivity",
                                                                                  "GR_bioactivity", "MC1R_bioactivity","ARE_bioactivity", "PPARa_bioactivity",
                                                                                  "PPARg_bioactivity", "PPRE_bioactivity", "PTGDR_bioactivity", "PTGER2_bioactivity",
                                                                                  "PTGIR_bioactivity", "RXRb_bioactivity", "CYP1A1_liver", 
                                                                                  "E2-EQ", "ERa_bioactivity", "ERE_bioactivity"))

wide_sdeffect2018 <- mean_sdeffect2018 %>% select("Site", "Measurement", "meanEffect") %>% spread(key = "Measurement", value = "meanEffect") %>%
  replace(is.na(.), 0) %>% column_to_rownames("Site") 
View(wide_sdeffect2018)

cluster_sdeffect2018 <- scale(wide_sdeffect2018, center = TRUE, scale = TRUE)
t_cluster_sdeffect2018 <- t(cluster_sdeffect2018)

set.seed(3)
pvclust_sdeffect2018 <- pvclust(t_cluster_sdeffect2018[,1:12], method.hclust = "average", method.dist = "correlation", nboot = 10000, parallel = FALSE, iseed = NULL)
clusterplot_sdeffect2018 <- plot(pvclust_sdeffect2018, hang = -1, cex = 0.5, print.pv = "au")

pvrect(pvclust_sdeffect2018, alpha = 0.80, pv = "au", max.only = TRUE)
finalcluster_sdeffect2018 <- pvpick(pvclust_sdeffect2018, alpha = 0.80, pv = "au", max.only = TRUE)
finalcluster_sdeffect2018

#2018

row_ha = rowAnnotation(`Effect Type` = c("Xenobiotic Response-Related", "Xenobiotic Response-Related", "Xenobiotic Response-Related",
                                       "Xenobiotic Response-Related","Xenobiotic Response-Related", "Xenobiotic Response-Related",
                                       "Xenobiotic Response-Related", "Xenobiotic Response-Related", "Xenobiotic Response-Related",
                                       "Xenobiotic Response-Related", "Xenobiotic Response-Related", "Xenobiotic Response-Related",
                                       "Xenobiotic Response-Related", "Xenobiotic Response-Related",
                                       "Endocrine-Related", "Endocrine-Related", "Endocrine-Related"), col = list(`Effect Type` = c("Xenobiotic Response-Related" = "black", "Endocrine-Related" = "darkred")))

View(heatmap_effect2018)
column_ha = HeatmapAnnotation(cluster = c("2", "NA", "4", "4", "3", "2", "4", "1", "1", "3", "NA", "4"), col = list(cluster = c("1" = "azure2", "2" = "azure3","3" = "azure4", "4" = "black", "NA" = "white")))

dend_effect_2018 <- as.dendrogram(pvclust_sdeffect2018)

heatmap_effect2018 <- as.matrix(t(cluster_sdeffect2018))
View(heatmap_effect2018)

htef2018 <- Heatmap(heatmap_effect2018, col = c("white", "gold", "darkblue"),
                    cluster_columns = dend_effect_2018, cluster_rows = FALSE, 
                    row_title = "Measured Effects", column_title = "Site",
                    column_title_side = "bottom", row_title_side = "right", cluster_row_slices = FALSE,
                    width = unit(15, "cm"), height = unit(7, "cm"), row_names_gp = grid::gpar(fontsize = 8), show_heatmap_legend = TRUE, 
                    name = " ", left_annotation = row_ha, top_annotation = column_ha)

draw(htef2018, heatmap_legend_side = "left")


