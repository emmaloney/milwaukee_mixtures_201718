#PCA test
library(dplyr)
library(tidyr)
library(tidyverse)
library(readxl)
library(writexl)
library(factoextra)
library(FactoMineR)
library(corrplot)
library(xlsx)
library(ggfortify)
library(ggplot2)
library(ggrepel)
library(gridExtra)

setwd("C:\\Users\\erinm\\OneDrive\\Desktop\\GLRI Project\\Milwauke\\Data for Papier\\Mixtures MS\\Effect Evaluation\\Effect Data")

effectdata <- read_excel("C:\\Users\\erinm\\OneDrive\\Desktop\\GLRI Project\\Milwauke\\Data for Papier\\Effect Data\\Effect_Data_for_all_Analyses.xlsx") 
names(effectdata)

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
  replace(is.na(.), 0)
View(wide_sdeffect2017)


###PCA
pca.effect_2017 <- prcomp(wide_sdeffect2017[,2:21], center = TRUE, scale = TRUE)
View(pca.effect_2017)
summary(pca.effect_2017)


#graph
#PCA 1 + 2 explains 68% of variance 
pca1_pca2 <- autoplot((pca.effect_2017), data = (wide_sdeffect2017),loadings = TRUE, 
                           loadings.colour = "red", loadings.label = TRUE, size = 3,
                           loadings.label.size = 4,loadings.label.repel = T,legend = F) + theme_classic() + geom_hline(yintercept = 0, linetype = "dotted") + geom_vline(xintercept = 0, linetype = "dotted") +
  theme(legend.title = element_text(size = 14), legend.text = element_text(size = 12)) + labs(title = "C") + theme(plot.title = element_text(size = 30)) +
  geom_text_repel(aes(label = Site), nudge_y = 0.01, nudge_x = 0.01, size = 5) + theme(text = element_text(size = 15))


#other evals####
#screeplot indicates that 1st 3 components explain the most variance

screeplot(pca.effect_2017, type= "lines", npcs = 14, main = "Screeplot Effect Data 2017")
abline(h = 1, col = "red", lty = 5)

cumpro <- cumsum(pca.effect_2017$sdev^2 / sum(pca.effect_2017$sdev^2))
plot(cumpro[0:7], xlab = "PC #", ylab = "Amount of explained variance", main = "Cumulative variance plot") 
abline(v = 3, col = "blue", lty = 5)
abline(h = 0.89, col = "blue", lty = 5)

sd <- pca.effect_2017$sdev
loadings <- pca.effect_2017$rotation
rownames(loadings) <- colnames(wide_sdeffect2017[2:21])
scores <- pca.effect_2017$x
View(scores)

#determine optimal # of principle components
var <- sd^2
varPercent <- var/sum(var) * 100

barplot(varPercent, xlab='PC', ylab='Percent Variance', names.arg=1:length(varPercent), las=1, ylim=c(0, max(varPercent)), col='gray')
abline(h = 1/ncol(wide_sdeffect2017)*100, col = 'red')

#because 1st 3 components explain more than 1 variable's worth of information (as indicated by red line) we should include all of these in our analysis#

#sum of squares of a loading for an individual principal component must sum to 1, so we can calculate loadings as if all variables contributed to component equally
#any variable with larger loading than equal value is contributes more than one variable's worth of information and is an important contributor to that principal component
View(loadings)
sqrt(1/ncol(wide_sdeffect2017[,2:21])) #cutoff for important loadings = 0.2236068

loadings_df <- as.data.frame(loadings)
loadings_df1 <- loadings_df %>% select((PC1:PC2))
loadings_df2 <- loadings_df1
loadings_df1$PC1 <- ifelse(abs(loadings_df1$PC1) < 0.2236068, "not important", loadings_df1$PC1)
loadings_df1$PC2 <- ifelse(abs(loadings_df1$PC2) < 0.2236068, "not important", loadings_df1$PC2)

View(loadings_df1)
write.xlsx(loadings_df2, "effect_loadings_2017.xlsx", row.names = TRUE)
write.xlsx(loadings_df1, "effect_loadings_2017_importance_annotation.xlsx", row.names = TRUE)

#2018 effect data####
effect2018  <- effect_data1 %>% filter(Year == "2018", !Site %in% c("BK", "MQ-BK", "Blank"))

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
  replace(is.na(.), 0)

names(wide_sdeffect2018)
pca.effect_2018 <- prcomp(wide_sdeffect2018[,2:18], center = TRUE, scale = TRUE)
summary(pca.effect_2018)

#graph
#PCA 1 + 2 explains ~77% of variance 
pca1_pca2_2018 <- autoplot(pca.effect_2018, data = wide_sdeffect2018, loadings = TRUE, 
                           loadings.colour = 'red', loadings.label = TRUE, size = 3, 
                           loadings.label.size = 4,loadings.label.repel = T, legend = F) + theme_classic() + geom_hline(yintercept = 0, linetype = "dotted") + geom_vline(xintercept = 0, linetype = "dotted") +
  theme(legend.title = element_text(size = 14), legend.text = element_text(size = 12)) + labs(title = "") + theme(plot.title = element_text(size = 30)) +
  geom_text_repel(aes(label = Site), nudge_y = 0.01, nudge_x = 0.01, size = 5)+ theme(text = element_text(size = 15))

#other evals#####
#screeplot indicates that 1st 2 components explain the most variance

screeplot(pca.effect_2018, type= "lines", npcs = 7, main = "Screeplot Effect Data 2018")
abline(h = 1, col = "red", lty = 5)

cumpro <- cumsum(pca.effect_2018$sdev^2 / sum(pca.effect_2018$sdev^2))
plot(cumpro[0:7], xlab = "PC #", ylab = "Amount of explained variance", main = "Cumulative variance plot") 
abline(v = 2, col = "blue", lty = 5)
abline(h = 0.88, col = "blue", lty = 5)

sd <- pca.effect_2018$sdev
loadings <- pca.effect_2018$rotation
rownames(loadings) <- colnames(wide_sdeffect2018[,2:18])
scores <- pca.effect_2018$x
View(scores)

#determine optimal # of principle components
var <- sd^2
varPercent <- var/sum(var) * 100

barplot(varPercent, xlab='PC', ylab='Percent Variance', names.arg=1:length(varPercent), las=1, ylim=c(0, max(varPercent)), col='gray')
abline(h = 1/ncol(wide_sdeffect2017)*100, col = 'red')

#because 1st 2 components explain more than 1 variable's worth of information (as indicated by red line) we should include all of these in our analysis#
varPercent[1:2]
sum(varPercent[1:2]) #explain 88.28959 % of variance

#sum of squares of a loading for an individual principal component must sum to 1, so we can calculate loadings as if all variables contributed to component equally
#any variable with larger loading than equal value is contributes more than one variable's worth of information and is an important contributor to that principal component
View(loadings)
sqrt(1/ncol(wide_sdeffect2018[,2:18])) #cutoff for important loadings = 0.2425356

loadings_df <- as.data.frame(loadings)
loadings_df1 <- loadings_df %>% select(PC1:PC2)
loadings_df2 <- loadings_df1
loadings_df2$PC1 <- ifelse(abs(loadings_df2$PC1) < 0.2425356, "not important", loadings_df2$PC1)
loadings_df2$PC2 <- ifelse(abs(loadings_df2$PC2) < 0.2425356, "not important", loadings_df2$PC2)

View(loadings_df1)

write.xlsx(loadings_df1, "effect_loadings_2018.xlsx", row.names = TRUE)
write.xlsx(loadings_df2, "effect_loadings_2018_annotated.xlsx", row.names = TRUE)

#graph effect PCA 2017 + 2018####

g1 <- grid.arrange(pca1_pca2, pca1_pca2_2018, nrow = 1)
ggsave("PCA_r1.jpeg", plot = g1, height = 20, width = 20)

#limit graphs to 'important contributors to PC#

