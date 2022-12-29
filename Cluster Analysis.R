

## cluster analysis looking at standardizing chemical concentrations across sites#####

library(xlsx)
library(dplyr)
library(ggplot2)
library(arsenal)
library(readxl)
library(writexl)
library(tidyr)
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

setwd("C:\\Users\\erinm\\OneDrive\\Desktop\\GLRI Project\\Milwauke\\Data for Papier\\Mixtures MS\\Chemistry Evaluation\\Chemistry Cluster Analysis")
chem.2017 <- read_excel("C:\\Users\\erinm\\OneDrive\\Desktop\\GLRI Project\\Milwauke\\Data for Papier\\Chemistry\\Compiled_Chemistry\\Chemistry Compiled_Final.xlsx", "2017") %>% filter(Media == "Water", !Site %in% c("MIM-BK", "UWM"))
names(chem.2017)
list(unique(chem.2017$Site))
list(unique(chem.2017$Detect.code))
list(unique(chem.2017$Chemical))

#WQ####
list(unique(chem.2017$Site))

chemistry_2017_1 <- chem.2017 %>% select(Site, Media, Year, Chemical, CAS, Detect.code, detect.con.ppb) %>% filter(Media == "Water", is.na(Detect.code) | Detect.code %in% c("E", "M"), !is.na(detect.con.ppb), !Site %in% c("MIM-BK", "UWM")) %>% rename("Concentration" = "detect.con.ppb")
chemistry_2017_1$Concentration <- as.numeric(chemistry_2017_1$Concentration)

Spec_class <- read_excel("C:\\Users\\erinm\\OneDrive\\Desktop\\GLRI Project\\Milwauke\\Data for Papier\\Mixtures MS\\Chemistry Evaluation\\Chemistry Cluster Analysis\\Chemistry_data.xlsx") %>% select(CAS, Class)
names(Spec_class)

chemistry_2017_2 <- left_join(chemistry_2017_1, Spec_class)
View(chemistry_2017_2)

#2017 ####
chemistry_2017 <- chemistry_2017_2 %>% filter(!Detect.code %in% c("<"), Chemical != "3-tert-Butyl-4-hydroxyanisole")
chemistry_2017 <-chemistry_2017 %>% arrange((Class)) %>% distinct()

#use below code for descriptive summary only - otherwise disregard
# chemistry_2017$Spec_Class <- as.factor(chemistry_2017$Spec_Class)
# summary((chemistry_2017 %>% select(CAS, Spec_Class) %>% distinct)$Spec_Class)

#cluster evaluation using average + correlation####
widechemistry_2017 <- chemistry_2017 %>%  select(Site, CAS, Concentration) %>% spread(CAS, Concentration) %>% replace(is.na(.), 0) %>% column_to_rownames("Site")

list_CAS_names <- list(names(widechemistry_2017))

cluster_chemistry2017<- scale(widechemistry_2017, center = TRUE, scale = TRUE)

tcluster_chemistry2017 <- t(cluster_chemistry2017)
View(tcluster_chemistry2017)

set.seed(3)
pvclust_chemistry2017 <- pvclust(tcluster_chemistry2017[,1:9], method.hclust = "average", method.dist = "corr", nboot = 10000, parallel = FALSE, iseed = NULL)
pvclust_chemistry2017_plot <- plot(pvclust_chemistry2017, hang = -1, cex = 0.5, print.pv = "au")

pvrect(pvclust_chemistry2017, alpha = 0.80, pv = "au", max.only = TRUE) #for now do 0.75

chemistry2017_clusters <- pvpick(pvclust_chemistry2017, alpha = 0.80, pv = "au", max.only = TRUE)
chemistry2017_clusters
beepr::beep()

#2018####
chemistry.2018 <- read_excel("C:\\Users\\erinm\\OneDrive\\Desktop\\GLRI Project\\Milwauke\\Data for Papier\\Chemistry\\Compiled_Chemistry\\Chemistry Compiled_Final.xlsx", "2018") %>% filter(Media == "Water", !Site %in% c("MIM-BK", "UWM"))
names(chemistry.2018)
list(unique(chemistry.2018$Site))
list(unique(chemistry.2018$Detect.code))
list(unique(chemistry.2018$Chemical))
View(chemistry.2018)

#WQ####
list(unique(chemistry.2018$Site))

chemistry_2018_1 <- chemistry.2018 %>% select(Site, Media, Year, Chemical, CAS, Detect.code, detect.con.ppb) %>% filter(Media == "Water", is.na(Detect.code) | Detect.code %in% c("E", "M"), !is.na(detect.con.ppb), !Site %in% c("UWM")) %>% rename("Concentration" = "detect.con.ppb")
chemistry_2018_1$Concentration <- as.numeric(chemistry_2018_1$Concentration)

chemistry.2018_MED <- chemistry.2018 %>% filter(Site == "MED") %>% select(Chemical, CAS, detect.con.ppb) %>%
  rename("Concentration" = "detect.con.ppb")
chemistry.2018_MED$Site <- "MED"
View(chemistry.2018_MED)

chemistry_2018_CAS_list <- chemistry_2018 %>% select(CAS, Chemical)
chemistry.2018_MED_1 <- left_join(chemistry_2018_CAS_list, chemistry.2018_MED) %>% distinct()
View(chemistry.2018_MED_1)

chemistry_2018_plusMED <- bind_rows(chemistry_2018_1, chemistry.2018_MED_1)
View(chemistry_2018_plusMED)  


chemistry_2018_2 <- left_join(chemistry_2018_plusMED, Spec_class)

chemistry_2018 <- chemistry_2018_2

chemistry_2018 <-chemistry_2018 %>% arrange((Class))%>% distinct()
View(chemistry_2018)

#use below code for descriptive summary only - otherwise disregard
chemistry_2018$Spec_Class <- as.factor(chemistry_2018$Spec_Class)
summary((chemistry_2018 %>% select(CAS, Spec_Class) %>% distinct)$Spec_Class)
View(chemistry_2018)

#continue###
View(widechemistry_2018)
widechemistry_2018<- chemistry_2018 %>% select(Site, CAS, Concentration) %>% distinct() %>%
  spread(CAS, Concentration) %>% column_to_rownames("Site") 

widechemistry_2018_1 <- widechemistry_2018 %>% replace(is.na(.), 0)
View(widechemistry_2018_1)
cluster_chemistry2018 <- scale(widechemistry_2018_1, center = TRUE, scale = TRUE)
View(cluster_chemistry2018)
list_CAS_names_2018 <- list(names(widechemistry_2018))

tcluster_chemistry2018 <- t(cluster_chemistry2018)
View(tcluster_chemistry2018)

set.seed(3)
?pvclust
pvclust_chemistry2018 <- pvclust(tcluster_chemistry2018[,1:12], method.hclust = "average", method.dist = "correlation", nboot = 10000, parallel = FALSE, iseed = NULL)
pvclust_chemistry2018_plot <- plot(pvclust_chemistry2018, hang = -1, cex = 0.5, print.pv = "au")

pvrect(pvclust_chemistry2018, alpha = 0.80, pv = "au", max.only = TRUE)
chemistry2018_clusters <- pvpick(pvclust_chemistry2018, alpha = 0.80, pv = "au", max.only = TRUE)
chemistry2018_clusters

#generate heat maps####
#2017 heatmap####
View(chemistry_2017)
dend_2017 <- as.dendrogram(pvclust_chemistry2017)
colour_dend_2017 <- color_branches(dend_2017, k = 2)

heatmap_chemistry2017 <- as.matrix(t(cluster_chemistry2017))
View(heatmap_chemistry2017)

list_CAS_names <- list(row.names(heatmap_chemistry2017))
list_CAS_names <- as.data.frame(list_CAS_names)
View(list_CAS_names)
names(list_CAS_names)

list_CAS_names <- list_CAS_names %>% rename("CAS" = "c..102.36.3....103.90.2....115.96.8....120.72.9....126.73.8...")

list(unique(chemistry_2017$Class))
class <- chemistry_2017 %>% select(CAS, Class) %>% distinct
class_2 <- left_join(list_CAS_names, class) %>% distinct()
class_3 <- class_2 %>% select(Class)
class_4 <- as.matrix(class_3)
View(class_4)
list(unique(class_3$Class))

row_ha = rowAnnotation(chemical_class = class_4, col = list(chemical_class = c("Industrial Chemicals" = "#0d0887",
                                                                               "Pesticides" = "#f0f921",
                                                                               "Fire Retardants" = "#fca636",
                                                                               "Pharmaceuticals and Personal Care Products" = "#e16462",
                                                                               "Fuels & PAHs" = "#b12a90",
                                                                               "Wastewater Indicators" = "#6a00a8")))
  
list(unique(chemistry_2017$CAS))
View(heatmap_chemistry2017)
column_ha = HeatmapAnnotation(cluster = c("NA", "2", "1", "1", "2", "1", "1", "1", "2"), col = list(cluster = c("1" = "darkgrey", "2" = "lightgrey","NA" = "white")))
htchem2017 <- Heatmap(heatmap_chemistry2017, col = c("white","gold", "darkblue"),
                      row_order = c("115-96-8", "126-73-8", "78-51-3",
                         "129-00-0","206-44-0", "50-32-8", "85-01-8", "90-12-0",
                         "136-85-6", "29385-43-1", "78-59-1", "80-05-7", "84-65-1", "86-74-8", "2315-67-5","120-72-9", 
                         "102-36-3", "134-62-3", "1912-24-9", "51218-45-2", "63-25-2",
                         "103-90-2", "137-58-6", "27203-92-5", "298-46-4", "396-01-0", "54-11-5", "57-53-4", "58-08-2", "59-05-2", "657-24-9", "723-46-6", "738-70-5", "83799-24-0", "93413-69-5","93413-62-8",
                         "486-56-6", "57-88-5", "83-46-5"),
        cluster_columns = dend_2017, cluster_rows = FALSE, row_title = "Detected Chemical", column_title = "Site",
        column_title_side = "bottom", row_title_side = "right", cluster_row_slices = FALSE,
        width = unit(10, "cm"), height = unit(13.5, "cm"), row_names_gp = grid::gpar(fontsize = 8), show_heatmap_legend = TRUE,
        name = " ", left_annotation = row_ha, top_annotation = column_ha)

draw(htchem2017, heatmap_legend_side = "left")

#2018 heatmap####
dend_2018 <- as.dendrogram(pvclust_chemistry2018)

heatmap_chemistry2018 <- as.matrix(tcluster_chemistry2018) %>% arrange()
View(heatmap_chemistry2018)

list_CAS_names <- list(row.names(heatmap_chemistry2018))
View(list_CAS_names)
list_CAS_names <- as.data.frame(list_CAS_names)
names(list_CAS_names)

CAS_names_2018 <- as.data.frame(list_CAS_names_2018) %>% rename("CAS" ="c..103.90.2....106.44.5....106.46.7....115.86.6....115.96.8...")
names(CAS_names_2018)
class_2018 <- chemistry_2018 %>% select(CAS, Class) %>% distinct() %>% filter(Class != "WQ") %>% filter(CAS != "5989-27-5")
D_lim <- class_2018 %>% filter(CAS == "5989-27-5", Class == "Industrial Chemicals")
class_2018 <- bind_rows(class_2018, D_lim)

class_2_2018 <- left_join(CAS_names_2018, class_2018) %>% distinct()
class_3_2018 <- class_2_2018 %>% select(Class)
class_3_2018$Class[is.na(class_3_2018$Class)] <- "Industrial Chemicals"
class_4_2018 <- as.matrix(class_3_2018)
View(class_3_2018)

row_ha = rowAnnotation(chemical_class = class_4_2018, col = list(chemical_class = c("Industrial Chemicals" = "#0d0887",
                                                                               "Pesticides" = "#f0f921",
                                                                               "Fire Retardants" = "#fca636",
                                                                               "Pharmaceuticals and Personal Care Products" = "#e16462",
                                                                               "Fuels & PAHs" = "#b12a90",
                                                                               "Wastewater Indicators" = "#6a00a8")))


heatmap_chemistry2018_1 <- as.data.frame(heatmap_chemistry2018)
names(heatmap_chemistry2018_1)
column_ha = HeatmapAnnotation(cluster = c("3", "2", "1", "1", "2", "2", "1", "3", "3", "1", "3", "2"), col = list(cluster = c("1" = "darkgrey", "2" = "lightgrey", "3" = "black", "NA" = "white")))

list(row.names(heatmap_chemistry2018))

htchem2018 <- Heatmap(heatmap_chemistry2018, col = c("white", "gold", "darkblue"),
        row_order = c("115-86-6", "126-73-8", "115-96-8", "13674-87-8", 
                      "120-12-7", "129-00-0", "206-44-0", "50-32-8", "85-01-8", "90-12-0",  "91-57-6",
                      "106-44-5",  "119-65-3","136-85-6", "140-66-9", "29385-43-1","75-25-2", "77-93-0", "84-65-1", "86-74-8","120-72-9", "1222-05-5","76-22-2","5989-27-5",
                      "106-46-7","134-62-3","1912-24-9", "51-03-6","51218-45-2","63-25-2", "87-86-5",
                      "103-90-2", "119-36-8", "122-11-2", "125-71-3", "137-58-6",
                      "148-79-8","27203-92-5","29122-68-7","298-46-4","3380-34-5", "34911-55-2",
                      "396-01-0",  "42399-41-7","486460-32-6", "51384-51-1","525-66-6",
                      "532-03-6","54-11-5","58-08-2","58-73-1","59277-89-3", "59729-33-8",
                      '657-24-9', "66357-35-5", "723-46-6", "738-70-5","486-56-6",
                      "76-42-6",  "76-99-3", "76824-35-6", "83799-24-0","89-78-1",  "90-82-4/299-42-3", "93413-69-5","93413-62-8",
                      "360-68-9", "83-46-5", "19466-47-8", "57-88-5"),
        cluster_columns = dend_2018, cluster_rows = FALSE, row_title = "Detected Chemical", column_title = "Site",
        column_title_side = "bottom", row_title_side = "right", cluster_row_slices = FALSE,
        width = unit(10, "cm"), height = unit(13.5, "cm"), row_names_gp = grid::gpar(fontsize = 8), show_heatmap_legend = TRUE,
        name = " ", left_annotation = row_ha, top_annotation = column_ha)

draw(htchem2018, heatmap_legend_side = "left")


#PCA on chemistry data? Ignore for now ####
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

#2017
View(chemistry_2017)
names(chemistry_2017)
spec_class_2017 <- chemistry_2017 %>% group_by(Spec_Class, Site) %>% summarize(sum_Con = sum(Concentration))

spec_class_2017_wide <- spec_class_2017 %>% spread(key = Spec_Class, value = sum_Con)
spec_class_2017_wide[is.na(spec_class_2017_wide)] <- 0

names(spec_class_2017_wide)

pca.chem_2017 <- prcomp(spec_class_2017_wide[,2:7], center = TRUE, scale = TRUE)
View(pca.chem_2017)
summary(pca.chem_2017)

#graph
#PCA 1 + 2 explains 70% of variance 
pca1_pca2_chem <- autoplot(pca.chem_2017, data = spec_class_2017_wide, loadings = TRUE, 
                      loadings.colour = 'red', loadings.label = TRUE, size = 2, 
                      loadings.label.size = 5,loadings.label.repel = T, legend = F) + theme_classic() + geom_hline(yintercept = 0, linetype = "dotted") + geom_vline(xintercept = 0, linetype = "dotted") +
  theme(legend.title = element_text(size = 14), legend.text = element_text(size = 12)) + labs(title = "C") + theme(plot.title = element_text(size = 30)) +
  geom_text_repel(aes(label = Site), nudge_y = 0.01, nudge_x = 0.01, size = 5) + xlim(-0.75, 0.75) + ylim(-0.75, 0.75)+ theme(text = element_text(size = 15))


#other evals####
#screeplot indicates that 1st 3 components explain the most variance

screeplot(pca.chem_2017, type= "lines", npcs = 14, main = "Screeplot Effect Data 2017")
abline(h = 1, col = "red", lty = 5)

cumpro <- cumsum(pca.chem_2017$sdev^2 / sum(pca.chem_2017$sdev^2))
plot(cumpro[0:7], xlab = "PC #", ylab = "Amount of explained variance", main = "Cumulative variance plot") 
abline(v = 3, col = "blue", lty = 5)
abline(h = 0.89, col = "blue", lty = 5)

sd <- pca.chem_2017$sdev
loadings <- pca.chem_2017$rotation
rownames(loadings) <- colnames(spec_class_2017_wide[2:7])
scores <- pca.chem_2017$x
View(scores)

#determine optimal # of principle components
var <- sd^2
varPercent <- var/sum(var) * 100

barplot(varPercent, xlab='PC', ylab='Percent Variance', names.arg=1:length(varPercent), las=1, ylim=c(0, max(varPercent)), col='gray')
abline(h = 1/ncol(spec_class_2017_wide)*100, col = 'red')

#because 1st 3 components explain more than 1 variable's worth of information (as indicated by red line) we should include all of these in our analysis#
varPercent[1:3]
sum(varPercent[1:3]) #explain 88.89603 % of variance

#sum of squares of a loading for an individual principal component must sum to 1, so we can calculate loadings as if all variables contributed to component equally
#any variable with larger loading than equal value is contributes more than one variable's worth of information and is an important contributor to that principal component
View(loadings)
sqrt(1/ncol(spec_class_2017_wide[,2:7])) #cutoff for important loadings = 0.4082483

loadings_df <- as.data.frame(loadings)
loadings_df1 <- loadings_df %>% select((PC1:PC2))
loadings_df2 <- loadings_df1
loadings_df1$PC1 <- ifelse(abs(loadings_df1$PC1) < 0.4082483, "not important", loadings_df1$PC1)
loadings_df1$PC2 <- ifelse(abs(loadings_df1$PC2) < 0.4082483, "not important", loadings_df1$PC2)

View(loadings_df1)
write.xlsx(loadings_df2, "chemistry_loadings_2017.xlsx", row.names = TRUE)
write.xlsx(loadings_df1, "chemistry_loadings_2017_importance_annotation.xlsx", row.names = TRUE)

#2018
View(chemistry_2018)
names(chemistry_2018)
spec_class_2018 <- chemistry_2018 %>% group_by(Spec_Class, Site) %>% summarize(sum_Con = sum(Concentration))

spec_class_2018_wide <- spec_class_2018 %>% spread(key = Spec_Class, value = sum_Con)
spec_class_2018_wide[is.na(spec_class_2018_wide)] <- 0

View(spec_class_2018_wide)

pca.chem_2018 <- prcomp(spec_class_2018_wide[,2:7], center = TRUE, scale = TRUE)
View(pca.chem_2018)
summary(pca.chem_2018) #81.6

#graph
#PCA 1 + 2 explains 81.6% of variance 
pca1_pca2_chem_2018 <- autoplot(pca.chem_2018, data = spec_class_2018_wide, loadings = TRUE, 
                           loadings.colour = 'red', loadings.label = TRUE, size = 2, 
                           loadings.label.size = 5,loadings.label.repel = T, legend = F) + theme_classic() + geom_hline(yintercept = 0, linetype = "dotted") + geom_vline(xintercept = 0, linetype = "dotted") +
  theme(legend.title = element_text(size = 14), legend.text = element_text(size = 12)) + labs(title = "D") + theme(plot.title = element_text(size = 30)) +
  geom_text_repel(aes(label = Site), nudge_y = 0.01, nudge_x = 0.01, size = 5) + xlim(-0.75, 0.75) + ylim(-0.75, 0.75)+ theme(text = element_text(size = 15))

pca1_pca2_chem_2018
#other evals####
#screeplot indicates that 1st 3 components explain the most variance

screeplot(pca.chem_2018, type= "lines", npcs = 14, main = "Screeplot Effect Data 2017")
abline(h = 1, col = "red", lty = 5)

cumpro <- cumsum(pca.chem_2018$sdev^2 / sum(pca.chem_2017$sdev^2))
plot(cumpro[0:7], xlab = "PC #", ylab = "Amount of explained variance", main = "Cumulative variance plot") 
abline(v = 3, col = "blue", lty = 5)
abline(h = 0.89, col = "blue", lty = 5)

sd <- pca.chem_2018$sdev
loadings <- pca.chem_2018$rotation
rownames(loadings) <- colnames(spec_class_2018_wide[2:7])
scores <- pca.chem_2018$x
View(scores)

#determine optimal # of principle components
var <- sd^2
varPercent <- var/sum(var) * 100

barplot(varPercent, xlab='PC', ylab='Percent Variance', names.arg=1:length(varPercent), las=1, ylim=c(0, max(varPercent)), col='gray')
abline(h = 1/ncol(spec_class_2017_wide)*100, col = 'red')

#because 1st 3 components explain more than 1 variable's worth of information (as indicated by red line) we should include all of these in our analysis#
varPercent[1:3]
sum(varPercent[1:3]) #explain 94.89653 % of variance

#sum of squares of a loading for an individual principal component must sum to 1, so we can calculate loadings as if all variables contributed to component equally
#any variable with larger loading than equal value is contributes more than one variable's worth of information and is an important contributor to that principal component
View(loadings)
sqrt(1/ncol(spec_class_2018_wide[,2:7])) #cutoff for important loadings = 0.4082483

loadings_df <- as.data.frame(loadings)
loadings_df1 <- loadings_df %>% select((PC1:PC2))
loadings_df2 <- loadings_df1
loadings_df1$PC1 <- ifelse(abs(loadings_df1$PC1) < 0.4082483, "not important", loadings_df1$PC1)
loadings_df1$PC2 <- ifelse(abs(loadings_df1$PC2) < 0.4082483, "not important", loadings_df1$PC2)

View(loadings_df1)
write.xlsx(loadings_df2, "chemistry_loadings_2018.xlsx", row.names = TRUE)
write.xlsx(loadings_df1, "chemistry_loadings_2018_importance_annotation.xlsx", row.names = TRUE)


