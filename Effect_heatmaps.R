#Effect Heat Maps
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


setwd("C:\\Users\\erinm\\OneDrive\\Desktop\\GLRI Project\\Milwauke\\Data for Papier\\Mixtures MS\\Effect Data")

effect_data <- read_excel("C:\\Users\\erinm\\OneDrive\\Desktop\\GLRI Project\\Milwauke\\Data for Papier\\Mixtures MS\\Effect Evaluation\\Effect Data\\Effect data - for R heatmap input.xlsx")
names(effect_data)

effect_data1 <- effect_data %>% select("Assay Name", "Assay Type", "Year", "Site", "Replicate", "Measurement", "Effect", "Units", "Sample Type") %>% filter(!Site %in% c("BK", "MQ-BK"))
View(effect_data1)


#split into study years####
effect2017  <- effect_data1 %>% filter(Year == "2017")
effect2018  <- effect_data1 %>% filter(Year == '2018')

View(effect2017)

#heatmaps
#2017
list(unique(effect2017$Measurement))
effect2017$Measurement <- ifelse(effect2017$Measurement == "17b_estradiol" & effect2017$`Sample Type` == "Male Fathead Minnow Serum",
                                  "17b_estradiol_M",
                                  ifelse(effect2017$Measurement == "17b_estradiol" & effect2017$`Sample Type` == "Female Fathead Minnow Serum",
                                         "17b_estradiol_F",
                                         effect2017$Measurement))

effect_2017 <- effect2017 %>% select(Site, Measurement, Effect) %>% filter(Measurement %in% c("E2-EQ", "CYP1A1_intestine", 
                                                                                                            "CYP1A1_liver", "CYP2AD6_intestine",
                                                                                                            "CYP2N13_intestine", "CYP3A_intestine",
                                                                                                            "CYP3A_liver", "UGT1A1_intestine",
                                                                                                            "UGT1A1_liver", "VTG",
                                                                                                            "17b_estradiol_M", "17b_estradiol_F", "Testosterone",
                                                                                                            "Body Weight", "Gonadosomatic Index",
                                                                                                            "Survival",
                                                                                                            "AhR", "PXR_cis","PXR_trans")) %>% 
  distinct() %>% group_by(Site, Measurement) %>% summarize(meanEffect = mean(Effect))

View(effect_2017)
effect_2017_wide <- effect_2017 %>% spread(key ="Measurement", value = "meanEffect")
names(effect_2017_wide)
effect_2017_wide[2:19] <- effect_2017_wide[2:19] %>% scale()

effect_2017_wide[is.na(effect_2017_wide)] <- -5

effect_2017_wide_1 <- effect_2017_wide  %>% relocate("17b_estradiol_M", .after = "VTG") %>% relocate("17b_estradiol_F", .after = "VTG") %>%
  relocate("Testosterone", .after = "17b_estradiol_M") %>% relocate("Body Weight", .before = "AhR") %>% relocate("Survival", .after = "Body Weight") %>% relocate("Gonadosomatic Index", .after = "Survival") %>%
  rename("GSI" = "Gonadosomatic Index")%>% relocate("PXR_cis", .after = "AhR") %>% relocate("PXR_trans", .after = "PXR_cis") %>%
  arrange(match(Site, c("MED", "MEF", "KKL", "MEC", "MET", "UCJ", "MIE", "MIM", "MIP"))) %>% column_to_rownames("Site") 
View(effect_2017_wide_1)

t_effect_2017_wide_1 <- t(effect_2017_wide_1)


superheat(X = t_effect_2017_wide_1,
          grid.hline.col = "white",
          grid.vline.col = "white",
          grid.hline = FALSE,
          grid.vline = FALSE,
          col.dendrogram = F,  #clustering
          row.dendrogram = F,  #clustering
          bottom.label.text.size = 4, 
          left.label.text.size = 4,
          left.label.col = 'white',
          bottom.label.col = 'white',
          left.label.text.col = 'black',
          heat.pal =c("#000004", "#57106e", "#bc3754", "#f98e09","#fcffa4"),
          clustering.method = 'hierarchical',
          bottom.label.text.angle = 0,
          legend = FALSE)

#copy to clipboard @ 850 X 610

#2018
list(unique(effect2018$Measurement))

effect_2018 <- effect2018 %>% select(Site, Measurement, Effect) %>% filter(Measurement %in% c("E2-EQ",
                                                                                              "CYP1A1_liver", "VTG",
                                                                                              "Body Weight", "Gonadosomatic Index",
                                                                                              "Survival",
                                                                                              "AhR", "PXR_cis","ERa", "ERE",
                                                                                              "PPARg" , "PTGER2")) %>% 
  distinct() %>% group_by(Site, Measurement) %>% summarize(meanEffect = mean(Effect))

effect_2018_wide <- effect_2018 %>% spread(key ="Measurement", value = "meanEffect")
names(effect_2018_wide)
effect_2018_wide[2:13] <- effect_2018_wide[2:13] %>% scale()

effect_2018_wide[is.na(effect_2018_wide)] <- -5

effect_2018_wide_1 <- effect_2018_wide  %>% relocate("ERa", .after = "VTG") %>%
  relocate("ERE", .after = "ERa") %>% relocate("Body Weight", .before = "AhR") %>% relocate("Survival", .after = "Body Weight") %>% relocate("Gonadosomatic Index", .after = "Survival") %>%
  rename("GSI" = "Gonadosomatic Index")%>% relocate("PXR_cis", .after = "AhR") %>% relocate("PTGER2", .after = "PXR_cis") %>% relocate("PPARg", .after = "PTGER2") %>%
  arrange(match(Site, c("MED", "MEF", "MIN", "CCM", "KKL", "MEC", "MET", "UCJ", "MIE", "MIM", "MIP", "JIP"))) %>% column_to_rownames("Site") 
View(effect_2018_wide_1)

t_effect_2018_wide_1 <- t(effect_2018_wide_1)


superheat(X = t_effect_2018_wide_1,
          grid.hline.col = "white",
          grid.vline.col = "white",
          grid.hline = FALSE,
          grid.vline = FALSE,
          col.dendrogram = F,  #clustering
          row.dendrogram = F,  #clustering
          bottom.label.text.size = 4, 
          left.label.text.size = 4,
          left.label.col = 'white',
          bottom.label.col = 'white',
          left.label.text.col = 'black',
          heat.pal =c("#000004", "#57106e", "#bc3754", "#f98e09","#fcffa4"),
          clustering.method = 'hierarchical',
          bottom.label.text.angle = 0,
          legend = FALSE)