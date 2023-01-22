#RF Wrangle

library(tidyverse)
library(superheat)
library(ggpubr)
library(readxl)
library(writexl)
library(viridisLite)

#read in RF_output

RF_output <- read_excel("C:\\Users\\erinm\\OneDrive\\Desktop\\GLRI Project\\Milwauke\\Data for Papier\\Mixtures MS\\Random_Forest\\Random_Forest\\RF_output.xlsx")
names(RF_output)

RF_effect_predictors <- RF_output %>% select(Var_list, Effect_Group, "Effect Type", Effect, "Year(s)") %>%
  rename("Predictor" = "Var_list")

RF_effect_predictors$Effect_Year <- paste(RF_effect_predictors$Effect, RF_effect_predictors$`Year(s)`)
RF_effect_predictors$Effect_Group_Type <- paste(RF_effect_predictors$Effect_Group, RF_effect_predictors$`Effect Type`)
list(unique(RF_effect_predictors$Effect_Group_Type))
RF_effect_predictors$Effect_Group_Type <- gsub("XM_related in_vitro", 0, RF_effect_predictors$Effect_Group_Type)
RF_effect_predictors$Effect_Group_Type <- gsub("XM_related in_vivo", 0.5, RF_effect_predictors$Effect_Group_Type)
RF_effect_predictors$Effect_Group_Type <- gsub("ER_related in_vitro", 1, RF_effect_predictors$Effect_Group_Type)
RF_effect_predictors$Effect_Group_Type <- gsub("ER_related in_vivo", 1.5, RF_effect_predictors$Effect_Group_Type)

RF <- RF_effect_predictors %>% select(-"Effect Type", -"Year(s)", -Effect, -Effect_Group)
RF$Val <- 1

RF_1 <- RF %>% distinct()  %>% 
  spread(key = "Predictor", value = "Val") %>% arrange(desc(Effect_Group_Type)) %>% select(-Effect_Group_Type)%>%
  column_to_rownames("Effect_Year") %>% 
  relocate(c("Acetaminophen", "Anthraquinone", "beta-Sitosterol", "Bupropion",
             "Caffeine", "Carbamazepine", "Cholesterol", "Cotinine",
             "DEET", "Desvenlafaxine", "Diethyl phthalate", "Fexofenadine",
             "Fluoranthene", "Gabapentin", "Isophorone", "Lidocaine", "Metformin", "Methyl-1H-benzotriazole",
             "Metolachlor", "Phenanthrene", "Pyrene", "Sulfamethoxazole",
             "Tributyl phosphate", "Tris(1,3-dichloro-2-propyl) phosphate", "Venlafaxine"), .before =  "11-deoxycortisol_Inhibition") %>%
  relocate(c("Benzophenone_Bupropion_Anthraquinone", "Carbazole_Indole", "Generally_Reactive", "Narcotic_Amine_PPCPs",
             "Narcotic_PPCPs_CNS_Stimulants_Antivirals", "Narcotic_PPCPs_Pesticides_Plasticizers", "PAHs_Fuels_Antioxidants",
             "Sterols_Menthol"), .before = "11-deoxycortisol_Inhibition") %>% 
  relocate(c("ADORA_HRH_TACR-Inhibition", "GRIN-Inhibition", "OPR-Activation", "SCN-Inhibition",
             "SLC-Inhibition", "TRPV1-Activation"), .after = "XBP1_Activation")

names(RF_1)

RF_1[is.na(RF_1)] <- 0
t_RF_1 <- t(RF_1)
View(t_RF_1)

superheat(X = t_RF_1,
          grid.hline.col = "black",
          grid.vline.col = "black",
          grid.hline = FALSE,
          grid.vline = FALSE,
          col.dendrogram = F,  #clustering
          row.dendrogram = F,  #clustering
          bottom.label.text.size = 4, 
          left.label.text.size = 4,
          left.label.col = 'white',
          bottom.label.col = 'white',
          left.label.text.col = 'black',
          heat.pal =viridis(n = 2, direction = -1, option = "magma"),
          clustering.method = 'hierarchical',
          bottom.label.text.angle = 90,
          legend = FALSE)


#save at 1500 X 1800
