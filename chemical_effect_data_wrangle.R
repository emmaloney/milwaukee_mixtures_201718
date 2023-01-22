#pre_RF data wrangle

library(dplyr)
library(tidyverse)
library(readxl)
library(writexl)

#grouping and binding of chemical data####
names(ECOTOX_groups)
ECOTOX_groups <- read_excel("C:\\Users\\erinm\\OneDrive\\Desktop\\GLRI Project\\Milwauke\\Data for Papier\\Mixtures MS\\Most_Recent_Draft\\Co-author & External Review\\SI_08_04_2022.xlsx", "Table S5", skip = 1) %>%
  select(CAS, "Chemical Name", "Final Mixture Group") %>%
  rename("Group" = "Final Mixture Group", "Chemical" = "Chemical Name") %>% filter(!is.na(CAS), Group != "NA - unclustered or < 2 co-occuring constituents") %>% select(-Chemical)
ECOTOX_groups$Grouping_type <- "ECOTOX"

Pharm_groups <- read_excel("C:\\Users\\erinm\\OneDrive\\Desktop\\GLRI Project\\Milwauke\\Data for Papier\\Mixtures MS\\Most_Recent_Draft\\Co-author & External Review\\SI_08_04_2022.xlsx", "Table S11", skip = 1)%>%
  select(CAS, Chemical, "Final Mixture Group for LoE Evaluation")%>% rename("Group" = "Final Mixture Group for LoE Evaluation") %>% filter(Group != "-")%>% select(-Chemical)
Pharm_groups$Grouping_type <- "Pharm"

ToxCast_groups <- read_excel("C:\\Users\\erinm\\OneDrive\\Desktop\\GLRI Project\\Milwauke\\Data for Papier\\Mixtures MS\\Most_Recent_Draft\\Co-author & External Review\\SI_08_04_2022.xlsx", "Table S9", skip = 3) %>%
  select("CAS...1", "Chemical Name...2", "Final Group") %>% rename("Group" = "Final Group", "CAS" = "CAS...1", "Chemical" = "Chemical Name...2") %>% filter(Group != "NA - ungrouped or <2 co-occuring constituents")%>% select(-Chemical)

ToxCast_groups$Grouping_type <- "ToxCast"

Final_Groupings <- bind_rows(ECOTOX_groups,Pharm_groups, ToxCast_groups)
View(Final_Groupings)

groups <- Final_Groupings %>% group_by(Group) %>% summarize(n = n_distinct(CAS))
groups1 <- left_join(groups, Final_Groupings) %>% distinct() %>% select(-n)

#read in chemical 2017 + 2018
chem_2017 <- read_excel("C:\\Users\\erinm\\OneDrive\\Desktop\\GLRI Project\\Milwauke\\Data for Papier\\Chemistry\\Compiled_Chemistry\\Chemistry Compiled_Final.xlsx", "2017") %>%
  filter(Detect == "Yes", !Site %in% c("MIM-BK", "UWM"), Media == "Water") %>% select(CAS, Chemical, detect.con.ppb, Site, Year)

chem_2018 <- read_excel("C:\\Users\\erinm\\OneDrive\\Desktop\\GLRI Project\\Milwauke\\Data for Papier\\Chemistry\\Compiled_Chemistry\\Chemistry Compiled_Final.xlsx", "2018")%>%
  filter(Detect == "Yes", !Site %in% c("MIM-BK", "UWM"), Media == "Water") %>% select(CAS, Chemical, detect.con.ppb, Site, Year)
View(chem_2018)
chemistry_compiled <- bind_rows((chem_2017 %>% select(CAS, Chemical)), (chem_2018 %>% select(CAS, Chemical))) %>% distinct()
chemistry_compiled$Group <- chemistry_compiled$Chemical
chemistry_compiled$Grouping_type <- "Single Chemical"
chemistry_compiled_1 <- chemistry_compiled %>% select(-Chemical)

#add to group 
groups4 <- bind_rows(groups1, chemistry_compiled_1)

#bind in concentration data
chem_2017_1 <- left_join(groups4, chem_2017) %>% filter(!is.na(Chemical))
chem_2017_1$Site_Year <- paste(chem_2017_1$Site, chem_2017_1$Year, sep = "_")
chem_2017_2 <- chem_2017_1 %>% select(-Site, -Year) %>% group_by(Group, Grouping_type, Site_Year) %>% summarize(Concentration = sum(detect.con.ppb))
chem_2017_3 <- chem_2017_1 %>% group_by(Site_Year, Group) %>% summarize(n_CAS = n_distinct(CAS))

chem_2017_4 <-left_join(chem_2017_2, chem_2017_3)
View(chem_2017_4)
chem_2017_4$Concentration <- ifelse(chem_2017_4$Grouping_type != "Single Chemical" & chem_2017_4$n_CAS == 1, 0, chem_2017_4$Concentration)
chem_2017_wide <- chem_2017_4 %>% ungroup() %>% select(-n_CAS, -Grouping_type) %>% spread(Group, Concentration)
chem_2017_wide[is.na(chem_2017_wide)] <- 0
chem_only_2017 <- chem_2017_wide[, colSums(chem_2017_wide !=0) > 0]

chem_2018_1 <- left_join(groups4, chem_2018) %>% filter(!is.na(Chemical)) %>% select(-Chemical)
chem_2018_1$Site_Year <- paste(chem_2018_1$Site, chem_2018_1$Year, sep = "_")
chem_2018_2 <- chem_2018_1 %>% select(-Site, -Year) %>% group_by(Group, Grouping_type, Site_Year) %>% summarize(Concentration = sum(detect.con.ppb))
chem_2018_3 <- chem_2018_1 %>% group_by(Site_Year, Group) %>% summarize(n_CAS = n_distinct(CAS))
chem_2018_4 <-left_join(chem_2018_2, chem_2018_3)
chem_2018_4$Concentration <- ifelse(chem_2018_4$Grouping_type != "Single Chemical" & chem_2018_4$n_CAS == 1, 0, chem_2018_4$Concentration)
chem_2018_wide <- chem_2018_4 %>% ungroup() %>% select(-n_CAS, -Grouping_type) %>% spread(Group, Concentration)
chem_2018_wide[is.na(chem_2018_wide)] <- 0
chem_only_2018 <- chem_2018_wide[, colSums(chem_2018_wide !=0) > 0]

chem_all_years <- bind_rows(chem_2017_4 %>% ungroup() %>% select(-n_CAS, - Grouping_type), chem_2018_4 %>% ungroup() %>%  select(-n_CAS, - Grouping_type))
chem_all_years_wide <- chem_all_years %>% spread(Group, Concentration)
chem_all_years_wide[is.na(chem_all_years_wide)] <- 0

#pull in site-specific data####
site_specific_effects <- read_excel("C:\\Users\\erinm\\OneDrive\\Desktop\\GLRI Project\\Milwauke\\Data for Papier\\Effect Data\\Effect_Data_for_all_Analyses.xlsx") %>% filter(!Site %in% c("UWM"))
list(unique(site_specific_effects$Measurement))
site_specific_effects  <- site_specific_effects %>% filter(!Site %in% c("BK", "MQ-BK", "Blank")) %>% filter(Measurement %in% c("E2-EQ", "17b_estradiol", "Testosterone", "Aryl hydrocarbon Receptor (AhR)",
                                                               "Estrogen Receptor alpha (ERa)", "Estrogen Response Element (ERE)",
                                                               "Antioxidant Response Element (ARE)-binding Nuclear factor (erythroid-derived 2)-like 2 (NRF2) (NRF2/ARE)",
                                                               "Peroxisome proliferator-activated receptor-gamma (PPARg)", "Peroxisome proliferator activating receptor (PPRE)",
                                                               "Peroxisome proliferator-activated receptor-a  (PPARa)",
                                                               "Pregnane-X-Receptor (PXR)", "CYP1A1_intestine", "CYP1A1_liver", "CYP2AD6_intestine",
                                                               "CYP2N13_intestine", "CYP3A_intestine", "CYP3A_liver", "UGT1A1_intestine", "UGT1A1_liver",
                                                               "Adrenoceptor Beta 1 (ADRB1)", "Glucocorticoid Receptor (GR)",
                                                               "Melanocortin 1 Receptor (MC1R)", "Prostaglandin D2 Receptor (PTGDR)",
                                                               "Prostaglandin E Receptor 2 (PTGER2)", "Prostaglandin I2 Receptor (PTGIR)", 
                                                               "Pregnane-X-Receptor (PXR)", 
                                                               "Retinoid X receptor-b (RXRb)"))


site_specific_effects$Site_Year <- paste(site_specific_effects$Site, site_specific_effects$Year, sep = "_")
site_specific_effects1 <- site_specific_effects %>% select(-Site, -Year)
names(site_specific_effects1)
list(unique(site_specific_effects1$`Assay Name`))

#format data####
#ATG####
ahr <- site_specific_effects1 %>% filter(Measurement == "Aryl hydrocarbon Receptor (AhR)") %>% select(Site_Year, Effect)
ahr <- left_join(ahr, chem_all_years_wide)
ahr[is.na(ahr)] <- 0

pxr_cis <- site_specific_effects %>% filter(Measurement == "Pregnane-X-Receptor (PXR)", `Assay Name` == "Attagene cis-Factorial") %>% select(Site_Year, Effect)
pxr_cis <- left_join(pxr_cis, chem_all_years_wide)
pxr_cis[is.na(pxr_cis)] <- 0

pxr_trans <- site_specific_effects1 %>%  filter(Measurement == "Pregnane-X-Receptor (PXR)", `Assay Name` == "Attagene trans-Factorial") %>% select(Site_Year, Effect)
pxr_trans <- left_join(pxr_trans, chem_all_years_wide)
pxr_trans[is.na(pxr_trans)] <- 0

ERa_trans <- site_specific_effects1 %>% filter(Measurement == "Estrogen Receptor alpha (ERa)") %>% select(Site_Year, Effect)
ERa_trans <- left_join(ERa_trans, chem_all_years_wide)
ERa_trans[is.na(ERa_trans)] <- 0

ERE_cis <- site_specific_effects %>% filter(Measurement == "Estrogen Response Element (ERE)")%>% select(Site_Year, Effect)
ERE_cis <- left_join(ERE_cis, chem_all_years_wide)
ERE_cis[is.na(ERE_cis)] <- 0

PPARg <- site_specific_effects %>% filter(Measurement ==  "Peroxisome proliferator-activated receptor-gamma (PPARg)") %>% select(Site_Year, Effect)
PPARg <- left_join(PPARg, chem_all_years_wide)
PPARg[is.na(PPARg)] <- 0

View(PPARg)

PPARa <- site_specific_effects %>% filter(Measurement ==  "Peroxisome proliferator-activated receptor-a  (PPARa)",) %>% select(Site_Year, Effect)
PPARa <- left_join(PPARa, chem_all_years_wide)
PPARa[is.na(PPARa)] <- 0
View(PPARg)

ARE <- site_specific_effects %>% filter(Measurement ==  "Antioxidant Response Element (ARE)-binding Nuclear factor (erythroid-derived 2)-like 2 (NRF2) (NRF2/ARE)",) %>% select(Site_Year, Effect)
ARE <- left_join(ARE, chem_all_years_wide)
ARE[is.na(ARE)] <- 0

PPRE <- site_specific_effects %>% filter(Measurement ==  "Peroxisome proliferator activating receptor (PPRE)") %>% select(Site_Year, Effect)
PPRE <- left_join(PPRE, chem_all_years_wide)
PPRE[is.na(PPRE)] <- 0

ADRB1 <- site_specific_effects %>% filter(Measurement ==  "Adrenoceptor Beta 1 (ADRB1)") %>% select(Site_Year, Effect)
ADRB1 <- left_join(ADRB1, chem_all_years_wide)
ADRB1[is.na(ADRB1)] <- 0

GR <- site_specific_effects %>% filter(Measurement ==  "Glucocorticoid Receptor (GR)") %>% select(Site_Year, Effect)
GR <- left_join(GR, chem_all_years_wide)
GR[is.na(GR)] <- 0

MC1R <- site_specific_effects %>% filter(Measurement ==  "Melanocortin 1 Receptor (MC1R)") %>% select(Site_Year, Effect)
MC1R <- left_join(MC1R, chem_all_years_wide)
MC1R[is.na(MC1R)] <- 0

PTGDR <- site_specific_effects %>% filter(Measurement ==  "Prostaglandin D2 Receptor (PTGDR)") %>% select(Site_Year, Effect)
PTGDR <- left_join(PTGDR, chem_all_years_wide)
PTGDR[is.na(PTGDR)] <- 0

PTGER <- site_specific_effects %>% filter(Measurement ==  "Prostaglandin E Receptor 2 (PTGER2)") %>% select(Site_Year, Effect)
PTGER <- left_join(PTGER, chem_all_years_wide)
PTGER[is.na(PTGER)] <- 0

PTGIR <- site_specific_effects %>% filter(Measurement ==  "Prostaglandin I2 Receptor (PTGIR)") %>% select(Site_Year, Effect)
PTGIR <- left_join(PTGIR, chem_all_years_wide)
PTGIR[is.na(PTGIR)] <- 0

RXRb <- site_specific_effects %>% filter(Measurement ==  "Retinoid X receptor-b (RXRb)") %>% select(Site_Year, Effect)
RXRb <- left_join(RXRb, chem_all_years_wide)
RXRb[is.na(RXRb)] <- 0

#t47D####
t47 <- site_specific_effects %>% filter(Measurement == "E2-EQ") %>% select(Site_Year, Effect)
t47 <- left_join(t47, chem_all_years_wide)
t47[is.na(t47)] <- 0

#RIA####
names(site_specific_effects1)
list(unique(site_specific_effects1$`Sample Type`))
ria <- site_specific_effects1 %>% filter(`Assay Name` == "RIA", Measurement == "17b_estradiol",`Sample Type`== "Male Fathead Minnow Serum") %>% select(Site_Year, Effect)
ria <- left_join(ria, chem_only_2017)

ria_F <-  site_specific_effects1 %>% filter(`Assay Name` == "RIA", Measurement == "17b_estradiol",`Sample Type`== "Female Fathead Minnow Serum") %>% select(Site_Year, Effect)
ria_F <- left_join(ria_F, chem_only_2017)

#format data for qPCR####
#read in chemical 2017 + 2018
chem_2017 <- read_excel("C:\\Users\\erinm\\OneDrive\\Desktop\\GLRI Project\\Milwauke\\Data for Papier\\Chemistry\\Compiled_Chemistry\\Chemistry Compiled_Final.xlsx", "2017") %>%
  filter(Detect == "Yes", !Site %in% c("MIM-BK", "UWM"), Media == "Water") %>% select(CAS, Chemical, detect.con.ppb, Site, Year)

chem_2018 <- read_excel("C:\\Users\\erinm\\OneDrive\\Desktop\\GLRI Project\\Milwauke\\Data for Papier\\Chemistry\\Compiled_Chemistry\\Chemistry Compiled_Final.xlsx", "2018")%>%
  filter(Detect == "Yes", !Site %in% c("MIM-BK", "UWM"), Media == "Water") %>% select(CAS, Chemical, detect.con.ppb, Site, Year)
View(chem_2018)
chemistry_compiled <- bind_rows((chem_2017 %>% select(CAS, Chemical)), (chem_2018 %>% select(CAS, Chemical))) %>% distinct()
chemistry_compiled$Group <- chemistry_compiled$Chemical
chemistry_compiled$Grouping_type <- "Single Chemical"
chemistry_compiled_1 <- chemistry_compiled %>% select(-Chemical)

#add to group 
groups4 <- bind_rows(groups1, chemistry_compiled_1)

#bind in concentration data
chem_2017_1 <- left_join(groups4, chem_2017) %>% filter(!is.na(Chemical))
chem_2017_1$Site_Year <- paste(chem_2017_1$Site, chem_2017_1$Year, sep = "_")
chem_2017_2 <- chem_2017_1 %>% select(-Site, -Year) %>% group_by(Group, Grouping_type, Site_Year) %>% summarize(Concentration = sum(detect.con.ppb))
chem_2017_3 <- chem_2017_1 %>% group_by(Site_Year, Group) %>% summarize(n_CAS = n_distinct(CAS))

chem_2017_4 <-left_join(chem_2017_2, chem_2017_3)
chem_2017_4$Concentration <- ifelse(chem_2017_4$Grouping_type != "Single Chemical" & chem_2017_4$n_CAS == 1, 0, chem_2017_4$Concentration)
chem_2017_wide <- chem_2017_4 %>% ungroup() %>% select(-n_CAS, -Grouping_type) %>% spread(Group, Concentration)
chem_2017_wide[is.na(chem_2017_wide)] <- 0
chem_only_2017 <- chem_2017_wide[, colSums(chem_2017_wide !=0) > 0]

chem_2018_1 <- left_join(groups4, chem_2018) %>% filter(!is.na(Chemical)) %>% select(-Chemical)
chem_2018_1$Site_Year <- paste(chem_2018_1$Site, chem_2018_1$Year, sep = "_")
chem_2018_2 <- chem_2018_1 %>% select(-Site, -Year) %>% group_by(Group, Grouping_type, Site_Year) %>% summarize(Concentration = sum(detect.con.ppb))
chem_2018_3 <- chem_2018_1 %>% group_by(Site_Year, Group) %>% summarize(n_CAS = n_distinct(CAS))
chem_2018_4 <-left_join(chem_2018_2, chem_2018_3)
chem_2018_4$Concentration <- ifelse(chem_2018_4$Grouping_type != "Single Chemical" & chem_2018_4$n_CAS == 1, 0, chem_2018_4$Concentration)
chem_2018_wide <- chem_2018_4 %>% ungroup() %>% select(-n_CAS, -Grouping_type) %>% spread(Group, Concentration)
chem_2018_wide[is.na(chem_2018_wide)] <- 0
chem_only_2018 <- chem_2018_wide[, colSums(chem_2018_wide !=0) > 0]

chem_all_years <- bind_rows(chem_2017_4 %>% ungroup() %>% select(-n_CAS, - Grouping_type), chem_2018_4 %>% ungroup() %>%  select(-n_CAS, - Grouping_type))
chem_all_years_wide <- chem_all_years %>% spread(Group, Concentration)
chem_all_years_wide[is.na(chem_all_years_wide)] <- 0


#2017
names(chem_only_2017)
chem_2017_long <- chem_only_2017 %>% gather(c("1-Methylnaphthalene":Venlafaxine), key = "Chemical", value = "Concentration")
chem_2017_con <- chem_2017_long %>% filter(Site_Year == "MED_2017") %>% rename("Control_con" = "Concentration") %>% select(-Site_Year)
chem_2017_nocon <- chem_2017_long
chem_2017_long1 <- left_join(chem_2017_nocon, chem_2017_con, by = "Chemical")
chem_2017_long1$adjusted_con <- (chem_2017_long1$Concentration + 1)/(chem_2017_long1$Control_con + 1)
chem_2017_long2 <- bind_rows(chem_2017_long1, (chem_2017_con %>% rename("adjusted_con" = "Control_con")))

qpcr_chem_2017_wide <- chem_2017_long2 %>% select(-Concentration, -Control_con) %>% spread(Chemical, adjusted_con)
View(qpcr_chem_2017_wide)

#2018
#no control detection - bind here 
names(chem_only_2018)
qpcr_chem_2018_wide <- chem_only_2018

#format qPCR data#### - filter out KKL for in vivo 2017 analyses 
names(site_specific_effects1)
list(unique(site_specific_effects1$Measurement))

#CYP1A1_intestine
cyp_1a1_intestine <- site_specific_effects1 %>% filter(`Assay Name` == "qPCR", Measurement == "CYP1A1_intestine") %>% select(Site_Year, Effect)
cyp_1a1_intestine <- left_join(cyp_1a1_intestine, qpcr_chem_2017_wide)
View(cyp_1a1_intestine)

#CYP1A1_liver_2017
cyp_1a1_liver_2017 <- site_specific_effects %>% filter(`Assay Name` == "qPCR", Measurement == "CYP1A1_liver", Year == "2017") %>% select(Site_Year, Effect)
cyp_1a1_liver_2017 <- left_join(cyp_1a1_liver_2017, qpcr_chem_2017_wide)

#CYP3A_liver
cyp3a_liver <- site_specific_effects %>% filter(`Assay Name` == "qPCR", Measurement == "CYP3A_liver", Year == "2017") %>% select(Site_Year, Effect)
cyp3a_liver <- left_join(cyp3a_liver, qpcr_chem_2017_wide)

#CYP3A_intestine
CYP3A_intestine <- site_specific_effects %>% filter(`Assay Name` == "qPCR", Measurement == "CYP3A_intestine", Year == "2017") %>% select(Site_Year, Effect)
CYP3A_intestine <- left_join(CYP3A_intestine, qpcr_chem_2017_wide)

#UGT_liver
UGT_liver <- site_specific_effects %>% filter(`Assay Name` == "qPCR", Measurement == "UGT1A1_liver", Year == "2017") %>% select(Site_Year, Effect)
UGT_liver <- left_join(UGT_liver, qpcr_chem_2017_wide)

#UGT_intestine
UGT_intestine <- site_specific_effects %>% filter(`Assay Name` == "qPCR", Measurement == "UGT1A1_intestine", Year == "2017") %>% select(Site_Year, Effect)
UGT_intestine <- left_join(UGT_intestine, qpcr_chem_2017_wide)

#CYP_2N13
CYP2N13 <- site_specific_effects %>% filter(`Assay Name` == "qPCR", Measurement == "CYP2N13_intestine", Year == "2017") %>% select(Site_Year, Effect)
CYP2N13 <- left_join(CYP2N13, qpcr_chem_2017_wide)

#CYP_2AD6
CYP_2AD6 <- site_specific_effects %>% filter(`Assay Name` == "qPCR", Measurement == "CYP2AD6_intestine", Year == "2017") %>% select(Site_Year, Effect)
CYP_2AD6 <- left_join(CYP_2AD6, qpcr_chem_2017_wide)

#CYP1A1_2018
cyp_1a1_2018 <- site_specific_effects %>% filter(`Assay Name` == "qPCR", Measurement == "CYP1A1_liver", Year == "2018") %>% select(Site_Year, Effect)
cyp_1a1_2018 <- left_join(cyp_1a1_2018, qpcr_chem_2018_wide)
cyp_1a1_2018[is.na(cyp_1a1_2018)] <- 0
View(cyp_1a1_2018)
