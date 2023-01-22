#Pharmaceutical Potential Mixture Analysis#

library(dplyr)
library(tidyverse)
library(ggpubr)
library(readxl)
library(writexl)
library(gridExtra)
library(viridisLite)

setwd("C:\\Users\\erinm\\OneDrive\\Desktop\\GLRI Project\\Milwauke\\Data for Papier\\Mixtures MS\\Pharmaceutical_Analysis")
#import + rename 
targets <- read_excel("DrugBank_Targets.xlsx") %>% rename("Uniprot_ID" = "UniprotIDs") %>% filter(Pharm_Activity == "Active")
targets$target <- targets$`DrugBank Targets`
names(targets)

SEQAPASS <- read_excel("susceptible_fish_final_annotated.xlsx") %>% rename("target" = "gene", "Level_1_SEQAPASS_Susceptible" = "Level_1", "Level_2_SEQAPASS_Susceptible" = "Level_2") %>% select(-general_taxa)

Table_s10 <- left_join(targets, SEQAPASS)
Table_s10$Susceptibility_Class <- ifelse(Table_s10$Level_1_SEQAPASS_Susceptible == "Susceptible" | Table_s10$Level_2_SEQAPASS_Susceptible == "Susceptible", "Fish_Susceptible", "Not_Fish_Susceptible")
Table_s10$Level_1_SEQAPASS_Susceptible[is.na(Table_s10$Level_1_SEQAPASS_Susceptible)] <- "Not_Susceptible"
Table_s10$Level_2_SEQAPASS_Susceptible[is.na(Table_s10$Level_2_SEQAPASS_Susceptible)] <- "Not_Susceptible"
Table_s10$LoE_for_Susceptibility[is.na(Table_s10$LoE_for_Susceptibility)] <- 0
Table_s10$Susceptibility_Class[is.na(Table_s10$Susceptibility_Class)] <- "Not_Fish_Susceptible"

write_xlsx(Table_s10, "Table_s10_preDAVID.xlsx")

DAVID <- read_excel("david_list_annotated_clusterHIGH_only.xlsx") %>% rename("target" = "Protein Target", "Cluster" = "Cluster_Group_n", "DAVID_Uniprot_ID" = "Uniprot_ID")

#bind + annotate
names(targets)
names(SEQAPASS)
names(DAVID)
names(cmax)

Table_s10_1 <- left_join(Table_s10, DAVID, by = "target")
Table_s10_1$Cluster[is.na(Table_s10_1$Cluster)] <- "-"
Table_s10_1$Enrichment_Score[is.na(Table_s10_1$Enrichment_Score)] <- "-"
Table_s10_1$`Final Classification`[is.na(Table_s10_1$`Final Classification`)] <- "-"
Table_s10_1$`Final Classification` <- ifelse(Table_s10_1$Susceptibility_Class == "Fish_Susceptible" & Table_s10_1$`Final Classification` == "-", "Unclustered",
                                             Table_s10_1$`Final Classification`)
View(Table_s10_1)

cmax <- read_excel("CMAX.xlsx") %>% select(-CAS) %>% rename(Chemical = "Chemical Name")
names(cmax)

Table_s10_2 <- left_join(Table_s10_1, cmax, by = "Chemical")
Table_s10_2$Class_Direction <- paste(Table_s10_2$`Final Classification`, Table_s10_2$Direction, sep = "-")

Table_s10_fin <- Table_s10_2 %>% relocate("Chemical", .before = "Uniprot_ID") %>% relocate("CAS_OG", .after = "Chemical") %>% relocate("Cmax":"units", .after = "Class_Direction")
write_xlsx(Table_s10_fin, "Table_s10.xlsx")

#pull together SI table for SEQAPASS output
SEQA_l1 <- read_excel("SEQAPASS_Input.xlsx", 1)
SEQA_l2 <- read_excel("SEQAPASS_Input.xlsx", 2)

SEQAPASS_SI_Table10 <- left_join(SEQA_l1, SEQA_l2)
write_xlsx(SEQAPASS_SI_Table10, "SEQAPASS_SI_Table10.xlsx")

#MCR analysis 
#pharmaceutical
pharm_for_MCR <- read_excel("Table_s10_EM_annotated_for_analysis.xlsx")
cmax <- read_excel("CMAX.xlsx")

#bind in cmax
pharm_cmax <- left_join(pharm_for_MCR, (cmax %>% select(-"Chemical Name", -CAS) %>% rename(CAS = CAS_OG)))
pharm_cmax$CAS <- ifelse(pharm_cmax$CAS == "299-42-3", "90-82-4/299-42-3", pharm_cmax$CAS)
pharm_cmax$CAS <- ifelse(pharm_cmax$CAS == "34841-39-9", "34911-55-2", pharm_cmax$CAS)

#chemistry
chem_2017 <- read_excel("C:\\Users\\erinm\\OneDrive\\Desktop\\GLRI Project\\Milwauke\\Data for Papier\\Chemistry\\Compiled_Chemistry\\Chemistry Compiled_Final.xlsx", "2017")%>%
  filter(!Site %in% c("MED", "UWM", "MIM-BK"), !Detect.code %in% c("<"), Media == "Water") %>% select(Site, Year, detect.con.ppb, CAS, Chemical) %>% 
  rename("Concentration_ugL" = "detect.con.ppb")
chem_2017$Site_Year <- paste(chem_2017$Site, chem_2017$Year)

chem_2018 <- read_excel("C:\\Users\\erinm\\OneDrive\\Desktop\\GLRI Project\\Milwauke\\Data for Papier\\Chemistry\\Compiled_Chemistry\\Chemistry Compiled_Final.xlsx", "2018") %>%
  filter(!Site %in% c("MED", "UWM"), !Detect.code %in% c("<"), Media == "Water") %>% select(Site, Year, detect.con.ppb, CAS, Chemical) %>% 
  rename("Concentration_ugL" = "detect.con.ppb")
chem_2018$Site_Year <- paste(chem_2018$Site, chem_2018$Year)

chem_2018$Year <- as.double(chem_2018$Year)
chemistry <- bind_rows(chem_2017, chem_2018) %>% select(-Site, -Year, -Chemical)
View(chemistry)

#bind pharm_chem
pharm_chemistry <- left_join(chemistry,pharm_cmax, by = "CAS") %>% relocate("CAS", .before = "Chemical") %>% filter(Fish_Susceptible =="Fish_Susceptible") %>%
  select(CAS, Chemical, Final_Mixture_Group, Cmax, Concentration_ugL, Site_Year) %>% distinct()
names(pharm_chemistry)
pharm_chemistry$Cmax <- as.numeric(pharm_chemistry$Cmax)
pharm_chemistry$Concentration_ugL <- as.numeric(pharm_chemistry$Concentration_ugL)
View(pharm_chemistry)

#calculate MCRs
#chemical_based
pharm_chemical <- pharm_chemistry %>% select(CAS, Chemical, Concentration_ugL, Cmax, Site_Year) %>% distinct()
pharm_chemical$TQ <- pharm_chemical$Concentration_ugL / pharm_chemical$Cmax

site_year_sumTQ <- pharm_chemical %>% group_by(Site_Year) %>% summarize(Site_Year_sumTQ = sum(TQ))
pharm_chemical1 <- left_join(pharm_chemical, site_year_sumTQ)

pharm_chemical1$MCR_chemical <- (pharm_chemical1$Site_Year_sumTQ)/ (pharm_chemical1$TQ)
pharm_chemical1$MCR_annotation_chemical <- ifelse(pharm_chemical1$MCR_chemical > 5, "> 5 (CR)",
                                             pharm_chemical1$Chemical)

#target_group based
names(pharm_chemistry)
pharm_chemistry_nocooccur <- pharm_chemistry %>% filter(Final_Mixture_Group != "NA - < 2 co-occuring constituents") %>% group_by(Final_Mixture_Group, Site_Year) %>% summarize(n_uniqueCAS = n_distinct(CAS))
View(pharm_chemistry_nocooccur)  

pharm_chemistry1 <- left_join(pharm_chemistry, pharm_chemistry_nocooccur)
pharm_chemistry1$n_uniqueCAS[is.na(pharm_chemistry1$n_uniqueCAS)] <- 0
pharm_chemistry1$Final_Mixture_Group <- ifelse(pharm_chemistry1$n_uniqueCAS == 1, "NA - < 2 co-occuring constituents", pharm_chemistry1$Final_Mixture_Group)

pharm_target <- pharm_chemistry1 %>% select(CAS, Chemical, Final_Mixture_Group, Cmax, Concentration_ugL, Site_Year, n_uniqueCAS) %>%
  filter(Final_Mixture_Group != "-", Final_Mixture_Group != "NA - < 2 co-occuring constituents") %>% distinct()

pharm_target$TQ <- pharm_target$Concentration_ugL / pharm_target$Cmax

site_year_group_sumTQ <- pharm_target %>% group_by(Site_Year, Final_Mixture_Group) %>% summarize(Site_Year_Group_sumTQ = sum(TQ)) %>% distinct()

pharm_target2 <- left_join(pharm_target, site_year_group_sumTQ)
pharm_target3 <- left_join(pharm_target2, site_year_sumTQ)

pharm_target3$MCR_target <- pharm_target3$Site_Year_sumTQ/(pharm_target3$Site_Year_Group_sumTQ)
pharm_target3$MCR_annotation_group <- ifelse(pharm_target3$MCR_target > 5, "> 5 (CR)",
                                         pharm_target3$Final_Mixture_Group)


View(pharm_target3)

#bind together for final SI spreadsheet
names(pharm_chemistry_4)
pharm_chemistry_final1 <- left_join(pharm_chemistry, pharm_chemical1)
pharm_chemistry_final2 <- left_join(pharm_chemistry_final1, pharm_target3) %>% select(-n_uniqueCAS)
pharm_chemistry_final2$Site_Year_Group_sumTQ[is.na(pharm_chemistry_final2$Site_Year_Group_sumTQ)] <- "-"
pharm_chemistry_final2$MCR_target[is.na(pharm_chemistry_final2$MCR_target)] <- "-"
pharm_chemistry_final2$MCR_annotation_group[is.na(pharm_chemistry_final2$MCR_annotation_group)] <- "> 5 (CR)"

write_xlsx(pharm_chemistry_final2, "pharm_LoE.xlsx")

#figures
UCJ_2017 <- data.frame(Site_Year = c("UCJ 2017"), Site_Year_sumTQ = c(0)) #add dataframe for UCJ because no pharmaceuticals detected in 2017
pharm_chemistry_final2 <- bind_rows(pharm_chemistry_final2, UCJ_2017)
pharm_chemical1 <- bind_rows(pharm_chemical1, UCJ_2017)

target_MCR_2017 <- pharm_chemistry_final2 %>% separate("Site_Year", into = c("site", "year")) %>% filter(year == "2017") %>%
  arrange((Site_Year_sumTQ)) %>% distinct() %>% rename("Target-Direction Group" = "MCR_annotation_group")%>%  
  ggbarplot("site", "TQ", orientation = "horiz", legend = "right", color = NA, title = "A", fill = "Target-Direction Group",
            palette = viridis(n = 6, direction = 1)) + 
  xlab("Site (2017)") + ylab("Pharmaceutical Potential (C/Cmax)") + ylim(0, 0.15) + rremove("legend")
target_MCR_2017


chem_MCR_2017 <- pharm_chemical1 %>% arrange((Site_Year_sumTQ)) %>%
  separate("Site_Year", into = c("site", "year")) %>% filter(year == "2017") %>%
  select(-Chemical) %>% distinct() %>% rename("Chemical" = "MCR_annotation_chemical")%>%
  ggbarplot("site", "TQ", orientation = "horiz", legend = "right", color = NA, title = "D", fill = "Chemical",
            palette = c("#000004", "#721f81", "#b73779", "#feb078")) + 
  xlab("Site (2017)") + ylab("Pharmaceutical Potential (C/Cmax)")+ ylim(0, 0.15)+ rremove("legend")

chem_MCR_2017

target_MCR_2018 <- pharm_chemistry_final2 %>% separate("Site_Year", into = c("site", "year")) %>% filter(year == "2018") %>%
  arrange((Site_Year_sumTQ)) %>% distinct() %>% rename("Target-Direction Group" = "MCR_annotation_group")%>%  
  ggbarplot("site", "TQ", orientation = "horiz", legend = "right", color = NA, title = "C", fill = "Target-Direction Group",
            palette = viridis(n = 6, direction = 1)) + 
  xlab("Site (2018)") + ylab("Pharmaceutical Potential (C/Cmax)")+ ylim(0, 0.15)+ rremove("legend")

target_MCR_2018

chem_MCR_2018 <- pharm_chemical1 %>% arrange((Site_Year_sumTQ)) %>%
  separate("Site_Year", into = c("site", "year")) %>% filter(year == "2018") %>%
  select(-Chemical) %>% distinct() %>% rename("Chemical" = "MCR_annotation_chemical")%>%
  ggbarplot("site", "TQ", orientation = "horiz", legend = "right", color = NA, title = "D", fill = "Chemical",
            palette = viridis(n = 6, option = "magma")) + 
  xlab("Site (2018)") + ylab("Pharmaceutical Potential (C/Cmax)")+ ylim(0, 0.15)+ rremove("legend")

chem_MCR_2018

#generate TQ vs. MCR plots
names(pharm_chemistry_final2)
list(unique(pharm_chemistry_final2$Final_Mixture_Group))

TQvMCR_group <- pharm_chemistry_final2  %>% select(Final_Mixture_Group, MCR_target, Site_Year_Group_sumTQ) %>%
  rename("Mixture" = "Final_Mixture_Group", "MCR" = "MCR_target", "TQ" = "Site_Year_Group_sumTQ") %>% distinct() %>% filter(TQ != "-")
TQvMCR_group$MCR <- as.double(TQvMCR_group$MCR)
TQvMCR_group$TQ <- as.double(TQvMCR_group$TQ)

TQvMCR_single <- pharm_chemistry_final2 %>% select(Chemical, MCR_chemical, TQ) %>%
  rename("Chemical" = "Chemical", "MCR" = "MCR_chemical", "TQ" = "TQ") %>% distinct()

TQvMCR_group_plot <-TQvMCR_group %>% filter(MCR < 10) %>% ggscatter("TQ", "MCR", color = "Mixture", 
                                                             palette = c("#3b528b"), title = "E", legend = "right", size = 3) + 
  geom_hline(yintercept = 5, linetype = "dashed") +
  geom_vline(xintercept = 0.01, linetype = "dashed") + 
  xlab("Toxicity Quotient") + ylab("Maximum Cumulative Ratio") + scale_x_log10() + ylim(0, 10)+ rremove("legend")

TQvMCR_group_plot


TQvMCR_chemical_plot <-TQvMCR_single %>% filter(MCR < 10) %>% ggscatter("TQ", "MCR", color = "Chemical", 
                                                                    palette = c("#3b0f70", "#8c2981", "#de4968", "#fe9f6d",
                                                                                "#fdda9c", "#fcfdbf"), title = "F", legend = "right", size = 3) + 
  geom_hline(yintercept = 5, linetype = "dashed") +
  geom_vline(xintercept = 0.01, linetype = "dashed") + 
  xlab("Toxicity Quotient") + ylab("Maximum Cumulative Ratio") + scale_x_log10()+ rremove("legend")

TQvMCR_chemical_plot

Pharm_MCR <- grid.arrange(target_MCR_2017, chem_MCR_2017, target_MCR_2018, chem_MCR_2018, TQvMCR_group_plot, TQvMCR_chemical_plot)
ggsave("Pharm_MCR.jpeg", Pharm_MCR, height = 15, width = 15)
