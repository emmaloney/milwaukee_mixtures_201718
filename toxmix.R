#toxeval-based mixture assessment 
library(tidyverse)
library(readxl)
library(writexl)
library(gridExtra)
library(ggpubr)
library(viridisLite)
library(stringr)

setwd("C:\\Users\\erinm\\OneDrive\\Desktop\\GLRI Project\\Milwauke\\Data for Papier\\Mixtures MS\\ToxCast_LoE")

#prepare data
chem_2017 <- read_excel("C:\\Users\\erinm\\OneDrive\\Desktop\\GLRI Project\\Milwauke\\Data for Papier\\Chemistry\\Compiled_Chemistry\\Chemistry Compiled_Final.xlsx", "2017") %>%
  filter(!Site %in% c("MED", "UWM", "MIM-BK"), !is.na(Site), Media == "Water")
list(unique(chem_2017$Site))

chem_2018 <- read_excel("C:\\Users\\erinm\\OneDrive\\Desktop\\GLRI Project\\Milwauke\\Data for Papier\\Chemistry\\Compiled_Chemistry\\Chemistry Compiled_final.xlsx", "2018") %>%
  filter(!Site %in% c("MED", "UWM"))
list(unique(chem_2018$Site))

#2017
names(chem_2017)
chem_2017_Chemicals <- chem_2017 %>% select(CAS, Chemical)
chem_2017_Chemicals$Class <- "No_Class"

chem_2017_Data <- chem_2017 %>% select(CAS, Site.code, Detect.code, detect.con.ppb, RL.ppb)
list(unique(chem_2017_Data$Detect.code))

chem_2017_Data$Value <- ifelse(is.na(chem_2017_Data$Detect.code), chem_2017_Data$detect.con.ppb,
                               ifelse(chem_2017_Data$Detect.code == "E", chem_2017_Data$detect.con.ppb,
                                      ifelse(chem_2017_Data$Detect.code == "M", (chem_2017_Data$RL.ppb/4), 
                                             "<")))
chem_2017_Data$Value <- as.numeric(chem_2017_Data$Value)
chem_2017_Data_1 <- chem_2017_Data %>% filter(!is.na(Value), !is.na(CAS)) %>% rename("SiteID" = "Site.code") %>% select(CAS, SiteID, Value)
chem_2017_Data_1$`Sample Date` <- 2017
View(chem_2017_Data_1)
chem_2017_Data_1$SiteID <- paste(chem_2017_Data_1$SiteID, chem_2017_Data_1$`Sample Date`, sep = "_")

View(chem_2017_Data)
chem_2017_Sites <- chem_2017 %>% select(Site.code, Site) %>% rename("SiteID" = "Site.code", "Short Name" = "Site")
chem_2017_Sites$Year <- 2017
chem_2017_Sites$SiteID <- paste(chem_2017_Sites$SiteID, chem_2017_Sites$Year, sep = "_")
chem_2017_Sites <- chem_2017_Sites %>% select(-Year)

#2018
names(chem_2018)
chem_2018_Chemicals <- chem_2018 %>% select(CAS, Chemical)
chem_2018_Chemicals$Class <- "No_Class"

chem_2018_Data <- chem_2018 %>% select(CAS, Site.code, Detect.code, detect.con.ppb, RL.ppb)
list(unique(chem_2017_Data$Detect.code))

chem_2018_Data$Value <- ifelse(is.na(chem_2018_Data$Detect.code), chem_2018_Data$detect.con.ppb,
                               ifelse(chem_2018_Data$Detect.code == "E", chem_2018_Data$detect.con.ppb,
                                      ifelse(chem_2018_Data$Detect.code == "M", (chem_2018_Data$RL.ppb/4), 
                                             "<")))
chem_2018_Data$Value <- as.numeric(chem_2018_Data$Value)
chem_2018_Data_1 <- chem_2018_Data %>% filter(!is.na(Value), !is.na(CAS)) %>% rename("SiteID" = "Site.code") %>% select(CAS, SiteID, Value)
chem_2018_Data_1$`Sample Date` <- 2018
chem_2018_Data_1$SiteID <- paste(chem_2018_Data_1$SiteID, chem_2018_Data_1$`Sample Date`, sep = "_")
View(chem_2018_Data_1)

chem_2018_Sites <- chem_2018 %>% select(Site.code, Site) %>% rename("SiteID" = "Site.code", "Short Name" = "Site")
chem_2018_Sites$Year <- 2018
chem_2018_Sites$SiteID <- paste(chem_2018_Sites$SiteID, chem_2018_Sites$Year, sep = "_")
chem_2018_Sites <- chem_2018_Sites %>% select(-Year)

#compile
Chemicals <- bind_rows(chem_2018_Chemicals, chem_2017_Chemicals) %>% distinct() %>% filter(!is.na(CAS))
Chemicals$CAS <- gsub("34911-55-2", "31677-93-7", Chemicals$CAS)
Chemicals$CAS <- gsub("59729-33-8", "59729-32-7", Chemicals$CAS)
Chemicals$CAS <- gsub("125-71-3", "125-69-9", Chemicals$CAS)
Chemicals$CAS <- gsub("83799-24-0", "153439-40-8", Chemicals$CAS)
Chemicals$CAS <- gsub("76-99-3", "1095-90-5", Chemicals$CAS)
Chemicals$CAS <- gsub("76-22-2", "464-48-2", Chemicals$CAS)
Chemicals$CAS <- gsub("58-73-1", "147-24-0", Chemicals$CAS)
Chemicals$CAS <- gsub("42399-41-7", "33286-225", Chemicals$CAS)
Chemicals$CAS <- gsub("27203-92-5", "36282-47-0", Chemicals$CAS)
Chemicals$CAS <- gsub("93413-69-5", "99300-78-4", Chemicals$CAS)
Chemicals$CAS <- gsub("90-82-4/299-42-3", "50-98-6", Chemicals$CAS)
Chemicals$CAS <- gsub("525-66-6", "318-98-9", Chemicals$CAS)

list(unique(Chemicals$CAS))

Data <- bind_rows(chem_2018_Data_1, chem_2017_Data_1) %>% filter(!is.na(CAS))
Data$CAS <- gsub("34911-55-2", "31677-93-7", Data$CAS)
Data$CAS <- gsub("59729-33-8", "59729-32-7", Data$CAS)
Data$CAS <- gsub("125-71-3", "125-69-9", Data$CAS)
Data$CAS <- gsub("83799-24-0", "153439-40-8", Data$CAS)
Data$CAS <- gsub("76-99-3", "1095-90-5", Data$CAS)
Data$CAS <- gsub("76-22-2", "464-48-2", Data$CAS)
Data$CAS <- gsub("58-73-1", "147-24-0", Data$CAS)
Data$CAS <- gsub("42399-41-7", "33286-225", Data$CAS)
Data$CAS <- gsub("27203-92-5", "36282-47-0", Data$CAS)
Data$CAS <- gsub("93413-69-5", "99300-78-4", Data$CAS)
Data$CAS <- gsub("90-82-4/299-42-3", "50-98-6", Data$CAS)
Data$CAS <- gsub("525-66-6", "318-98-9", Data$CAS)

list(unique(Data$CAS))

Sites <- bind_rows(chem_2017_Sites, chem_2018_Sites) %>% distinct
names(Sites)
list(unique(Sites$`Short Name`))

#write to xlsx
write_xlsx(Chemicals, "toxEval_chemicals.xlsx")
write_xlsx(Data, "toxEval_data.xlsx")
write_xlsx(Sites, "toxEval_sites.xlsx")


#start with toxEval analysis
library(toxEval)

tox_list <- create_toxEval("toxmix_toxEval.xlsx")
ACClong <- get_ACC(tox_list$chem_info$CAS)
ACClong <- remove_flags(ACClong)
list(unique(ACClong$flags))

cleaned_ep <- clean_endPoint_info(end_point_info)
filtered_ep <- filter_groups(cleaned_ep, 
                             groupCol = 'intended_target_family',
                             assays = c("ATG", "NVS", "OT", "TOX21", "CEETOX", "CLD", "TANGUAY", "NHEERL_PADILLA", "NCCT_SIMMONS", "ACEA"),
                             remove_groups = c('Background Measurement','Undefined', 'Cell Cycle', 'Cell Morphology'))

chemicalSummary <- get_chemical_summary(tox_list, 
                                        ACClong, 
                                        filtered_ep)

View(chemicalSummary)

#annotate Site_Year
names(chemicalSummary)
list(unique(chemicalSummary$endPoint))
chemicalSummary1 <- chemicalSummary
chemicalSummary1$site_year <- paste(chemicalSummary1$shortName, chemicalSummary1$date, sep = "_")

write_xlsx(chemicalSummary1, "ToxCast_grouping.xlsx")

#generate TQ's 
chemicalSummary2 <- read_excel("ToxCast_grouping.xlsx")
Sum_Chem_Site <- chemicalSummary2 %>% group_by(chnm, site_year) %>% summarize(sumEAR_chemsite = sum(EAR))
Sum_Site_Year <- chemicalSummary2 %>% group_by(site_year) %>% summarize(sumEAR_siteyear = sum(EAR))
chemicalSummary3 <- left_join(chemicalSummary2, Sum_Chem_Site)
chemicalSummary4 <- left_join(chemicalSummary3, Sum_Site_Year)

#read in Dan's file
endpoint_grouping <- read_csv("C:\\Users\\erinm\\OneDrive\\Desktop\\GLRI Project\\Milwauke\\Data for Papier\\Mixtures MS\\ToxCast_LoE\\ToxCast_MIE-DESKTOP-EL682HP\\genes_endpoints_for_review_DLV_corrected 11-09-2020.csv")
endpoint_grouping$gene_group <- paste(endpoint_grouping$gene_symbol_DLV, endpoint_grouping$`Activation/inhibition/binding`, sep = "_")
names(endpoint_grouping)
View(endpoint_grouping)
endpoint_grouping1 <- endpoint_grouping %>% select(aenm, gene_group) %>% rename("endPoint" = "aenm")

#join annotation to chemSummary file
chemicalSummary5 <- left_join(chemicalSummary4,endpoint_grouping1)

#remove all of the capitalization-related mistakes

chemicalSummary6 <- chemicalSummary5 %>% filter(!is.na(gene_group))
write_xlsx(chemicalSummary6, "toxcast_file_for_genegroup_annotation.xlsx")
list(unique(chemicalSummary6$gene_group))
View(chemicalSummary6)

chemicalSummary6$gene_group <- gsub("Ache_Inhibition", "ACHE_Inhibition", chemicalSummary6$gene_group)
chemicalSummary6$gene_group <- gsub("Adra2a_Binding", "ADRA2A_Binding", chemicalSummary6$gene_group)
chemicalSummary6$gene_group <- gsub("Ar_Binding", "AR_Binding", chemicalSummary6$gene_group)
chemicalSummary6$gene_group <- gsub("Chrm3_Binding", "CHRM3_Binding", chemicalSummary6$gene_group)
chemicalSummary6$gene_group <- gsub("Gabra1_Inhibition", "GABRA1_Inhibition", chemicalSummary6$gene_group)
chemicalSummary6$gene_group <- gsub("Oprm1_Binding", "OPRM1_Binding", chemicalSummary6$gene_group)
chemicalSummary6$gene_group <- gsub("Slc6a2_Binding", "SLC6A2_Binding", chemicalSummary6$gene_group)
chemicalSummary6$gene_group <- gsub("Slc6a3_Binding", "SLC6A3_Binding", chemicalSummary6$gene_group)
chemicalSummary6$gene_group <- gsub("Slc6a4_Binding", "SLC6A4_Binding", chemicalSummary6$gene_group)
chemicalSummary6$gene_group <- gsub("Thrb_Inhibition", "THRB_Inhibition", chemicalSummary6$gene_group)
chemicalSummary6$gene_group <- gsub("Tspo_Inhibition", "TSPO_Inhibition", chemicalSummary6$gene_group)

list(unique(chemicalSummary6$gene_group))

#regroup - eliminating any 'groups' that only consist of one single chemical & any group_site_years that only consist of a single chemical
gene_group_sum <- chemicalSummary6 %>% group_by(gene_group) %>% summarize(n_distinct_chem = n_distinct(chnm))
chemicalSummary7 <- left_join(chemicalSummary6, gene_group_sum)
chemicalSummary7$gene_group <- ifelse(chemicalSummary7$n_distinct_chem == 1, chemicalSummary7$chnm, chemicalSummary7$gene_group)

write_xlsx(chemicalSummary7, "toxcast_mixture_grouping_r2.xlsx")

chemicalSummary8 <- chemicalSummary7 %>% select(-n_distinct_chem)
gene_group_siteyear <- chemicalSummary8 %>% group_by(gene_group, site_year) %>% summarize(n_distinct_chem = n_distinct(chnm))
chemicalSummary9 <- left_join(chemicalSummary8, gene_group_siteyear)
chemicalSummary9$gene_group <- ifelse(chemicalSummary9$n_distinct_chem == 1, chemicalSummary9$chnm, chemicalSummary9$gene_group)

chemicalSummary10 <- chemicalSummary9 %>% select(-n_distinct_chem) %>% filter(!is.na(gene_group))
chemicalSummary10$gene_group_type <- ifelse(chemicalSummary10$chnm == chemicalSummary10$gene_group, "Single Compound", "Group")
list(unique(chemicalSummary10$gene_group_type))
View(chemicalSummary10)

#generate MCRs
Sum_gene_site <- chemicalSummary10 %>% group_by(gene_group, site_year) %>% summarize(sumEAR_genesite = sum(EAR))
chemicalSummary11 <- left_join(chemicalSummary10, Sum_gene_site)

chemicalSummary11$MCR_chem <- (chemicalSummary11$sumEAR_siteyear)/(chemicalSummary11$sumEAR_chemsite)
chemicalSummary11$MCR_gene <-  (chemicalSummary11$sumEAR_siteyear)/(chemicalSummary11$sumEAR_genesite)

#pare down columns
names(chemicalSummary11)
chemicalSummary12 <- chemicalSummary11 %>% select(-site, -date, -Bio_category, -shortName, -Class)

#Annotate MCRs 
chemicalSummary12$MCR_annotation_chemical <- ifelse(chemicalSummary12$MCR_chem > 5, "> 5 (CR)",
                                                    chemicalSummary12$chnm)

chemicalSummary12$MCR_annotation_group <- ifelse(chemicalSummary12$MCR_gene > 5, "> 5 (CR)",
                                                 chemicalSummary12$gene_group)


write_xlsx(chemicalSummary12, "ToxCast_output.xlsx")

#Generate MCR plot
names(chemicalSummary12)
chemicalSummary12 <- read_excel("ToxCast_output.xlsx")
chemicalSummary12$MCR_annotation_group<- gsub("4-(1,1,3,3-Tetramethylbutyl)phenol", "4-tert-octylphenol", chemicalSummary12$MCR_annotation_group)
chemicalSummary12$MCR_annotation_group<- gsub("Diphenhydramine hydrochloride", "Diphenhydramine", chemicalSummary12$MCR_annotation_group)
View(chemicalSummary12)
list(unique(chemicalSummary12$site_year))

ToxCast_MCR_plot_2017 <- chemicalSummary12 %>% filter(gene_group_type == "Group", site_year %in% c("MIE_2017", "MIM_2017", "MIP_2017",
                                                                                                   "MEC_2017", "KKL_2017", "MET_2017", "UCJ_2017", "MEF_2017")) %>%
  arrange((sumEAR_siteyear)) %>% select(site_year, sumEAR_genesite, MCR_annotation_group) %>% rename("Target-Direction Group" = "MCR_annotation_group") %>% 
  separate("site_year", into = c("site", "year")) %>%
  distinct()%>%
  ggbarplot("site", "sumEAR_genesite", fill = "Target-Direction Group", orientation = "horiz", legend = "right", color = NA, title = "A", 
            palette = c("#440154", "#fde725", "#c8e020", "#90d743", "#21918c", "#287c8e", "#31688e", "#443983")) + 
  xlab("Site (2017)") + ylab("Exposure Activity Ratio (EAR)") + geom_vline(xintercept = 0.01) +rremove("legend")
  

ToxCast_MCR_chem_plot_2017 <- chemicalSummary12 %>% filter(site_year %in% c("MIE_2017", "MIM_2017", "MIP_2017",
                                                                                                   "MEC_2017", "KKL_2017", "MET_2017", "UCJ_2017", "MEF_2017")) %>% arrange((sumEAR_siteyear)) %>%
  separate("site_year", into = c("site", "year")) %>%
  select(site, sumEAR_chemsite, MCR_annotation_chemical) %>% rename("Chemical" = "MCR_annotation_chemical") %>%
  distinct()%>%
  ggbarplot("site", "sumEAR_chemsite", fill = "Chemical", orientation = "horiz", legend = "right", color = NA, title = "B", 
            palette = c("#000004", "#febf84", "#fe9f6d", "#f1605d", "#a8327d", "#8c2981", "#0b0924")) + 
  xlab("Site (2017)") + ylab("Exposure Activity Ratio (EAR)") + geom_vline(xintercept = 0.01) +rremove("legend")

ToxCast_MCR_plot_2018 <- chemicalSummary12 %>% filter(gene_group_type == "Group", site_year %in% c("CCM_2018", "JIP_2018", "MIN_2018", "KKL_2018", "MEC_2018", "MIE_2018", "MIM_2018", "MIP_2018", "UCJ_2018", "MEF_2018")) %>%
  separate("site_year", into = c("site", "year")) %>%
  arrange((sumEAR_siteyear)) %>% select(site, sumEAR_genesite, MCR_annotation_group) %>% rename("Target-Direction Group" = "MCR_annotation_group") %>% 
  distinct()%>%
  ggbarplot("site", "sumEAR_genesite", fill = "Target-Direction Group", orientation = "horiz", legend = "right", color = NA, title = "C", 
            palette = c("#440154", "#90d743", "#20a486", "#443983")) + 
  xlab("Site (2018)") + ylab("Exposure Activity Ratio (EAR)") + geom_vline(xintercept = 0.01) +rremove("legend") +ylim(0, 0.2)


ToxCast_MCR_chem_plot_2018 <- chemicalSummary12 %>% filter(site_year %in% c("CCM_2018", "JIP_2018", "MIN_2018", "KKL_2018", "MEC_2018", "MIE_2018", "MIM_2018", "MIP_2018", "UCJ_2018", "MEF_2018")) %>% arrange((sumEAR_siteyear)) %>%
  select(site_year, sumEAR_chemsite, MCR_annotation_chemical) %>% rename("Chemical" = "MCR_annotation_chemical") %>%
  separate("site_year", into = c("site", "year")) %>%
  distinct()%>%
  ggbarplot("site", "sumEAR_chemsite", fill = "Chemical", orientation = "horiz", legend = "right", color = NA, title = "D", 
            palette = c("#000004", "#fcfdbf", "#fddea0", "#febf84", "#de4968", "#c43c75", "#8c2981", "#3b0f70")) + 
  xlab("Site (2018)") + ylab("Exposure Activity Ratio (EAR)") + geom_vline(xintercept = 0.01) +rremove("legend")+ylim(0, 0.2)


#generate MCR vs. EAR plot
names(chemicalSummary12)
EARvMCR_group <- chemicalSummary12 %>% filter(gene_group_type != "Single Compound") %>% select(gene_group, MCR_gene, sumEAR_genesite, site_year) %>%
  rename("Mixture" = "gene_group", "MCR" = "MCR_gene", "EAR" = "sumEAR_genesite") %>% distinct()

EARvMCR_single <- chemicalSummary12 %>% select(chnm, MCR_chem, sumEAR_chemsite, site_year) %>%
  rename("Chemical" = "chnm", "MCR" = "MCR_chem", "EAR" = "sumEAR_chemsite") %>% distinct()
EARvMCR_single$Chemical <- gsub("4-(1,1,3,3-Tetramethylbutyl)phenol", "4-tert-octylphenol", EARvMCR_single$Chemical )
EARvMCR_single$Chemical <- gsub("Diphenhydramine hydrochloride", "Diphenhydramine", EARvMCR_single$Chemical)

EARvMCRgroup_plot <-EARvMCR_group %>% filter(MCR < 10) %>% ggscatter("EAR", "MCR", color = "Mixture", 
                                                             palette = viridis(n = 13, direction = -1), title = "E", legend = "right", size = 3) + scale_x_log10() +
  geom_vline(xintercept = 0.01, linetype = "dashed")  + geom_hline(yintercept = 5, linetype = "dashed")+
  xlab("Toxicity Quotient") + ylab("Cumulative Ratio (min)")+ rremove("legend")
EARvMCRgroup_plot

ggsave("toxCast_LoE_C.jpeg",EARvMCRgroup_plot, height = 5, width = 7.5)

EARvMCRsingle_plot <-EARvMCR_single %>% filter(MCR < 10) %>% ggscatter("EAR", "MCR", color = "Chemical", 
                                                                     palette = viridis(n = 16, option = "magma", direction = -1), title = "F", legend = "right", size = 3) + scale_x_log10() +
  geom_vline(xintercept = 0.01, linetype = "dashed")  + geom_hline(yintercept = 5, linetype = "dashed")+
  xlab("Toxicity Quotient") + ylab("Cumulative Ratio (min)") + rremove("legend")
EARvMCRsingle_plot

ggsave("toxCast_LoE_D.jpeg",EARvMCRsingle_plot, height = 5, width = 7.5)

ToxCast_plot <- grid.arrange(ToxCast_MCR_plot_2017, ToxCast_MCR_chem_plot_2017, ToxCast_MCR_plot_2018,ToxCast_MCR_chem_plot_2018, EARvMCRgroup_plot, EARvMCRsingle_plot)
ToxCast_plot
ggsave("ToxCast_MCR.jpeg", ToxCast_plot, height = 15, width = 15)
