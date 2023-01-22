#MCR Wrangle#
library(tidyverse)
library(readxl)
library(writexl)
library(ggpubr)
library(gridExtra)
library(viridisLite)

setwd("C:\\Users\\erinm\\OneDrive\\Desktop\\GLRI Project\\Milwauke\\Data for Papier\\Mixtures MS\\ECOTOX Mixtures")

#read in chemistry
chem_2017 <- read_excel("C:\\Users\\erinm\\OneDrive\\Desktop\\GLRI Project\\Milwauke\\Data for Papier\\Chemistry\\Compiled_Chemistry\\Chemistry Compiled_Final.xlsx", "2017") %>% filter(Media == "Water", Detect == "Yes", !Site %in% c("UWM", "MIM-BK", "MED")) %>%
  select(Site, Year, Chemical, CAS, detect.con.ppb)
chem_2018 <- read_excel("C:\\Users\\erinm\\OneDrive\\Desktop\\GLRI Project\\Milwauke\\Data for Papier\\Chemistry\\Compiled_Chemistry\\Chemistry Compiled_Final.xlsx", "2018") %>% filter(Media == "Water", Detect == "Yes", !Site %in% c("UWM", "MIM-BK", "MED")) %>%
  select(Site, Year, Chemical, CAS, detect.con.ppb)
chem_2018$Year <- as.double(chem_2018$Year)

chemistry_compiled <- bind_rows(chem_2017, chem_2018)
list(unique(chemistry_compiled$Chemical))

chemistry_compiled$Site_Year <- paste(chemistry_compiled$Site, chemistry_compiled$Year, sep = "_")
chemistry_compiled$Concentration <- chemistry_compiled$detect.con.ppb / 10^3
View(chemistry_compiled)

single_chemical_list <- chemistry_compiled %>% select(CAS, Chemical) %>% distinct()
single_chemical_list$Group <- "Single Chemical"

#read in EC
names(EC_Group_ECOTOX)
EC_Group_ECOTOX <- read_excel("Final_Mixtures.xlsx") %>% select(CAS, "Chemical Name", "Final Mixture Group") %>% rename("Group" = "Final Mixture Group") %>%
  rename("Chemical" = "Chemical Name") %>% filter(!Group %in% c(NA, "NA - Unclustered or < no co-occurrence of mixture group"))
Final_ECOTOX_groupchem <- bind_rows(EC_Group_ECOTOX, single_chemical_list)

list(unique(EC_Group_ECOTOX$Group))
EC_for_MCR <- read_excel("milwaukee_EC_for_MCR_needsannotation_EM_annotated.xlsx") %>% select(-Chemical) 
names(EC_for_MCR)

ECOTOX_MCR_data <- left_join(Final_ECOTOX_groupchem, EC_for_MCR) %>% select(-Chemical)
View(ECOTOX_MCR_data)

#create MCR document + calculate TQ
MCR_doc <- left_join(chemistry_compiled, ECOTOX_MCR_data) %>% select(-Site, -Year, -detect.con.ppb)
MCR_doc$EC_estimate <- as.numeric(MCR_doc$EC_estimate)
MCR_doc$TQ <- MCR_doc$Concentration / MCR_doc$EC_estimate
View(MCR_doc)

TQ_site_year <- MCR_doc %>% filter(Group == "Single Chemical") %>% group_by(Site_Year) %>% summarize(sumTQ = sum(TQ, na.rm = TRUE))
View(TQ_site_year)

#create MCR document
MCR_doc1 <- left_join(MCR_doc, TQ_site_year, by = "Site_Year")

#regroup - if ndistinct chemical in group/site/year, replace with chemical (n> 2 in group = group)
n_distinct_chem <- MCR_doc1 %>% filter(Group != "Single Chemical") %>% group_by(Site_Year, Group) %>% summarize(n_distinct_CAS = n_distinct(CAS))
MCR_doc2 <- left_join(MCR_doc1, n_distinct_chem)
MCR_doc2$n_distinct_CAS[is.na(MCR_doc2$n_distinct_CAS)] <- 1
MCR_doc2$Group <- ifelse(MCR_doc2$n_distinct_CAS == 1, "Single Chemical", MCR_doc2$Group)

#carry out TQ for group
TQ_group <- MCR_doc2 %>% filter(Group != "Single Chemical") %>% group_by(Site_Year, Group) %>% summarize(sumTQ_group = sum(TQ))

MCR_doc3 <- left_join(MCR_doc2, TQ_group) %>% relocate("Group", .before = "CAS") %>% select(-n_distinct_CAS)

MCR_doc3$MCR_group <-MCR_doc3$sumTQ/ MCR_doc3$sumTQ_group 
MCR_doc3$MCR_chem <-MCR_doc3$sumTQ/ MCR_doc3$TQ 

MCR_doc3$MCR_Group_Annotation <- ifelse(MCR_doc3$MCR_group > 5, "> 5 (CR)",  MCR_doc3$Group) 
MCR_doc3$MCR_Chemical_Annotation <- ifelse(MCR_doc3$MCR_chem > 5, "> 5 (CR)",  MCR_doc3$Chemical) 

View(MCR_doc3)
write_xlsx(MCR_doc3, "ECOTOX_MCR_output.xlsx")

#create graph for ECOTOX only data
MCR_doc3$MCR_Group_Annotation[is.na(MCR_doc3$MCR_Group_Annotation)] <- "> 5 (CR)" 
MCR_doc3 <- MCR_doc3 %>% separate("Site_Year", into = c("Site", "Year"), sep = "_")

MCR_2017_graph <- MCR_doc3 %>% filter(Year == "2017")  %>% arrange((sumTQ)) %>% distinct() %>% rename("Structure-Mode of Action Group" = "MCR_Group_Annotation") %>% 
  ggbarplot("Site", "TQ", orientation = "horiz", legend = "right", color = "NA", title = "A", fill = "Structure-Mode of Action Group",
            palette =  c("#440154", "#3b528b")) + 
  xlab("Site (2017)") + ylab("Ecotoxicological Potential (C/96-h LC50)") +rremove("legend") + ylim(0, 0.12)

MCR_2018_graph <- MCR_doc3 %>% filter(Year == "2018")  %>% arrange((sumTQ)) %>% distinct() %>% rename("Structure-Mode of Action Group" = "MCR_Group_Annotation") %>% 
  ggbarplot("Site", "TQ", orientation = "horiz", legend = "right", color = "NA", title = "C", fill = "Structure-Mode of Action Group",
            palette =  viridis(n = 4)) + 
  xlab("Site (2018)") + ylab("Ecotoxicological Potential (C/96-h LC50)") +rremove("legend")

MCR_chem_2017 <- MCR_doc3 %>% filter(Year == "2017") %>% arrange((sumTQ)) %>% distinct() %>% select(-Chemical) %>% rename("Chemical" = "MCR_Chemical_Annotation") %>% 
  ggbarplot("Site", "TQ", orientation = "horiz", legend = "right", color = "NA", title = "B", fill = "Chemical",
            palette =  c("black", "#b73779", "#fc8961")) + 
  xlab("Site (2017)") + ylab("Ecotoxicological Potential (C/96-h LC50)") +rremove("legend") + ylim(0, 0.12)

MCR_chem_2018 <- MCR_doc3 %>% filter(Year == "2018") %>% arrange((sumTQ)) %>% distinct() %>% select(-Chemical) %>% rename("Chemical" = "MCR_Chemical_Annotation") %>% 
  ggbarplot("Site", "TQ", orientation = "horiz", legend = "right", color = "NA", title = "D", fill = "Chemical",
            palette =  viridis(n = 5, option = "magma")) + 
  xlab("Site (2018)") + ylab("Ecotoxicological Potential (C/96-h LC50)") +rremove("legend")

#generate TQ vs. MCR plot
TQvMCR_group <- MCR_doc3 %>% filter(Group != "Single Compound") %>% select(Group, MCR_group, sumTQ_group) %>%
  rename("Chemical/Mixture" = "Group", "MCR" = "MCR_group", "TQ" = "sumTQ_group") %>% distinct()

TQvMCR_single <- MCR_doc3 %>% select(Chemical, MCR_chem, TQ) %>%
  rename("Chemical" = "Chemical", "MCR" = "MCR_chem", "TQ" = "TQ") %>% distinct()

TQvMCR_data <- bind_rows(TQvMCR_group, TQvMCR_single) %>% group_by(`Chemical/Mixture`) %>% summarize(MCR = min(MCR), TQ = max(TQ))

TQvMCR_group_plot <-TQvMCR_group %>% filter(MCR < 5) %>% ggscatter("TQ", "MCR", color = "Chemical/Mixture", 
                                                             palette = c("#440154", "#3b528b"), title = "E", legend = "right", size = 3) + scale_x_log10() + 
 geom_hline(yintercept = 5, linetype = "dashed") +
  geom_vline(xintercept = 0.01, linetype = "dashed") + ylim(0, 10) + 
  xlab("Toxic Units (max)") + ylab("Cumulative Ratio (min)") +rremove("legend")

ggsave("ECOTOX_LoE_C.jpeg",TQvMCR_group_plot, height = 5, width = 7.5)


TQvMCR_chem_plot <-TQvMCR_single %>% filter(MCR < 5) %>% ggscatter("TQ", "MCR", color = "Chemical", 
                                                                   palette = c("#51127c", "#b73779", "#fc8961"), title = "F", legend = "right", size = 3) + scale_x_log10() + 
  geom_hline(yintercept = 5, linetype = "dashed") +
  geom_vline(xintercept = 0.01, linetype = "dashed") + ylim(0, 10) + 
  xlab("Toxic Units (max)") + ylab("Cumulative Ratio (min)") +rremove("legend")

ggsave("ECOTOX_LoE_D.jpeg",TQvMCR_chem_plot, height = 5, width = 7.5)

MCR <- grid.arrange(MCR_2017_graph, MCR_chem_2017, MCR_2018_graph, MCR_chem_2018, TQvMCR_group_plot, TQvMCR_chem_plot, nrow = 3)
MCR
ggsave("ECOTOX_MCR.jpeg", MCR, height = 15, width = 15)

