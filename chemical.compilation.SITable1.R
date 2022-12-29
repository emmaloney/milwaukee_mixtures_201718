#collating long list of detected concentrations for SI
library(dplyr)
library(tidyr)
library(writexl)
library(readxl)
getwd()
setwd("C:\\Users\\erinm\\OneDrive\\Desktop\\GLRI Project\\Milwauke\\Data for Papier\\Chemistry")
chem.2017 <- read_excel("C:\\Users\\erinm\\OneDrive\\Desktop\\GLRI Project\\Milwauke\\Data for Papier\\Chemistry\\Compiled_Chemistry\\Chemistry Compiled_Final.xlsx", "2017")
chem.2018 <- read_excel("C:\\Users\\erinm\\OneDrive\\Desktop\\GLRI Project\\Milwauke\\Data for Papier\\Chemistry\\Compiled_Chemistry\\Chemistry Compiled_Final.xlsx", "2018")
View(chem.2018)
names(chem.2017)
chem.2017.1 <- chem.2017 %>% filter(!Site %in% c("UWM", "MIM-BK")) %>% select(Site, detect.con.ppb, Chemical, CAS, Media, Detect.code, "RL.ppb") %>% rename("Chemical Name" = "Chemical", "Environmental Compartment" = "Media", "Concentration" = "detect.con.ppb", "Study Site" = "Site", "Detect Code" = "Detect.code", "Reporting Limit" = "RL.ppb") %>%
  relocate(CAS, .after = `Study Site`) %>% relocate(Concentration, .after = `Environmental Compartment`) %>% relocate("Detect Code", .before = "Concentration")
View(chem.2017.1)
chem.2017.1$`Study Year` <- "2017"
chem.2017.1$`Detection Units` <- "ug/L"

chem.2017.2 <- chem.2017.1 %>% relocate(`Study Year`, .before = `Study Site`)
View(chem.2017.2)

names(chem.2018)
chem.2018.1 <- chem.2018 %>% select(Site, detect.con.ppb, Chemical, CAS, Media, Detect.code, "RL.ppb") %>% rename("Chemical Name" = "Chemical", "Environmental Compartment" = "Media", "Concentration" = "detect.con.ppb", "Study Site" = "Site", "Detect Code" = "Detect.code", "Reporting Limit" = "RL.ppb") %>%
  relocate(CAS, .after = `Study Site`) %>% relocate(Concentration, .after = `Environmental Compartment`) %>% relocate("Detect Code", .before = "Concentration")
View(chem.2018.1)
chem.2018.1$`Study Year` <- "2018"
chem.2018.1$`Detection Units` <- "ug/L"

chem.2018.2 <- chem.2018.1 %>% relocate(`Study Year`, .before = `Study Site`)
View(chem.2018.2)

all.chem <- bind_rows(chem.2017.2, chem.2018.2) %>% filter(`Study Site` != "UWM")
View(all.chem)

#bind in site names (long)
SI_table1 <- read_excel("C:\\Users\\erinm\\OneDrive\\Desktop\\GLRI Project\\Milwauke\\Data for Papier\\Chemistry\\site_longnam.xlsx")
names(SI_table1)
site_name <- SI_table1 %>% select("Site Code", "USGS Station ID (STAID)", "Site Name") %>% rename("Study Site" = "Site Code")

all_chem_semifin <- left_join(all.chem, site_name) %>% relocate(c("USGS Station ID (STAID)":"Site Name"), .before = "Study Site")

names(all_chem_semifin)

all_chem_semifin$`Sample Type` <- ifelse(all_chem_semifin$`Study Year` == "2018" & all_chem_semifin$`Study Site` %in% c("MIE", "JIP"), "Grab Sample (end of study)", "Composite Sample")
all_chem_semifin$`Sample Notes` <- ifelse(all_chem_semifin$`Study Year` == "2017" & all_chem_semifin$`Environmental Compartment` == "Water" & all_chem_semifin$`Study Site` %in% c("KKL", "MET", "MIM", "MIP"), "Replacement Sample used Due to Shipping Issues", "")

View(all_chem_semifin)

write_xlsx(all_chem_semifin, "SI_Table_2_Milwaukee_Prioritization.xlsx")

#Add extra for SI for Milwaukee Mixtures MS

CAS <- all_chem_semifin %>% select(CAS)
write_xlsx(CAS, "CAS_for_Mixtures_SI.xlsx")

CAS_MW <- read_excel("CAS_MW.xlsx")

all_chem_semifin_1 <- left_join(all_chem_semifin, CAS_MW) %>% filter(`Environmental Compartment` == "Water")
names(all_chem_semifin_1)
list(unique(all_chem_semifin_1$`Detection Units`))

all_chem_semifin_1$`Concentration (mol/L)` <- (all_chem_semifin_1$Concentration /10^6)/(all_chem_semifin_1$MW)
all_chem_semifin_1$`Concentration (mol/L)`  <- as.numeric(all_chem_semifin_1$`Concentration (mol/L)` )
View(all_chem_semifin_1)

all_chem_semifin_1$MW <- ifelse(is.na(all_chem_semifin_1$MW), "-", all_chem_semifin_1$MW)
all_chem_semifin_1$`Concentration (mol/L)` <- ifelse(is.na(all_chem_semifin_1$`Concentration (mol/L)`), 0, all_chem_semifin_1$`Concentration (mol/L)`)

all_chem_semifin_1_site_con <- all_chem_semifin_1 %>% group_by(`Study Site`, `Study Year`) %>% summarize(`Cumulative Concentration (Site; mol/L)` = sum(`Concentration (mol/L)`, na.rm = TRUE))
View(all_chem_semifin_1_site_con)
names(all_chem_semifin_1_site_con)

all_chem_semifin_2 <- left_join(all_chem_semifin_1, all_chem_semifin_1_site_con)
all_chem_semifin_2$Percent_contribution <- (all_chem_semifin_2$`Concentration (mol/L)`/all_chem_semifin_2$`Cumulative Concentration (Site; mol/L)`) * 100
all_chem_semifin_2$Concentration_Contribution_Note <- ifelse(all_chem_semifin_2$Percent_contribution >= 5, "Important Contributor", "")

#add in reporting limit for each chemical!!!####


write_xlsx(all_chem_semifin_2, "SI_Table2_Milwaukee_Mixtures.xlsx")

#summary for text
names(all_chem_semifin_1)
list(unique(all_chem_semifin_summary$`Detect Code`))
all_chem_semifin_summary <- all_chem_semifin_1 %>% group_by(`Study Year`, `Study Site`) %>%
  summarize(n_measured_chem = n_distinct(CAS))
all_chem_detect_summary <- all_chem_semifin_1 %>% group_by(`Study Year`, `Study Site`) %>% filter(is.na(`Detect Code`) | `Detect Code` %in% c("E", "M")) %>%
  summarize(n_measured_chem = n_distinct(CAS))
all_chem_detect_year <- all_chem_semifin_1 %>% group_by(`Study Year`) %>% filter(is.na(`Detect Code`) | `Detect Code` %in% c("E", "M")) %>%
  summarize(n_measured_chem = n_distinct(CAS))

View(all_chem_semifin_summary)
View(all_chem_detect_summary %>% filter(`Study Year` == "2017"))

View(all_chem_semifin_summary)

all_chem_detect_summary_avg <- all_chem_detect_summary %>% group_by(`Study Year`) %>% summarize(mean_n = mean(n_measured_chem), sum_n = sum(n_measured_chem))
all_chem_detect_summary_avg
