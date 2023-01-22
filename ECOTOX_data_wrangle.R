#ECOTOX data wrangle

library(readxl)
library(writexl)
library(tidyverse)
library(gridExtra)

setwd("C:\\Users\\erinm\\OneDrive\\Desktop\\GLRI Project\\Milwauke\\Data for Papier\\Mixtures MS\\ECOTOX Mixtures")
ecotox_data <- read_excel("ECOTOX_export_for_readacross.xlsx")
CAS <- read_excel("96hLC50_for_Aster.xls") %>% select(-DTXSID)
names(CAS)

ecotox_v1 <- left_join(ecotox_data, CAS) %>% select(-c("BCF 1 Value Op":"BCF 3 Unit",
                                                    "CAS Number", "Chemical Name"))%>%
  rename("Chemical" = "PREFERRED_NAME", "CAS" = "CASRN") %>% relocate("Chemical", .before = "Chemical Grade") %>%
  relocate("CAS", .before = "Chemical")
n_chem <- ecotox_v1 %>% summarize(n_distinct(CAS)) #1282 distinct chemicals

ecotox_v1$`Chemical Purity` <- as.numeric(ecotox_v1$`Chemical Purity`)
names(ecotox_v1)
View(ecotox_v1)

#based on OECD standards for 96-h fish assessment#
#filter to exclude unmeasured LC50s - double check to see how much overlap there is - include nominal to expand data availability
#filter to include only static, renewal, or flow-through exposures
#filter to include only fresh water exposures
#filter to have >= 5 doses
#filter to have chemical purity > 70 %
#exclude formulated products to avoid unintentionally including mixture assessments 
#excluded egg measurements (included > juvenile stages b/c only juvenile was too data limiting)


ecotox_v2 <- ecotox_v1 %>% filter(`Chemical Purity` > 70 | is.na(`Chemical Purity`), 
                                  !`Chemical Analysis` %in% c("Unmeasured", "Unmeasured values (some measured values reported in article)"),
                                  `Media Type` == "Fresh water", `Conc 1 Type (Standardized)` != "Formulation",
                                  !`Organism Lifestage` %in% c("Egg"), !`Number of Doses` %in% c("1", "2","3", "4"),
                                  !`Conc 1 Mean Op (Standardized)` %in% c(">","<","~"))
n_chem <- ecotox_v2 %>% summarize(n_distinct(CAS)) #767 distinct chemicals 

#parameters
list(unique(ecotox_v2$`Chemical Purity`)) # NA or > 70 % purity - will add 'formulated' only for compounds that are not pesticides (e.g. PAHs)
list(unique(ecotox_v2$`Chemical Analysis`)) #measured
list(unique(ecotox_v2$`Media Type`)) #fresh water
list(unique(ecotox_v2$`Conc 1 Type (Standardized)`)) #fresh water
list(unique(ecotox_v2$`Conc 1 Units (Standardized)`)) # AI mg/L
list(unique(ecotox_v2$`Organism Lifestage`)) # remove egg
list(unique(ecotox_v2$`Number of Doses`)) # remove < 5 doses 
list(unique(ecotox_v2$`Conc 1 Mean Op (Standardized)`)) # remove operations
list(unique(ecotox_v2$`Summary of Additional Parameters`)) # remove operations


#pull out chemicals missing in unformulated data####
unformulated <- ecotox_v1 %>% filter(`Chemical Purity` > 70 | is.na(`Chemical Purity`), 
                                     !`Chemical Analysis` %in% c("Unmeasured", "Unmeasured values (some measured values reported in article)"),
                                     `Media Type` == "Fresh water", `Conc 1 Type (Standardized)` != "Formulation",
                                     !`Organism Lifestage` %in% c("Egg"), !`Number of Doses` %in% c("1", "2","3", "4"),
                                     !`Conc 1 Mean Op (Standardized)` %in% c(">","<","~"))

formulated <- ecotox_v1 %>% filter(`Chemical Purity` > 70 | is.na(`Chemical Purity`), 
                                   !`Chemical Analysis` %in% c("Unmeasured", "Unmeasured values (some measured values reported in article)"),
                                   `Media Type` == "Fresh water", `Conc 1 Type (Standardized)` == "Formulation",
                                   !`Organism Lifestage` %in% c("Egg"), !`Number of Doses` %in% c("1", "2","3", "4"),
                                   !`Conc 1 Mean Op (Standardized)` %in% c(">","<","~"))

missing_in_unformulated <- anti_join(formulated, unformulated, by = "CAS") #nada 
View(missing_in_unformulated)

#add in formulated missing chemicals that meet other specifications####
ecotox_v3 <- bind_rows(ecotox_v2, missing_in_unformulated)
n_chem <- ecotox_v3 %>% summarize(n_distinct(CAS)) #771 distinct chemicals 


#pull out chemicals missing in nominal data#### 
measured <- ecotox_v1 %>% filter(`Chemical Purity` > 70 | is.na(`Chemical Purity`), 
                                     !`Chemical Analysis` %in% c("Unmeasured", "Unmeasured values (some measured values reported in article)"),
                                     `Media Type` == "Fresh water", `Conc 1 Type (Standardized)` != "Formulation",
                                     !`Organism Lifestage` %in% c("Egg"), !`Number of Doses` %in% c("1", "2","3", "4"),
                                     !`Conc 1 Mean Op (Standardized)` %in% c(">","<","~"))

nominal <- ecotox_v1 %>% filter(`Chemical Purity` > 70 | is.na(`Chemical Purity`), 
                                   `Chemical Analysis` %in% c("Unmeasured", "Unmeasured values (some measured values reported in article)"),
                                   `Media Type` == "Fresh water", `Conc 1 Type (Standardized)` != "Formulation",
                                   !`Organism Lifestage` %in% c("Egg"), !`Number of Doses` %in% c("1", "2","3", "4"),
                                   !`Conc 1 Mean Op (Standardized)` %in% c(">","<","~"))

missing_in_measured <- anti_join(nominal, measured, by = "CAS") 
View(missing_in_measured)

#add in unmeasured chemicals that meet other specifications####
ecotox_v4 <- bind_rows(ecotox_v3, missing_in_measured) %>% filter(!is.na(`Conc 1 Mean (Standardized)`))
n_chem <- ecotox_v4 %>% summarize(n_distinct(CAS)) #873 distinct chemicals 

#QAQC data####
#generate boxplots to identify outliers####
names(ecotox_v4)
ecotox_v4$`Conc 1 Mean (Standardized)` <- as.numeric(ecotox_v4$`Conc 1 Mean (Standardized)`)

boxplot <- boxplot(ecotox_v4$`Conc 1 Mean (Standardized)` ~ ecotox_v4$CAS)
boxplot$out

outliers <- subset(ecotox_v4, ecotox_v4$`Conc 1 Mean (Standardized)` %in% boxplot$out)
View(outliers)

boxplot_info <- ecotox_v4 %>% group_by(CAS) %>% summarise(Q1 = quantile(`Conc 1 Mean (Standardized)`, 0.25), Q3 = quantile(`Conc 1 Mean (Standardized)`, 0.75))

outliers1 <- left_join(outliers, boxplot_info, by = "CAS")
View(outliers1)

outliers1$outlier_class <- ifelse(outliers1$`Conc 1 Mean (Standardized)` < outliers1$Q1, "potential_lower_outlier", "potential_upper_outlier")

outliers2 <- outliers1 %>% filter(outlier_class == "potential_lower_outlier")
View(outliers2)

##identify lower outliers for data limited compounds####
summary <- ecotox_v4 %>% group_by(CAS) %>% summarize(min2_conc = nth(`Conc 1 Mean (Standardized)`, 2, order_by = `Conc 1 Mean (Standardized)`), min_conc = min(`Conc 1 Mean (Standardized)`), n_EC = length(`Conc 1 Mean (Standardized)`), n_Study = length(unique(Title)))

min_outliers <- left_join(ecotox_v4, summary)
min_outliers$min_min2 <- min_outliers$min2_conc / min_outliers$min_conc

min_outliers$n_EC <- as.numeric(min_outliers$n_EC)
min_outliers$min_min2 <- as.numeric(min_outliers$min_min2)
min_outliers$flag_EC <- ifelse(min_outliers$n_EC < 5 & min_outliers$min_min2 >= 10 | min_outliers$n_EC < 5 & is.na(min_outliers$min2_conc), "Potential Lower Outlier", "")
min_outliers$flag_Study <- ifelse(min_outliers$n_Study < 5 & min_outliers$min_min2 >= 10 | min_outliers$n_Study < 5 & is.na(min_outliers$min2_conc), "Potential Lower Outlier", "")

min_outliers1 <- min_outliers %>% filter(flag_EC == "Potential Lower Outlier" | flag_Study == "Potential Lower Outlier")

#bind 
outliers3 <- bind_rows(outliers2, min_outliers1) %>% distinct()
View(outliers3)

write_xlsx(outliers3, "ECOTOX_outliers.xlsx")

#write out potential upper outliers
outliers_upper <- outliers1 %>% filter(outlier_class == "potential_upper_outlier")
write_xlsx(outliers_upper, "ECOTOX_upper_outliers.xlsx")


#read in outlier chemicals ####

ECOTOX_outliers <- read_excel("ECOTOX_outliers_EM_annotated.xlsx", 1)
ECOTOX_outliers2 <- read_excel("ECOTOX_upper_outliers_EM_annotated.xlsx", 1)
ECOTOX_outliers2$`Chemical Purity` <- as.character(ECOTOX_outliers2$`Chemical Purity`)
ECOTOX_outliers <- bind_rows(ECOTOX_outliers, ECOTOX_outliers2)

names(ECOTOX_outliers)
exclude <- ECOTOX_outliers %>% filter(Action == "Exclude") %>% select(-c("Q1":"...56"))
exclude$`Chemical Purity` <- as.double(exclude$`Chemical Purity`)
names(exclude)
ecotox_v5 <- anti_join(ecotox_v4, exclude)

replace <- read_excel("ECOTOX_outliers_EM_annotated.xlsx", 2) %>% select(-c("Q1": "Action"))
replace2 <- read_excel("ECOTOX_upper_outliers_EM_annotated.xlsx", 2)%>% select(-c("Q1": "Action"))
replace <- bind_rows(replace, replace2)
names(replace)
replace$`Chemical Purity` <- as.double(replace$`Chemical Purity`)

ecotox_v6 <- bind_rows(ecotox_v5, replace)
names(ecotox_v6)

gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
ecotox_v7 <- ecotox_v6 %>% group_by(Chemical, CAS) %>% summarize(Concentration = gm_mean(`Conc 1 Mean (Standardized)`)) %>% filter(!is.na(CAS))
View(ecotox_v7)

write_xlsx(ecotox_v7, "ECOTOX_EC.xlsx")

