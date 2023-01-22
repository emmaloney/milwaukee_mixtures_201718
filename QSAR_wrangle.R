#QSAR MOA data wrangle
library(readxl)
library(writexl)
library(tidyverse)
library(gridExtra)

ecotox_EC <- read_excel("ECOTOX_EC.xlsx")

#Aster
aster <- read_excel("ASTER_output.xlsx") %>% select(Input, MOA) %>% rename("CAS" = "Input")
names(aster)
View(aster)

ecotox_CAS <- ecotox_EC %>% select(CAS)

QSAR <- left_join(ecotox_CAS, aster, by = "CAS") %>% rename("ASTER_MOA" = "MOA")
QSAR[is.na(QSAR)] <- "N.A."
View(QSAR)

#vega
vega <- read_excel("VEGA_Class.xlsx") %>% select(-MOA_Verhaar_Prediction, -MOA_Verhaar_Assessment)
names(vega)

QSAR1 <- left_join(QSAR, vega)

#oasis
oasis.v2 <- read_excel("OASIS_output.v2.csv.xlsx")
names(oasis.v2)

QSAR2 <- left_join(QSAR1, oasis.v2) %>% select(-c(MOA_TEST_Prediction, MOA_IRFMN_Prediction, MOA_OASIS))
names(QSAR2)

#double check that VEGA and OASIS Verhaar classifications are the same

#wrangle
list(unique(QSAR2$MOA_TEST_Assessment))
QSAR2$MOA_TEST_Assessment <- ifelse(QSAR2$MOA_TEST_Assessment %in% c("Narcosis (good reliability)",
                                                                     "Narcosis (low reliability)",
                                                                     "Narcosis (moderate reliability)",
                                                                     "Narcosis (EXPERIMENTAL value)"),
                                    "Narcosis", QSAR2$MOA_TEST_Assessment)
QSAR2$MOA_TEST_Assessment <- ifelse(QSAR2$MOA_TEST_Assessment %in% c("Reactivity (EXPERIMENTAL value)",
                                                                     "Reactivity (low reliability)",
                                                                     "Reactivity (moderate reliability)",
                                                                     "Reactivity (good reliability)"),
                                    "Reactive", QSAR2$MOA_TEST_Assessment)
QSAR2$MOA_TEST_Assessment <- ifelse(QSAR2$MOA_TEST_Assessment %in% c("N/A", "NA"), "N.A.",
                                    QSAR2$MOA_TEST_Assessment)
QSAR2$MOA_TEST_Assessment <- ifelse(QSAR2$MOA_TEST_Assessment %in% c("Neurotoxicity (EXPERIMENTAL value)",
                                                                     "Neurotoxicity (good reliability)",
                                                                     "Neurotoxicity (low reliability)"), "Neurotoxic",
                                    QSAR2$MOA_TEST_Assessment)
QSAR2$MOA_TEST_Assessment <- ifelse(QSAR2$MOA_TEST_Assessment %in% c("Uncoupler (low reliability)",
                                                                     "Uncoupler (moderate reliability)",
                                                                     "Uncoupler (EXPERIMENTAL value)",
                                                                     "Uncoupler (good reliability)"), "Uncoupler",
                                    QSAR2$MOA_TEST_Assessment)
QSAR2$MOA_TEST_Assessment <- ifelse(QSAR2$MOA_TEST_Assessment %in% c("AChE Inhibition (low reliability)",
                                                                     "AChE Inhibition (EXPERIMENTAL value)",
                                                                     "AChE Inhibition (moderate reliability)",
                                                                     "AChE Inhibition (good reliability)"),
                                                                     "AChE Inhibitor",
                                    QSAR2$MOA_TEST_Assessment)
QSAR2$MOA_TEST_Assessment <- ifelse(QSAR2$MOA_TEST_Assessment %in% c("nACHr Antagonism (low reliability)"),
                                    "nAChr antagonist",
                                    QSAR2$MOA_TEST_Assessment)
QSAR2$MOA_TEST_Assessment <- ifelse(QSAR2$MOA_TEST_Assessment %in% c("Anticoagulation (low reliability)"),
                                    "Anticoagulant",
                                    QSAR2$MOA_TEST_Assessment)

names(QSAR2)

#add qualifiers
list(unique(QSAR2$Verhaar_Class))
QSAR2$ASTER_Qual <- ifelse(QSAR2$ASTER_MOA %in% c("Nonpolar narcosis", "Polar Narcosis", "Ester narcosis"), "N",
                           ifelse(QSAR2$ASTER_MOA %in% c("N.A."), "U",
                                  "S"))
QSAR2$TEST_Qual <- ifelse(QSAR2$MOA_TEST_Assessment %in% c("Narcosis"), "N",
                                    ifelse(QSAR2$MOA_TEST_Assessment %in% c("N.A."), "U",
                                           "S"))
QSAR2$OASIS_Qual <- ifelse(QSAR2$OASIS_MOA %in% c("Basesurface narcotics", "Narcotic Amine",
                                                             "Phenols and Anilines", "Alpha-, Beta-unsaturated alcochols",
                                                             "Esters"), "N",
                           ifelse(QSAR2$OASIS_MOA %in% c("Aldehydes", "Reactive unspecified"), "S",
                                  "U")) 
QSAR2$Verhaar_Qual <- ifelse(QSAR2$Verhaar_Class %in% c("Class 1 (narcosis or baseline toxicity)", 
                                                                 "Class 2 (less inert compounds)"), "N",
                             ifelse(QSAR2$Verhaar_Class %in% c("Class 3 (unspecific reactivity)",
                                                                        "Class 4 (compounds and groups of compounds acting by a specific mechanism)"),
                                    "S", "U"))

QSAR2$CONSENSUS_MOA <- paste(QSAR2$ASTER_Qual, QSAR2$TEST_Qual, QSAR2$OASIS_Qual, QSAR2$Verhaar_Qual, sep = "")
View(QSAR2)

list(unique(QSAR2$CONSENSUS_MOA))

QSAR2$MOA_Classification <- ifelse(QSAR2$CONSENSUS_MOA %in% c("SSUU", "SSNS", "SSSN", "SSSS", "SSSU", "USSS", "SNSS",
                                                              "NSSS", "USSU"), "S",
                                   ifelse(QSAR2$CONSENSUS_MOA %in% c("UNNU", "NNNU", "NNNN", "UNNN", "NNUU",
                                                                     "NSNN", "NNSN", "NNNS", "NNUN",
                                                                     "NUNN", "NUNU"), "N",
                                          "U"))
QSAR2$Confidence_Score <- ifelse(QSAR2$CONSENSUS_MOA %in% c("NNNN", "SSSS"), 3,
                                 ifelse(QSAR2$CONSENSUS_MOA %in% c("SSNS", "NNNU", "UNNN", "SSSN", "NSNN", "SSSU",
                                                                   "USSS", "SNSS", "NNSN", "NNNS", "NNUN",
                                                                   "NUNN", "SNNN"), 2,
                                        ifelse(QSAR2$CONSENSUS_MOA %in% c("SSUU", "UNNU", "NNUU", "SUSU",
                                                                          "USSU", "NUNU"), 1,
                                               "0")))

#pull in Milwaukee chemicals
milwaukee_consensus <- read_excel("Milwaukee_chem_QSAR_Consensus.xlsx")

missing_milwaukee <- anti_join((milwaukee_consensus), (QSAR2), by = "CAS")
View(missing_milwaukee)

#rename and bind
names(QSAR3)
names(missing_milwaukee)

QSAR3 <- QSAR2 %>% rename("TEST_MOA" = "MOA_TEST_Assessment", "IRFMN_MOA" = "MOA_IRFMN_Assessment", "Verhaar_MOA" = "Verhaar_Class")
missing_milwaukee1 <- missing_milwaukee %>% select(-c("Chemical_Name", "QSAR_Classification", "QSAR_Class_detail")) %>%
  rename("ASTER_MOA" = "ASTER", "ASTER_Qual" = "ASTER_code", "TEST_MOA" = "TEST", "TEST_Qual" = "TEST_code",
         "OASIS_MOA" = "OASIS", "OASIS_Qual" = "OASIS_code", "Verhaar_MOA" = "Verhaar", "Verhaar_Qual" = "Verhaar_code",
         "CONSENSUS_MOA" = "Compiled_Code", "MOA_Classification" = "Consensus_MOA", "Confidence_Score" = "MOA_Confidence_Score")

QSAR3$Confidence_Score <- as.double(QSAR3$Confidence_Score)

whole_chemical_set <- bind_rows(QSAR3, missing_milwaukee1)
View(whole_chemical_set)

write_xlsx(whole_chemical_set, "Final_Consensus_MOA.xlsx")
