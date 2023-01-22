#SEQAPASS meddling#

library(dplyr)
library(tidyverse)
library(readxl)
library(writexl)
library(ggplot2)
library(ggpubr)
library(stringr)

#read in L1 & modify file
list.files(path = "SEQAPASS_Outputs\\L1", pattern = ".csv")
raw.files <- data_frame(filename = list.files(path = "SEQAPASS_Outputs\\L1", pattern = ".csv"))
raw.file.paths <- raw.files %>% mutate(filepath = paste0("SEQAPASS_Outputs\\L1\\", filename))
raw.file.paths %>% head(3)

read.csv.and.add.filename <- function(filepath){
  read_csv(filepath, col_types = cols_only(`Taxonomic Group` = "c", `Scientific Name` = "c", `Common Name` = "c", `Susceptibility Prediction` = "c")) %>%
    mutate(filepath = filepath)}

raw.data <- raw.file.paths %>% rowwise() %>% do(., read.csv.and.add.filename(file = .$filepath)) %>% 
  separate(filepath, c("SEQAPASS", "Outputs", "Level", "Primary", "Report", "gene"), sep = "_") %>% select(-c("SEQAPASS", "Outputs", "Level", "Primary", "Report"))
raw.data$gene <- gsub(".csv", "", raw.data$gene)
raw.data$Level <- "L1"
View(raw.data)

list(unique(raw.data$gene))

#read in L2 & modify file 
list.files(path = "SEQAPASS_Outputs\\L2", pattern = ".csv")
raw.files.L2 <- data_frame(filename = list.files(path = "SEQAPASS_Outputs\\L2", pattern = ".csv"))
raw.file.paths.L2 <- raw.files.L2 %>% mutate(filepath = paste0("SEQAPASS_Outputs\\L2\\", filename))
raw.file.paths.L2 %>% head(3)

read.csv.and.add.filename <- function(filepath){
  read_csv(filepath, col_types = cols_only(`Taxonomic Group` = "c", `Scientific Name` = "c", `Common Name` = "c", `Susceptibility Prediction` = "c")) %>%
    mutate(filepath = filepath)}

raw.data.L2 <- raw.file.paths.L2 %>% rowwise() %>% do(., read.csv.and.add.filename(file = .$filepath)) %>% 
  separate(filepath, c("SEQAPASS", "Outputs", "Level", "Primary", "Report", "gene"), sep = "_") %>% select(-c("SEQAPASS", "Outputs", "Level", "Primary", "Report"))
raw.data.L2$gene <- gsub(".csv", "", raw.data.L2$gene)
raw.data.L2$Level <- "L2"
View(raw.data.L2)
list(unique(raw.data.L2$gene))

#bind together levels 
full_pharm_seqapass <- bind_rows(raw.data, raw.data.L2)
names(full_pharm_seqapass)
list(unique(full_pharm_seqapass$`Susceptibility Prediction`))

susceptible <- full_pharm_seqapass %>% filter(`Susceptibility Prediction` == "Y") %>% rename("Tax_group" = "Taxonomic Group")
susceptible$general_taxa <- ifelse(susceptible$Tax_group %in% c("Mammalia", "Aves" , "Lepidosauria", "Testudines", "Crocodylia", "Leptocardii", "Ascidiacea", "Aconoidasida", "Trichoplacidae", "Rotosphaerida", "Choanoflagellata", "Tentaculata", "Hexactinellida",  "Amphibia", "Appendicularia", "Rotaliida", "Heterolobosea sp. OSA",
                                                                "Heterolobosea sp. BA", "Telonemida", "Cristamonadida", "Hypotrichomonadida", "Creneidae", "Eukaryota sp. ATCC 50646", "Ancyromonadida", "Acrasidae", "Hicanonectes", "Placozoa sp. H4", "Dicyemidae", "Retortamonadidae", "Hicanonectes", "Percolomonas", "Psalteriomonadidae"), "Other_Eukaryotes",
                                   ifelse(susceptible$Tax_group %in% c("Actinopteri", "Chondrichthyes", "Cladistia", "Coelacanthimorpha", "Hyperoartia", "Myxini", "Dipnoi"), "Fish",
                                          ifelse(susceptible$Tax_group %in% c("Bivalvia", "Scyphozoa", "Asteroidea", "Crinoidea", "Holothuroidea", "Echinoidea", "Insecta", "Collembola", "Merostomata", "Monogononta", "Malacostraca", "Gastropoda", "Arachnida", "Cephalopoda", "Anthozoa", "Clitellata", 
                                                                              "Hexanauplia", "Chromadorea", "Polychaeta", "Hydrozoa", "Trematoda", "Cestoda", "Branchiopoda", "Enoplea", "Rhabditophora", "Solenogastres", 
                                                                              "Caudofoveata", "Demospongiae", "Myxozoa", "Diplopoda", "Pauropoda", "Chilopoda", "Pycnogonida", "Diplura", "Symphyla", "Monogenea", "Bdelloidea", "Rotaliida", "Caudofoveata", "Solenogastres", "Protura", "Enteropneusta",
                                                                              "Eutardigrada", "Palaeonemertea", "Rhopaluridae", "Priapulimorpha", "Enopla", "Udeonychophora", "Palaeonemertea", "Pilidiophora", "Acoela", "Nemertodermatida", "Rhynchonellata", "Dimorpha"), "Invertebrates", 
                                                 ifelse(susceptible$Tax_group %in% c("Magnoliopsida", "Sphagnopsida", "Marchantiopsida", "Bryopsida", "Anthocerotopsida", "Leiosporocerotopsida", "Polypodiopsida", "Polytrichopsida","Andreaeobryopsida",
                                                                                     "Jungermanniopsida", "Lycopodiopsida", "Andreaeopsida","Haplomitriopsida", "Ginkgoopsida", "Pinopsida", "Gnetopsida", "Cycadopsida", "Picocystophyceae", "Tetraphidopsida", 
                                                                                     "Takakiopsida", "Oedipodiopsida", "Chlorophyta sp. W-c", "Prasinophyte sp. SL-175", "Prasinophyte sp. RS-11"), "Other_Plants", 
                                                        ifelse(susceptible$Tax_group %in% c("Lingulata","Basidiobolomycetes", "Mucoromycetes", "Endogonomycetes", "Glomeromycetes", "Chytridiomycetes", "Zoopagomycetes",
                                                                                            "Saccharomycetes", "Mortierellomycetes", "Monoblepharidomycetes", "Tremellomycetes", "Agaricomycetes", "Harpellomycetes", "Taphrinomycetes", "Sordariomycetes", "Dacrymycetes",
                                                                                            "Eurotiomycetes", "Kickxellomycetes", "Leotiomycetes", "Pucciniomycetes", "Dothideomycetes", "Ustilaginomycetes", "Malasseziomycetes", "Fungal sp. No.11243", "Exobasidiomycetes", "Pezizomycetes",
                                                                                            "Neocallimastigomycetes", "Entomophthoromycetes", "Dimargaritomycetes", "Schizosaccharomycetes", "Saitoella", "Orbiliomycetes", "Wallemiomycetes", "Fungal sp. No.14919", 
                                                                                            "Neolectomycetes", "Pneumocystidomycetes", "Rozella", "Microbotryomycetes", "Paramicrosporidium", "Blastocladiomycetes", "Mixiomycetes", "Lecanoromycetes", "Xylonomycetes", "Nematocida",
                                                                                            "Mitosporidium", "Pansporoblastina", "Tubulinosematoidea", "Apansporoblastina", "Amphiacanthidae", "Pseudoloma", "Ordosporidae", "Culicosporidae", "Mucoromycetes", "Endogonomycetes",
                                                                                            "Zoopagomycetes", "Mortierellomycetes", "Harpellomycetes", "Taphrinomycetes", "Dacrymycetes", "Malasseziomycetes", "Fungal sp. No.11243", "Exobasidiomycetes", "Schizosaccharomycetes", "Saitoella",
                                                                                            "Orbiliomycetes", "Wallemiomycetes", "Fungal sp. No.14919", "Neolectomycetes", "Pneumocystidomycetes", "Rotosphaerida", "Microbotryomycetes", "Paramicrosporidium", "Endoreticulatus", "Asellariomycetes", "Hormodochis",
                                                                                            "Ramicandelaberales", "Humicolopsis", "Enteropsectra", "Lichinomycetes", "Uncultured fungus", "Arthoniomycetes", "Calcarea"), "Fungi/Metazoa",
                                                               ifelse(susceptible$Tax_group %in% c("Diplonemea", "Kipferlia", "Ichthyosporea", "Longamoebia", "Euglyphida", "Variosea", "Eumycetozoa", "Jakobida", "Oxymonadida", "Diplomonadida", "Tritrichomonadida", "Trichomonadida", "Euglyphida",
                                                                                                   "Vahlkampfiidae", "Filasterea", "Euglenida", "Kinetoplastea", "Apusomonadidae", "Mastigamoebida", "Heterotrichea", "Carpediemonas", "Tsukubamonadida", "Trimastigidae", "Dysnectes", "Cercomonadida", 
                                                                                                   "Malawimonadidae", "Trichonymphida", "Spirotrichonymphida", "Antonospora", "Ascetosporea", "Barbatosporales", "Thaumatomonadida", "Spongomonadida", "Flabellinia", "Thecofilosea", "Allogromia", 
                                                                                                   "Gymnophrys", "Astrorhizida", "Glissomonadida", "Rzehakinidae", "Elardia", "Stephanopogonidae", "Centroplasthelida", "Corycida", "Corallomyxa", "Ergobibamus", "Chromerida sp. RM11", "Colpodea",
                                                                                                   "Colponemida", "Breviatea", "Retortamonadidae", "Corallochytrium", "Leukarachnion sp. ATCC PRA-24", "Echinamoebida"), "Protozoa",
                                                                      ifelse(susceptible$Tax_group %in% c("Chlorarachniophyceae", "Florideophyceae", "Trebouxiophyceae", "Bangiophyceae", "Compsopogonophyceae", "Stylonematophyceae", "Zygnemophyceae",
                                                                                                          "Pyramimonadophyceae", "Chlorophyceae", "Ulvophyceae", "Coleochaetophyceae", "Klebsormidiophyceae", "Chlorokybophyceae", "Charophyceae", "Palmophyllophyceae", 
                                                                                                          "Pedinophyceae", "Chlorodendrophyceae", "Pycnococcaceae", "Mesostigmatophyceae", "Mamiellophyceae", "Chloropicophyceae", "Nephroselmidophyceae", "Rhodellophyceae", "Phaeophyceae", "Pelagophyceae",
                                                                                                          "Dinophyceae", "Coscinodiscophyceae", "Glaucocystophyceae", "Vitrellaceae", "Chrysomerophyceae", "Phaeothamniophyceae", "Chrysomerophyceae",
                                                                                                          "Aurearenophyceae", "Pinguiophyceae", "Aurearenophyceae", "Stylonematophyceae", "Glaucocystophyceae", "Cryptophyceae", "Compsopogonophyceae", "Eustigmatophyceae", "Bacillariophyceae", "Bolidophyceae", "Xanthophyceae",
                                                                                                          "Chrysoparadoxa", "Mediophyceae", "Fragilariophyceae", "Chrysophyceae", "Raphidophyceae", "Synurophyceae", "Dictyochophyceae", "Conoidasida", "Chromeraceae", "Schizocladia"), "Algae/Diatoms",
                                                                             ifelse(susceptible$Tax_group %in% c("Oligohymenophorea", "Litostomatea", "Spirotrichea", "Perkinsida", "Oomycota","Bigyra", "Haptophyta", "Plasmodiophorida", "Reticulomyxidae", "Palpitomonas", "Polycystinea", "Armophorea", "Leannia",
                                                                                                                 "Miliolida", "Crithionina",  "Phyllopharyngea", "Bacillariophyta sp. SL64/78c",  "Phyllopharyngea","Phyllopharyngea", "Nassophorea", "Plagiopylea", "Karyorelictea", "Prostomatea", "Sticholonchida"), "other_Chromista",
                                                                                    "Prokaryotes/Archaea"))))))))

list(unique(susceptible$general_taxa))
options(max.print=1000000)
list(unique((susceptible %>% filter(general_taxa == "Prokaryotes/Archaea"))$`Common Name`))

#fish
susceptible_fish <- susceptible %>% filter(general_taxa == "Fish") %>% group_by(gene, general_taxa, Level) %>% summarize(sus_level = n_distinct(`Susceptibility Prediction`))
susceptible_fish_gene <- susceptible_fish %>% select(gene) %>% distinct()
susceptible_fish_L1 <- susceptible_fish %>% filter(Level == "L1") %>% rename("Level_1" = "Level")%>% select(-"sus_level")
susceptible_fish_L1$Level_1 <- gsub("L1", "Susceptible", susceptible_fish_L1$Level_1) 
susceptible_fish_L2 <- susceptible_fish %>% filter(Level == "L2")%>% rename("Level_2" = "Level")%>% select(-"sus_level")
susceptible_fish_L2$Level_2 <- gsub("L2", "Susceptible", susceptible_fish_L2$Level_2)

susceptible_fish_final <- left_join(susceptible_fish_gene, susceptible_fish_L1)
susceptible_fish_final <- left_join(susceptible_fish_final, susceptible_fish_L2)

View(susceptible_fish_final)

write_xlsx(susceptible_fish_final, "susceptible_fish_final.xlsx")

#bind in accession IDs for DAVID
target_IDS <- read_excel("SEQAPASS_Input.xlsx") %>% select(Uniprot_ID, Protein_Target) %>% rename("gene" = "Protein_Target")

david_list <- left_join(susceptible_fish_final, target_IDS) %>% ungroup() %>% select(gene, Uniprot_ID)

write_xlsx(david_list, "david_list.xlsx")
