#Structure Wrangle#

library(tidyverse)
library(readxl)
library(writexl)
library(ggpubr)
library(pvclust)
library(vegan)
library(rdist)
library(beepr)
library(boot)
library(parallel)
library(cluster)
library(clValid)
library(ade4)

length(cl <- makeCluster(detectCores()))
binary.sqrt <- function(x){
  jac.dist <- dist(t(x), method="binary")
  sqrt.jac.dist <- sqrt(jac.dist)
  return(sqrt.jac.dist)
}

binary6 <- function(x){
  dist <- dist.binary(t(x), method=6) # change method here - see ?dist.binary
  return(dist)
}
binary5 <- function(x){
  dist <- dist.binary(t(x), method=5) # change method here - see ?dist.binary
  return(dist)
}
binary4 <- function(x){
  dist <- dist.binary(t(x), method=4) # change method here - see ?dist.binary
  return(dist)
}
binary3 <- function(x){
  dist <- dist.binary(t(x), method=3) # change method here - see ?dist.binary
  return(dist)
}

binary2 <- function(x){
  dist <- dist.binary(t(x), method=2) # change method here - see ?dist.binary
  return(dist)
}

binary7 <- function(x){
  dist <- dist.binary(t(x), method=7) # change method here - see ?dist.binary
  return(dist)
}
binary8 <- function(x){
  dist <- dist.binary(t(x), method=8) # change method here - see ?dist.binary
  return(dist)
}

binary1 <- function(x){
  dist <- dist.binary(t(x), method=1) # change method here - see ?dist.binary
  return(dist)
}
binary9 <- function(x){
  dist <- dist.binary(t(x), method=9) # change method here - see ?dist.binary
  return(dist)
}

binary10 <- function(x){
  dist <- dist.binary(t(x), method=10) # change method here - see ?dist.binary
  return(dist)
}

#read in consensus MOA & ECOTOX output documents & bind
consensus_moa <- read_excel("Final_Consensus_MOA_EM_annotated.xlsx")
ecotox_ec <- read_excel("ECOTOX_EC.xlsx")

all_data_so_far <- left_join(ecotox_ec, consensus_moa, by = "CAS")

#add in missing data
missing_data <- anti_join(consensus_moa, ecotox_ec, by = "CAS")
names(missing_data)
all_data_so_far1 <- bind_rows(all_data_so_far, missing_data)
View(all_data_so_far1)

#write out CAS for toxprint data pull 
write_xlsx((all_data_so_far1 %>% select(CAS) %>% distinct()), "CAS_for_toxprint.xlsx")
write_xlsx(all_data_so_far1, "data_for_cluster_annotation.xlsx")

#read in toxprint data
tox_print <- read_csv("tox_print.csv") %>% select(-c("FOUND_BY", DTXSID, PREFERRED_NAME)) %>% rename("CAS" = "INPUT")
names(tox_print)

#bind in tox_print & group chemicals
data_plus_toxprint <- left_join(all_data_so_far1, tox_print, by = "CAS")

#milwaukee data####
milwaukee_CAS <- read_excel("Milwaukee_chem_QSAR_Consensus.xlsx") %>% select(CAS)

milwaukee_data_so_far <- left_join(milwaukee_CAS, all_data_so_far1)
milwaukee_data_toxprint <- left_join((milwaukee_data_so_far %>% select(-Chemical)), tox_print, by = "CAS")
View(milwaukee_data_toxprint)

#split into QSAR groups
narcosis_mil <- milwaukee_data_toxprint %>% filter(QSAR_class_detail %in% c("Narcosis", "Unknown"))
g_reactive_mil <- milwaukee_data_toxprint %>% filter(QSAR_class_detail %in% c("Reactive", "Unknown"))

#narcosis_unknown chemicals
narcosis_mil_cluster <- narcosis_mil %>% select("CAS","atom:element_main_group":"ring:polycycle_tricyclo_benzvalene") %>% drop_na() %>% 
  filter(!CAS %in% c("29385-43-1", "34911-55-2")) %>%
  mutate_at(vars("atom:element_main_group":"ring:polycycle_tricyclo_benzvalene"), as.numeric)%>% 
  column_to_rownames("CAS")

narcosis_mil_cluster <- as.matrix(narcosis_mil_cluster)
t_narcosis_mil_cluster<- t(narcosis_mil_cluster)
t_narcosis_mil_cluster <- as.matrix(t_narcosis_mil_cluster)
View(t_narcosis_mil_cluster)
#pvclust
set.seed(3)
pvclust_narcosis_unknown_mil <- pvclust(t_narcosis_mil_cluster, method.hclust = "average", method.dist = binary6,
                            nboot = 10000, parallel = FALSE, iseed = NULL)
beep(5)
pvclust_narcosis_plot <- plot(pvclust_narcosis_unknown_mil, hang = -1, cex = 0.5)

pvrect(pvclust_narcosis_unknown_mil, alpha = 0.80, pv = "au", max.only = FALSE) 
print(pvclust_narcosis_unknown_mil)

pvclust_narcosis_clusters <- pvpick(pvclust_narcosis_unknown_mil, alpha = 0.80, pv = "au", max.only = FALSE)

View(pvclust_narcosis_clusters$clusters)
pvclust_narcosis_clusters$clusters %>% select(Value)

list <- do.call(rbind, pvclust_narcosis_clusters$clusters)
t_list <- t(list)
t_list <- as.data.frame(t_list)
names(t_list)

t_list_1 <- t_list %>% gather(c("V1":"V35"), key = "Narcosis_Cluster", value = "CAS") %>% distinct()
t_list_1$Narcosis_Cluster <- gsub("V", "", t_list_1$Narcosis_Cluster)

narcosis_unknown_clusters <- t_list_1
narcosis_unknown_clusters_for_elim <- narcosis_unknown_clusters %>% group_by(Narcosis_Cluster) %>%
  summarize(n_distinct_CAS = n_distinct(CAS))
narcosis_unknown_clusters_2 <- left_join(narcosis_unknown_clusters, narcosis_unknown_clusters_for_elim) %>%
  filter(n_distinct_CAS < 10) 

narcosis_unknown_clusters_3 <- narcosis_unknown_clusters_2 %>% group_by(CAS) %>% summarize(max_CAS = max(n_distinct_CAS)) %>%
  rename(n_distinct_CAS = max_CAS)
View(narcosis_unknown_clusters_3)

narcosis_unknown_clusters_4 <- left_join(narcosis_unknown_clusters_3, narcosis_unknown_clusters_2) %>%
  rename(n_Chem = n_distinct_CAS)
View(narcosis_unknown_clusters_4)

#bind to data so far
milwaukee_data_so_far1 <- left_join(milwaukee_data_so_far, narcosis_unknown_clusters_4, by = "CAS")
View(milwaukee_data_so_far1)

#generallyreactive_unknown chemicals
names(g_reactive_mil)
g_reactive_unknown_cluster <- g_reactive_mil %>% select("CAS","atom:element_main_group":"ring:polycycle_tricyclo_benzvalene") %>% drop_na() %>%
  filter(CAS != "29385-43-1") %>%
  mutate_at(vars("atom:element_main_group":"ring:polycycle_tricyclo_benzvalene"), as.numeric)%>% 
  column_to_rownames("CAS")
View(g_reactive_unknown_cluster)
g_reactive_unknown_cluster <- as.matrix(g_reactive_unknown_cluster)
t_g_reactive_unknown_cluster<- t(g_reactive_unknown_cluster)
t_g_reactive_unknown_cluster <- as.matrix(t_g_reactive_unknown_cluster)
View(t_g_reactive_unknown_cluster)

#pvclust
set.seed(3)
pvclust_greactive_unknown <- pvclust(t_g_reactive_unknown_cluster, method.hclust = "average", method.dist = binary6,
                                    nboot = 10000, parallel = FALSE, iseed = NULL)

beep()
greactive_unknown_plot <- plot(pvclust_greactive_unknown, hang = -1, cex = 0.5)

pvrect(pvclust_greactive_unknown, alpha = 0.80, pv = "au", max.only = FALSE) 
print(pvclust_greactive_unknown)
ggsave("greactive_unknown_plot.jpeg", pvclust_greactive_unknown)

pvclust_greactive_clusters <- pvpick(pvclust_greactive_unknown, alpha = 0.80, pv = "au", max.only = FALSE)

list <- do.call(rbind, pvclust_greactive_clusters$clusters)
t_list <- t(list)
t_list <- as.data.frame(t_list)
names(t_list)

t_list_1 <- t_list %>% gather(c("V1":"V23"), key = "Reactive_Cluster", value = "CAS") %>% distinct()
t_list_1$Reactive_Cluster <- gsub("V", "", t_list_1$Reactive_Cluster)

generallyreactive_unknown_clusters <- t_list_1
generallyreactive_unknown_clusters_for_elim <- generallyreactive_unknown_clusters %>% group_by(Reactive_Cluster) %>%
  summarize(n_distinct_CAS = n_distinct(CAS))
generallyreactive_unknown_clusters_2 <- left_join(generallyreactive_unknown_clusters, generallyreactive_unknown_clusters_for_elim) %>%
  filter(n_distinct_CAS < 10) 

generallyreactive_unknown_clusters_3 <- generallyreactive_unknown_clusters_2 %>% group_by(CAS) %>% summarize(max_CAS = max(n_distinct_CAS)) %>%
  rename(n_distinct_CAS = max_CAS)

generallyreactive_unknown_clusters_4 <- left_join(generallyreactive_unknown_clusters_3, generallyreactive_unknown_clusters_2) %>%
  rename(n_Chem_reactive = n_distinct_CAS)

View(generallyreactive_unknown_clusters_4)

#bind to data so far
milwaukee_data_so_far3 <- left_join(milwaukee_data_so_far1, generallyreactive_unknown_clusters_4, by = "CAS")
View(milwaukee_data_so_far3)

#bind in chnm
milwaukee_chnm <- read_excel("Milwaukee_chem_QSAR_Consensus.xlsx") %>% select(Chemical_Name, CAS)

milwaukee_data_so_far4 <- left_join(milwaukee_data_so_far3, milwaukee_chnm)
milwaukee_data_so_far4$Narcosis_Cluster <- paste("N", milwaukee_data_so_far4$Narcosis_Cluster, sep = "")
milwaukee_data_so_far4$Reactive_Cluster <- paste("R", milwaukee_data_so_far4$Reactive_Cluster, sep = "")

#bind in class from GLRI to name groups - if possible 
chemical_class <- read_excel("CAS_Class_GLRI.xlsx") %>% select(CAS, Class)

milwaukee_data_so_far5 <- left_join(milwaukee_data_so_far4, chemical_class)
write_xlsx(milwaukee_data_so_far5, "milwaukee_data_so_far6.xlsx")

#all data fo read-across####
#split into QSAR groups
list(unique(data_plus_toxprint$QSAR_class_detail))

narcosis_unknown <- data_plus_toxprint %>% filter(QSAR_class_detail %in% c("Narcosis", "Unknown")) %>%
  filter(`atom:element_main_group` != "N/A") %>% distinct()%>% filter(!CAS %in% c("7446-18-6", "95-76-1"))

reactive_unknown <- data_plus_toxprint %>% filter(QSAR_class_detail%in% c("Reactive","Unknown")) %>%
  filter(`atom:element_main_group` != "N/A") %>% distinct()%>% filter(CAS != "7446-18-6")

#reactive cluster
reactive_cluster <- reactive_unknown %>% select("CAS","atom:element_main_group":"ring:polycycle_tricyclo_benzvalene") %>% drop_na() %>% 
  mutate_at(vars("atom:element_main_group":"ring:polycycle_tricyclo_benzvalene"), as.numeric)%>% 
  column_to_rownames("CAS") %>% drop_na()
View(reactive_cluster)

reactive_cluster <- as.matrix(reactive_cluster)
t_reactive_cluster <- t(reactive_cluster)
t_reactive_cluster <- as.matrix(t_reactive_cluster)

#pvclust
set.seed(3)
pvclust_reactive <- pvclust(t_reactive_cluster, method.hclust = "complete", method.dist = binary6,
                            nboot = 10000, parallel = FALSE, iseed = NULL)

pvclust_reactive_plot <- plot(pvclust_reactive, hang = -1, cex = 0.5)
beep()
pvrect(pvclust_reactive, alpha = 0.90, pv = "au", max.only = TRUE) 
print(pvclust_reactive)

pvclust_reactive_clusters <- pvpick(pvclust_reactive, alpha = 0.90, pv = "au", max.only = TRUE)

list <- do.call(rbind, pvclust_reactive_clusters$clusters)
t_list <- t(list)
t_list <- as.data.frame(t_list)
names(t_list)

t_list_1 <- t_list %>% gather(c("V1":"V81"), key = "Cluster", value = "CAS")
t_list_1$Cluster <- gsub("V", "", t_list_1$Cluster)
View(t_list_1)

reactive_clusters <- t_list_1
reactive_clusters$Cluster <- paste("R", reactive_clusters$Cluster, sep = "")
write_xlsx(reactive_clusters, "reactive_clusters_RA.xlsx")

###continue code####
#narcotic_unknown chemicals
names(narcosis_unknown)
narcosis_unknown_cluster <- narcosis_unknown %>% select("CAS","atom:element_main_group":"ring:polycycle_tricyclo_benzvalene") %>% 
  mutate_at(vars("atom:element_main_group":"ring:polycycle_tricyclo_benzvalene"), as.numeric)%>% 
  column_to_rownames("CAS") %>% drop_na() 

narcosis_unknown_cluster <- as.matrix(narcosis_unknown_cluster)
t_narcosis_unknown_cluster <- t(narcosis_unknown_cluster)
t_narcosis_unknown_cluster <- as.matrix(t_narcosis_unknown_cluster)

#pvclust
set.seed(3)
CLUST.narcosis <- pvclust(t_narcosis_unknown_cluster, method.hclust = "complete", method.dist = binary6,
                            nboot = 10000, parallel = FALSE, iseed = NULL)
beep()
pvclust_narcosis_plot <- plot(CLUST.narcosis, hang = -1, cex = 0.5)

pvrect(CLUST.narcosis, alpha = 0.90, pv = "au", max.only = TRUE) 
print(pvclust_narcosis)

pvclust_narcosis_clusters <- pvpick(CLUST.narcosis, alpha = 0.90, pv = "au", max.only = TRUE)

View(pvclust_narcosis_clusters$clusters)
pvclust_narcosis_clusters$clusters %>% select(Value)

list <- do.call(rbind, pvclust_narcosis_clusters$clusters)
t_list <- t(list)
t_list <- as.data.frame(t_list)
names(t_list)
t_list_1 <- t_list %>% gather(c("V1":"V115"), key = "Cluster", value = "CAS") %>% distinct()
t_list_1$Cluster <- gsub("V", "", t_list_1$Cluster)
View(t_list_1)

narcosis_clusters <- t_list_1

write_xlsx(narcosis_clusters, "narcosis_clusters.xlsx")

#pull in cluster data
reactive_clusters <- read_excel("reactive_clusters_RA.xlsx")
reactive_clusters$Cluster <- paste("R", reactive_clusters$Cluster, sep = "")
narcosis_clusters <- read_excel("narcosis_clusters.xlsx")
narcosis_clusters$Cluster <- paste("N", narcosis_clusters$Cluster, sep = "")

RA_clusters <- bind_rows(reactive_clusters, narcosis_clusters)
names(RA_clusters)

all_RA_data <- left_join((data_plus_toxprint %>% select(CAS:QSAR_class_detail)),RA_clusters) 
all_RA_data$Cluster <- ifelse(is.na(all_RA_data$Cluster), "unclustered", all_RA_data$Cluster)

#write CAS for logP derivation
write_xlsx((all_RA_data %>% select(CAS)), "CAS_for_logP_pull.xlsx")

#just use log P for discrimination####
phys_chem <- read_excel("phys_chem_output.xlsx", 2) %>% select("CAS", "LogP") %>%
  rename("Log_P" = "LogP")
names(phys_chem)

#join to dataset
all_RA_data1 <- left_join(all_RA_data, phys_chem) %>% distinct()

#write dataset & check it out
write_xlsx(all_RA_data1, "data_all_with_physchem.xlsx")

#annotate dataset & read back in
annotated_data_w_physchem <- read_excel("data_all_with_physchem_EM_annotated.xlsx")
View(annotated_data_w_physchem)

#milwaukee data####
milwaukee_chem_list <- read_excel("milwaukee_data_so_far6_EM_annotated_for_analysis.xlsx") %>% select(Chemical, CAS)
View(milwaukee_chem_list)
ECOTOX_EC <- read_excel("ECOTOX_EC.xlsx")

milwaukee_with_EC <- left_join(milwaukee_chem_list, (ECOTOX_EC %>% select(-Chemical))) %>% rename("EC_estimate" = "Concentration")
milwaukee_with_EC$EC_Source <- ifelse(is.na(milwaukee_with_EC$EC_estimate), "", "ECOTOX Knowledgebase")
milwaukee_with_EC1 <- milwaukee_with_EC %>% filter(!is.na(EC_estimate))
View(milwaukee_with_EC1)

#direct RA
list(unique(annotated_data_w_physchem$Notes))
direct_RA <- annotated_data_w_physchem %>% filter(Notes == "Direct_RA")
View(direct_RA)
noEC <- direct_RA %>% filter(is.na(Concentration))
EC <- direct_RA %>% filter(!is.na(Concentration)) %>% select(CAS, Cluster, Concentration) %>% rename("CAS_2" = "CAS", "EC_estimate" = "Concentration")
direct_RA <- left_join(noEC, EC) %>% select(-Concentration, -c(ASTER_MOA:Confidence_Score), -Notes)
direct_RA$EC_Source <- "Direct Read-across "
View(direct_RA)

milwaukee_with_EC2 <- left_join(milwaukee_chem_list, direct_RA) %>% filter(!is.na(EC_estimate))
View(milwaukee_with_EC2)

milwaukee_with_EC3 <- bind_rows(milwaukee_with_EC1, milwaukee_with_EC2)
View(milwaukee_with_EC3)

#readacross based on logP - 10 nearest, then calculate mean weighted
names(annotated_data_w_physchem)
list(unique(annotated_data_w_physchem$Notes))
data_for_eval <- annotated_data_w_physchem %>% filter(Notes == "logP_RA")
data_for_eval$Log_P <- as.numeric(data_for_eval$Log_P)
mean_logP <- data_for_eval %>% group_by(Cluster) %>% summarize(mean_logP = mean(Log_P))
View(mean_logP)

data_for_eval1 <- left_join(data_for_eval, mean_logP)

#split apart cluster groups into matrices to evaluate distance
list(unique(data_for_eval1$Cluster))


#N113####
N113<- data_for_eval1 %>% filter(Cluster == "N113")
N113$CAS[which(is.na(N113$Concentration))] #"76-22-2" 

N113 <- N113 %>% select(CAS, Log_P) %>% column_to_rownames("CAS")
(N113)
N113 <- as.matrix(N113)
N113_dist <- pdist(N113, metric = "euclidean", p = 2)

t_N113_dist <- t(N113_dist)
(t_N113_dist)

N113_1 <- as.data.frame(cbind(N113, t_N113_dist)) %>% rownames_to_column() %>% rename("CAS" = "rowname") %>%
  gather(c("V2":"V15"), key = "CAS_2", value = "distance")
View(N113_1)

N113_1$CAS_2 <- gsub("V3", "464-48-2", N113_1$CAS_2)
N113_1$CAS_2 <- gsub("V4", "700-58-3", N113_1$CAS_2)
N113_1$CAS_2 <- gsub("V5", "126-81-8", N113_1$CAS_2)
N113_1$CAS_2 <- gsub("V6", "514-10-3", N113_1$CAS_2)
N113_1$CAS_2 <- gsub("V7", "79-77-6", N113_1$CAS_2)
N113_1$CAS_2 <- gsub("V8", "108-94-1", N113_1$CAS_2)
N113_1$CAS_2 <- gsub("V9", "5989-27-5", N113_1$CAS_2)
N113_1$CAS_2 <- gsub("V10", "77-73-6", N113_1$CAS_2)
N113_1$CAS_2 <- gsub("V11", "78-59-1", N113_1$CAS_2)
N113_1$CAS_2 <- gsub("V12", "5835-26-7", N113_1$CAS_2)
N113_1$CAS_2 <- gsub("V13", "471-77-2", N113_1$CAS_2)
N113_1$CAS_2 <- gsub("V14", "498-66-8", N113_1$CAS_2)
N113_1$CAS_2 <- gsub("V15", "76-22-2", N113_1$CAS_2)
N113_1$CAS_2 <- gsub("V2", "10293-06-8", N113_1$CAS_2)

N113_1_estimates <- N113_1 %>% filter(CAS %in% c("76-22-2"))
N113_1_estimates$distance <- N113_1_estimates$distance + 1
N113_1_estimates$weight <- 1/N113_1_estimates$distance
sum(N113_1_estimates1$weight)
View(N113_1_estimates)

#bind in EC estimates
N113_1_ECOTOX <- data_for_eval1 %>% filter(Cluster == "N113") %>% select(CAS, Concentration) %>% rename("CAS_2" = "CAS")
View(N113_1_ECOTOX)
N113_2 <- left_join(N113_1_estimates, N113_1_ECOTOX) %>% filter(CAS_2 != "76-22-2") %>% slice_max(weight, n = 10)
N113_2$log_Con <- log10(N113_2$Concentration)

N113_2$Con_dist_wt <- N113_2$log_Con * N113_2$weight
N113_3 <- N113_2 %>% group_by(CAS) %>% summarize(sum_con_dist = sum(Con_dist_wt), sum_wt = sum(weight))
N113_3$log_estimate <- N113_3$sum_con_dist / N113_3$sum_wt
N113_3$EC_estimate <- 10^(N113_3$log_estimate)
View(N113_3)

#bind together all
N113_final <- left_join(N113_2, N113_3)
View(N113_final)

shortform_EC <- N113_final %>% select(c(CAS, EC_estimate))
longform_EC <- N113_final

#N53####
N53<- data_for_eval1 %>% filter(Cluster == "N53")
N53$CAS[which(is.na(N53$Concentration))] #"89-78-1"

N53 <- N53 %>% select(CAS, Log_P) %>% column_to_rownames("CAS")
N53 <- as.matrix(N53)
N53_dist <- pdist(N53, metric = "euclidean", p = 2)

t_N53_dist<- t(N53_dist)

N53_1 <- as.data.frame(cbind(N53, N53_dist)) %>% rownames_to_column() %>% rename("CAS" = "rowname") %>%
  gather(c("V2":"V4"), key = "CAS_2", value = "distance")

N53_1$CAS_2 <- gsub("V2", "2216-51-5", N53_1$CAS_2)
N53_1$CAS_2 <- gsub("V3", "108-93-0", N53_1$CAS_2)
N53_1$CAS_2 <- gsub("V4", "89-78-1", N53_1$CAS_2)

N53_1_estimates <- N53_1 %>% filter(CAS %in% c("89-78-1"))
N53_1_estimates$distance <- N53_1_estimates$distance + 1
N53_1_estimates$weight <- 1/N53_1_estimates$distance

#bind in EC estimates
N53 <- data_for_eval1 %>% filter(Cluster == "N53") %>% select(CAS, Concentration) %>% rename("CAS_2" = "CAS")
N53_estimates <- left_join(N53_1_estimates, N53) %>% filter(CAS_2 !="89-78-1") 
N53_estimates$logCon <- log10(N53_estimates$Concentration)
N53_estimates$Con_dist_wt <- N53_estimates$logCon * N53_estimates$weight
N53_estimates_1 <- N53_estimates %>% group_by(CAS) %>% slice_max(weight, n = 10) %>% summarize(sum_con_dist = sum(Con_dist_wt), sum_wt = sum(weight))
N53_estimates_1$log_estimate <- N53_estimates_1$sum_con_dist/ N53_estimates_1$sum_wt
N53_estimates_1$EC_estimate <- 10^(N53_estimates_1$log_estimate)
View(N53_estimates_1)

#bind together all
N53_final <- left_join((N53_estimates %>% group_by(CAS) %>% slice_max(Con_dist_wt, n = 10)), N53_estimates_1)
View(N53_final)

shortform_EC <- bind_rows(shortform_EC, (N53_final %>% select(c(CAS, EC_estimate))))
longform_EC <- bind_rows(longform_EC, N53_final)


#N110####
N110<- data_for_eval1 %>% filter(Cluster == "N110")
N110$CAS[which(is.na(N110$Concentration))] #  "93413-62-8" "27203-92-5" "93413-69-5"

N110 <- N110 %>% select(CAS, Log_P) %>% column_to_rownames("CAS")
N110 <- as.matrix(N110)
N110_dist <- pdist(N110, metric = "euclidean", p = 2)

t_N110_dist <- t(N110_dist)

N110
t_N110_dist

N110_1 <- as.data.frame(cbind(N110, t_N110_dist)) %>% rownames_to_column() %>% rename("CAS" = "rowname") %>%
  gather(c("V2":"V7"), key = "CAS_2", value = "distance")
N110_1$CAS_2 <- gsub("V2", "3923-52-2", N110_1$CAS_2)
N110_1$CAS_2 <- gsub("V3", "127-66-2", N110_1$CAS_2)
N110_1$CAS_2 <- gsub("V4", "80-06-8", N110_1$CAS_2)
N110_1$CAS_2 <- gsub("V5", "93413-62-8", N110_1$CAS_2)
N110_1$CAS_2 <- gsub("V6", "27203-92-5", N110_1$CAS_2)
N110_1$CAS_2 <- gsub("V7", "93413-69-5", N110_1$CAS_2)


N110_1_estimates <- N110_1 %>% filter(CAS %in% c("93413-62-8", "27203-92-5", "93413-69-5"))
N110_1_estimates$distance <- N110_1_estimates$distance + 1
N110_1_estimates$weight <- 1/N110_1_estimates$distance

#bind in EC estimates
N110 <- data_for_eval1 %>% filter(Cluster == "N110") %>% select(CAS, Concentration) %>% rename("CAS_2" = "CAS")
N110_estimates <- left_join(N110_1_estimates, N110) %>% filter(!CAS_2 %in% c("93413-62-8", "27203-92-5", "93413-69-5"))
N110_estimates$log_Con <- log10(N110_estimates$Concentration)
N110_estimates$Con_dist_wt <- N110_estimates$log_Con * N110_estimates$weight

N110_estimates_1 <- N110_estimates %>% group_by(CAS) %>% slice_max(weight, n = 10) %>% summarize(sum_con_dist = sum(Con_dist_wt), sum_wt = sum(weight))
N110_estimates_1$log_estimate <- N110_estimates_1$sum_con_dist/ N110_estimates_1$sum_wt
N110_estimates_1$EC_estimate <- 10^(N110_estimates_1$log_estimate)

#bind together all
N110_final <- left_join((N110_estimates %>% group_by(CAS) %>% slice_max(weight, n = 10)), N110_estimates_1)
View(N110_final)
shortform_EC <- bind_rows(shortform_EC, (N110_final %>% select(c(CAS, EC_estimate))))
longform_EC <- bind_rows(longform_EC, N110_final)

#N56####
N56<- data_for_eval1 %>% filter(Cluster == "N56")
N56$CAS[which(is.na(N56$Concentration))] #"75-25-2"

N56 <- N56 %>% select(CAS, Log_P) %>% column_to_rownames("CAS")
N56 <- as.matrix(N56)
N56_dist <- pdist(N56, metric = "euclidean", p = 2)

t_N56_dist <- t(N56_dist)

N56
t_N56_dist

N56_1 <- as.data.frame(cbind(N56, t_N56_dist)) %>% rownames_to_column() %>% rename("CAS" = "rowname") %>%
  gather(c("V2":"V7"), key = "CAS_2", value = "distance")
N56_1$CAS_2 <- gsub("V2", "71-55-6", N56_1$CAS_2)
N56_1$CAS_2 <- gsub("V3", "56-23-5", N56_1$CAS_2)
N56_1$CAS_2 <- gsub("V4", "67-66-3", N56_1$CAS_2)
N56_1$CAS_2 <- gsub("V5", "75-09-2", N56_1$CAS_2)
N56_1$CAS_2 <- gsub("V6", "75-47-8", N56_1$CAS_2)
N56_1$CAS_2 <- gsub("V7", "75-25-2", N56_1$CAS_2)


N56_1_estimates <- N56_1 %>% filter(CAS %in% c("75-25-2"))
N56_1_estimates$distance <- 1 + N56_1_estimates$distance
N56_1_estimates$weight <- 1/N56_1_estimates$distance

#bind in EC estimates
N56 <- data_for_eval1 %>% filter(Cluster == "N56") %>% select(CAS, Concentration) %>% rename("CAS_2" = "CAS")
N56_estimates <- left_join(N56_1_estimates, N56) %>% filter(!CAS_2 %in% c("75-25-2"))
N56_estimates$logCon <- log10(N56_estimates$Concentration)
N56_estimates$Con_dist_wt <- N56_estimates$logCon * N56_estimates$weight
N56_estimates_1 <- N56_estimates %>% group_by(CAS) %>% slice_max(weight, n = 10) %>% summarize(sum_con_dist = sum(Con_dist_wt), sum_wt = sum(weight))
N56_estimates_1$log_estimate <- N56_estimates_1$sum_con_dist/ N56_estimates_1$sum_wt
N56_estimates_1$EC_estimate <- 10^(N56_estimates_1$log_estimate)

View(N56_estimates_1)

#bind together all
N56_final <- left_join((N56_estimates %>% group_by(CAS) %>% slice_max(weight, n = 10)), N56_estimates_1)
View(N56_final)

shortform_EC <- bind_rows(shortform_EC, (N56_final %>% select(c(CAS, EC_estimate))))
longform_EC <- bind_rows(longform_EC, N56_final)


#N93####
N93<- data_for_eval1 %>% filter(Cluster == "N93")
N93$CAS[which(is.na(N93$Concentration))] #"90-12-0" "91-57-6"

N93 <- N93 %>% select(CAS, Log_P) %>% column_to_rownames("CAS")
N93 <- as.matrix(N93)
N93_dist <- pdist(N93, metric = "euclidean", p = 2)

t_N93_dist <- t(N93_dist)

N93
t_N93_dist

N93_1 <- as.data.frame(cbind(N93, t_N93_dist)) %>% rownames_to_column() %>% rename("CAS" = "rowname") %>%
  gather(c("V2":"V13"), key = "CAS_2", value = "distance")
N93_1$CAS_2 <- gsub("V2", "95-63-6", N93_1$CAS_2)
N93_1$CAS_2 <- gsub("V3", "2245-38-7", N93_1$CAS_2)
N93_1$CAS_2 <- gsub("V4", "121-14-2", N93_1$CAS_2)
N93_1$CAS_2 <- gsub("V5", "5673-07-4", N93_1$CAS_2)
N93_1$CAS_2 <- gsub("V6", "99-08-1", N93_1$CAS_2)
N93_1$CAS_2 <- gsub("V7", "621-08-9", N93_1$CAS_2)
N93_1$CAS_2 <- gsub("V8", "108-38-3", N93_1$CAS_2)
N93_1$CAS_2 <- gsub("V9", "95-47-6", N93_1$CAS_2)
N93_1$CAS_2 <- gsub("V10", "106-42-3", N93_1$CAS_2)
N93_1$CAS_2 <- gsub("V11", "108-88-3", N93_1$CAS_2)
N93_1$CAS_2 <- gsub("V12", "90-12-0", N93_1$CAS_2)
N93_1$CAS_2 <- gsub("V13", "91-57-6", N93_1$CAS_2)

N93_1_estimates <- N93_1 %>% filter(CAS %in% c("90-12-0", "91-57-6"))
N93_1_estimates$distance <- N93_1_estimates$distance + 1
N93_1_estimates$weight <- 1/N93_1_estimates$distance

#bind in EC estimates
N93 <- data_for_eval1 %>% filter(Cluster == "N93") %>% select(CAS, Concentration) %>% rename("CAS_2" = "CAS")
N93_estimates <- left_join(N93_1_estimates, N93) %>% filter(!CAS_2 %in% c("90-12-0", "91-57-6"))
N93_estimates$logCon <- log10(N93_estimates$Concentration)
N93_estimates$Con_dist_wt <- N93_estimates$logCon * N93_estimates$weight

N93_estimates_1 <- N93_estimates %>% group_by(CAS) %>% slice_max(weight, n = 10) %>% summarize(sum_con_dist = sum(Con_dist_wt), sum_wt = sum(weight))
N93_estimates_1$log_estimate <- N93_estimates_1$sum_con_dist/ N93_estimates_1$sum_wt
N93_estimates_1$EC_estimate <- 10^(N93_estimates_1$log_estimate)

View(N93_estimates_1)

#bind together all
N93_final <- left_join((N93_estimates %>% group_by(CAS) %>% slice_max(Con_dist_wt, n = 10)), N93_estimates_1)
View(N93_final)

shortform_EC <- bind_rows(shortform_EC, (N93_final %>% select(c(CAS, EC_estimate))))
longform_EC <- bind_rows(longform_EC, N93_final)

#RR71####
RR71<- data_for_eval1 %>% filter(Cluster == "RR71")
RR71$CAS[which(is.na(RR71$Concentration))] #"13674-87-8" "115-96-8"

RR71 <- RR71 %>% select(CAS, Log_P) %>% column_to_rownames("CAS")
RR71 <- as.matrix(RR71)
RR71_dist <- pdist(RR71, metric = "euclidean", p = 2)

t_RR71_dist <- t(RR71_dist)

RR71
t_RR71_dist

RR71_1 <- as.data.frame(cbind(RR71, t_RR71_dist)) %>% rownames_to_column() %>% rename("CAS" = "rowname") %>%
  gather(c("V2":"V8"), key = "CAS_2", value = "distance")
RR71_1$CAS_2 <- gsub("V2", "109-64-8", RR71_1$CAS_2)
RR71_1$CAS_2 <- gsub("V3", "31502-57-5", RR71_1$CAS_2)
RR71_1$CAS_2 <- gsub("V4", "111-91-1", RR71_1$CAS_2)
RR71_1$CAS_2 <- gsub("V5", "5407-04-5", RR71_1$CAS_2)
RR71_1$CAS_2 <- gsub("V6", "106-89-8", RR71_1$CAS_2)
RR71_1$CAS_2 <- gsub("V7", "13674-87-8", RR71_1$CAS_2)
RR71_1$CAS_2 <- gsub("V8", "115-96-8", RR71_1$CAS_2)


RR71_1_estimates <- RR71_1 %>% filter(CAS %in% c("13674-87-8", "115-96-8"))
RR71_1_estimates$distance <- RR71_1_estimates$distance + 1
RR71_1_estimates$weight <- 1/RR71_1_estimates$distance

#bind in EC estimates
RR71 <- data_for_eval1 %>% filter(Cluster == "RR71") %>% select(CAS, Concentration) %>% rename("CAS_2" = "CAS")
RR71_estimates <- left_join(RR71_1_estimates, RR71)%>%  filter(!CAS_2 %in% c("13674-87-8", "115-96-8"))
RR71_estimates$logCon <- log10(RR71_estimates$Concentration)
RR71_estimates$Con_dist_wt <- RR71_estimates$logCon * RR71_estimates$weight
RR71_estimates_1 <- RR71_estimates %>% group_by(CAS) %>% slice_max(weight, n = 10) %>% summarize(sum_con_dist = sum(Con_dist_wt), sum_wt = sum(weight))
RR71_estimates_1$log_estimate <- RR71_estimates_1$sum_con_dist/ RR71_estimates_1$sum_wt
RR71_estimates_1$EC_estimate <- 10^(RR71_estimates_1$log_estimate)
View(RR71_estimates_1)

#bind together all
RR71_final <- left_join((RR71_estimates %>% group_by(CAS) %>% slice_max(weight, n = 10)), RR71_estimates_1)
View(RR71_final)

shortform_EC <- bind_rows(shortform_EC, (RR71_final %>% select(c(CAS, EC_estimate))))
longform_EC <- bind_rows(longform_EC, RR71_final)


#N103####
N103<- data_for_eval1 %>% filter(Cluster == "N103")
N103$CAS[which(is.na(N103$Concentration))] # "86-74-8"  "120-72-9"

N103 <- N103 %>% select(CAS, Log_P) %>% column_to_rownames("CAS")
N103 <- as.matrix(N103)
N103_dist <- pdist(N103, metric = "euclidean", p = 2)

t_N103_dist <- t(N103_dist)

N103
t_N103_dist

N103_1 <- as.data.frame(cbind(N103, t_N103_dist)) %>% rownames_to_column() %>% rename("CAS" = "rowname") %>%
  gather(c("V2":"V8"), key = "CAS_2", value = "distance")
N103_1$CAS_2 <- gsub("V2", "271-89-6", N103_1$CAS_2)
N103_1$CAS_2 <- gsub("V3", "496-16-2", N103_1$CAS_2)
N103_1$CAS_2 <- gsub("V4", "83-34-1", N103_1$CAS_2)
N103_1$CAS_2 <- gsub("V5", "132-64-9", N103_1$CAS_2)
N103_1$CAS_2 <- gsub("V6", "109-97-7", N103_1$CAS_2)
N103_1$CAS_2 <- gsub("V7", "86-74-8", N103_1$CAS_2)
N103_1$CAS_2 <- gsub("V8", "120-72-9", N103_1$CAS_2)

N103_1_estimates <- N103_1 %>% filter(CAS %in% c("86-74-8", "120-72-9"))
N103_1_estimates$distance <- N103_1_estimates$distance + 1
N103_1_estimates$weight <- 1/N103_1_estimates$distance

#bind in EC estimates
N103 <- data_for_eval1 %>% filter(Cluster == "N103") %>% select(CAS, Concentration) %>% rename("CAS_2" = "CAS")
N103_estimates <- left_join(N103_1_estimates, N103) %>% filter(!CAS_2 %in% c("86-74-8", "120-72-9"))
N103_estimates$logCon <- log10(N103_estimates$Concentration)
N103_estimates$Con_dist_wt <- N103_estimates$logCon * N103_estimates$weight

N103_estimates_1 <- N103_estimates %>% group_by(CAS) %>% slice_max(weight, n = 10) %>% summarize(sum_con_dist = sum(Con_dist_wt), sum_wt = sum(weight))
N103_estimates_1$log_estimate <- N103_estimates_1$sum_con_dist/ N103_estimates_1$sum_wt
N103_estimates_1$EC_estimate <- 10^N103_estimates_1$log_estimate


#bind together all
N103_final <- left_join((N103_estimates %>% group_by(CAS) %>% slice_max(weight, n = 10)), N103_estimates_1)
View(N103_final)

shortform_EC <- bind_rows(shortform_EC, (N103_final %>% select(c(CAS, EC_estimate))))
longform_EC <- bind_rows(longform_EC, N103_final)

#N40####
N40<- data_for_eval1 %>% filter(Cluster == "N40")
N40$CAS[which(is.na(N40$Concentration))] # "140-66-9" 

N40 <- N40 %>% select(CAS, Log_P) %>% column_to_rownames("CAS")
N40 <- as.matrix(N40)
N40_dist <- pdist(N40, metric = "euclidean", p = 2)

t_N40_dist <- t(N40_dist)

N40
t_N40_dist

N40_1 <- as.data.frame(cbind(N40, t_N40_dist)) %>% rownames_to_column() %>% rename("CAS" = "rowname") %>%
  gather(c("V2":"V4"), key = "CAS_2", value = "distance")
N40_1$CAS_2 <- gsub("V2", "732-26-3", N40_1$CAS_2)
N40_1$CAS_2 <- gsub("V3", "98-54-4", N40_1$CAS_2)
N40_1$CAS_2 <- gsub("V4", "140-66-9", N40_1$CAS_2)

N40_1_estimates <- N40_1 %>% filter(CAS %in% c("140-66-9"))
N40_1_estimates$distance <- N40_1_estimates$distance + 1
N40_1_estimates$weight <- 1 /N40_1_estimates$distance

#bind in EC estimates
N40 <- data_for_eval1 %>% filter(Cluster == "N40") %>% select(CAS, Concentration) %>% rename("CAS_2" = "CAS")
N40_estimates <- left_join(N40_1_estimates, N40)%>% filter(!CAS_2 %in% c("140-66-9"))
N40_estimates$logCon <- log10(N40_estimates$Concentration)
N40_estimates$Con_dist_wt <- N40_estimates$logCon * N40_estimates$weight
N40_estimates_1 <- N40_estimates %>% group_by(CAS) %>% slice_max(weight, n = 10) %>% summarize(sum_con_dist = sum(Con_dist_wt), sum_wt = sum(weight))

N40_estimates_1$log_estimate <- N40_estimates_1$sum_con_dist/ N40_estimates_1$sum_wt
N40_estimates_1$EC_estimate <- 10^(N40_estimates_1$log_estimate)
View(N40_estimates)

#bind together all
N40_final <- left_join((N40_estimates %>% group_by(CAS) %>% slice_max(weight, n = 10)), N40_estimates_1)
View(N40_final)

shortform_EC <- bind_rows(shortform_EC, (N40_final %>% select(c(CAS, EC_estimate))))
longform_EC <- bind_rows(longform_EC, N40_final)

#RR60####
RR60<- data_for_eval1 %>% filter(Cluster == "RR60")
RR60$CAS[which(is.na(RR60$Concentration))] #"57-53-4"

RR60 <- RR60 %>% select(CAS, Log_P) %>% column_to_rownames("CAS")
RR60 <- as.matrix(RR60)
RR60_dist <- pdist(RR60, metric = "euclidean", p = 2)

t_RR60_dist <- t(RR60_dist)

RR60
t_RR60_dist

RR60_1 <- as.data.frame(cbind(RR60, t_RR60_dist)) %>% rownames_to_column() %>% rename("CAS" = "rowname") %>%
  gather(c("V2":"V5"), key = "CAS_2", value = "distance")
RR60_1$CAS_2 <- gsub("V2", "55406-53-6", RR60_1$CAS_2)
RR60_1$CAS_2 <- gsub("V3", "17804-35-2", RR60_1$CAS_2)
RR60_1$CAS_2 <- gsub("V4", "761-65-9", RR60_1$CAS_2)
RR60_1$CAS_2 <- gsub("V5", "57-53-4", RR60_1$CAS_2)

RR60_1_estimates <- RR60_1 %>% filter(CAS %in% c("57-53-4"))
RR60_1_estimates$distance <- RR60_1_estimates$distance + 1
RR60_1_estimates$weight <- 1/RR60_1_estimates$distance

#bind in EC estimates
RR60 <- data_for_eval1 %>% filter(Cluster == "RR60") %>% select(CAS, Concentration) %>% rename("CAS_2" = "CAS")
RR60_estimates <- left_join(RR60_1_estimates, RR60) %>% filter(!CAS_2 %in% c("57-53-4"))
RR60_estimates$logCon <- log10(RR60_estimates$Concentration)
RR60_estimates$Con_dist_wt <- RR60_estimates$logCon * RR60_estimates$weight
RR60_estimates_1 <- RR60_estimates %>% group_by(CAS) %>% slice_max(weight, n = 10) %>% summarize(sum_con_dist = sum(Con_dist_wt), sum_wt = sum(weight))

RR60_estimates_1$log_estimate <- RR60_estimates_1$sum_con_dist/ RR60_estimates_1$sum_wt
RR60_estimates_1$EC_estimate <- 10^(RR60_estimates_1$log_estimate)

#bind together all
RR60_final <- left_join((RR60_estimates %>% group_by(CAS) %>% slice_max(weight, n = 10)), RR60_estimates_1)
View(RR60_final)

shortform_EC <- bind_rows(shortform_EC, (RR60_final %>% select(c(CAS, EC_estimate))))
longform_EC <- bind_rows(longform_EC, RR60_final)


#RR81####
RR81<- data_for_eval1 %>% filter(Cluster == "RR81")
RR81$CAS[which(is.na(RR81$Concentration))] #"83-46-5" "57-88-5"

RR81 <- RR81 %>% select(CAS, Log_P) %>% column_to_rownames("CAS")
RR81 <- as.matrix(RR81)
RR81_dist <- pdist(RR81, metric = "euclidean", p = 2)

t_RR81_dist <- t(RR81_dist)

RR81
t_RR81_dist

RR81_1 <- as.data.frame(cbind(RR81, t_RR81_dist)) %>% rownames_to_column() %>% rename("CAS" = "rowname") %>%
  gather(c("V2":"V13"), key = "CAS_2", value = "distance")
RR81_1$CAS_2 <- gsub("V2", "514-10-3", RR81_1$CAS_2)
RR81_1$CAS_2 <- gsub("V3", "584-79-2", RR81_1$CAS_2)
RR81_1$CAS_2 <- gsub("V4", "79-77-6", RR81_1$CAS_2)
RR81_1$CAS_2 <- gsub("V5", "1740-19-8", RR81_1$CAS_2)
RR81_1$CAS_2 <- gsub("V6", "78-59-1", RR81_1$CAS_2)
RR81_1$CAS_2 <- gsub("V7", "5835-26-7", RR81_1$CAS_2)
RR81_1$CAS_2 <- gsub("V8", "87-72-9", RR81_1$CAS_2)
RR81_1$CAS_2 <- gsub("V9", "51596-11-3", RR81_1$CAS_2)
RR81_1$CAS_2 <- gsub("V10", "471-77-2", RR81_1$CAS_2)
RR81_1$CAS_2 <- gsub("V11", "28434-00-6", RR81_1$CAS_2)
RR81_1$CAS_2 <- gsub("V12", "83-46-5", RR81_1$CAS_2)
RR81_1$CAS_2 <- gsub("V13", "57-88-5", RR81_1$CAS_2)

RR81_1_estimates <- RR81_1 %>% filter(CAS %in% c("83-46-5","57-88-5"))
RR81_1_estimates$distance <- RR81_1_estimates$distance + 1
RR81_1_estimates$weight <- 1 / RR81_1_estimates$distance

#bind in EC estimates
RR81 <- data_for_eval1 %>% filter(Cluster == "RR81") %>% select(CAS, Concentration) %>% rename("CAS_2" = "CAS")
RR81_estimates <- left_join(RR81_1_estimates, RR81) %>% filter(!CAS_2 %in% c("83-46-5","57-88-5"))
RR81_estimates$logCon <- log10(RR81_estimates$Concentration)
RR81_estimates$Con_dist_wt <- RR81_estimates$logCon * RR81_estimates$weight

RR81_estimates_1 <- RR81_estimates %>% group_by(CAS) %>% slice_max(weight, n = 10) %>% summarize(sum_con_dist = sum(Con_dist_wt), sum_wt = sum(weight))
RR81_estimates_1$log_estimate <- RR81_estimates_1$sum_con_dist/ RR81_estimates_1$sum_wt
RR81_estimates_1$EC_estimate <- 10^(RR81_estimates_1$log_estimate)

#bind together all
RR81_final <- left_join((RR81_estimates %>% group_by(CAS) %>% slice_max(weight, n = 10)), RR81_estimates_1)
View(RR81_final)

shortform_EC <- bind_rows(shortform_EC, (RR81_final %>% select(c(CAS, EC_estimate))))
longform_EC <- bind_rows(longform_EC, RR81_final)


#RR63####
RR63<- data_for_eval1 %>% filter(Cluster == "RR63")
RR63$CAS[which(is.na(RR63$Concentration))] #"532-03-6"

RR63 <- RR63 %>% select(CAS, Log_P) %>% column_to_rownames("CAS")
RR63 <- as.matrix(RR63)
RR63_dist <- pdist(RR63, metric = "euclidean", p = 2)

t_RR63_dist <- t(RR63_dist)

RR63
t_RR63_dist

RR63_1 <- as.data.frame(cbind(RR63, t_RR63_dist)) %>% rownames_to_column() %>% rename("CAS" = "rowname") %>%
  gather(c("V2":"V13"), key = "CAS_2", value = "distance")
RR63_1$CAS_2 <- gsub("V2", "116-06-3", RR63_1$CAS_2)
RR63_1$CAS_2 <- gsub("V3", "1111-78-0", RR63_1$CAS_2)
RR63_1$CAS_2 <- gsub("V4", "63-25-2", RR63_1$CAS_2)
RR63_1$CAS_2 <- gsub("V5", "1563-66-2", RR63_1$CAS_2)
RR63_1$CAS_2 <- gsub("V6", "16752-77-5", RR63_1$CAS_2)
RR63_1$CAS_2 <- gsub("V7", "68-12-2", RR63_1$CAS_2)
RR63_1$CAS_2 <- gsub("V8", "607-00-1", RR63_1$CAS_2)
RR63_1$CAS_2 <- gsub("V9", "23135-22-0", RR63_1$CAS_2)
RR63_1$CAS_2 <- gsub("V10", "114-26-1", RR63_1$CAS_2)
RR63_1$CAS_2 <- gsub("V11", "3206-31-3", RR63_1$CAS_2)
RR63_1$CAS_2 <- gsub("V12", "51-79-6", RR63_1$CAS_2)
RR63_1$CAS_2 <- gsub("V13", "532-03-6", RR63_1$CAS_2)

RR63_1_estimates <- RR63_1 %>% filter(CAS %in% c("532-03-6"))
RR63_1_estimates$distance <- RR63_1_estimates$distance + 1
RR63_1_estimates$weight <- 1/RR63_1_estimates$distance

#bind in EC estimates
RR63 <- data_for_eval1 %>% filter(Cluster == "RR63") %>% select(CAS, Concentration) %>% rename("CAS_2" = "CAS")
RR63_estimates <- left_join(RR63_1_estimates, RR63)%>% filter(!CAS_2 %in% c("532-03-6")) 
RR63_estimates$logCon <- log10(RR63_estimates$Concentration)
RR63_estimates$Con_dist_wt <- RR63_estimates$logCon * RR63_estimates$weight

RR63_estimates_1 <- RR63_estimates %>% group_by(CAS) %>% slice_max(weight, n = 10) %>% summarize(sum_con_dist = sum(Con_dist_wt), sum_wt = sum(weight))
RR63_estimates_1$log_estimate <- RR63_estimates_1$sum_con_dist/ RR63_estimates_1$sum_wt
RR63_estimates_1$EC_estimate <- 10^(RR63_estimates_1$log_estimate)

#bind together all
RR63_final <- left_join((RR63_estimates %>% group_by(CAS) %>% slice_max(weight, n = 10)), RR63_estimates_1)
View(RR63_final)

shortform_EC <- bind_rows(shortform_EC, (RR63_final %>% select(c(CAS, EC_estimate))))
longform_EC <- bind_rows(longform_EC, RR63_final)


#N100####
N100<- data_for_eval1 %>% filter(Cluster == "N100")
N100$CAS[which(is.na(N100$Concentration))] #"1912-24-9"

N100 <- N100 %>% select(CAS, Log_P) %>% column_to_rownames("CAS")
N100 <- as.matrix(N100)
N100_dist <- pdist(N100, metric = "euclidean", p = 2)

t_N100_dist <- t(N100_dist)

N100
t_N100_dist

N100_1 <- as.data.frame(cbind(N100, t_N100_dist)) %>% rownames_to_column() %>% rename("CAS" = "rowname") %>%
  gather(c("V2":"V6"), key = "CAS_2", value = "distance")
N100_1$CAS_2 <- gsub("V2", "834-12-8", N100_1$CAS_2)
N100_1$CAS_2 <- gsub("V3", "21725-46-2", N100_1$CAS_2)
N100_1$CAS_2 <- gsub("V4", "22936-86-3", N100_1$CAS_2)
N100_1$CAS_2 <- gsub("V5", "122-34-9", N100_1$CAS_2)
N100_1$CAS_2 <- gsub("V6", "1912-24-9", N100_1$CAS_2)

N100_1_estimates <- N100_1 %>% filter(CAS %in% c("1912-24-9"))
N100_1_estimates$distance <- N100_1_estimates$distance + 1
N100_1_estimates$weight <- 1/N100_1_estimates$distance

#bind in EC estimates
N100 <- data_for_eval1 %>% filter(Cluster == "N100") %>% select(CAS, Concentration) %>% rename("CAS_2" = "CAS")
N100_estimates <- left_join(N100_1_estimates, N100) %>% filter(!CAS_2 %in% c("1912-24-9"))
N100_estimates$logCon <- log10(N100_estimates$Concentration)
N100_estimates$Con_dist_wt <- N100_estimates$logCon * N100_estimates$weight

N100_estimates_1 <- N100_estimates %>% group_by(CAS) %>% slice_max(weight, n = 10) %>% summarize(sum_con_dist = sum(Con_dist_wt), sum_wt = sum(weight))
N100_estimates_1$log_estimate <- N100_estimates_1$sum_con_dist/ N100_estimates_1$sum_wt
N100_estimates_1$EC_estimate <- 10^(N100_estimates_1$log_estimate)

#bind together all
N100_final <- left_join((N100_estimates %>% group_by(CAS) %>% slice_max(weight, n = 10)), N100_estimates_1)
View(N100_final)

shortform_EC <- bind_rows(shortform_EC, (N100_final %>% select(c(CAS, EC_estimate))))
longform_EC <- bind_rows(longform_EC, N100_final)


#N108####
N108<- data_for_eval1 %>% filter(Cluster == "N108")
N108$CAS[which(is.na(N108$Concentration))] # "119-36-8"

N108 <- N108 %>% select(CAS, Log_P) %>% column_to_rownames("CAS")
N108 <- as.matrix(N108)
N108_dist <- pdist(N108, metric = "euclidean", p = 2)

t_N108_dist <- t(N108_dist)

N108
t_N108_dist

N108_1 <- as.data.frame(cbind(N108, t_N108_dist)) %>% rownames_to_column() %>% rename("CAS" = "rowname") %>%
  gather(c("V2":"V19"), key = "CAS_2", value = "distance")
N108_1$CAS_2 <- gsub("V2", "94-09-7", N108_1$CAS_2)
N108_1$CAS_2 <- gsub("V3", "84-66-2", N108_1$CAS_2)
N108_1$CAS_2 <- gsub("V4", "5372-81-6", N108_1$CAS_2)
N108_1$CAS_2 <- gsub("V5", "5292-45-5", N108_1$CAS_2)
N108_1$CAS_2 <- gsub("V6", "131-11-3", N108_1$CAS_2)
N108_1$CAS_2 <- gsub("V7", "84-62-8", N108_1$CAS_2)
N108_1$CAS_2 <- gsub("V8", "93-89-0", N108_1$CAS_2)
N108_1$CAS_2 <- gsub("V9", "118-61-6", N108_1$CAS_2)
N108_1$CAS_2 <- gsub("V10", "2150-47-2", N108_1$CAS_2)
N108_1$CAS_2 <- gsub("V11", "2905-69-3", N108_1$CAS_2)
N108_1$CAS_2 <- gsub("V12", "42087-80-9", N108_1$CAS_2)
N108_1$CAS_2 <- gsub("V13", "1126-46-1", N108_1$CAS_2)
N108_1$CAS_2 <- gsub("V14", "619-50-1", N108_1$CAS_2)
N108_1$CAS_2 <- gsub("V15", "133-11-9", N108_1$CAS_2)
N108_1$CAS_2 <- gsub("V16", "118-55-8", N108_1$CAS_2)
N108_1$CAS_2 <- gsub("V17", "532-32-1", N108_1$CAS_2)
N108_1$CAS_2 <- gsub("V18", "54-21-7", N108_1$CAS_2)
N108_1$CAS_2 <- gsub("V19", "119-36-8", N108_1$CAS_2)

N108_1_estimates <- N108_1 %>% filter(CAS %in% c("119-36-8"))
N108_1_estimates$distance <- N108_1_estimates$distance + 1
N108_1_estimates$weight <- 1/N108_1_estimates$distance

#bind in EC estimates
N108 <- data_for_eval1 %>% filter(Cluster == "N108") %>% select(CAS, Concentration) %>% rename("CAS_2" = "CAS")
N108_estimates <- left_join(N108_1_estimates, N108) %>% filter(!CAS_2 %in% c("119-36-8"))
N108_estimates$logCon <- log10(N108_estimates$Concentration)
N108_estimates$Con_dist_wt <- N108_estimates$logCon * N108_estimates$weight

N108_estimates_1 <- N108_estimates %>% group_by(CAS) %>% slice_max(weight, n = 10) %>% summarize(sum_con_dist = sum(Con_dist_wt), sum_wt = sum(weight))
N108_estimates_1$log_estimate <- N108_estimates_1$sum_con_dist/ N108_estimates_1$sum_wt
N108_estimates_1$EC_estimate <- 10^N108_estimates_1$log_estimate

#bind together all
N108_final <- left_join((N108_estimates %>% group_by(CAS) %>% slice_max(weight, n = 10)), N108_estimates_1)
View(N108_final)

shortform_EC <- bind_rows(shortform_EC, (N108_final %>% select(c(CAS, EC_estimate))))
longform_EC <- bind_rows(longform_EC, N108_final)


#RR80####
RR80<- data_for_eval1 %>% filter(Cluster == "RR80")
RR80$CAS[which(is.na(RR80$Concentration))] # "84-65-1"

RR80 <- RR80 %>% select(CAS, Log_P) %>% column_to_rownames("CAS")
RR80 <- as.matrix(RR80)
RR80_dist <- pdist(RR80, metric = "euclidean", p = 2)

t_RR80_dist <- t(RR80_dist)

RR80
t_RR80_dist

RR80_1 <- as.data.frame(cbind(RR80, t_RR80_dist)) %>% rownames_to_column() %>% rename("CAS" = "rowname") %>%
  gather(c("V2":"V6"), key = "CAS_2", value = "distance")
RR80_1$CAS_2 <- gsub("V2", "66-76-2", RR80_1$CAS_2)
RR80_1$CAS_2 <- gsub("V3", "3347-22-6", RR80_1$CAS_2)
RR80_1$CAS_2 <- gsub("V4", "58-27-5", RR80_1$CAS_2)
RR80_1$CAS_2 <- gsub("V5", "83-79-4", RR80_1$CAS_2)
RR80_1$CAS_2 <- gsub("V6", "84-65-1", RR80_1$CAS_2)

RR80_1_estimates <- RR80_1 %>% filter(CAS %in% c("84-65-1"))
RR80_1_estimates$distance <- 1+RR80_1_estimates$distance
RR80_1_estimates$weight <- 1/RR80_1_estimates$distance

#bind in EC estimates
RR80 <- data_for_eval1 %>% filter(Cluster == "RR80") %>% select(CAS, Concentration) %>% rename("CAS_2" = "CAS")
RR80_estimates <- left_join(RR80_1_estimates, RR80) %>% filter(!CAS_2 %in% c("84-65-1"))
RR80_estimates$logCon <- log10(RR80_estimates$Concentration)
RR80_estimates$Con_dist_wt <- RR80_estimates$logCon * RR80_estimates$weight
RR80_estimates_1 <- RR80_estimates %>% group_by(CAS) %>% slice_max(weight, n = 10) %>% summarize(sum_con_dist = sum(Con_dist_wt), sum_wt = sum(weight))
RR80_estimates_1$log_estimate <- RR80_estimates_1$sum_con_dist/ RR80_estimates_1$sum_wt
RR80_estimates_1$EC_estimate <- 10^(RR80_estimates_1$log_estimate)

#bind together all
RR80_final <- left_join((RR80_estimates %>% group_by(CAS) %>% slice_max(weight, n = 10)), RR80_estimates_1)
View(RR80_final)

shortform_EC <- bind_rows(shortform_EC, (RR80_final %>% select(c(CAS, EC_estimate))))
longform_EC <- bind_rows(longform_EC, RR80_final)

#RR65####
RR65<- data_for_eval1 %>% filter(Cluster == "RR65")
RR65$CAS[which(is.na(RR65$Concentration))] # "77-93-0"

RR65 <- RR65 %>% select(CAS, Log_P) %>% column_to_rownames("CAS")
RR65 <- as.matrix(RR65)
RR65_dist <- pdist(RR65, metric = "euclidean", p = 2)

t_RR65_dist <- t(RR65_dist)

RR65
t_RR65_dist

RR65_1 <- as.data.frame(cbind(RR65, t_RR65_dist)) %>% rownames_to_column() %>% rename("CAS" = "rowname") %>%
  gather(c("V2":"V6"), key = "CAS_2", value = "distance")
RR65_1$CAS_2 <- gsub("V2", "14064-10-9", RR65_1$CAS_2)
RR65_1$CAS_2 <- gsub("V3", "105-53-3", RR65_1$CAS_2)
RR65_1$CAS_2 <- gsub("V4", "383-63-1", RR65_1$CAS_2)
RR65_1$CAS_2 <- gsub("V5", "121-75-5", RR65_1$CAS_2)
RR65_1$CAS_2 <- gsub("V6", "77-93-0", RR65_1$CAS_2)

RR65_1_estimates <- RR65_1 %>% filter(CAS %in% c("77-93-0"))
RR65_1_estimates$distance <- RR65_1_estimates$distance +1
RR65_1_estimates$weight <- 1/RR65_1_estimates$distance

#bind in EC estimates
RR65 <- data_for_eval1 %>% filter(Cluster == "RR65") %>% select(CAS, Concentration) %>% rename("CAS_2" = "CAS")
RR65_estimates <- left_join(RR65_1_estimates, RR65) %>% filter(!CAS_2 %in% c("77-93-0"))
RR65_estimates$logCon <- log10(RR65_estimates$Concentration)

RR65_estimates$Con_dist_wt <- RR65_estimates$logCon * RR65_estimates$weight
RR65_estimates_1 <- RR65_estimates %>% group_by(CAS) %>% slice_max(weight, n = 10) %>% summarize(sum_con_dist = sum(Con_dist_wt), sum_wt = sum(weight))
RR65_estimates_1$log_estimate <- RR65_estimates_1$sum_con_dist/ RR65_estimates_1$sum_wt
RR65_estimates_1$EC_estimate <- 10^(RR65_estimates_1$log_estimate)

#bind together all
RR65_final <- left_join((RR65_estimates %>% group_by(CAS) %>% slice_max(weight, n = 10)), RR65_estimates_1)
View(RR65_final)

shortform_EC <- bind_rows(shortform_EC, (RR65_final %>% select(c(CAS, EC_estimate))))
longform_EC <- bind_rows(longform_EC, RR65_final)

write_xlsx(shortform_EC, "RA_short.xlsx")
write_xlsx(longform_EC, "RA_long.xlsx")


##add in RA data####
milwaukee_with_EC4 <- left_join(milwaukee_chem_list, longform_EC)%>% filter(!is.na(EC_estimate))
milwaukee_with_EC4$EC_Source <- "LogP Read-across"
  
milwaukee_with_EC3$Log_P <- as.double(milwaukee_with_EC3$Log_P)
  
milwaukee_with_EC5 <- bind_rows(milwaukee_with_EC3, milwaukee_with_EC4)
View(milwaukee_with_EC5)

#QSAR####
#fill in data gaps with QSAR
list(unique(annotated_data_w_physchem$Notes))
data_for_qsar <- annotated_data_w_physchem %>% filter(Notes == "QSAR")
View(data_for_qsar)

ECOSAR <- read_excel("TOX_QSAR\\ECOSAR_Output.xlsx") %>% filter(INSIDE_AD == "Yes") %>% select(CAS, ECOSAR_pred_mgL)
TEST <- read_excel("TOX_QSAR\\TEST_output.xlsx") %>% select(CAS, TEST_pred_mgL)
VEGA_FHM_KNN <- read_excel("TOX_QSAR\\VEGA_Estimates.xlsx") %>% 
  filter(FHM_KNN_ACF > 0.75) %>% select(CAS, FHM_KNN_mgL)
VEGA_KNN <- read_excel("TOX_QSAR\\VEGA_Estimates.xlsx") %>% 
  filter(KNN_ACF > 0.75) %>% select(CAS, KNN_mgL)
VEGA_NIC <- read_excel("TOX_QSAR\\VEGA_Estimates.xlsx") %>% 
  filter(NIC_ACF > 0.75) %>% select(CAS, NIC_mgL)

data_for_qsar_1 <- left_join(data_for_qsar, ECOSAR)
data_for_qsar_2 <- left_join(data_for_qsar_1, TEST)
data_for_qsar_3 <- left_join(data_for_qsar_2, VEGA_FHM_KNN)
data_for_qsar_4 <- left_join(data_for_qsar_3, VEGA_KNN)
data_for_qsar_5 <- left_join(data_for_qsar_4, VEGA_NIC)

gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}

data_for_qsar_EC_est <- data_for_qsar_5 %>% gather(ECOSAR_pred_mgL:NIC_mgL, key = QSAR_type, value = EC_est) 
data_for_qsar_EC_est$EC_est <- as.numeric(data_for_qsar_EC_est$EC_est)
data_for_qsar_EC_est1 <- data_for_qsar_EC_est %>% group_by(CAS) %>%
  summarize("EC_estimate" = gm_mean(EC_est))

View(data_for_qsar_EC_est1)

QSAR_output <- left_join(data_for_qsar_5, data_for_qsar_EC_est1) %>% select(-c(Concentration:Notes))
names(QSAR_output)

#bind to dataset

milwaukee_with_EC6 <- left_join(milwaukee_chem_list, QSAR_output)%>% filter(!is.na(EC_estimate))
View(milwaukee_with_EC6)
milwaukee_with_EC6$EC_Source <- "QSAR"

milwaukee_with_EC7 <- bind_rows(milwaukee_with_EC5, milwaukee_with_EC6)
View(milwaukee_with_EC7)
list(unique(milwaukee_with_EC7$CAS))

write_xlsx(milwaukee_with_EC7, "milwaukee_EC_forSI.xlsx")

milwaukee_with_EC8 <- milwaukee_with_EC7 %>% select(CAS, Chemical, EC_estimate) %>% distinct()
View(milwaukee_with_EC8)

write_xlsx(milwaukee_with_EC8, "milwaukee_EC_for_MCR_needsannotation.xlsx")
