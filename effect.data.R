library(dplyr)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(factoextra)
library(gridExtra)
library(readxl)
library(FSA)
library(car)
library(FSA)
library(DescTools)
library(rcompanion)
library(multcompView)
library(fuzzyjoin)
library(PMCMRplus)

#Site-specific Effect Data
effectdata <- read_excel("C:\\Users\\erinm\\OneDrive\\Desktop\\GLRI Project\\Milwauke\\Data for Papier\\Effect Data\\Effect_Data_for_all_Analyses.xlsx") 
effectdata$Site <- factor(effectdata$Site, levels = c("KKL", "MEC", "MEF", "MET", "MIE", "MIM", "MIP", "UCJ","CCM", "JIP", "MIN","MED"))
names(effectdata)
getwd()
setwd("C:\\Users\\erinm\\OneDrive\\Desktop\\GLRI Project\\Milwauke\\Data for Papier\\Effect Data")

################ qPCR #############################
#2017 

qPCR.2017 <- effectdata %>% filter(`Assay Name`== "qPCR", Year == "2017")

CYP1A1.2017 <- qPCR.2017 %>% filter(Measurement %in% c("CYP1A1_intestine", "CYP1A1_liver"))
CYP1A1.2017$Organ <- ifelse(CYP1A1.2017$Measurement == "CYP1A1_intestine", "Intestine",
                            "Liver")
CYP1A1.2017$Measurement <- "CYP1A1"
View(CYP1A1.2017)

fig.1 <- CYP1A1.2017%>% ggbarplot(x = "Site", y = "Effect", fill = "Organ", color = "black", add = "mean_se", palette = c("grey", "black"), xlab = "", ylab = "mRNA (relative # of copies)", 
                                       title = "A", position = position_dodge(0.9), x.text.angle = 45)+ font("title", size = 20) + font("legend.text", size = 15) + font("legend.title", size = 20)

CYP2AD6 <- qPCR.2017 %>% filter(Measurement == "CYP2AD6_intestine")

fig.2 <- CYP2AD6 %>% ggbarplot(x = "Site", y = "Effect", fill = "grey", color = "black", add = "mean_se", xlab = "", ylab = "mRNA (relative # of copies)", 
                                  title = "C", position = position_dodge(0.9), x.text.angle = 45)+ font("title", size = 20)

CYP2N13 <-  qPCR.2017 %>% filter(Measurement == "CYP2N13_intestine")

fig.3 <- CYP2N13 %>% ggbarplot(x = "Site", y = "Effect", fill = "grey", color = "black", add = "mean_se", xlab = "", ylab = "", 
                               title = "D", position = position_dodge(0.9), x.text.angle = 45)+ font("title", size = 20)

CYP3A.2017 <- qPCR.2017 %>% filter(Measurement %in% c("CYP3A_liver", "CYP3A_intestine"))
CYP3A.2017$Organ <- ifelse(CYP3A.2017$Measurement == "CYP3A_liver", "Liver",
                           "Intestine")

fig.4 <- CYP3A.2017%>% ggbarplot(x = "Site", y = "Effect", fill = "Organ", color = "black", add = "mean_se", palette = c("grey", "black"), xlab = "", ylab = "", 
                                  title = "B", position = position_dodge(0.9), x.text.angle = 45) + font("title", size = 20)+ font("legend.text", size = 15) + font("legend.title", size = 20)

UGT.2017 <- qPCR.2017 %>% filter(Measurement %in% c("UGT1A1_liver", "UGT1A1_intestine"))
UGT.2017$Organ <- ifelse(UGT.2017$Measurement == "UGT1A1_liver", "Liver", "Intestine")

fig.5 <- UGT.2017 %>% ggbarplot(x = "Site", y = "Effect", fill = "Organ", color = "black", add = "mean_se", palette = c("grey", "black"), xlab = "Site", ylab = "mRNA (relative # of copies)", 
                                 title = "E", position = position_dodge(0.9), legend = NULL, x.text.angle = 45) + rremove("legend")+ font("title", size = 20)

VTG.2017 <- qPCR.2017 %>% filter(Measurement == "VTG")

fig.6 <- VTG.2017 %>% ggbarplot(x = "Site", y = "Effect", fill = "black", color = "black", add = "mean_se", xlab = "Site", ylab = "", 
                                title = "A", position = position_dodge(0.9), x.text.angle = 45)+ font("title", size = 20)

qPCR_2017 <- grid.arrange(fig.1, fig.4, fig.2, fig.3, fig.5)
ggsave("qPCR_figures.jpeg", qPCR_2017, height = 15, width = 10)

#2018 
qPCR.2018 <- effectdata %>% filter(`Assay Name` == "qPCR", Year == "2018")
qPCR.2018$Organ <- "liver"
View(qPCR.2018)

CYP1A1.2018 <- qPCR.2018 %>% filter(Measurement == "CYP1A1_liver")
fig.1 <- CYP1A1.2018 %>% ggbarplot(x = "Site", y = "Effect", fill = "black", color = "black", add = "mean_se", xlab = "Site", ylab = "mRNA (relative # of copies)", 
                                title = "A", position = position_dodge(0.9), x.text.angle = 45) + font("title", size = 17)

VTG.2018 <- qPCR.2018 %>% filter(Measurement == "VTG")
fig.2 <- VTG.2018 %>% ggbarplot(x = "Site", y = "Effect", fill = "black", color = "black", add = "mean_se", xlab = "Site", ylab = "mRNA (relative # of copies)", 
                                   title = "B", position = position_dodge(0.9), x.text.angle = 45)+ font("title", size = 17)

qPCR_2018 <- grid.arrange(fig.1, nrow = 1)

ggsave("qPCR_2018.jpeg", qPCR_2018, height = 5, width = 10)

#VTG Fig

vtg <- grid.arrange(fig.6, fig.2, nrow = 1)
ggsave("vtg.jpeg", vtg, height = 5, width = 10)

#identify statistically significant effects in QPCR####
#qPCR_2017####
names(qPCR.2017)
qPCR.2017$Site <- factor(qPCR.2017$Site, levels = c("MED", "MEC", "MEF", "MET", "MIE", "MIM", "MIP", "UCJ"))

cyp1a1_liver <- qPCR.2017 %>% filter(Measurement == "CYP1A1_liver") %>% select(Site, Effect)
cyp1a1_intestine <- qPCR.2017 %>% filter(Measurement == "CYP1A1_intestine") %>% select(Site, Effect)
cyp2n13 <- qPCR.2017 %>% filter(Measurement == "CYP2N13_intestine") %>% select(Site, Effect)
cyp2ad6 <- qPCR.2017 %>% filter(Measurement == "CYP2AD6_intestine") %>% select(Site, Effect)
cyp3a_liver <- qPCR.2017 %>% filter(Measurement == "CYP3A_liver") %>% select(Site, Effect)
cyp3a_intestine <- qPCR.2017 %>% filter(Measurement == "CYP3A_intestine") %>% select(Site, Effect)
ugt1a1_liver <- qPCR.2017 %>% filter(Measurement == "UGT1A1_liver") %>% select(Site, Effect)
ugt1a1_intestine <- qPCR.2017 %>% filter(Measurement == "UGT1A1_intestine") %>% select(Site, Effect)
vtg <- qPCR.2017 %>% filter(Measurement == "VTG") %>% select(Site, Effect)

#anova
#cyp1a1_liver
View(cyp1a1_liver)
res.aov <- aov(Effect ~ Site, data = cyp1a1_liver)
summary(res.aov) #sd
plot(res.aov,1)
leveneTest(Effect ~ Site, data = cyp1a1_liver) #variance not homogenous
plot(res.aov, 2)
aov_residuals <- residuals(object = res.aov)
shapiro.test(x=aov_residuals)#not normal 

cyp1a1_liver$log_Effect <- log10(cyp1a1_liver$Effect) # try with log transformation
res.aov <- aov(log_Effect ~ Site, data = cyp1a1_liver)
summary(res.aov) #sd
plot(res.aov,1)
leveneTest(log_Effect ~ Site, data = cyp1a1_liver) #variance not homogenous
plot(res.aov, 2)
aov_residuals <- residuals(object = res.aov)
shapiro.test(x=aov_residuals)#normal 

kruskalTest(Effect ~ Site, data = cyp1a1_liver) # sigdif
cyp1a1_liver_kw <- kwManyOneDunnTest(Effect~Site, data = cyp1a1_liver,
         p.adjust.method="bonferroni")
summary(cyp1a1_liver_kw)

#cyp1a1_intestine
cyp1a1_intestine$log_Effect <- log10(cyp1a1_intestine$Effect) # try with log transformation
res.aov <- aov(Effect ~ Site, data = cyp1a1_intestine)
summary(res.aov) #sd
plot(res.aov,1)
leveneTest(Effect ~ Site, data = cyp1a1_intestine) #variance not homogenous
plot(res.aov, 2)
aov_residuals <- residuals(object = res.aov)
shapiro.test(x=aov_residuals)#not normal 

res.aov <- aov(log_Effect ~ Site, data = cyp1a1_intestine) #try with log_Effect
summary(res.aov) #sd
plot(res.aov,1)
leveneTest(log_Effect ~ Site, data = cyp1a1_intestine) #variance not homogenous
plot(res.aov, 2)
aov_residuals <- residuals(object = res.aov)
shapiro.test(x=aov_residuals)#normal 
DunnettTest(x=cyp1a1_intestine$log_Effect, g=cyp1a1_intestine$Site, control = "MED")

#cyp2ad6
res.aov <- aov(Effect ~ Site, data = cyp2ad6)
summary(res.aov) #sd
plot(res.aov,1)
leveneTest(Effect ~ Site, data = cyp2ad6) #variance not homogenous
plot(res.aov, 2)
aov_residuals <- residuals(object = res.aov)
shapiro.test(x=aov_residuals)#not normal 

cyp2ad6$log_Effect <- log10(cyp2ad6$Effect) # try with log transformation
res.aov <- aov(log_Effect ~ Site, data = cyp2ad6)
summary(res.aov) #sd
plot(res.aov,1)
leveneTest(log_Effect ~ Site, data = cyp2ad6) #variance not homogenous
plot(res.aov, 2)
aov_residuals <- residuals(object = res.aov)
shapiro.test(x=aov_residuals)#not normal 

kruskalTest(Effect ~ Site, data = cyp2ad6) # sigdif
cyp2ad6_liver_kw <- kwManyOneDunnTest(Effect~Site, data = cyp2ad6,
                                     p.adjust.method="bonferroni")
summary(cyp2ad6_liver_kw)

#cyp2n13
res.aov <- aov(Effect ~ Site, data = cyp2n13)
summary(res.aov) #sd
plot(res.aov,1)
leveneTest(Effect ~ Site, data = cyp2n13) #variance not homogenous
plot(res.aov, 2)
aov_residuals <- residuals(object = res.aov)
shapiro.test(x=aov_residuals)#not normal 

cyp2n13$log_Effect <- log10(cyp2n13$Effect) # try with log transformation
res.aov <- aov(log_Effect ~ Site, data = cyp2n13)
summary(res.aov) #sd
plot(res.aov,1)
leveneTest(log_Effect ~ Site, data = cyp2n13) #variance homogenous
plot(res.aov, 2)
aov_residuals <- residuals(object = res.aov)
shapiro.test(x=aov_residuals)#not normal 

kruskalTest(Effect ~ Site, data = cyp2n13) # sigdif
cyp2n13_kw <- kwManyOneDunnTest(Effect~Site, data = cyp2n13,
                                      p.adjust.method="bonferroni")
summary(cyp2n13_kw)


#cyp3a_liver
res.aov <- aov(Effect ~ Site, data = cyp3a_liver)
summary(res.aov) #notsd
plot(res.aov,1)
leveneTest(Effect ~ Site, data = cyp3a_liver) #variance homogenous
plot(res.aov, 2)
aov_residuals <- residuals(object = res.aov)
shapiro.test(x=aov_residuals)#not normal 

cyp3a_liver$log_Effect <- log10(cyp3a_liver$Effect) # try with log transformation
res.aov <- aov(log_Effect ~ Site, data = cyp3a_liver)
summary(res.aov) #not sd
plot(res.aov,1)
leveneTest(log_Effect ~ Site, data = cyp3a_liver) #variance homogenous
plot(res.aov, 2)
aov_residuals <- residuals(object = res.aov)
shapiro.test(x=aov_residuals)#not normal 
DunnettTest(x=cyp3a_liver$log_Effect, g=cyp3a_liver$Site, control = "MED")

#cyp3a_intestine
res.aov <- aov(Effect ~ Site, data = cyp3a_intestine)
summary(res.aov) #sd
plot(res.aov,1)
leveneTest(Effect ~ Site, data = cyp3a_intestine) #variance not homogenous
plot(res.aov, 2)
aov_residuals <- residuals(object = res.aov)
shapiro.test(x=aov_residuals)#not normal 

cyp3a_intestine$log_Effect <- log10(cyp3a_intestine$Effect) # try with log transformation
res.aov <- aov(log_Effect ~ Site, data = cyp3a_intestine)
summary(res.aov) #sd
plot(res.aov,1)
leveneTest(log_Effect ~ Site, data = cyp3a_intestine) #variance homogenous
plot(res.aov, 2)
aov_residuals <- residuals(object = res.aov)
shapiro.test(x=aov_residuals)#not normal 

kruskalTest(Effect ~ Site, data = cyp3a_intestine) # sigdif
cyp3a_int_kw <- kwManyOneDunnTest(Effect~Site, data = cyp3a_intestine,
                                      p.adjust.method="bonferroni")
summary(cyp3a_int_kw)

#ugt1a1_liver
res.aov <- aov(Effect ~ Site, data = ugt1a1_liver)
summary(res.aov) #sd
plot(res.aov,1)
leveneTest(Effect ~ Site, data = ugt1a1_liver) #variance not homogenous
plot(res.aov, 2)
aov_residuals <- residuals(object = res.aov)
shapiro.test(x=aov_residuals)#not normal 

ugt1a1_liver$log_Effect <- log10(ugt1a1_liver$Effect) # try with log transformation
res.aov <- aov(log_Effect ~ Site, data = ugt1a1_liver)
summary(res.aov) #sd
plot(res.aov,1)
leveneTest(log_Effect ~ Site, data = ugt1a1_liver) #variance homogenous
plot(res.aov, 2)
aov_residuals <- residuals(object = res.aov)
shapiro.test(x=aov_residuals)# normal 
DunnettTest(x=ugt1a1_liver$log_Effect, g=ugt1a1_liver$Site, control = "MED")

#ugt1a1_intestine
res.aov <- aov(Effect ~ Site, data = ugt1a1_intestine)
summary(res.aov) #sd
plot(res.aov,1)
leveneTest(Effect ~ Site, data = ugt1a1_intestine) #variance not homogenous
plot(res.aov, 2)
aov_residuals <- residuals(object = res.aov)
shapiro.test(x=aov_residuals)#not normal 

ugt1a1_intestine$log_Effect <- log10(ugt1a1_intestine$Effect) # try with log transformation
res.aov <- aov(log_Effect ~ Site, data = ugt1a1_intestine)
summary(res.aov) #sd
plot(res.aov,1)
leveneTest(log_Effect ~ Site, data = ugt1a1_intestine) #variance homogenous
plot(res.aov, 2)
aov_residuals <- residuals(object = res.aov)
shapiro.test(x=aov_residuals)#not normal 

kruskalTest(Effect ~ Site, data = ugt1a1_intestine) # sigdif
ugt1a1_intestine_kw <- kwManyOneDunnTest(Effect~Site, data = ugt1a1_intestine,
                                  p.adjust.method="bonferroni")
summary(ugt1a1_intestine_kw)


#vtg
res.aov <- aov(Effect ~ Site, data = vtg)
summary(res.aov) #sd
plot(res.aov,1)
leveneTest(Effect ~ Site, data = vtg) #variance not homogenous
plot(res.aov, 2)
aov_residuals <- residuals(object = res.aov)
shapiro.test(x=aov_residuals)#not normal 

vtg$log_Effect <- log10(vtg$Effect) # try with log transformation
res.aov <- aov(log_Effect ~ Site, data = vtg)
summary(res.aov) #sd
plot(res.aov,1)
leveneTest(log_Effect ~ Site, data = vtg) #variance homogenous
plot(res.aov, 2)
aov_residuals <- residuals(object = res.aov)
shapiro.test(x=aov_residuals)#normal 

DunnettTest(x=vtg$log_Effect, g=vtg$Site, control = "MED")

#2018####
list(unique(qPCR.2018$Site))
qPCR.2018$Site <- factor(qPCR.2018$Site, levels = c("MED", "CCM", "JIP", "MEC", "MEF", "MIE", "MIM", "MIN", "MIP", "UCJ"))

cyp1a1_liver <- qPCR.2018 %>% filter(Measurement == "CYP1A1_liver") %>% select(Site, Effect)
vtg <- qPCR.2018 %>% filter(Measurement == "VTG") %>% select(Site, Effect)

#cyp1a1_liver
res.aov <- aov(Effect ~ Site, data = cyp1a1_liver)
summary(res.aov) #sd
plot(res.aov,1)
leveneTest(Effect ~ Site, data = cyp1a1_liver) #variance not homogenous
plot(res.aov, 2)
aov_residuals <- residuals(object = res.aov)
shapiro.test(x=aov_residuals)#not normal 

cyp1a1_liver$log_Effect <- log10(cyp1a1_liver$Effect) # try with log transformation
res.aov <- aov(log_Effect ~ Site, data = cyp1a1_liver)
summary(res.aov) #sd
plot(res.aov,1)
leveneTest(log_Effect ~ Site, data = cyp1a1_liver) #variance nothomogenous
plot(res.aov, 2)
aov_residuals <- residuals(object = res.aov)
shapiro.test(x=aov_residuals)#not normal 

kruskalTest(Effect ~ Site, data = cyp1a1_liver) # sigdif
cyp1a1_liver_kw <- kwManyOneDunnTest(Effect~Site, data = cyp1a1_liver,
                                   p.adjust.method="bonferroni")

summary(cyp1a1_liver_kw)

#vtg
res.aov <- aov(Effect ~ Site, data = vtg)
summary(res.aov) #sd
plot(res.aov,1)
leveneTest(Effect ~ Site, data = vtg) #variance not homogenous
plot(res.aov, 2)
aov_residuals <- residuals(object = res.aov)
shapiro.test(x=aov_residuals)#not normal 

vtg$log_Effect <- log10(vtg$Effect) # try with log transformation
res.aov <- aov(log_Effect ~ Site, data = vtg)
summary(res.aov) #sd
plot(res.aov,1)
leveneTest(log_Effect ~ Site, data = vtg) #variance homogenous
plot(res.aov, 2)
aov_residuals <- residuals(object = res.aov)
shapiro.test(x=aov_residuals)#not normal 

kruskalTest(Effect ~ Site, data = vtg) # sigdif
vtg_2018_kw <- kwManyOneDunnTest(Effect~Site, data = vtg,
                                     p.adjust.method="bonferroni")

summary(vtg_2018_kw)

################## T47 ##################################################### 
effectdata <- read_excel("C:\\Users\\erinm\\OneDrive\\Desktop\\GLRI Project\\Milwauke\\Data for Papier\\Effect Data\\Effect_Data_for_all_Analyses.xlsx") 
T47<- effectdata %>% filter(Measurement == "E2-EQ")
View(T47)
list(unique(T47$Site))
T47$Effect <- as.numeric(T47$Effect)
T47$Year <- as.factor(T47$Year)

#2017
T47_2017 <- T47 %>% filter(Year == "2017")
list(unique(T47_2017$Site))

#compute mean MQ-BK and 3 SD for cut-off
BK_2017 <- T47_2017 %>% filter(Site == "MQ-BK") %>% summarize(mean_Effect = mean(Effect), sd = sd(Effect))
BK_2017$E2_cutoff <- BK_2017$mean_Effect + (3*BK_2017$sd)
View(BK_2017)

T47_2017_1 <- bind_cols((T47_2017 %>% filter(Site != "MQ-BK") %>% group_by(Site) %>% summarize(mean_Effect = mean(Effect))), (BK_2017 %>% select(E2_cutoff)))
T47_2017_1$significant <- ifelse(T47_2017_1$mean_Effect > T47_2017_1$E2_cutoff, "Yes", "No")
View(T47_2017_1)

t47.2017 <- T47 %>% filter(Year == "2017", !Site %in% c("MQ-BK", "UWM"))

fig.1 <- t47.2017 %>% ggbarplot(x = "Site", y = "Effect", fill = "black", color = "black", add = "mean_sd", xlab = "Site", ylab = "E2-Equivalents (ng/L)",
                                x.text.angle = 90, title = "A") + font("title", size = 30) + font("xlab", size = 15) + font("ylab", size = 15) + font("xy.text", size = 12) + font("legend.text", size = 15) + font("legend.title", size = 15) + ylim(0, 2)+
  geom_hline(yintercept = 0.3308273, colour = "red", size = 1, linetype = "dashed")


##2018
T47_2018 <- T47 %>% filter(Year == "2018")
BK_2018 <- T47_2018 %>% filter(Site == "BK") %>% summarize(mean_Effect = mean(Effect), sd = sd(Effect))
BK_2018$E2_cutoff <- BK_2018$mean_Effect + (3*BK_2018$sd)
View(BK_2018)

T47_2018_1 <- bind_cols((T47_2018 %>% filter(Site != "BK") %>% group_by(Site) %>% summarize(mean_Effect = mean(Effect))), (BK_2018 %>% select(E2_cutoff)))
T47_2018_1$significant <- ifelse(T47_2018_1$mean_Effect > T47_2018_1$E2_cutoff, "Yes", "No")
View(T47_2018_1)

t47.2018 <- T47 %>% filter(Year == "2018", Site != "BK")
t47.2018$Site <- as.factor(t47.2018$Site)
View(t47.2018)

#make fig
fig.2 <- t47.2018 %>% ggbarplot(x = "Site", y = "Effect", fill = "black", color = "black", add = "mean_sd", xlab = "Site", ylab = "E2-Equivalents (ng/L)", x.text.angle = 90, lab.size = 10, lab.vjust = -4.5, title = "B") + font("title", size = 30) + font("xlab", size = 15) + font("ylab", size = 15) +
  font("xy.text", size = 12) + font("legend.text", size = 15) + font("legend.title", size = 15) + coord_cartesian(ylim = c(0, 2.0)) + geom_hline(yintercept = 0.4130249, colour = "red", size = 1, linetype = "dashed")

T47_graph <- grid.arrange(fig.1, fig.2, nrow = 1)
ggsave("T47_graph.jpeg", T47_graph, height = 5, width = 10)

##################### RIA #################################################### 
RIA <- effectdata %>% filter(`Assay Name` == "RIA")
View(RIA)
RIA$Sex <- ifelse(RIA$`Sample Type` == "Male Fathead Minnow Serum", "Male", "Female")

male <- RIA %>% filter(Sex == "Male")
male.E <- male %>% filter(Measurement == "17b_estradiol")
male.T <- male %>% filter(Measurement == "Testosterone")

female <- RIA %>% filter(Sex == "Female")
female.E <- female %>% filter(Measurement == "17b_estradiol")

E <- RIA %>% filter(Measurement == "17b_estradiol")

fig.1 <- E %>% ggbarplot(x = "Site", y = "Effect", fill = "Sex", add = "mean_sd", error.plot = "upper_errorbar", xlab = "Site", 
                              ylab = "Estrogen Concentration (ng/mL)", title = "A", palette = c("grey", "black"), position = position_dodge(0.9))

fig.2 <- male.T %>% ggbarplot(x = "Site", y = "Effect", fill = "black", color = "black", add = "mean_sd", error.plot = "upper_errorbar", xlab = "Site", ylab = "Testosterone Concentration (ng/mL)", title = "B") 

RIA_figs <- grid.arrange(fig.1, fig.2, nrow = 1)
ggsave("RIA_figures.jpeg", RIA_figs, height = 5, width = 10)

#female E2
list(unique(female.E$Site))
female.E$Site <- factor(female.E$Site, levels = c("MED", "MEC", "MEF", "MET", "MIE", "MIM", "MIP", "UCJ"))

res.aov <- aov(Effect ~ Site, data = female.E)
summary(res.aov) #sd
plot(res.aov,1)
leveneTest(Effect ~ Site, data = female.E) #variance not homogenous
plot(res.aov, 2)
aov_residuals <- residuals(object = res.aov)
shapiro.test(x=aov_residuals)#not normal 

female.E$log_Effect <- log10((female.E$Effect)) # try with log transformation
res.aov <- aov(log_Effect ~ Site, data = female.E)
summary(res.aov) #sd
plot(res.aov,1)
leveneTest(log_Effect ~ Site, data = female.E) #variance homogenous
plot(res.aov, 2)
aov_residuals <- residuals(object = res.aov)
shapiro.test(x=aov_residuals)#not normal 

kruskalTest(Effect ~ Site, data = female.E) # sigdif
female_E_kw <- kwManyOneDunnTest(Effect~Site, data = female.E,
                                     p.adjust.method="bonferroni")

summary(female_E_kw)


#male E2
male.E$Site <- factor(male.E$Site, levels = c("MED", "MEC", "MEF", "MET", "MIE", "MIM", "MIP", "UCJ"))

res.aov <- aov(Effect ~ Site, data = male.E)
summary(res.aov) #sd
plot(res.aov,1)
leveneTest(Effect ~ Site, data = male.E) #variance not homogenous
plot(res.aov, 2)
aov_residuals <- residuals(object = res.aov)
shapiro.test(x=aov_residuals)#not normal 

male.E$log_Effect <- log10((1+ male.E$Effect)) # try with log transformation
res.aov <- aov(log_Effect ~ Site, data = male.E)
summary(res.aov) #sd
plot(res.aov,1)
leveneTest(log_Effect ~ Site, data = male.E) #variance homogenous
plot(res.aov, 2)
aov_residuals <- residuals(object = res.aov)
shapiro.test(x=aov_residuals)#not normal 

kruskalTest(Effect ~ Site, data = male.E) # sigdif
male_E_kw <- kwManyOneDunnTest(Effect~Site, data = male.E,
                                 p.adjust.method="bonferroni")
summary(male_E_kw)

#T
male.T$Site <- factor(male.T$Site, levels = c("MED", "MEC", "MEF", "MET", "MIE", "MIM", "MIP", "UCJ"))
list(unique(male.T$Site ))
male.T <- male.T %>% filter(!is.na(Effect))
male.T <- male.T %>% filter(!is.na(Site))

res.aov <- aov(Effect ~ Site, data = male.T)
summary(res.aov) #sd
plot(res.aov,1)
leveneTest(Effect ~ Site, data = male.T) #variance not homogenous
plot(res.aov, 2)
aov_residuals <- residuals(object = res.aov)
shapiro.test(x=aov_residuals)#not normal 

male.T$log_Effect <- log10((1 + male.T$Effect)) # try with log transformation
res.aov <- aov(log_Effect ~ Site, data = male.T)
summary(res.aov) #sd
plot(res.aov,1)
leveneTest(log_Effect ~ Site, data = male.T) #variance homogenous
plot(res.aov, 2)
aov_residuals <- residuals(object = res.aov)
shapiro.test(x=aov_residuals)#not normal 
DunnettTest(x=male.T$log_Effect, g=male.T$Site, control = "MED")


############## body condition + GSI #####################################

###2017 Survival###
survival.2017 <- effectdata %>% filter(Measurement == "Survival", Year == "2017", Site != "KKL")
survival.2017$Site <- factor(survival.2017$Site, levels = c("MED", "MEC", "MEF", "MET", "MIE", "MIM", "MIP", "UCJ"))
res.aov <- aov(Effect ~ Site, data = survival.2017)
summary(res.aov) #not sd
plot(res.aov,1)
leveneTest(Effect ~ Site, data = survival.2017) #variance not homogenous
plot(res.aov, 2)
aov_residuals <- residuals(object = res.aov)
shapiro.test(x=aov_residuals)#not normal 

survival.2017$log_Effect <- log10((survival.2017$Effect)) # try with log transformation
res.aov <- aov(log_Effect ~ Site, data = survival.2017)
summary(res.aov) #sd
plot(res.aov,1)
leveneTest(log_Effect ~ Site, data = survival.2017) #variance homogenous
plot(res.aov, 2)
aov_residuals <- residuals(object = res.aov)
shapiro.test(x=aov_residuals)#not normal 

kruskalTest(Effect ~ Site, data = survival.2017) # sigdif


#2017 Weight 
Weight.2017 <- effectdata %>% filter(Measurement == "Body Weight", Year == "2017", Site !="KKL", Effect != 0)
Weight.2017$Sex <- ifelse(Weight.2017$`Sample Type` == "Male Fathead Minnow (whole)", "Male", "Female")

f.W.2017 <- Weight.2017 %>% filter(`Sample Type` == "Female Fathead Minnow (whole)")
res.aov <- aov(Effect ~ Site, data = f.W.2017)
summary(res.aov) #not sd
plot(res.aov,1)
leveneTest(Effect ~ Site, data = f.W.2017) #variance homogenous
plot(res.aov, 2)
aov_residuals <- residuals(object = res.aov)
shapiro.test(x=aov_residuals)#normal 

f.W.2017$log_Effect <- log10((f.W.2017$Effect)) # try with log transformation
res.aov <- aov(log_Effect ~ Site, data = f.W.2017)
summary(res.aov) #sd
plot(res.aov,1)
leveneTest(log_Effect ~ Site, data = f.W.2017) #variance homogenous
plot(res.aov, 2)
aov_residuals <- residuals(object = res.aov)
shapiro.test(x=aov_residuals)#not normal 


m.W.2017 <- Weight.2017 %>% filter(`Sample Type` == "Male Fathead Minnow (whole)")
res.aov <- aov(Effect ~ Site, data = m.W.2017)
summary(res.aov) #not sd
plot(res.aov,1)
leveneTest(Effect ~ Site, data = m.W.2017) #variance homogenous
plot(res.aov, 2)
aov_residuals <- residuals(object = res.aov)
shapiro.test(x=aov_residuals)#normal 


#2017 GSI
GSI.2017 <- effectdata %>% filter(Measurement == "Gonadosomatic Index", Year == "2017")
GSI.2017$Sex <- ifelse(GSI.2017$`Sample Type` == "Male Fathead Minnow (whole)", "Male", "Female")
View(GSI.2017)

f.GSI.2017 <- GSI.2017 %>% filter(`Sample Type` == "Female Fathead Minnow (whole)")
res.aov <- aov(Effect ~ Site, data = f.GSI.2017)
summary(res.aov) #not sd
plot(res.aov,1)
leveneTest(Effect ~ Site, data = f.GSI.2017) #variance homogenous
plot(res.aov, 2)
aov_residuals <- residuals(object = res.aov)
shapiro.test(x=aov_residuals)#normal 

m.GSI.2017 <- GSI.2017 %>% filter(`Sample Type` == "Male Fathead Minnow (whole)")
res.aov <- aov(Effect ~ Site, data = m.GSI.2017)
summary(res.aov) #not sd
plot(res.aov,1)
leveneTest(Effect ~ Site, data = m.GSI.2017) #variance homogenous
plot(res.aov, 2)
aov_residuals <- residuals(object = res.aov)
shapiro.test(x=aov_residuals)#normal 


#2018 Survival
survival.2018 <- effectdata %>% filter(Measurement == "Survival", Year == "2018")
View(survival.2018)
res.aov <- aov(Effect ~ Site, data = survival.2018)
summary(res.aov) #not sd
plot(res.aov,1)
leveneTest(Effect ~ Site, data = survival.2018) #variance homogenous
plot(res.aov, 2)
aov_residuals <- residuals(object = res.aov)
shapiro.test(x=aov_residuals)#normal 

survival.2018$log_Effect <- log10((survival.2018$Effect)) # try with log transformation
res.aov <- aov(log_Effect ~ Site, data = survival.2018)
summary(res.aov) #sd
plot(res.aov,1)
leveneTest(log_Effect ~ Site, data = survival.2018) #variance homogenous
plot(res.aov, 2)
aov_residuals <- residuals(object = res.aov)
shapiro.test(x=aov_residuals)#not normal 

shapiro.test(x=aov_residuals)#normal 

kruskalTest(Effect ~ Site, data = survival.2018)


#2018 Weight
Weight.2018 <- effectdata %>% filter(Year == "2018", Measurement == "Body Weight", Effect != 0)
Weight.2018$Sex <- ifelse(Weight.2018$`Sample Type` == "Male Fathead Minnow (whole)", "Male", "Female")
View(Weight.2018)

f.W.2018 <- Weight.2018 %>% filter(`Sample Type` == "Female Fathead Minnow (whole)")
res.aov <- aov(Effect ~ Site, data = f.W.2018)
summary(res.aov) #not sd
plot(res.aov,1)
leveneTest(Effect ~ Site, data = f.W.2018) #variance homogenous
plot(res.aov, 2)
aov_residuals <- residuals(object = res.aov)
shapiro.test(x=aov_residuals)#normal 

f.W.2018$log_Effect <- log10((f.W.2018$Effect)) # try with log transformation
res.aov <- aov(log_Effect ~ Site, data = f.W.2018)
summary(res.aov) #sd
plot(res.aov,1)
leveneTest(log_Effect ~ Site, data = f.W.2018) #variance homogenous
plot(res.aov, 2)
aov_residuals <- residuals(object = res.aov)
shapiro.test(x=aov_residuals)#not normal 

m.W.2018 <- Weight.2018 %>% filter(`Sample Type` == "Male Fathead Minnow (whole)")
res.aov <- aov(Effect ~ Site, data = m.W.2018)
summary(res.aov) #not sd
plot(res.aov,1)
leveneTest(Effect ~ Site, data = m.W.2018) #variance homogenous
plot(res.aov, 2)
aov_residuals <- residuals(object = res.aov)
shapiro.test(x=aov_residuals)#normal 


#2018 GSI
GSI.2018 <- effectdata %>% filter(Year == "2018", Measurement == "Gonadosomatic Index", Effect != 0)
GSI.2018$Sex <- ifelse(GSI.2018$`Sample Type` == "Male Fathead Minnow (whole)", "Male", "Female")
View(GSI.2018)

f.GSI.2018 <- GSI.2018 %>% filter(`Sample Type` == "Female Fathead Minnow (whole)")
res.aov <- aov(Effect ~ Site, data = f.GSI.2018)
summary(res.aov) #not sd
plot(res.aov,1)
leveneTest(Effect ~ Site, data = f.GSI.2018) #variance homogenous
plot(res.aov, 2)
aov_residuals <- residuals(object = res.aov)
shapiro.test(x=aov_residuals)#normal 

m.GSI.2018 <- GSI.2018 %>% filter(`Sample Type` == "Male Fathead Minnow (whole)")
res.aov <- aov(Effect ~ Site, data = m.GSI.2018)
summary(res.aov) #not sd
plot(res.aov,1)
leveneTest(Effect ~ Site, data = m.GSI.2018) #variance homogenous
plot(res.aov, 2)
aov_residuals <- residuals(object = res.aov)
shapiro.test(x=aov_residuals)#normal 

#2017 deformities####
list(unique(effectdata$Measurement))
deformities_2017 <- effectdata %>%  filter(Year == "2017", Site !="KKL", Measurement == "Deformity/Injury")
res.aov <- aov(Effect ~ Site, data = deformities_2017)
summary(res.aov) #not sd
plot(res.aov,1)
leveneTest(Effect ~ Site, data = deformities_2017) #variance homogenous
plot(res.aov, 2)
aov_residuals <- residuals(object = res.aov)
shapiro.test(x=aov_residuals)#normal 

deformities_2017$logEffect <- log10((deformities_2017$Effect + 1))
res.aov <- aov(logEffect ~ Site, data = deformities_2017)
summary(res.aov) #not sd
plot(res.aov,1)
leveneTest(logEffect ~ Site, data = deformities_2017) #variance homogenous
plot(res.aov, 2)
aov_residuals <- residuals(object = res.aov)
shapiro.test(x=aov_residuals)#normal 
kruskalTest(Effect ~ Site, data = deformities_2017) # notsd

#2018 deformities####
list(unique(effectdata$Measurement))
deformities_2018 <- effectdata %>%  filter(Year == "2018", Site !="UWM", Measurement == "Deformity/Injury")
res.aov <- aov(Effect ~ Site, data = deformities_2018)
summary(res.aov) #not sd
plot(res.aov,1)
leveneTest(Effect ~ Site, data = deformities_2018) #variance homogenous
plot(res.aov, 2)
aov_residuals <- residuals(object = res.aov)
shapiro.test(x=aov_residuals)#normal 

deformities_2018$logEffect <- log10((deformities_2018$Effect + 1))
res.aov <- aov(logEffect ~ Site, data = deformities_2018)
summary(res.aov) #not sd
plot(res.aov,1)
leveneTest(logEffect ~ Site, data = deformities_2018) #variance homogenous
plot(res.aov, 2)
aov_residuals <- residuals(object = res.aov)
shapiro.test(x=aov_residuals)#normal 

kruskalTest(Effect ~ Site, data = deformities_2018) # notsd

##generate figures 

fig.A <- survival.2017 %>% ggbarplot(x = "Site", y = "Effect", color = "black", 
                                      fill = c("gray"), add = "mean_sd",
                                      xlab = "Site", ylab = "Survival (%)", title = "A", x.text.angle = 45) + font("title", size = 25) + font("xlab", size = 18) + font("ylab", size = 18) + font("legend.text", size = 18)+ font("legend.title", size = 18) +font("xy", size = 15)+ rremove("legend")

fig.B <- survival.2018 %>% ggbarplot(x = "Site", y = "Effect", color = "black", 
                                     fill = c("grey"), add = "mean_sd",
                                     xlab = "Site", ylab = "Survival (%)", title = "B", x.text.angle = 45) + font("title", size = 25) + font("xlab", size = 18) + font("ylab", size = 18) + font("legend.text", size = 18)+ font("legend.title", size = 18) +font("xy", size = 15)+ rremove("legend")

fig.C <- deformities_2017 %>% ggbarplot(x = "Site", y = "Effect", color = "black", 
                                     fill = c("gray"), add = "mean_sd",
                                     xlab = "Site", ylab = "Deformities/Injuries (%)", title = "C", x.text.angle = 45) + font("title", size = 25) + font("xlab", size = 18) + font("ylab", size = 18) + font("legend.text", size = 18)+ font("legend.title", size = 18) +font("xy", size = 15)+ rremove("legend") +
  ylim(0, 100)

fig.D <- deformities_2018 %>% ggbarplot(x = "Site", y = "Effect", color = "black", 
                                     fill = c("grey"), add = "mean_sd",
                                     xlab = "Site", ylab = "Deformities/Injuries (%)", title = "D", x.text.angle = 45) + font("title", size = 25) + font("xlab", size = 18) + font("ylab", size = 18) + font("legend.text", size = 18)+ font("legend.title", size = 18) +font("xy", size = 15)+ rremove("legend") + ylim(0, 100)


fig.E <- Weight.2017 %>% ggbarplot(x = "Site", y = "Effect", fill = "Sample Type", color = "black", 
                                   palette = c("white", "black"), add = "mean_sd",
                                   xlab = "Site", ylab = "Wet Weight (g)", title = "E", position = position_dodge(0.9), x.text.angle = 45) + font("title", size = 25) + font("xlab", size = 18) + font("ylab", size = 18) + font("legend.text", size = 18)+ font("legend.title", size = 18) +font("xy", size = 15)

fig.F <- Weight.2018 %>% ggbarplot(x = "Site", y = "Effect", fill = "Sample Type", color = "black", 
                                   palette = c("white", "black"), add = "mean_sd",
                                   xlab = "Site", ylab = "Wet Weight (g)", title = "F", position = position_dodge(0.9), x.text.angle = 45) + font("title", size = 25) + font("xlab", size = 18) + font("ylab", size = 18) + font("legend.text", size = 18) + font("legend.title", size = 18) +font("xy", size = 15)

fig.G <- GSI.2017 %>% ggbarplot(x = "Site", y = "Effect", fill = "Sample Type", color = "black", 
                                   palette = c("white", "black"), add = "mean_sd",
                                   xlab = "Site", ylab = "Gonadosomatic Index (%)", title = "G", position = position_dodge(0.9), x.text.angle = 45) + font("title", size = 25) + font("xlab", size = 18) + font("ylab", size = 18) + font("legend.text", size = 18) + font("legend.title", size = 18)+font("xy", size = 15)


fig.H <- GSI.2018 %>% ggbarplot(x = "Site", y = "Effect", fill = "Sample Type", color = "black", 
                                palette = c("white", "black"), add = "mean_sd",
                                xlab = "Site", ylab = "Gonadosomatic Index (%)", title = "H", position = position_dodge(0.9), x.text.angle = 45) + font("title", size = 25) + font("xlab", size = 18) + font("ylab", size = 18) + font("legend.text", size = 18) + font("legend.title", size = 18) +font("xy", size = 15)

body_condition <- grid.arrange(fig.A, fig.B, fig.C, fig.D, fig.E, fig.F,fig.G, fig.H, nrow = 4)
ggsave("body_condition.jpeg", body_condition, height = 25, width = 20)

