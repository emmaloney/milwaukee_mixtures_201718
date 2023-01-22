#Random Forest Mixture Groups#

library(tidyverse)
library(dplyr)
library(ggplot2)
library(randomForest) 
library(permimp)
library(party)
library(readxl)
library(writexl)
library(data.table) # don't think I need this after all but uncomment if something doesn't work
library(pdp)
library(factoextra)
library(gridExtra)
library(standardize)
library(ggpubr)
library(beepr)
getwd()
setwd("C:\\Users\\erinm\\OneDrive\\Desktop\\GLRI Project\\Milwauke\\Data for Papier\\Mixtures MS\\Random_Forest\\Random_Forest")

#ahR ATG####
names(ahr)
ahr_sub <- ahr[,c(2:222)]
ahr_sub[,1] <- round(ahr_sub[,1],2)
View(ahr_sub)

p <- ncol(ahr_sub) - 1

#cforest - all variables
set.seed(3)
cf.all <- cforest(Effect~., data =ahr_sub, control = cforest_unbiased(mtry =ceiling(p/3), ntree = 5000))
cf.pred <- predict(cf.all, OOB = T)
sqrt(mean((ahr_sub$Effect-cf.pred)^2)) #1.223567
mean(abs(ahr_sub$Effect - cf.pred)) #0.8962163

#plot predicted vs. observed
plot(y=cf.pred, x= ahr_sub$Effect)
abline(a=0, b=1)

# Unconditional variable importance (threshold=1)
set.seed(3) 
cf.upi <- permimp(cf.all, conditional=TRUE, threshold=1, thresholdDiagnostics = TRUE)
as.matrix(sort(cf.upi$values, decreasing=T)) # CF, unconditional

# Conditional variable importance (threshold=0.9)
set.seed(3) 
cf.cpi <- permimp(cf.all, conditional=TRUE, threshold=0.90, thresholdDiagnostics = TRUE)
as.matrix(sort(cf.cpi$values, decreasing=T)) # CF, conditional

# cf.cpi$perTree # VI values per tree
plot(cf.cpi, type="bar")

### Loop to recursively eliminate least informative variables

# Pull out response and predictors - columns 5, 7-53
ahr_preds <- ahr[, c(2:222)] # original predictors
ahr_preds[,1] <- round(ahr_preds[,1], 2)
ahr_preds_loop <- ahr_preds 
View(ahr_preds_loop)

p <- ncol(ahr_preds)-1

# Save OOB error, variable importance
VI.list <- NULL
var.list <- rep(NA,p)

OOB.rmse <- rep(NA,p)
OOB.mae <- rep(NA,p)

train.rmse <- rep(NA,p)
train.mae <- rep(NA,p)

#run RF - conditional
set.seed(3) # Set seed for reproducibility

for(i in 1:p){
  print(i)
  
  nump <- ncol(ahr_preds_loop)-1 # Number predictors
  
  set.seed(3) 
  cf <- cforest(Effect~., data=ahr_preds_loop, control=cforest_unbiased(mtry=ceiling(nump/3), ntree=5000))
  cf.pred <- predict(cf, OOB=T)
  cf.pred.tr <- predict(cf, OOB=F)
  
  # OOB error
  OOB.rmse[i] <- sqrt(mean((ahr_preds$Effect-cf.pred)^2)) # OOB rmse
  OOB.mae[i] <- mean(abs(ahr_preds$Effect-cf.pred)) # oob mae
  
  # Training error
  train.rmse[i] <- sqrt(mean((ahr_preds$Effect-cf.pred.tr)^2)) # training rmse
  train.mae[i] <- mean(abs(ahr_preds$Effect-cf.pred.tr)) #  training mae
  
  
  # Variable importance
  set.seed(3) 
  cf.cpi <- permimp(cf, conditional=TRUE, threshold=0.9) # set threshold to 1 for unconditional
  
  VI.list[[i]] <- as.matrix(sort(cf.cpi$values, decreasing=F))
  var.list[i] <- row.names(VI.list[[i]])[1] # Least informative variable to remove
  
  
  ahr_preds_loop <- ahr_preds_loop[,!colnames(ahr_preds_loop)==var.list[i]] # remove variable
  
}
beep()

# Plot OOB error - minimum is the best subset of predictors
plot(OOB.rmse, xlab="Iteration", type='l')
points(OOB.rmse, pch=16)
plot(OOB.mae, xlab="Iteration", type='l')
points(OOB.mae, pch=16)

# training error (not to be used for optimizing number predictors)
plot(train.rmse, xlab="Iteration", type='l')
points(train.rmse, pch=16)

which(OOB.rmse==min(OOB.rmse)) #215
which(OOB.mae==min(OOB.mae)) #215

var.list[which(OOB.mae==min(OOB.mae)):length(OOB.mae)]

#double check with rounded error#
OOB.rmse_1 <- round(OOB.rmse, 2)
OOB.mae_1 <- round(OOB.mae, 2)

which(OOB.rmse_1==min(OOB.rmse_1))
which(OOB.mae_1==min(OOB.mae_1)) 

var.list[which(OOB.mae_1==min(OOB.mae_1)):length(OOB.mae_1)] # same

# Save output - name files in a way that lets you know what differs between different runs of the algorithm

saveRDS(VI.list, "cf_VI_at_iter_thresh90_ahr.rds") # R object (list) with the variable importances at each iteration of the algorithm

RFE_info <- data.frame(Worst_Var=var.list, OOB_rmse=OOB.rmse, OOB_mae=OOB.mae, Train_rmse=train.rmse, Train_mae=train.mae)

write.csv(RFE_info, "cf_RFE_info_thresh90_ahr.csv", row.names = FALSE)

#set mtry = # of predictors left in model / 3 (nump/3) - and round it up to solid number 
# Refit final model
ahr_final <- ahr %>% dplyr::select(Effect, Acetaminophen,
                                   AR_Activation, androstenedione_Inhibition, Anthraquinone,
                                   Cotinine, Caffeine)

set.seed(3)
cf.final <- cforest(Effect~., data=ahr_final, control=cforest_unbiased(mtry=3, ntree=5000))

# OOB predictions
ahr_final.pred.oob <- data.frame(predict(cf.final, OOB=T)) %>% rename(Predicted_OOB=Effect)
ahr_final_final_wPred <- cbind(ahr_final, ahr_final.pred.oob)

# OOB error
mean(abs(ahr_final$Effect-ahr_final.pred.oob$Predicted_OOB)) #3.213091
sqrt(mean((ahr_final$Effect-ahr_final.pred.oob$Predicted_OOB)^2)) #3.931485

# Predicted (OOB) vs. observed
plot(Predicted_OOB~Effect, data=ahr_final_final_wPred, pch=16)
abline(a=0, b=1)

#PXR_Cis ATG####
names(pxr_cis)
pxr_cis_sub <- pxr_cis[,c(2:222)]
pxr_cis_sub[,1] <- round(pxr_cis_sub[,1],2)
View(pxr_cis_sub)

p <- ncol(pxr_cis_sub) - 1

#cforest - all variables
set.seed(3)
cf.all <- cforest(Effect~., data =pxr_cis_sub, control = cforest_unbiased(mtry =ceiling(p/3), ntree = 5000))
cf.pred <- predict(cf.all, OOB = T)
sqrt(mean((pxr_cis_sub$Effect-cf.pred)^2)) #1.617086
mean(abs(pxr_cis_sub$Effect - cf.pred)) #1.28555

#plot predicted vs. observed
plot(y=cf.pred, x= pxr_cis_sub$Effect)
abline(a=0, b=1)

# Unconditional variable importance (threshold=1)
set.seed(3) 
cf.upi <- permimp(cf.all, conditional=TRUE, threshold=1, thresholdDiagnostics = TRUE)
as.matrix(sort(cf.upi$values, decreasing=T)) # CF, unconditional

# Conditional variable importance (threshold=0.9)
set.seed(3) 
cf.cpi <- permimp(cf.all, conditional=TRUE, threshold=0.90, thresholdDiagnostics = TRUE)
as.matrix(sort(cf.cpi$values, decreasing=T)) # CF, conditional

# cf.cpi$perTree # VI values per tree
plot(cf.cpi, type="bar")

### Loop to recursively eliminate least informative variables

# Pull out response and predictors - columns 5, 7-53
pxr_cis_preds <- pxr_cis[, c(2:222)] # original predictors
pxr_cis_preds[,1] <- round(pxr_cis_preds[,1], 2)
pxr_cis_preds_loop <- pxr_cis_preds 
View(pxr_cis_preds_loop)

p <- ncol(pxr_cis_preds)-1

# Save OOB error, variable importance
VI.list <- NULL
var.list <- rep(NA,p)

OOB.rmse <- rep(NA,p)
OOB.mae <- rep(NA,p)

train.rmse <- rep(NA,p)
train.mae <- rep(NA,p)

#run RF - conditional
set.seed(3) # Set seed for reproducibility

for(i in 1:p){
  print(i)
  
  nump <- ncol(pxr_cis_preds_loop)-1 # Number predictors
  
  set.seed(3) 
  cf <- cforest(Effect~., data=pxr_cis_preds_loop, control=cforest_unbiased(mtry=ceiling(nump/3), ntree=5000))
  cf.pred <- predict(cf, OOB=T)
  cf.pred.tr <- predict(cf, OOB=F)
  
  # OOB error
  OOB.rmse[i] <- sqrt(mean((pxr_cis_preds$Effect-cf.pred)^2)) # OOB rmse
  OOB.mae[i] <- mean(abs(pxr_cis_preds$Effect-cf.pred)) # oob mae
  
  # Training error
  train.rmse[i] <- sqrt(mean((pxr_cis_preds$Effect-cf.pred.tr)^2)) # training rmse
  train.mae[i] <- mean(abs(pxr_cis_preds$Effect-cf.pred.tr)) #  training mae
  
  
  # Variable importance
  set.seed(3) 
  cf.cpi <- permimp(cf, conditional=TRUE, threshold=0.9) # set threshold to 1 for unconditional
  
  VI.list[[i]] <- as.matrix(sort(cf.cpi$values, decreasing=F))
  var.list[i] <- row.names(VI.list[[i]])[1] # Least informative variable to remove
  
  
  pxr_cis_preds_loop <- pxr_cis_preds_loop[,!colnames(pxr_cis_preds_loop)==var.list[i]] # remove variable
  
}
beep()

# Plot OOB error - minimum is the best subset of predictors
plot(OOB.rmse, xlab="Iteration", type='l')
points(OOB.rmse, pch=16)
plot(OOB.mae, xlab="Iteration", type='l')
points(OOB.mae, pch=16)

# training error (not to be used for optimizing number predictors)
plot(train.rmse, xlab="Iteration", type='l')
points(train.rmse, pch=16)

which(OOB.rmse==min(OOB.rmse)) #218
which(OOB.mae==min(OOB.mae)) #216

var.list[which(OOB.mae==min(OOB.mae)):length(OOB.mae)]
# [1] "beta-Sitosterol" 

OOB.rmse_1 <- round(OOB.rmse, 2)
OOB.mae_1 <- round(OOB.mae, 2)

which(OOB.rmse_1==min(OOB.rmse_1)) #218
which(OOB.mae_1==min(OOB.mae_1)) #218

var.list[which(OOB.mae_1==min(OOB.mae_1)):length(OOB.mae_1)] # same

# Save output - name files in a way that lets you know what differs between different runs of the algorithm

saveRDS(VI.list, "cf_VI_at_iter_thresh90_pxr_cis.rds") # R object (list) with the variable importances at each iteration of the algorithm

RFE_info <- data.frame(Worst_Var=var.list, OOB_rmse=OOB.rmse, OOB_mae=OOB.mae, Train_rmse=train.rmse, Train_mae=train.mae)

write.csv(RFE_info, "cf_RFE_info_thresh90_pxr_cis.csv", row.names = FALSE)

#set mtry = # of predictors left in model / 3 (nump/3) - and round it up to solid number 
# Refit final model
pxr_cis_final <- pxr_cis %>% dplyr::select(Effect,"beta-Sitosterol")

set.seed(3)
cf.final <- cforest(Effect~., data=pxr_cis_final, control=cforest_unbiased(mtry=1, ntree=5000))

# OOB predictions
pxr_cis_final.pred.oob <- data.frame(predict(cf.final, OOB=T)) %>% rename(Predicted_OOB=Effect)
pxr_cis_final_wPred <- cbind(pxr_cis_final, pxr_cis_final.pred.oob)

# OOB error
mean(abs(pxr_cis_final$Effect-pxr_cis_final.pred.oob$Predicted_OOB)) #0.5229215
sqrt(mean((pxr_cis_final$Effect-pxr_cis_final.pred.oob$Predicted_OOB)^2)) #0.8773672

# Predicted (OOB) vs. observed
plot(Predicted_OOB~Effect, data=pxr_cis_final_wPred, pch=16)
abline(a=0, b=1)

#pxr_trans ATG####
names(pxr_trans)
pxr_trans_sub <- pxr_trans[,c(2:222)]
pxr_trans_sub[,1] <- round(pxr_trans_sub[,1],2)
View(pxr_trans_sub)

p <- ncol(pxr_trans_sub) - 1

#cforest - all variables
set.seed(3)
cf.all <- cforest(Effect~., data =pxr_trans_sub, control = cforest_unbiased(mtry =ceiling(p/3), ntree = 5000))
cf.pred <- predict(cf.all, OOB = T)
sqrt(mean((pxr_trans_sub$Effect-cf.pred)^2)) #0.5643476
mean(abs(pxr_trans_sub$Effect - cf.pred)) #0.4455936

#plot predicted vs. observed
plot(y=cf.pred, x= pxr_trans_sub$Effect)
abline(a=0, b=1)

# Unconditional variable importance (threshold=1)
set.seed(3) 
cf.upi <- permimp(cf.all, conditional=TRUE, threshold=1, thresholdDiagnostics = TRUE)
as.matrix(sort(cf.upi$values, decreasing=T)) # CF, unconditional

# Conditional variable importance (threshold=0.9)
set.seed(3) 
cf.cpi <- permimp(cf.all, conditional=TRUE, threshold=0.90, thresholdDiagnostics = TRUE)
as.matrix(sort(cf.cpi$values, decreasing=T)) # CF, conditional

# cf.cpi$perTree # VI values per tree
plot(cf.cpi, type="bar")

### Loop to recursively eliminate least informative variables

# Pull out response and predictors - columns 5, 7-53
pxr_trans_preds <- pxr_trans[, c(2:222)] # original predictors
pxr_trans_preds[,1] <- round(pxr_trans_preds[,1], 2)
pxr_trans_preds_loop <- pxr_trans_preds 
View(pxr_trans_preds_loop)

p <- ncol(pxr_trans_preds)-1

# Save OOB error, variable importance
VI.list <- NULL
var.list <- rep(NA,p)

OOB.rmse <- rep(NA,p)
OOB.mae <- rep(NA,p)

train.rmse <- rep(NA,p)
train.mae <- rep(NA,p)

#run RF - conditional
set.seed(3) # Set seed for reproducibility

for(i in 1:p){
  print(i)
  
  nump <- ncol(pxr_trans_preds_loop)-1 # Number predictors
  
  set.seed(3) 
  cf <- cforest(Effect~., data=pxr_trans_preds_loop, control=cforest_unbiased(mtry=ceiling(nump/3), ntree=5000))
  cf.pred <- predict(cf, OOB=T)
  cf.pred.tr <- predict(cf, OOB=F)
  
  # OOB error
  OOB.rmse[i] <- sqrt(mean((pxr_trans_preds$Effect-cf.pred)^2)) # OOB rmse
  OOB.mae[i] <- mean(abs(pxr_trans_preds$Effect-cf.pred)) # oob mae
  
  # Training error
  train.rmse[i] <- sqrt(mean((pxr_trans_preds$Effect-cf.pred.tr)^2)) # training rmse
  train.mae[i] <- mean(abs(pxr_trans_preds$Effect-cf.pred.tr)) #  training mae
  
  
  # Variable importance
  set.seed(3) 
  cf.cpi <- permimp(cf, conditional=TRUE, threshold=0.9) # set threshold to 1 for unconditional
  
  VI.list[[i]] <- as.matrix(sort(cf.cpi$values, decreasing=F))
  var.list[i] <- row.names(VI.list[[i]])[1] # Least informative variable to remove
  
  
  pxr_trans_preds_loop <- pxr_trans_preds_loop[,!colnames(pxr_trans_preds_loop)==var.list[i]] # remove variable
  
}
beep()

# Plot OOB error - minimum is the best subset of predictors
plot(OOB.rmse, xlab="Iteration", type='l')
points(OOB.rmse, pch=16)
plot(OOB.mae, xlab="Iteration", type='l')
points(OOB.mae, pch=16)

# training error (not to be used for optimizing number predictors)
plot(train.rmse, xlab="Iteration", type='l')
points(train.rmse, pch=16)

which(OOB.rmse==min(OOB.rmse)) #207
which(OOB.mae==min(OOB.mae)) #218

var.list[which(OOB.mae==min(OOB.mae)):length(OOB.mae)]
# [1] "N,N-diethyl-meta-toluamide" "Sterols"                    "beta-Sitosterol"     

#double check with rounded error#
OOB.rmse_1 <- round(OOB.rmse, 2)
OOB.mae_1 <- round(OOB.mae, 2)

which(OOB.rmse_1==min(OOB.rmse_1)) #215
which(OOB.mae_1==min(OOB.mae_1)) #220

var.list[which(OOB.mae_1==min(OOB.mae_1)):length(OOB.mae_1)] # same

# Save output - name files in a way that lets you know what differs between different runs of the algorithm

saveRDS(VI.list, "cf_VI_at_iter_thresh90_pxr_trans.rds") # R object (list) with the variable importances at each iteration of the algorithm

RFE_info <- data.frame(Worst_Var=var.list, OOB_rmse=OOB.rmse, OOB_mae=OOB.mae, Train_rmse=train.rmse, Train_mae=train.mae)

write.csv(RFE_info, "cf_RFE_info_thresh90_pxr_trans.csv", row.names = FALSE)

#set mtry = # of predictors left in model / 3 (nump/3) - and round it up to solid number 
# Refit final model
pxr_trans_final <- pxr_trans %>% dplyr::select(Effect,
                                               "beta-Sitosterol")

set.seed(3)
cf.final <- cforest(Effect~., data=pxr_trans_final, control=cforest_unbiased(mtry=1, ntree=5000))

# OOB predictions
pxr_trans_final.pred.oob <- data.frame(predict(cf.final, OOB=T)) %>% rename(Predicted_OOB=Effect)
pxr_trans_final_wPred <- cbind(pxr_trans_final, pxr_trans_final.pred.oob)

# OOB error
mean(abs(pxr_trans_final$Effect-pxr_trans_final.pred.oob$Predicted_OOB)) #0.3405522
sqrt(mean((pxr_trans_final$Effect-pxr_trans_final.pred.oob$Predicted_OOB)^2)) #0.5289651

# Predicted (OOB) vs. observed
plot(Predicted_OOB~Effect, data=pxr_trans_final_wPred, pch=16)
abline(a=0, b=1)

#ERa_trans####
names(ERa_trans)
ERa_trans_sub <- ERa_trans[,c(2:222)]
ERa_trans_sub[,1] <- round(ERa_trans_sub[,1],2)
View(ERa_trans_sub)

p <- ncol(ERa_trans_sub) - 1

#cforest - all variables
set.seed(3)
cf.all <- cforest(Effect~., data =ERa_trans_sub, control = cforest_unbiased(mtry =ceiling(p/3), ntree = 5000))
cf.pred <- predict(cf.all, OOB = T)
sqrt(mean((ERa_trans_sub$Effect-cf.pred)^2)) #0.3777912
mean(abs(ERa_trans_sub$Effect - cf.pred)) #0.2889635

#plot predicted vs. observed
plot(y=cf.pred, x= ERa_trans_sub$Effect)
abline(a=0, b=1)

# Unconditional variable importance (threshold=1)
set.seed(3) 
cf.upi <- permimp(cf.all, conditional=TRUE, threshold=1, thresholdDiagnostics = TRUE)
as.matrix(sort(cf.upi$values, decreasing=T)) # CF, unconditional

# Conditional variable importance (threshold=0.9)
set.seed(3) 
cf.cpi <- permimp(cf.all, conditional=TRUE, threshold=0.90, thresholdDiagnostics = TRUE)
as.matrix(sort(cf.cpi$values, decreasing=T)) # CF, conditional

# cf.cpi$perTree # VI values per tree
plot(cf.cpi, type="bar")

### Loop to recursively eliminate least informative variables #

# Pull out response and predictors - columns 5, 7-53
ERa_trans_preds <- ERa_trans[, c(2:222)] # original predictors
ERa_trans_preds[,1] <- round(ERa_trans_preds[,1], 2)
ERa_trans_preds_loop <- ERa_trans_preds 
View(ERa_trans_preds_loop)

p <- ncol(ERa_trans_preds)-1

# Save OOB error, variable importance
VI.list <- NULL
var.list <- rep(NA,p)

OOB.rmse <- rep(NA,p)
OOB.mae <- rep(NA,p)

train.rmse <- rep(NA,p)
train.mae <- rep(NA,p)

#run RF - conditional
set.seed(3) # Set seed for reproducibility

for(i in 1:p){
  print(i)
  
  nump <- ncol(ERa_trans_preds_loop)-1 # Number predictors
  
  set.seed(3) 
  cf <- cforest(Effect~., data=ERa_trans_preds_loop, control=cforest_unbiased(mtry=ceiling(nump/3), ntree=5000))
  cf.pred <- predict(cf, OOB=T)
  cf.pred.tr <- predict(cf, OOB=F)
  
  # OOB error
  OOB.rmse[i] <- sqrt(mean((ERa_trans_preds$Effect-cf.pred)^2)) # OOB rmse
  OOB.mae[i] <- mean(abs(ERa_trans_preds$Effect-cf.pred)) # oob mae
  
  # Training error
  train.rmse[i] <- sqrt(mean((ERa_trans_preds$Effect-cf.pred.tr)^2)) # training rmse
  train.mae[i] <- mean(abs(ERa_trans_preds$Effect-cf.pred.tr)) #  training mae
  
  
  # Variable importance
  set.seed(3) 
  cf.cpi <- permimp(cf, conditional=TRUE, threshold=0.9) # set threshold to 1 for unconditional
  
  VI.list[[i]] <- as.matrix(sort(cf.cpi$values, decreasing=F))
  var.list[i] <- row.names(VI.list[[i]])[1] # Least informative variable to remove
  
  
  ERa_trans_preds_loop <- ERa_trans_preds_loop[,!colnames(ERa_trans_preds_loop)==var.list[i]] # remove variable
  
}
beep()

# Plot OOB error - minimum is the best subset of predictors
plot(OOB.rmse, xlab="Iteration", type='l')
points(OOB.rmse, pch=16)
plot(OOB.mae, xlab="Iteration", type='l')
points(OOB.mae, pch=16)

# training error (not to be used for optimizing number predictors)
plot(train.rmse, xlab="Iteration", type='l')
points(train.rmse, pch=16)

which(OOB.rmse==min(OOB.rmse)) #219
which(OOB.mae==min(OOB.mae)) #219

var.list[which(OOB.mae==min(OOB.mae)):length(OOB.mae)]


#double check with rounded error#
OOB.rmse_1 <- round(OOB.rmse, 2)
OOB.mae_1 <- round(OOB.mae, 2)
which(OOB.rmse_1==min(OOB.rmse_1)) #219
which(OOB.mae_1==min(OOB.mae_1)) #219

var.list[which(OOB.mae_1==min(OOB.mae_1)):length(OOB.mae_1)]

#set mtry = # of predictors left in model / 3 (nump/3) - and round it up to solid number 
# Refit final model
era_trans_final <- ERa_trans %>% dplyr::select(Effect, "Anthraquinone", Industrial_PPCP_Narcotics)

set.seed(3)
cf.final <- cforest(Effect~., data=era_trans_final, control=cforest_unbiased(mtry=1, ntree=5000))

# OOB predictions
era_trans_final.pred.oob <- data.frame(predict(cf.final, OOB=T)) %>% rename(Predicted_OOB=Effect)
era_trans_final_wPred <- cbind(era_trans_final, era_trans_final.pred.oob)

# OOB error
mean(abs(era_trans_final$Effect-era_trans_final.pred.oob$Predicted_OOB)) #0.2679142
sqrt(mean((era_trans_final$Effect-era_trans_final.pred.oob$Predicted_OOB)^2)) #0.3637479

# Predicted (OOB) vs. observed
plot(Predicted_OOB~Effect, data=era_trans_final_wPred, pch=16)
abline(a=0, b=1)


#ERE_cis####
names(ERE_cis)
ERE_cis_sub <- ERE_cis[,c(2:222)]
ERE_cis_sub[,1] <- round(ERE_cis_sub[,1],2)
View(ERE_cis_sub)

p <- ncol(ERE_cis_sub) - 1

#cforest - all variables
set.seed(3)
cf.all <- cforest(Effect~., data =ERE_cis_sub, control = cforest_unbiased(mtry =ceiling(p/3), ntree = 5000))
cf.pred <- predict(cf.all, OOB = T)
sqrt(mean((ERE_cis_sub$Effect-cf.pred)^2)) #0.4530935
mean(abs(ERE_cis_sub$Effect - cf.pred)) #0.340722

#plot predicted vs. observed
plot(y=cf.pred, x= ERE_cis_sub$Effect)
abline(a=0, b=1)

# Unconditional variable importance (threshold=1)
set.seed(3) 
cf.upi <- permimp(cf.all, conditional=TRUE, threshold=1, thresholdDiagnostics = TRUE)
as.matrix(sort(cf.upi$values, decreasing=T)) # CF, unconditional

# Conditional variable importance (threshold=0.9)
set.seed(3) 
cf.cpi <- permimp(cf.all, conditional=TRUE, threshold=0.90, thresholdDiagnostics = TRUE)
as.matrix(sort(cf.cpi$values, decreasing=T)) # CF, conditional

# cf.cpi$perTree # VI values per tree
plot(cf.cpi, type="bar")

### Loop to recursively eliminate least informative variables ##

# Pull out response and predictors - columns 5, 7-53
ERE_cis_preds <- ERE_cis[, c(2:222)] # original predictors
ERE_cis_preds[,1] <- round(ERE_cis_preds[,1], 2)
ERE_cis_preds_loop <- ERE_cis_preds 
View(ERE_cis_preds_loop)

p <- ncol(ERE_cis_preds)-1

# Save OOB error, variable importance
VI.list <- NULL
var.list <- rep(NA,p)

OOB.rmse <- rep(NA,p)
OOB.mae <- rep(NA,p)

train.rmse <- rep(NA,p)
train.mae <- rep(NA,p)

#run RF - conditional
set.seed(3) # Set seed for reproducibility

for(i in 1:p){
  print(i)
  
  nump <- ncol(ERE_cis_preds_loop)-1 # Number predictors
  
  set.seed(3) 
  cf <- cforest(Effect~., data=ERE_cis_preds_loop, control=cforest_unbiased(mtry=ceiling(nump/3), ntree=5000))
  cf.pred <- predict(cf, OOB=T)
  cf.pred.tr <- predict(cf, OOB=F)
  
  # OOB error
  OOB.rmse[i] <- sqrt(mean((ERE_cis_preds$Effect-cf.pred)^2)) # OOB rmse
  OOB.mae[i] <- mean(abs(ERE_cis_preds$Effect-cf.pred)) # oob mae
  
  # Training error
  train.rmse[i] <- sqrt(mean((ERE_cis_preds$Effect-cf.pred.tr)^2)) # training rmse
  train.mae[i] <- mean(abs(ERE_cis_preds$Effect-cf.pred.tr)) #  training mae
  
  
  # Variable importance
  set.seed(3) 
  cf.cpi <- permimp(cf, conditional=TRUE, threshold=0.9) # set threshold to 1 for unconditional
  
  VI.list[[i]] <- as.matrix(sort(cf.cpi$values, decreasing=F))
  var.list[i] <- row.names(VI.list[[i]])[1] # Least informative variable to remove
  
  
  ERE_cis_preds_loop <- ERE_cis_preds_loop[,!colnames(ERE_cis_preds_loop)==var.list[i]] # remove variable
  
}
beep()

# Plot OOB error - minimum is the best subset of predictors
plot(OOB.rmse, xlab="Iteration", type='l')
points(OOB.rmse, pch=16)
plot(OOB.mae, xlab="Iteration", type='l')
points(OOB.mae, pch=16)

# training error (not to be used for optimizing number predictors)
plot(train.rmse, xlab="Iteration", type='l')
points(train.rmse, pch=16)

which(OOB.rmse==min(OOB.rmse)) #218
which(OOB.mae==min(OOB.mae)) #218

var.list[which(OOB.mae==min(OOB.mae)):length(OOB.mae)]

#Sterols 

#double check with rounded error#
OOB.rmse_1 <- round(OOB.rmse, 2)
OOB.mae_1 <- round(OOB.mae, 2)
which(OOB.rmse_1==min(OOB.rmse_1)) #218
which(OOB.mae_1==min(OOB.mae_1)) #218

# Save output - name files in a way that lets you know what differs between different runs of the algorithm

saveRDS(VI.list, "cf_VI_at_iter_thresh90_ere_cis.rds") # R object (list) with the variable importances at each iteration of the algorithm

RFE_info <- data.frame(Worst_Var=var.list, OOB_rmse=OOB.rmse, OOB_mae=OOB.mae, Train_rmse=train.rmse, Train_mae=train.mae)

write.csv(RFE_info, "cf_RFE_info_thresh90_ere_cis.csv", row.names = FALSE)

#set mtry = # of predictors left in model / 3 (nump/3) - and round it up to solid number 
# Refit final model
ere_cis_final <- ERE_cis %>% dplyr::select(Effect, "Sterols","beta-Sitosterol", Industrial_PPCP_Narcotics)

set.seed(3)
cf.final <- cforest(Effect~., data=ere_cis_final, control=cforest_unbiased(mtry=1, ntree=5000))

# OOB predictions
ere_cis_final.pred.oob <- data.frame(predict(cf.final, OOB=T)) %>% rename(Predicted_OOB=Effect)
ere_cis_final_wPred <- cbind(ere_cis_final, ere_cis_final.pred.oob)

# OOB error
mean(abs(ere_cis_final$Effect-ere_cis_final.pred.oob$Predicted_OOB)) #0.295328
sqrt(mean((ere_cis_final$Effect-ere_cis_final.pred.oob$Predicted_OOB)^2)) #0.3919638

# Predicted (OOB) vs. observed
plot(Predicted_OOB~Effect, data=ere_cis_final_wPred, pch=16)
abline(a=0, b=1)

#PPARg####
names(PPARg)
PPARg_sub <- PPARg[,c(2:222)]
PPARg_sub[,1] <- round(PPARg_sub[,1],2)

p <- ncol(PPARg_sub) - 1

#cforest - all variables
set.seed(3)
cf.all <- cforest(Effect~., data =PPARg_sub, control = cforest_unbiased(mtry =ceiling(p/3), ntree = 5000))
cf.pred <- predict(cf.all, OOB = T)
sqrt(mean((PPARg_sub$Effect-cf.pred)^2)) #0.5857959
mean(abs(PPARg_sub$Effect - cf.pred)) #0.3668146

#plot predicted vs. observed
plot(y=cf.pred, x= PPARg_sub$Effect)
abline(a=0, b=1)

# Unconditional variable importance (threshold=1)
set.seed(3) 
cf.upi <- permimp(cf.all, conditional=TRUE, threshold=1, thresholdDiagnostics = TRUE)
as.matrix(sort(cf.upi$values, decreasing=T)) # CF, unconditional

# Conditional variable importance (threshold=0.9)
set.seed(3) 
cf.cpi <- permimp(cf.all, conditional=TRUE, threshold=0.90, thresholdDiagnostics = TRUE)
as.matrix(sort(cf.cpi$values, decreasing=T)) # CF, conditional

# cf.cpi$perTree # VI values per tree
plot(cf.cpi, type="bar")

### Loop to recursively eliminate least informative variables ##

# Pull out response and predictors - columns 5, 7-53
PPARg_preds <- PPARg[, c(2:222)] # original predictors
PPARg_preds[,1] <- round(PPARg_preds[,1], 2)
PPARg_preds_loop <- PPARg_preds 
View(PPARg_preds_loop)

p <- ncol(PPARg_preds)-1

# Save OOB error, variable importance
VI.list <- NULL
var.list <- rep(NA,p)

OOB.rmse <- rep(NA,p)
OOB.mae <- rep(NA,p)

train.rmse <- rep(NA,p)
train.mae <- rep(NA,p)

#run RF - conditional
set.seed(3) # Set seed for reproducibility

for(i in 1:p){
  print(i)
  
  nump <- ncol(PPARg_preds_loop)-1 # Number predictors
  
  set.seed(3) 
  cf <- cforest(Effect~., data=PPARg_preds_loop, control=cforest_unbiased(mtry=ceiling(nump/3), ntree=5000))
  cf.pred <- predict(cf, OOB=T)
  cf.pred.tr <- predict(cf, OOB=F)
  
  # OOB error
  OOB.rmse[i] <- sqrt(mean((PPARg_preds$Effect-cf.pred)^2)) # OOB rmse
  OOB.mae[i] <- mean(abs(PPARg_preds$Effect-cf.pred)) # oob mae
  
  # Training error
  train.rmse[i] <- sqrt(mean((PPARg_preds$Effect-cf.pred.tr)^2)) # training rmse
  train.mae[i] <- mean(abs(PPARg_preds$Effect-cf.pred.tr)) #  training mae
  
  
  # Variable importance
  set.seed(3) 
  cf.cpi <- permimp(cf, conditional=TRUE, threshold=0.9) # set threshold to 1 for unconditional
  
  VI.list[[i]] <- as.matrix(sort(cf.cpi$values, decreasing=F))
  var.list[i] <- row.names(VI.list[[i]])[1] # Least informative variable to remove
  
  
  PPARg_preds_loop <- PPARg_preds_loop[,!colnames(PPARg_preds_loop)==var.list[i]] # remove variable
  
}
beep()

# Plot OOB error - minimum is the best subset of predictors
plot(OOB.rmse, xlab="Iteration", type='l')
points(OOB.rmse, pch=16)
plot(OOB.mae, xlab="Iteration", type='l')
points(OOB.mae, pch=16)

# training error (not to be used for optimizing number predictors)
plot(train.rmse, xlab="Iteration", type='l')
points(train.rmse, pch=16)

which(OOB.rmse==min(OOB.rmse)) #219
which(OOB.mae==min(OOB.mae)) #219

var.list[which(OOB.rmse==min(OOB.rmse)):length(OOB.rmse)] # switched to rmse because oob.mae min = 219
# 
# [1] "Cotinine"      "Anthraquinone"

#double check with rounded error#
OOB.rmse_1 <- round(OOB.rmse, 2)
OOB.mae_1 <- round(OOB.mae, 2)
which(OOB.rmse_1==min(OOB.rmse_1)) #216
which(OOB.mae_1==min(OOB.mae_1)) #216

var.list[which(OOB.mae_1==min(OOB.mae_1)):length(OOB.mae_1)] #same

# Save output - name files in a way that lets you know what differs between different runs of the algorithm

saveRDS(VI.list, "cf_VI_at_iter_thresh90_pparg.rds") # R object (list) with the variable importances at each iteration of the algorithm

RFE_info <- data.frame(Worst_Var=var.list, OOB_rmse=OOB.rmse, OOB_mae=OOB.mae, Train_rmse=train.rmse, Train_mae=train.mae)

write.csv(RFE_info, "cf_RFE_info_thresh90_pparg.csv", row.names = FALSE)

#set mtry = # of predictors left in model / 3 (nump/3) - and round it up to solid number 
# Refit final model
pparg_final <- PPARg %>% dplyr::select(Effect, Anthraquinone, Cotinine)
                                       
set.seed(3)
cf.final <- cforest(Effect~., data=pparg_final, control=cforest_unbiased(mtry=1, ntree=5000))

# OOB predictions
pparg_final.pred.oob <- data.frame(predict(cf.final, OOB=T)) %>% rename(Predicted_OOB=Effect)
pparg_final_wPred <- cbind(pparg_final, pparg_final.pred.oob)

# OOB error
mean(abs(pparg_final$Effect-pparg_final.pred.oob$Predicted_OOB)) #0.3138849
sqrt(mean((pparg_final$Effect-pparg_final.pred.oob$Predicted_OOB)^2)) #0.6535772

# Predicted (OOB) vs. observed
plot(Predicted_OOB~Effect, data=pparg_final_wPred, pch=16)
abline(a=0, b=1)

#PPARa####

names(PPARa)
PPARa_sub <- PPARa[,c(2:222)]
PPARa_sub[,1] <- round(PPARa_sub[,1],2)
View(PPARa_sub)

p <- ncol(PPARa_sub) - 1

#cforest - all variables
set.seed(3)
cf.all <- cforest(Effect~., data =PPARa_sub, control = cforest_unbiased(mtry =ceiling(p/3), ntree = 5000))
cf.pred <- predict(cf.all, OOB = T)
sqrt(mean((PPARa_sub$Effect-cf.pred)^2)) #0.5977465
mean(abs(PPARa_sub$Effect - cf.pred)) #0.3174237

#plot predicted vs. observed
plot(y=cf.pred, x= PPARa_sub$Effect)
abline(a=0, b=1)

# Unconditional variable importance (threshold=1)
set.seed(3) 
cf.upi <- permimp(cf.all, conditional=TRUE, threshold=1, thresholdDiagnostics = TRUE)
as.matrix(sort(cf.upi$values, decreasing=T)) # CF, unconditional

# Conditional variable importance (threshold=0.9)
set.seed(3) 
cf.cpi <- permimp(cf.all, conditional=TRUE, threshold=0.90, thresholdDiagnostics = TRUE)
as.matrix(sort(cf.cpi$values, decreasing=T)) # CF, conditional

# cf.cpi$perTree # VI values per tree
plot(cf.cpi, type="bar")

### Loop to recursively eliminate least informative variables ##

# Pull out response and predictors - columns 5, 7-53
PPARa_preds <- PPARa[, c(2:222)] # original predictors
PPARa_preds[,1] <- round(PPARa_preds[,1], 2)
PPARa_preds_loop <- PPARa_preds 
View(PPARa_preds_loop)

p <- ncol(PPARa_preds)-1

# Save OOB error, variable importance
VI.list <- NULL
var.list <- rep(NA,p)

OOB.rmse <- rep(NA,p)
OOB.mae <- rep(NA,p)

train.rmse <- rep(NA,p)
train.mae <- rep(NA,p)

#run RF - conditional
set.seed(3) # Set seed for reproducibility

for(i in 1:p){
  print(i)
  
  nump <- ncol(PPARa_preds_loop)-1 # Number predictors
  
  set.seed(3) 
  cf <- cforest(Effect~., data=PPARa_preds_loop, control=cforest_unbiased(mtry=ceiling(nump/3), ntree=5000))
  cf.pred <- predict(cf, OOB=T)
  cf.pred.tr <- predict(cf, OOB=F)
  
  # OOB error
  OOB.rmse[i] <- sqrt(mean((PPARa_preds$Effect-cf.pred)^2)) # OOB rmse
  OOB.mae[i] <- mean(abs(PPARa_preds$Effect-cf.pred)) # oob mae
  
  # Training error
  train.rmse[i] <- sqrt(mean((PPARa_preds$Effect-cf.pred.tr)^2)) # training rmse
  train.mae[i] <- mean(abs(PPARa_preds$Effect-cf.pred.tr)) #  training mae
  
  
  # Variable importance
  set.seed(3) 
  cf.cpi <- permimp(cf, conditional=TRUE, threshold=0.9) # set threshold to 1 for unconditional
  
  VI.list[[i]] <- as.matrix(sort(cf.cpi$values, decreasing=F))
  var.list[i] <- row.names(VI.list[[i]])[1] # Least informative variable to remove
  
  
  PPARa_preds_loop <- PPARa_preds_loop[,!colnames(PPARa_preds_loop)==var.list[i]] # remove variable
  
}
beep()

# Plot OOB error - minimum is the best subset of predictors
plot(OOB.rmse, xlab="Iteration", type='l')
points(OOB.rmse, pch=16)
plot(OOB.mae, xlab="Iteration", type='l')
points(OOB.mae, pch=16)

# training error (not to be used for optimizing number predictors)
plot(train.rmse, xlab="Iteration", type='l')
points(train.rmse, pch=16)

which(OOB.rmse==min(OOB.rmse)) #204
which(OOB.mae==min(OOB.mae)) #206

var.list[which(OOB.rmse==min(OOB.rmse)):length(OOB.rmse)] # switched to rmse because oob.mae min = 219

#double check with rounded error#
OOB.rmse_1 <- round(OOB.rmse, 2)
OOB.mae_1 <- round(OOB.mae, 2)
which(OOB.rmse_1==min(OOB.rmse_1)) #220
which(OOB.mae_1==min(OOB.mae_1)) #218

var.list[which(OOB.mae_1==min(OOB.mae_1)):length(OOB.mae_1)] #same

# Save output - name files in a way that lets you know what differs between different runs of the algorithm

saveRDS(VI.list, "cf_VI_at_iter_thresh90_ppara.rds") # R object (list) with the variable importances at each iteration of the algorithm

RFE_info <- data.frame(Worst_Var=var.list, OOB_rmse=OOB.rmse, OOB_mae=OOB.mae, Train_rmse=train.rmse, Train_mae=train.mae)

write.csv(RFE_info, "cf_RFE_info_thresh90_ppara.csv", row.names = FALSE)

#set mtry = # of predictors left in model / 3 (nump/3) - and round it up to solid number 
# Refit final model
ppara_final <- PPARa %>% dplyr::select(Effect, Camphor, Anthraquinone, "Industrial_PPCP_Narcotics")

set.seed(3)
cf.final <- cforest(Effect~., data=ppara_final, control=cforest_unbiased(mtry=1, ntree=5000))

# OOB predictions
ppara_final.pred.oob <- data.frame(predict(cf.final, OOB=T)) %>% rename(Predicted_OOB=Effect)
ppara_final_wPred <- cbind(ppara_final, ppara_final.pred.oob)

# OOB error
mean(abs(ppara_final$Effect-ppara_final.pred.oob$Predicted_OOB)) #0.1003553
sqrt(mean((ppara_final$Effect-ppara_final.pred.oob$Predicted_OOB)^2)) #0.1354775

# Predicted (OOB) vs. observed
plot(Predicted_OOB~Effect, data=ppara_final_wPred, pch=16)
abline(a=0, b=1)

#ARE####
names(ARE)
ARE_sub <- ARE[,c(2:222)]
ARE_sub[,1] <- round(ARE_sub[,1],2)

p <- ncol(ARE_sub) - 1

#cforest - all variables
set.seed(3)
cf.all <- cforest(Effect~., data =ARE_sub, control = cforest_unbiased(mtry =ceiling(p/3), ntree = 5000))
cf.pred <- predict(cf.all, OOB = T)
sqrt(mean((ARE_sub$Effect-cf.pred)^2)) #0.2313008
mean(abs(ARE_sub$Effect - cf.pred)) #0.1881474

#plot predicted vs. observed
plot(y=cf.pred, x= ARE_sub$Effect)
abline(a=0, b=1)

# Unconditional variable importance (threshold=1)
set.seed(3) 
cf.upi <- permimp(cf.all, conditional=TRUE, threshold=1, thresholdDiagnostics = TRUE)
as.matrix(sort(cf.upi$values, decreasing=T)) # CF, unconditional

# Conditional variable importance (threshold=0.9)
set.seed(3) 
cf.cpi <- permimp(cf.all, conditional=TRUE, threshold=0.90, thresholdDiagnostics = TRUE)
as.matrix(sort(cf.cpi$values, decreasing=T)) # CF, conditional

# cf.cpi$perTree # VI values per tree
plot(cf.cpi, type="bar")

### Loop to recursively eliminate least informative variables ##

# Pull out response and predictors - columns 5, 7-53
ARE_preds <- ARE[, c(2:222)] # original predictors
ARE_preds[,1] <- round(ARE_preds[,1], 2)
ARE_preds_loop <- ARE_preds 
View(ARE_preds_loop)

p <- ncol(ARE_preds)-1

# Save OOB error, variable importance
VI.list <- NULL
var.list <- rep(NA,p)

OOB.rmse <- rep(NA,p)
OOB.mae <- rep(NA,p)

train.rmse <- rep(NA,p)
train.mae <- rep(NA,p)

#run RF - conditional
set.seed(3) # Set seed for reproducibility

for(i in 1:p){
  print(i)
  
  nump <- ncol(ARE_preds_loop)-1 # Number predictors
  
  set.seed(3) 
  cf <- cforest(Effect~., data=ARE_preds_loop, control=cforest_unbiased(mtry=ceiling(nump/3), ntree=5000))
  cf.pred <- predict(cf, OOB=T)
  cf.pred.tr <- predict(cf, OOB=F)
  
  # OOB error
  OOB.rmse[i] <- sqrt(mean((ARE_preds$Effect-cf.pred)^2)) # OOB rmse
  OOB.mae[i] <- mean(abs(ARE_preds$Effect-cf.pred)) # oob mae
  
  # Training error
  train.rmse[i] <- sqrt(mean((ARE_preds$Effect-cf.pred.tr)^2)) # training rmse
  train.mae[i] <- mean(abs(ARE_preds$Effect-cf.pred.tr)) #  training mae
  
  
  # Variable importance
  set.seed(3) 
  cf.cpi <- permimp(cf, conditional=TRUE, threshold=0.9) # set threshold to 1 for unconditional
  
  VI.list[[i]] <- as.matrix(sort(cf.cpi$values, decreasing=F))
  var.list[i] <- row.names(VI.list[[i]])[1] # Least informative variable to remove
  
  
  ARE_preds_loop <- ARE_preds_loop[,!colnames(ARE_preds_loop)==var.list[i]] # remove variable
  
}
beep()

# Plot OOB error - minimum is the best subset of predictors
plot(OOB.rmse, xlab="Iteration", type='l')
points(OOB.rmse, pch=16)
plot(OOB.mae, xlab="Iteration", type='l')
points(OOB.mae, pch=16)

# training error (not to be used for optimizing number predictors)
plot(train.rmse, xlab="Iteration", type='l')
points(train.rmse, pch=16)

which(OOB.rmse==min(OOB.rmse)) #219
which(OOB.mae==min(OOB.mae)) #220

var.list[which(OOB.mae==min(OOB.mae)):length(OOB.mae)] # switched to rmse because oob.mae min = 219

#double check with rounded error#
OOB.rmse_1 <- round(OOB.rmse, 2)
OOB.mae_1 <- round(OOB.mae, 2)
which(OOB.rmse_1==min(OOB.rmse_1)) #221
which(OOB.mae_1==min(OOB.mae_1)) #221

var.list[which(OOB.mae_1==min(OOB.mae_1)):length(OOB.mae_1)] #same

# Save output - name files in a way that lets you know what differs between different runs of the algorithm

saveRDS(VI.list, "cf_VI_at_iter_thresh90_are.rds") # R object (list) with the variable importances at each iteration of the algorithm

RFE_info <- data.frame(Worst_Var=var.list, OOB_rmse=OOB.rmse, OOB_mae=OOB.mae, Train_rmse=train.rmse, Train_mae=train.mae)

write.csv(RFE_info, "cf_RFE_info_thresh90_are.csv", row.names = FALSE)

#set mtry = # of predictors left in model / 3 (nump/3) - and round it up to solid number 
# Refit final model
are_final <- ARE %>% dplyr::select(Effect,
                                   Caffeine)


set.seed(3)
cf.final <- cforest(Effect~., data=are_final, control=cforest_unbiased(mtry=1, ntree=5000))

# OOB predictions
are_final.pred.oob <- data.frame(predict(cf.final, OOB=T)) %>% rename(Predicted_OOB=Effect)
are_final_wPred <- cbind(are_final, are_final.pred.oob)

# OOB error
mean(abs(are_final$Effect-are_final.pred.oob$Predicted_OOB)) #0.1635839
sqrt(mean((are_final$Effect-are_final.pred.oob$Predicted_OOB)^2)) #0.2221913

# Predicted (OOB) vs. observed
plot(Predicted_OOB~Effect, data=are_final_wPred, pch=16)
abline(a=0, b=1)


#PPRE####
names(PPRE)
PPRE_sub <- PPRE[,c(2:222)]
PPRE_sub[,1] <- round(PPRE_sub[,1],2)

p <- ncol(PPRE_sub) - 1

#cforest - all variables
set.seed(3)
cf.all <- cforest(Effect~., data =PPRE_sub, control = cforest_unbiased(mtry =ceiling(p/3), ntree = 5000))
cf.pred <- predict(cf.all, OOB = T)
sqrt(mean((PPRE_sub$Effect-cf.pred)^2)) #0.1486993
mean(abs(PPRE_sub$Effect - cf.pred)) #0.1129969

#plot predicted vs. observed
plot(y=cf.pred, x= PPRE_sub$Effect)
abline(a=0, b=1)

# Unconditional variable importance (threshold=1)
set.seed(3) 
cf.upi <- permimp(cf.all, conditional=TRUE, threshold=1, thresholdDiagnostics = TRUE)
as.matrix(sort(cf.upi$values, decreasing=T)) # CF, unconditional

# Conditional variable importance (threshold=0.9)
set.seed(3) 
cf.cpi <- permimp(cf.all, conditional=TRUE, threshold=0.90, thresholdDiagnostics = TRUE)
as.matrix(sort(cf.cpi$values, decreasing=T)) # CF, conditional

# cf.cpi$perTree # VI values per tree
plot(cf.cpi, type="bar")

### Loop to recursively eliminate least informative variables ##

# Pull out response and predictors - columns 5, 7-53
PPRE_preds <- PPRE[, c(2:222)] # original predictors
PPRE_preds[,1] <- round(PPRE_preds[,1], 2)
PPRE_preds_loop <- PPRE_preds 
View(PPRE_preds_loop)

p <- ncol(PPRE_preds)-1

# Save OOB error, variable importance
VI.list <- NULL
var.list <- rep(NA,p)

OOB.rmse <- rep(NA,p)
OOB.mae <- rep(NA,p)

train.rmse <- rep(NA,p)
train.mae <- rep(NA,p)

#run RF - conditional
set.seed(3) # Set seed for reproducibility

for(i in 1:p){
  print(i)
  
  nump <- ncol(PPRE_preds_loop)-1 # Number predictors
  
  set.seed(3) 
  cf <- cforest(Effect~., data=PPRE_preds_loop, control=cforest_unbiased(mtry=ceiling(nump/3), ntree=5000))
  cf.pred <- predict(cf, OOB=T)
  cf.pred.tr <- predict(cf, OOB=F)
  
  # OOB error
  OOB.rmse[i] <- sqrt(mean((PPRE_preds$Effect-cf.pred)^2)) # OOB rmse
  OOB.mae[i] <- mean(abs(PPRE_preds$Effect-cf.pred)) # oob mae
  
  # Training error
  train.rmse[i] <- sqrt(mean((PPRE_preds$Effect-cf.pred.tr)^2)) # training rmse
  train.mae[i] <- mean(abs(PPRE_preds$Effect-cf.pred.tr)) #  training mae
  
  
  # Variable importance
  set.seed(3) 
  cf.cpi <- permimp(cf, conditional=TRUE, threshold=0.9) # set threshold to 1 for unconditional
  
  VI.list[[i]] <- as.matrix(sort(cf.cpi$values, decreasing=F))
  var.list[i] <- row.names(VI.list[[i]])[1] # Least informative variable to remove
  
  
  PPRE_preds_loop <- PPRE_preds_loop[,!colnames(PPRE_preds_loop)==var.list[i]] # remove variable
  
}
beep()

# Plot OOB error - minimum is the best subset of predictors
plot(OOB.rmse, xlab="Iteration", type='l')
points(OOB.rmse, pch=16)
plot(OOB.mae, xlab="Iteration", type='l')
points(OOB.mae, pch=16)

# training error (not to be used for optimizing number predictors)
plot(train.rmse, xlab="Iteration", type='l')
points(train.rmse, pch=16)

which(OOB.rmse==min(OOB.rmse)) #218
which(OOB.mae==min(OOB.mae)) #219

var.list[which(OOB.rmse==min(OOB.rmse)):length(OOB.rmse)] # switched to rmse because oob.mae min = 219

#double check with rounded error#
OOB.rmse_1 <- round(OOB.rmse, 2)
OOB.mae_1 <- round(OOB.mae, 2)
which(OOB.rmse_1==min(OOB.rmse_1)) #219
which(OOB.mae_1==min(OOB.mae_1)) #219

var.list[which(OOB.mae_1==min(OOB.mae_1)):length(OOB.mae_1)] #same

# Save output - name files in a way that lets you know what differs between different runs of the algorithm

saveRDS(VI.list, "cf_VI_at_iter_thresh90_ppre.rds") # R object (list) with the variable importances at each iteration of the algorithm

RFE_info <- data.frame(Worst_Var=var.list, OOB_rmse=OOB.rmse, OOB_mae=OOB.mae, Train_rmse=train.rmse, Train_mae=train.mae)

write.csv(RFE_info, "cf_RFE_info_thresh90_ppre.csv", row.names = FALSE)

#set mtry = # of predictors left in model / 3 (nump/3) - and round it up to solid number 
# Refit final model
ppre_final <- PPRE %>% dplyr::select(Effect,AR_Activation,Cotinine)

set.seed(3)
cf.final <- cforest(Effect~., data=ppre_final, control=cforest_unbiased(mtry =1, ntree=5000))

# OOB predictions
ppre_final.pred.oob <- data.frame(predict(cf.final, OOB=T)) %>% rename(Predicted_OOB=Effect)
ppre_final_wPred <- cbind(ppre_final, ppre_final.pred.oob)

# OOB error
mean(abs(ppre_final$Effect-ppre_final.pred.oob$Predicted_OOB)) #0.1143592
sqrt(mean((ppre_final$Effect-ppre_final.pred.oob$Predicted_OOB)^2)) #0.1437593

# Predicted (OOB) vs. observed
plot(Predicted_OOB~Effect, data=ppre_final_wPred, pch=16)
abline(a=0, b=1)


#ADRB1####
names(ADRB1)
ADRB1_sub <- ADRB1[,c(2:222)]
ADRB1_sub[,1] <- round(ADRB1_sub[,1],2)
View(ADRB1_sub)

p <- ncol(ADRB1_sub) - 1

#cforest - all variables
set.seed(3)
cf.all <- cforest(Effect~., data =ADRB1_sub, control = cforest_unbiased(mtry =ceiling(p/3), ntree = 5000))
cf.pred <- predict(cf.all, OOB = T)
sqrt(mean((ADRB1_sub$Effect-cf.pred)^2)) #0.1450891
mean(abs(ADRB1_sub$Effect - cf.pred)) #0.117667

#plot predicted vs. observed
plot(y=cf.pred, x= ADRB1_sub$Effect)
abline(a=0, b=1)

# Unconditional variable importance (threshold=1)
set.seed(3) 
cf.upi <- permimp(cf.all, conditional=TRUE, threshold=1, thresholdDiagnostics = TRUE)
as.matrix(sort(cf.upi$values, decreasing=T)) # CF, unconditional

# Conditional variable importance (threshold=0.9)
set.seed(3) 
cf.cpi <- permimp(cf.all, conditional=TRUE, threshold=0.90, thresholdDiagnostics = TRUE)
as.matrix(sort(cf.cpi$values, decreasing=T)) # CF, conditional

# cf.cpi$perTree # VI values per tree
plot(cf.cpi, type="bar")

### Loop to recursively eliminate least informative variables ##

# Pull out response and predictors - columns 5, 7-53
ADRB1_preds <- ADRB1[, c(2:223)] # original predictors
ADRB1_preds[,1] <- round(ADRB1_preds[,1], 2)
ADRB1_preds_loop <- ADRB1_preds 
View(ADRB1_preds_loop)

p <- ncol(ADRB1_preds)-1

# Save OOB error, variable importance
VI.list <- NULL
var.list <- rep(NA,p)

OOB.rmse <- rep(NA,p)
OOB.mae <- rep(NA,p)

train.rmse <- rep(NA,p)
train.mae <- rep(NA,p)

#run RF - conditional
set.seed(3) # Set seed for reproducibility

for(i in 1:p){
  print(i)
  
  nump <- ncol(ADRB1_preds_loop)-1 # Number predictors
  
  set.seed(3) 
  cf <- cforest(Effect~., data=ADRB1_preds_loop, control=cforest_unbiased(mtry=ceiling(nump/3), ntree=5000))
  cf.pred <- predict(cf, OOB=T)
  cf.pred.tr <- predict(cf, OOB=F)
  
  # OOB error
  OOB.rmse[i] <- sqrt(mean((ADRB1_preds$Effect-cf.pred)^2)) # OOB rmse
  OOB.mae[i] <- mean(abs(ADRB1_preds$Effect-cf.pred)) # oob mae
  
  # Training error
  train.rmse[i] <- sqrt(mean((ADRB1_preds$Effect-cf.pred.tr)^2)) # training rmse
  train.mae[i] <- mean(abs(ADRB1_preds$Effect-cf.pred.tr)) #  training mae
  
  
  # Variable importance
  set.seed(3) 
  cf.cpi <- permimp(cf, conditional=TRUE, threshold=0.9) # set threshold to 1 for unconditional
  
  VI.list[[i]] <- as.matrix(sort(cf.cpi$values, decreasing=F))
  var.list[i] <- row.names(VI.list[[i]])[1] # Least informative variable to remove
  
  
  ADRB1_preds_loop <- ADRB1_preds_loop[,!colnames(ADRB1_preds_loop)==var.list[i]] # remove variable
  
}
beep()

# Plot OOB error - minimum is the best subset of predictors
plot(OOB.rmse, xlab="Iteration", type='l')
points(OOB.rmse, pch=16)
plot(OOB.mae, xlab="Iteration", type='l')
points(OOB.mae, pch=16)

# training error (not to be used for optimizing number predictors)
plot(train.rmse, xlab="Iteration", type='l')
points(train.rmse, pch=16)

which(OOB.rmse==min(OOB.rmse))#NA
which(OOB.mae==min(OOB.mae)) #NA

#model failure



#GR####
names(GR)
GR_sub <- GR[,c(2:222)]
GR_sub[,1] <- round(GR_sub[,1],2)
View(GR_sub)

p <- ncol(GR_sub) - 1

#cforest - all variables
set.seed(3)
cf.all <- cforest(Effect~., data =GR_sub, control = cforest_unbiased(mtry =ceiling(p/3), ntree = 5000))
cf.pred <- predict(cf.all, OOB = T)
sqrt(mean((GR_sub$Effect-cf.pred)^2)) #0.1738014
mean(abs(GR_sub$Effect - cf.pred)) #0.0991572

#plot predicted vs. observed
plot(y=cf.pred, x= GR_sub$Effect)
abline(a=0, b=1)

# Unconditional variable importance (threshold=1)
set.seed(3) 
cf.upi <- permimp(cf.all, conditional=TRUE, threshold=1, thresholdDiagnostics = TRUE)
as.matrix(sort(cf.upi$values, decreasing=T)) # CF, unconditional

# Conditional variable importance (threshold=0.9)
set.seed(3) 
cf.cpi <- permimp(cf.all, conditional=TRUE, threshold=0.90, thresholdDiagnostics = TRUE)
as.matrix(sort(cf.cpi$values, decreasing=T)) # CF, conditional

# cf.cpi$perTree # VI values per tree
plot(cf.cpi, type="bar")

### Loop to recursively eliminate least informative variables ##

# Pull out response and predictors - columns 5, 7-53
GR_preds <- GR[, c(2:222)] # original predictors
GR_preds[,1] <- round(GR_preds[,1], 2)
GR_preds_loop <- GR_preds 

p <- ncol(GR_preds)-1

# Save OOB error, variable importance
VI.list <- NULL
var.list <- rep(NA,p)

OOB.rmse <- rep(NA,p)
OOB.mae <- rep(NA,p)

train.rmse <- rep(NA,p)
train.mae <- rep(NA,p)

#run RF - conditional
set.seed(3) # Set seed for reproducibility

for(i in 1:p){
  print(i)
  
  nump <- ncol(GR_preds_loop)-1 # Number predictors
  
  set.seed(3) 
  cf <- cforest(Effect~., data=GR_preds_loop, control=cforest_unbiased(mtry=ceiling(nump/3), ntree=5000))
  cf.pred <- predict(cf, OOB=T)
  cf.pred.tr <- predict(cf, OOB=F)
  
  # OOB error
  OOB.rmse[i] <- sqrt(mean((GR_preds$Effect-cf.pred)^2)) # OOB rmse
  OOB.mae[i] <- mean(abs(GR_preds$Effect-cf.pred)) # oob mae
  
  # Training error
  train.rmse[i] <- sqrt(mean((GR_preds$Effect-cf.pred.tr)^2)) # training rmse
  train.mae[i] <- mean(abs(GR_preds$Effect-cf.pred.tr)) #  training mae
  
  
  # Variable importance
  set.seed(3) 
  cf.cpi <- permimp(cf, conditional=TRUE, threshold=0.9) # set threshold to 1 for unconditional
  
  VI.list[[i]] <- as.matrix(sort(cf.cpi$values, decreasing=F))
  var.list[i] <- row.names(VI.list[[i]])[1] # Least informative variable to remove
  
  
  GR_preds_loop <- GR_preds_loop[,!colnames(GR_preds_loop)==var.list[i]] # remove variable
  
}
beep()

# Plot OOB error - minimum is the best subset of predictors
plot(OOB.rmse, xlab="Iteration", type='l')
points(OOB.rmse, pch=16)
plot(OOB.mae, xlab="Iteration", type='l')
points(OOB.mae, pch=16)

# training error (not to be used for optimizing number predictors)
plot(train.rmse, xlab="Iteration", type='l')
points(train.rmse, pch=16)

which(OOB.rmse==min(OOB.rmse)) #218
which(OOB.mae==min(OOB.mae)) #218

var.list[which(OOB.mae==min(OOB.mae)):length(OOB.mae)] # switched to rmse because oob.mae min = 219

#double check with rounded error#
OOB.rmse_1 <- round(OOB.rmse, 2)
OOB.mae_1 <- round(OOB.mae, 2)
which(OOB.rmse_1==min(OOB.rmse_1)) #218
which(OOB.mae_1==min(OOB.mae_1)) #218

var.list[which(OOB.mae_1==min(OOB.mae_1)):length(OOB.mae_1)] #same

# Save output - name files in a way that lets you know what differs between different runs of the algorithm

saveRDS(VI.list, "cf_VI_at_iter_thresh90_gr.rds") # R object (list) with the variable importances at each iteration of the algorithm

RFE_info <- data.frame(Worst_Var=var.list, OOB_rmse=OOB.rmse, OOB_mae=OOB.mae, Train_rmse=train.rmse, Train_mae=train.mae)

write.csv(RFE_info, "cf_RFE_info_thresh90_gr.csv", row.names = FALSE)

#set mtry = # of predictors left in model / 3 (nump/3) - and round it up to solid number 
# Refit final model
gr_final <- GR %>% dplyr::select(Effect,Triamterene, Metformin,
                                 PPCPs_Pesticides_Plasticizers)

set.seed(3)
cf.final <- cforest(Effect~., data=gr_final, control=cforest_unbiased(mtry=1, ntree=5000))

# OOB predictions
gr_final.pred.oob <- data.frame(predict(cf.final, OOB=T)) %>% rename(Predicted_OOB=Effect)
gr_final_wPred <- cbind(gr_final, gr_final.pred.oob)

# OOB error
mean(abs(gr_final$Effect-gr_final.pred.oob$Predicted_OOB)) #0.09439614
sqrt(mean((gr_final$Effect-gr_final.pred.oob$Predicted_OOB)^2)) #0.1648086

# Predicted (OOB) vs. observed
plot(Predicted_OOB~Effect, data=gr_final_wPred, pch=16)
abline(a=0, b=1)

#MC1R####
names(MC1R)
MC1R_sub <- MC1R[,c(2:222)]
MC1R_sub[,1] <- round(MC1R_sub[,1],2)

p <- ncol(MC1R_sub) - 1

#cforest - all variables
set.seed(3)
cf.all <- cforest(Effect~., data =MC1R_sub, control = cforest_unbiased(mtry =ceiling(p/3), ntree = 5000))
cf.pred <- predict(cf.all, OOB = T)
sqrt(mean((MC1R_sub$Effect-cf.pred)^2)) #0.1953119
mean(abs(MC1R_sub$Effect - cf.pred)) #0.1482523

#plot predicted vs. observed
plot(y=cf.pred, x= MC1R_sub$Effect)
abline(a=0, b=1)

# Unconditional variable importance (threshold=1)
set.seed(3) 
cf.upi <- permimp(cf.all, conditional=TRUE, threshold=1, thresholdDiagnostics = TRUE)
as.matrix(sort(cf.upi$values, decreasing=T)) # CF, unconditional

# Conditional variable importance (threshold=0.9)
set.seed(3) 
cf.cpi <- permimp(cf.all, conditional=TRUE, threshold=0.90, thresholdDiagnostics = TRUE)
as.matrix(sort(cf.cpi$values, decreasing=T)) # CF, conditional

# cf.cpi$perTree # VI values per tree
plot(cf.cpi, type="bar")

### Loop to recursively eliminate least informative variables ##

# Pull out response and predictors - columns 5, 7-53
MC1R_preds <- MC1R[, c(2:223)] # original predictors
MC1R_preds[,1] <- round(MC1R_preds[,1], 2)
MC1R_preds_loop <- MC1R_preds 
View(MC1R_preds_loop)

p <- ncol(MC1R_preds)-1

# Save OOB error, variable importance
VI.list <- NULL
var.list <- rep(NA,p)

OOB.rmse <- rep(NA,p)
OOB.mae <- rep(NA,p)

train.rmse <- rep(NA,p)
train.mae <- rep(NA,p)

#run RF - conditional
set.seed(3) # Set seed for reproducibility

for(i in 1:p){
  print(i)
  
  nump <- ncol(MC1R_preds_loop)-1 # Number predictors
  
  set.seed(3) 
  cf <- cforest(Effect~., data=MC1R_preds_loop, control=cforest_unbiased(mtry=ceiling(nump/3), ntree=5000))
  cf.pred <- predict(cf, OOB=T)
  cf.pred.tr <- predict(cf, OOB=F)
  
  # OOB error
  OOB.rmse[i] <- sqrt(mean((MC1R_preds$Effect-cf.pred)^2)) # OOB rmse
  OOB.mae[i] <- mean(abs(MC1R_preds$Effect-cf.pred)) # oob mae
  
  # Training error
  train.rmse[i] <- sqrt(mean((MC1R_preds$Effect-cf.pred.tr)^2)) # training rmse
  train.mae[i] <- mean(abs(MC1R_preds$Effect-cf.pred.tr)) #  training mae
  
  
  # Variable importance
  set.seed(3) 
  cf.cpi <- permimp(cf, conditional=TRUE, threshold=0.9) # set threshold to 1 for unconditional
  
  VI.list[[i]] <- as.matrix(sort(cf.cpi$values, decreasing=F))
  var.list[i] <- row.names(VI.list[[i]])[1] # Least informative variable to remove
  
  
  MC1R_preds_loop <- MC1R_preds_loop[,!colnames(MC1R_preds_loop)==var.list[i]] # remove variable
  
}
beep()

# Plot OOB error - minimum is the best subset of predictors
plot(OOB.rmse, xlab="Iteration", type='l')
points(OOB.rmse, pch=16)
plot(OOB.mae, xlab="Iteration", type='l')
points(OOB.mae, pch=16)

# training error (not to be used for optimizing number predictors)
plot(train.rmse, xlab="Iteration", type='l')
points(train.rmse, pch=16)

#model failure##


#PTGDR####
names(PTGDR)
PTGDR_sub <- PTGDR[,c(2:222)]
PTGDR_sub[,1] <- round(PTGDR_sub[,1],2)
View(PTGDR_sub)

p <- ncol(PTGDR_sub) - 1

#cforest - all variables
set.seed(3)
cf.all <- cforest(Effect~., data =PTGDR_sub, control = cforest_unbiased(mtry =ceiling(p/3), ntree = 5000))
cf.pred <- predict(cf.all, OOB = T)
sqrt(mean((PTGDR_sub$Effect-cf.pred)^2)) #0.3834848
mean(abs(PTGDR_sub$Effect - cf.pred)) #0.3083478

#plot predicted vs. observed
plot(y=cf.pred, x= PTGDR_sub$Effect)
abline(a=0, b=1)

# Unconditional variable importance (threshold=1)
set.seed(3) 
cf.upi <- permimp(cf.all, conditional=TRUE, threshold=1, thresholdDiagnostics = TRUE)
as.matrix(sort(cf.upi$values, decreasing=T)) # CF, unconditional

# Conditional variable importance (threshold=0.9)
set.seed(3) 
cf.cpi <- permimp(cf.all, conditional=TRUE, threshold=0.90, thresholdDiagnostics = TRUE)
as.matrix(sort(cf.cpi$values, decreasing=T)) # CF, conditional

# cf.cpi$perTree # VI values per tree
plot(cf.cpi, type="bar")

### Loop to recursively eliminate least informative variables ##

# Pull out response and predictors - columns 5, 7-53
PTGDR_preds <- PTGDR[, c(2:223)] # original predictors
PTGDR_preds[,1] <- round(PTGDR_preds[,1], 2)
PTGDR_preds_loop <- PTGDR_preds 

p <- ncol(PTGDR_preds)-1

# Save OOB error, variable importance
VI.list <- NULL
var.list <- rep(NA,p)

OOB.rmse <- rep(NA,p)
OOB.mae <- rep(NA,p)

train.rmse <- rep(NA,p)
train.mae <- rep(NA,p)

#run RF - conditional
set.seed(3) # Set seed for reproducibility

for(i in 1:p){
  print(i)
  
  nump <- ncol(PTGDR_preds_loop)-1 # Number predictors
  
  set.seed(3) 
  cf <- cforest(Effect~., data=PTGDR_preds_loop, control=cforest_unbiased(mtry=ceiling(nump/3), ntree=5000))
  cf.pred <- predict(cf, OOB=T)
  cf.pred.tr <- predict(cf, OOB=F)
  
  # OOB error
  OOB.rmse[i] <- sqrt(mean((PTGDR_preds$Effect-cf.pred)^2)) # OOB rmse
  OOB.mae[i] <- mean(abs(PTGDR_preds$Effect-cf.pred)) # oob mae
  
  # Training error
  train.rmse[i] <- sqrt(mean((PTGDR_preds$Effect-cf.pred.tr)^2)) # training rmse
  train.mae[i] <- mean(abs(PTGDR_preds$Effect-cf.pred.tr)) #  training mae
  
  
  # Variable importance
  set.seed(3) 
  cf.cpi <- permimp(cf, conditional=TRUE, threshold=0.9) # set threshold to 1 for unconditional
  
  VI.list[[i]] <- as.matrix(sort(cf.cpi$values, decreasing=F))
  var.list[i] <- row.names(VI.list[[i]])[1] # Least informative variable to remove
  
  
  PTGDR_preds_loop <- PTGDR_preds_loop[,!colnames(PTGDR_preds_loop)==var.list[i]] # remove variable
  
}
beep()

# Plot OOB error - minimum is the best subset of predictors
plot(OOB.rmse, xlab="Iteration", type='l')
points(OOB.rmse, pch=16)
plot(OOB.mae, xlab="Iteration", type='l')
points(OOB.mae, pch=16)

# training error (not to be used for optimizing number predictors)
plot(train.rmse, xlab="Iteration", type='l')
points(train.rmse, pch=16)

which(OOB.rmse==min(OOB.rmse)) 
which(OOB.mae==min(OOB.mae)) 

#model failure###


#PTGER####
names(PTGER)
PTGER_sub <- PTGER[,c(2:222)]
PTGER_sub[,1] <- round(PTGER_sub[,1],2)

p <- ncol(PTGER_sub) - 1

#cforest - all variables
set.seed(3)
cf.all <- cforest(Effect~., data =PTGER_sub, control = cforest_unbiased(mtry =ceiling(p/3), ntree = 5000))
cf.pred <- predict(cf.all, OOB = T)
sqrt(mean((PTGER_sub$Effect-cf.pred)^2)) #0.4443157
mean(abs(PTGER_sub$Effect - cf.pred)) #0.3113312

#plot predicted vs. observed
plot(y=cf.pred, x= PTGER_sub$Effect)
abline(a=0, b=1)

# Unconditional variable importance (threshold=1)
set.seed(3) 
cf.upi <- permimp(cf.all, conditional=TRUE, threshold=1, thresholdDiagnostics = TRUE)
as.matrix(sort(cf.upi$values, decreasing=T)) # CF, unconditional

# Conditional variable importance (threshold=0.9)
set.seed(3) 
cf.cpi <- permimp(cf.all, conditional=TRUE, threshold=0.90, thresholdDiagnostics = TRUE)
as.matrix(sort(cf.cpi$values, decreasing=T)) # CF, conditional

# cf.cpi$perTree # VI values per tree
plot(cf.cpi, type="bar")

### Loop to recursively eliminate least informative variables ##

# Pull out response and predictors - columns 5, 7-53
PTGER_preds <- PTGER[, c(2:223)] # original predictors
PTGER_preds[,1] <- round(PTGER_preds[,1], 2)
PTGER_preds_loop <- PTGER_preds 
View(PTGER_preds_loop)

p <- ncol(PTGER_preds)-1

# Save OOB error, variable importance
VI.list <- NULL
var.list <- rep(NA,p)

OOB.rmse <- rep(NA,p)
OOB.mae <- rep(NA,p)

train.rmse <- rep(NA,p)
train.mae <- rep(NA,p)

#run RF - conditional
set.seed(3) # Set seed for reproducibility

for(i in 1:p){
  print(i)
  
  nump <- ncol(PTGER_preds_loop)-1 # Number predictors
  
  set.seed(3) 
  cf <- cforest(Effect~., data=PTGER_preds_loop, control=cforest_unbiased(mtry=ceiling(nump/3), ntree=5000))
  cf.pred <- predict(cf, OOB=T)
  cf.pred.tr <- predict(cf, OOB=F)
  
  # OOB error
  OOB.rmse[i] <- sqrt(mean((PTGER_preds$Effect-cf.pred)^2)) # OOB rmse
  OOB.mae[i] <- mean(abs(PTGER_preds$Effect-cf.pred)) # oob mae
  
  # Training error
  train.rmse[i] <- sqrt(mean((PTGER_preds$Effect-cf.pred.tr)^2)) # training rmse
  train.mae[i] <- mean(abs(PTGER_preds$Effect-cf.pred.tr)) #  training mae
  
  
  # Variable importance
  set.seed(3) 
  cf.cpi <- permimp(cf, conditional=TRUE, threshold=0.9) # set threshold to 1 for unconditional
  
  VI.list[[i]] <- as.matrix(sort(cf.cpi$values, decreasing=F))
  var.list[i] <- row.names(VI.list[[i]])[1] # Least informative variable to remove
  
  
  PTGER_preds_loop <- PTGER_preds_loop[,!colnames(PTGER_preds_loop)==var.list[i]] # remove variable
  
}
beep()

# Plot OOB error - minimum is the best subset of predictors
plot(OOB.rmse, xlab="Iteration", type='l')
points(OOB.rmse, pch=16)
plot(OOB.mae, xlab="Iteration", type='l')
points(OOB.mae, pch=16)

# training error (not to be used for optimizing number predictors)
plot(train.rmse, xlab="Iteration", type='l')
points(train.rmse, pch=16)

which(OOB.rmse==min(OOB.rmse)) #219
which(OOB.mae==min(OOB.mae)) #213

var.list[which(OOB.rmse==min(OOB.rmse)):length(OOB.rmse)] # switched to rmse because oob.mae min = 219

#model failure###


#PTGIR####
names(PTGIR)
PTGIR_sub <- PTGIR[,c(2:222)]
PTGIR_sub[,1] <- round(PTGIR_sub[,1],2)
View(PTGIR_sub)

p <- ncol(PTGIR_sub) - 1

#cforest - all variables
set.seed(3)
cf.all <- cforest(Effect~., data =PTGIR_sub, control = cforest_unbiased(mtry =ceiling(p/3), ntree = 5000))
cf.pred <- predict(cf.all, OOB = T)
sqrt(mean((PTGIR_sub$Effect-cf.pred)^2)) #0.2115806
mean(abs(PTGIR_sub$Effect - cf.pred)) #0.1608252

#plot predicted vs. observed
plot(y=cf.pred, x= PTGIR_sub$Effect)
abline(a=0, b=1)

# Unconditional variable importance (threshold=1)
set.seed(3) 
cf.upi <- permimp(cf.all, conditional=TRUE, threshold=1, thresholdDiagnostics = TRUE)
as.matrix(sort(cf.upi$values, decreasing=T)) # CF, unconditional

# Conditional variable importance (threshold=0.9)
set.seed(3) 
cf.cpi <- permimp(cf.all, conditional=TRUE, threshold=0.90, thresholdDiagnostics = TRUE)
as.matrix(sort(cf.cpi$values, decreasing=T)) # CF, conditional

# cf.cpi$perTree # VI values per tree
plot(cf.cpi, type="bar")

### Loop to recursively eliminate least informative variables ##

# Pull out response and predictors - columns 5, 7-53
PTGIR_preds <- PTGIR[, c(2:223)] # original predictors
PTGIR_preds[,1] <- round(PTGIR_preds[,1], 2)
PTGIR_preds_loop <- PTGIR_preds 
View(PTGIR_preds_loop)

p <- ncol(PTGIR_preds)-1

# Save OOB error, variable importance
VI.list <- NULL
var.list <- rep(NA,p)

OOB.rmse <- rep(NA,p)
OOB.mae <- rep(NA,p)

train.rmse <- rep(NA,p)
train.mae <- rep(NA,p)

#run RF - conditional
set.seed(3) # Set seed for reproducibility

for(i in 1:p){
  print(i)
  
  nump <- ncol(PTGIR_preds_loop)-1 # Number predictors
  
  set.seed(3) 
  cf <- cforest(Effect~., data=PTGIR_preds_loop, control=cforest_unbiased(mtry=ceiling(nump/3), ntree=5000))
  cf.pred <- predict(cf, OOB=T)
  cf.pred.tr <- predict(cf, OOB=F)
  
  # OOB error
  OOB.rmse[i] <- sqrt(mean((PTGIR_preds$Effect-cf.pred)^2)) # OOB rmse
  OOB.mae[i] <- mean(abs(PTGIR_preds$Effect-cf.pred)) # oob mae
  
  # Training error
  train.rmse[i] <- sqrt(mean((PTGIR_preds$Effect-cf.pred.tr)^2)) # training rmse
  train.mae[i] <- mean(abs(PTGIR_preds$Effect-cf.pred.tr)) #  training mae
  
  
  # Variable importance
  set.seed(3) 
  cf.cpi <- permimp(cf, conditional=TRUE, threshold=0.9) # set threshold to 1 for unconditional
  
  VI.list[[i]] <- as.matrix(sort(cf.cpi$values, decreasing=F))
  var.list[i] <- row.names(VI.list[[i]])[1] # Least informative variable to remove
  
  
  PTGIR_preds_loop <- PTGIR_preds_loop[,!colnames(PTGIR_preds_loop)==var.list[i]] # remove variable
  
}
beep()

# Plot OOB error - minimum is the best subset of predictors
plot(OOB.rmse, xlab="Iteration", type='l')
points(OOB.rmse, pch=16)
plot(OOB.mae, xlab="Iteration", type='l')
points(OOB.mae, pch=16)

# training error (not to be used for optimizing number predictors)
plot(train.rmse, xlab="Iteration", type='l')
points(train.rmse, pch=16)

which(OOB.rmse==min(OOB.rmse)) #219
which(OOB.mae==min(OOB.mae)) #213

var.list[which(OOB.rmse==min(OOB.rmse)):length(OOB.rmse)] # switched to rmse because oob.mae min = 219

#model failure


#RXRb####
names(RXRb)
RXRb_sub <- RXRb[,c(2:222)]
RXRb_sub[,1] <- round(RXRb_sub[,1],2)
View(RXRb_sub)

p <- ncol(RXRb_sub) - 1

#cforest - all variables
set.seed(3)
cf.all <- cforest(Effect~., data =RXRb_sub, control = cforest_unbiased(mtry =ceiling(p/3), ntree = 5000))
cf.pred <- predict(cf.all, OOB = T)
sqrt(mean((RXRb_sub$Effect-cf.pred)^2)) #0.2049511
mean(abs(RXRb_sub$Effect - cf.pred)) #0.1501685

#plot predicted vs. observed
plot(y=cf.pred, x= RXRb_sub$Effect)
abline(a=0, b=1)

# Unconditional variable importance (threshold=1)
set.seed(3) 
cf.upi <- permimp(cf.all, conditional=TRUE, threshold=1, thresholdDiagnostics = TRUE)
as.matrix(sort(cf.upi$values, decreasing=T)) # CF, unconditional

# Conditional variable importance (threshold=0.9)
set.seed(3) 
cf.cpi <- permimp(cf.all, conditional=TRUE, threshold=0.90, thresholdDiagnostics = TRUE)
as.matrix(sort(cf.cpi$values, decreasing=T)) # CF, conditional

# cf.cpi$perTree # VI values per tree
plot(cf.cpi, type="bar")

### Loop to recursively eliminate least informative variables ##

# Pull out response and predictors - columns 5, 7-53
RXRb_preds <- RXRb[, c(2:222)] # original predictors
RXRb_preds[,1] <- round(RXRb_preds[,1], 2)
RXRb_preds_loop <- RXRb_preds 
View(RXRb_preds_loop)

p <- ncol(RXRb_preds)-1

# Save OOB error, variable importance
VI.list <- NULL
var.list <- rep(NA,p)

OOB.rmse <- rep(NA,p)
OOB.mae <- rep(NA,p)

train.rmse <- rep(NA,p)
train.mae <- rep(NA,p)

#run RF - conditional
set.seed(3) # Set seed for reproducibility

for(i in 1:p){
  print(i)
  
  nump <- ncol(RXRb_preds_loop)-1 # Number predictors
  
  set.seed(3) 
  cf <- cforest(Effect~., data=RXRb_preds_loop, control=cforest_unbiased(mtry=ceiling(nump/3), ntree=5000))
  cf.pred <- predict(cf, OOB=T)
  cf.pred.tr <- predict(cf, OOB=F)
  
  # OOB error
  OOB.rmse[i] <- sqrt(mean((RXRb_preds$Effect-cf.pred)^2)) # OOB rmse
  OOB.mae[i] <- mean(abs(RXRb_preds$Effect-cf.pred)) # oob mae
  
  # Training error
  train.rmse[i] <- sqrt(mean((RXRb_preds$Effect-cf.pred.tr)^2)) # training rmse
  train.mae[i] <- mean(abs(RXRb_preds$Effect-cf.pred.tr)) #  training mae
  
  
  # Variable importance
  set.seed(3) 
  cf.cpi <- permimp(cf, conditional=TRUE, threshold=0.9) # set threshold to 1 for unconditional
  
  VI.list[[i]] <- as.matrix(sort(cf.cpi$values, decreasing=F))
  var.list[i] <- row.names(VI.list[[i]])[1] # Least informative variable to remove
  
  
  RXRb_preds_loop <- RXRb_preds_loop[,!colnames(RXRb_preds_loop)==var.list[i]] # remove variable
  
}
beep()

# Plot OOB error - minimum is the best subset of predictors
plot(OOB.rmse, xlab="Iteration", type='l')
points(OOB.rmse, pch=16)
plot(OOB.mae, xlab="Iteration", type='l')
points(OOB.mae, pch=16)

# training error (not to be used for optimizing number predictors)
plot(train.rmse, xlab="Iteration", type='l')
points(train.rmse, pch=16)

which(OOB.rmse==min(OOB.rmse)) #221
which(OOB.mae==min(OOB.mae)) #221

var.list[which(OOB.rmse==min(OOB.rmse)):length(OOB.rmse)] # switched to rmse because oob.mae min = 219

#double check with rounded error#
OOB.rmse_1 <- round(OOB.rmse, 2)
OOB.mae_1 <- round(OOB.mae, 2)
which(OOB.rmse_1==min(OOB.rmse_1)) #220
which(OOB.mae_1==min(OOB.mae_1)) #220

var.list[which(OOB.mae_1==min(OOB.mae_1)):length(OOB.mae_1)] #same

# Save output - name files in a way that lets you know what differs between different runs of the algorithm

saveRDS(VI.list, "cf_VI_at_iter_thresh90_rxrb.rds") # R object (list) with the variable importances at each iteration of the algorithm

RFE_info <- data.frame(Worst_Var=var.list, OOB_rmse=OOB.rmse, OOB_mae=OOB.mae, Train_rmse=train.rmse, Train_mae=train.mae)

write.csv(RFE_info, "cf_RFE_info_thresh90_rxrb.csv", row.names = FALSE)

#set mtry = # of predictors left in model / 3 (nump/3) - and round it up to solid number 
# Refit final model
rxrb_final <- RXRb %>% dplyr::select(Effect, Acetaminophen, Anthraquinone, RXRB_Activation, Industrial_PPCP_Narcotics)

set.seed(3)
cf.final <- cforest(Effect~., data=rxrb_final, control=cforest_unbiased(mtry=1, ntree=5000))

# OOB predictions
rxrb_final.pred.oob <- data.frame(predict(cf.final, OOB=T)) %>% rename(Predicted_OOB=Effect)
rxrb_final_wPred <- cbind(rxrb_final, rxrb_final.pred.oob)

# OOB error
mean(abs(rxrb_final$Effect-rxrb_final.pred.oob$Predicted_OOB)) #0.1374372
sqrt(mean((rxrb_final$Effect-rxrb_final.pred.oob$Predicted_OOB)^2)) #0.1904174

# Predicted (OOB) vs. observed
plot(Predicted_OOB~Effect, data=rxrb_final_wPred, pch=16)
abline(a=0, b=1)

#t47D####
names(t47)
t47_sub <- t47[,c(2:222)]
t47_sub[,1] <- round(t47_sub[,1],2)
View(t47_sub)

p <- ncol(t47_sub) - 1

#cforest - all variables
set.seed(3)
cf.all <- cforest(Effect~., data =t47_sub, control = cforest_unbiased(mtry =ceiling(p/3), ntree = 5000))
cf.pred <- predict(cf.all, OOB = T)
sqrt(mean((t47_sub$Effect-cf.pred)^2)) #0.165846
mean(abs(t47_sub$Effect - cf.pred)) #0.1235471

#plot predicted vs. observed
plot(y=cf.pred, x= t47_sub$Effect)
abline(a=0, b=1)

# Unconditional variable importance (threshold=1)
set.seed(3) 
cf.upi <- permimp(cf.all, conditional=TRUE, threshold=1, thresholdDiagnostics = TRUE)
as.matrix(sort(cf.upi$values, decreasing=T)) # CF, unconditional

# Conditional variable importance (threshold=0.9)
set.seed(3) 
cf.cpi <- permimp(cf.all, conditional=TRUE, threshold=0.90, thresholdDiagnostics = TRUE)
as.matrix(sort(cf.cpi$values, decreasing=T)) # CF, conditional

# cf.cpi$perTree # VI values per tree
plot(cf.cpi, type="bar")

### Loop to recursively eliminate least informative variables ##

# Pull out response and predictors - columns 5, 7-53
t47_preds <- t47[, c(2:222)] # original predictors
t47_preds[,1] <- round(t47_preds[,1], 2)
t47_preds_loop <- t47_preds 
View(t47_preds_loop)

p <- ncol(t47_preds)-1

# Save OOB error, variable importance
VI.list <- NULL
var.list <- rep(NA,p)

OOB.rmse <- rep(NA,p)
OOB.mae <- rep(NA,p)

train.rmse <- rep(NA,p)
train.mae <- rep(NA,p)

#run RF - conditional
set.seed(3) # Set seed for reproducibility

for(i in 1:p){
  print(i)
  
  nump <- ncol(t47_preds_loop)-1 # Number predictors
  
  set.seed(3) 
  cf <- cforest(Effect~., data=t47_preds_loop, control=cforest_unbiased(mtry=ceiling(nump/3), ntree=5000))
  cf.pred <- predict(cf, OOB=T)
  cf.pred.tr <- predict(cf, OOB=F)
  
  # OOB error
  OOB.rmse[i] <- sqrt(mean((t47_preds$Effect-cf.pred)^2)) # OOB rmse
  OOB.mae[i] <- mean(abs(t47_preds$Effect-cf.pred)) # oob mae
  
  # Training error
  train.rmse[i] <- sqrt(mean((t47_preds$Effect-cf.pred.tr)^2)) # training rmse
  train.mae[i] <- mean(abs(t47_preds$Effect-cf.pred.tr)) #  training mae
  
  
  # Variable importance
  set.seed(3) 
  cf.cpi <- permimp(cf, conditional=TRUE, threshold=0.9) # set threshold to 1 for unconditional
  
  VI.list[[i]] <- as.matrix(sort(cf.cpi$values, decreasing=F))
  var.list[i] <- row.names(VI.list[[i]])[1] # Least informative variable to remove
  
  
  t47_preds_loop <- t47_preds_loop[,!colnames(t47_preds_loop)==var.list[i]] # remove variable
  
}
beep()

# Plot OOB error - minimum is the best subset of predictors
plot(OOB.rmse, xlab="Iteration", type='l')
points(OOB.rmse, pch=16)
plot(OOB.mae, xlab="Iteration", type='l')
points(OOB.mae, pch=16)

# training error (not to be used for optimizing number predictors)
plot(train.rmse, xlab="Iteration", type='l')
points(train.rmse, pch=16)

which(OOB.rmse==min(OOB.rmse)) #217
which(OOB.mae==min(OOB.mae)) #215

var.list[which(OOB.mae==min(OOB.mae)):length(OOB.mae)]

# [1] "Cholesterol"               "Dextromethorphan"          "Methyl-1H-benzotriazole"   "Industrial_PPCP_Narcotics"
# [5] "beta-Sitosterol"           "Sterols"   

#double check with rounded error#
OOB.rmse_1 <- round(OOB.rmse, 2)
OOB.mae_1 <- round(OOB.mae, 2)
which(OOB.rmse_1==min(OOB.rmse_1)) #218
which(OOB.mae_1==min(OOB.mae_1)) #213

var.list[which(OOB.mae_1==min(OOB.mae_1)):length(OOB.mae_1)] #same

# Save output - name files in a way that lets you know what differs between different runs of the algorithm

saveRDS(VI.list, "cf_VI_at_iter_thresh90_t47.rds") # R object (list) with the variable importances at each iteration of the algorithm

RFE_info <- data.frame(Worst_Var=var.list, OOB_rmse=OOB.rmse, OOB_mae=OOB.mae, Train_rmse=train.rmse, Train_mae=train.mae)

write.csv(RFE_info, "cf_RFE_info_thresh90_t47.csv", row.names = FALSE)

#set mtry = # of predictors left in model / 3 (nump/3) - and round it up to solid number 
# Refit final model
t47_final <- t47 %>% dplyr::select(Effect, "Methyl-1H-benzotriazole", Cholesterol, 
                                   Sterols, Dextromethorphan, "beta-Sitosterol", Industrial_PPCP_Narcotics)
set.seed(3)
cf.final <- cforest(Effect~., data=t47_final, control=cforest_unbiased(mtry=2, ntree=5000))

# OOB predictions
t47_final.pred.oob <- data.frame(predict(cf.final, OOB=T)) %>% rename(Predicted_OOB=Effect)
t47_final_wPred <- cbind(t47_final, t47_final.pred.oob)

# OOB error
mean(abs(t47_final$Effect-t47_final.pred.oob$Predicted_OOB)) #0.1219773
sqrt(mean((t47_final$Effect-t47_final.pred.oob$Predicted_OOB)^2)) #0.1624588

# Predicted (OOB) vs. observed
plot(Predicted_OOB~Effect, data=t47_final_wPred, pch=16)
abline(a=0, b=1)

#RIA####
#ria M### code needs to modified if it is re-run
names(ria)
ria_sub <- ria[,c(2:111)]
ria_sub[,1] <- round(ria_sub[,1],2)

p <- ncol(ria_sub) - 1

#cforest - all variables
set.seed(3)
cf.all <- cforest(Effect~., data =ria_sub, control = cforest_unbiased(mtry= ceiling(p/3), ntree = 5000))
cf.pred <- predict(cf.all, OOB = T)
sqrt(mean((ria_sub$Effect-cf.pred)^2)) #1.534635  
mean(abs(ria_sub$Effect - cf.pred)) #0.4158677

#plot predicted vs. observed
plot(y=cf.pred, x= ria_sub$Effect)
abline(a=0, b=1)

# Unconditional variable importance (threshold=1)
set.seed(3) 
cf.upi <- permimp(cf.all, conditional=TRUE, threshold=1, thresholdDiagnostics = TRUE)
as.matrix(sort(cf.upi$values, decreasing=T)) # CF, unconditional

# Conditional variable importance (threshold=0.9)
set.seed(3) 
cf.cpi <- permimp(cf.all, conditional=TRUE, threshold=0.90, thresholdDiagnostics = TRUE)
as.matrix(sort(cf.cpi$values, decreasing=T)) # CF, conditional

# cf.cpi$perTree # VI values per tree
plot(cf.cpi, type="bar")

### Loop to recursively eliminate least informative variables ##

# Pull out response and predictors - columns 5, 7-53
ria_preds <- ria[, c(2:111)] # original predictors
ria_preds[,1] <- round(ria_preds[,1], 2)
ria_preds_loop <- ria_preds 
View(ria_preds_loop)

p <- ncol(ria_preds)-1

# Save OOB error, variable importance
VI.list <- NULL
var.list <- rep(NA,p)

OOB.rmse <- rep(NA,p)
OOB.mae <- rep(NA,p)

train.rmse <- rep(NA,p)
train.mae <- rep(NA,p)

#run RF - conditional
set.seed(3) # Set seed for reproducibility

for(i in 1:p){
  print(i)
  
  nump <- ncol(ria_preds_loop)-1 # Number predictors
  
  set.seed(3) 
  cf <- cforest(Effect~., data=ria_preds_loop, control=cforest_unbiased(mtry=ceiling(nump/3), ntree=5000))
  cf.pred <- predict(cf, OOB=T)
  cf.pred.tr <- predict(cf, OOB=F)
  
  # OOB error
  OOB.rmse[i] <- sqrt(mean((ria_preds$Effect-cf.pred)^2)) # OOB rmse
  OOB.mae[i] <- mean(abs(ria_preds$Effect-cf.pred)) # oob mae
  
  # Training error
  train.rmse[i] <- sqrt(mean((ria_preds$Effect-cf.pred.tr)^2)) # training rmse
  train.mae[i] <- mean(abs(ria_preds$Effect-cf.pred.tr)) #  training mae
  
  
  # Variable importance
  set.seed(3) 
  cf.cpi <- permimp(cf, conditional=TRUE, threshold=0.9) # set threshold to 1 for unconditional
  
  VI.list[[i]] <- as.matrix(sort(cf.cpi$values, decreasing=F))
  var.list[i] <- row.names(VI.list[[i]])[1] # Least informative variable to remove
  
  
  ria_preds_loop <- ria_preds_loop[,!colnames(ria_preds_loop)==var.list[i]] # remove variable
  
}
beep()

# Plot OOB error - minimum is the best subset of predictors
plot(OOB.rmse, xlab="Iteration", type='l')
points(OOB.rmse, pch=16)
plot(OOB.mae, xlab="Iteration", type='l')
points(OOB.mae, pch=16)

# training error (not to be used for optimizing number predictors)
plot(train.rmse, xlab="Iteration", type='l')
points(train.rmse, pch=16)

which(OOB.rmse==min(OOB.rmse)) #99
which(OOB.mae==min(OOB.mae)) #106

var.list[which(OOB.mae==min(OOB.mae)):length(OOB.mae)]

#double check with rounded error#
OOB.rmse_1 <- round(OOB.rmse, 2)
OOB.mae_1 <- round(OOB.mae, 2)
which(OOB.rmse_1==min(OOB.rmse_1)) #106
which(OOB.mae_1==min(OOB.mae_1)) #106

# Save output - name files in a way that lets you know what differs between different runs of the algorithm

saveRDS(VI.list, "cf_VI_at_iter_thresh90_ria.rds") # R object (list) with the variable importances at each iteration of the algorithm

RFE_info <- data.frame(Worst_Var=var.list, OOB_rmse=OOB.rmse, OOB_mae=OOB.mae, Train_rmse=train.rmse, Train_mae=train.mae)

write.csv(RFE_info, "cf_RFE_info_thresh90_ria.csv", row.names = FALSE)

#set mtry = # of predictors left in model / 3 (nump/3) - and round it up to solid number 
# Refit final model
ria_final <- ria %>% dplyr::select(Effect, NR1I2_Activation,
                                   Desvenlafaxine, "SLC-Inhibition", Venlafaxine)
                                   
set.seed(3)
cf.final <- cforest(Effect~., data=ria_final, control=cforest_unbiased(mtry=ceiling(1), ntree=5000))

# OOB predictions
ria_final.pred.oob <- data.frame(predict(cf.final, OOB=T)) %>% rename(Predicted_OOB=Effect)
ria_final_wPred <- cbind(ria_final, ria_final.pred.oob)

# OOB error
mean(abs(ria_final$Effect-ria_final.pred.oob$Predicted_OOB)) #0.408791
sqrt(mean((ria_final$Effect-ria_final.pred.oob$Predicted_OOB)^2)) #1.536085

# Predicted (OOB) vs. observed
plot(Predicted_OOB~Effect, data=ria_final_wPred, pch=16)
abline(a=0, b=1)

#RIA F####
names(ria_F)
ria_F_sub <- ria_F[,c(2:111)]
ria_F_sub[,1] <- round(ria_F_sub[,1],2)
View(ria_F_sub)

p <- ncol(ria_F_sub) - 1

#cforest - all variables
set.seed(3)
cf.all <- cforest(Effect~., data =ria_F_sub, control = cforest_unbiased(mtry =ceiling(p/3), ntree = 5000))
cf.pred <- predict(cf.all, OOB = T)
sqrt(mean((ria_F_sub$Effect-cf.pred)^2)) #2.748557
mean(abs(ria_F_sub$Effect - cf.pred)) #1.935697

# Unconditional variable importance (threshold=1)
set.seed(3) 
cf.upi <- permimp(cf.all, conditional=TRUE, threshold=1, thresholdDiagnostics = TRUE)
as.matrix(sort(cf.upi$values, decreasing=T)) # CF, unconditional

# Conditional variable importance (threshold=0.9)
set.seed(3) 
cf.cpi <- permimp(cf.all, conditional=TRUE, threshold=0.90, thresholdDiagnostics = TRUE)
as.matrix(sort(cf.cpi$values, decreasing=T)) # CF, conditional

# cf.cpi$perTree # VI values per tree
plot(cf.cpi, type="bar")

### Loop to recursively eliminate least informative variables ##

# Pull out response and predictors - columns 5, 7-53
ria_F_preds <- ria_F[, c(2:111)] # original predictors
ria_F_preds[,1] <- round(ria_F_preds[,1], 2)
ria_F_preds_loop <- ria_F_preds 
View(ria_F_preds_loop)

p <- ncol(ria_F_preds)-1

# Save OOB error, variable importance
VI.list <- NULL
var.list <- rep(NA,p)

OOB.rmse <- rep(NA,p)
OOB.mae <- rep(NA,p)

train.rmse <- rep(NA,p)
train.mae <- rep(NA,p)

#run RF - conditional
set.seed(3) # Set seed for reproducibility

for(i in 1:p){
  print(i)
  
  nump <- ncol(ria_F_preds_loop)-1 # Number predictors
  
  set.seed(3) 
  cf <- cforest(Effect~., data=ria_F_preds_loop, control=cforest_unbiased(mtry=ceiling(nump/3), ntree=5000))
  cf.pred <- predict(cf, OOB=T)
  cf.pred.tr <- predict(cf, OOB=F)
  
  # OOB error
  OOB.rmse[i] <- sqrt(mean((ria_F_preds$Effect-cf.pred)^2)) # OOB rmse
  OOB.mae[i] <- mean(abs(ria_F_preds$Effect-cf.pred)) # oob mae
  
  # Training error
  train.rmse[i] <- sqrt(mean((ria_F_preds$Effect-cf.pred.tr)^2)) # training rmse
  train.mae[i] <- mean(abs(ria_F_preds$Effect-cf.pred.tr)) #  training mae
  
  
  # Variable importance
  set.seed(3) 
  cf.cpi <- permimp(cf, conditional=TRUE, threshold=0.9) # set threshold to 1 for unconditional
  
  VI.list[[i]] <- as.matrix(sort(cf.cpi$values, decreasing=F))
  var.list[i] <- row.names(VI.list[[i]])[1] # Least informative variable to remove
  
  
  ria_F_preds_loop <- ria_F_preds_loop[,!colnames(ria_F_preds_loop)==var.list[i]] # remove variable
  
}
beep()

# Plot OOB error - minimum is the best subset of predictors
plot(OOB.rmse, xlab="Iteration", type='l')
points(OOB.rmse, pch=16)
plot(OOB.mae, xlab="Iteration", type='l')
points(OOB.mae, pch=16)

# training error (not to be used for optimizing number predictors)
plot(train.rmse, xlab="Iteration", type='l')
points(train.rmse, pch=16)

which(OOB.rmse==min(OOB.rmse)) #107
which(OOB.mae==min(OOB.mae)) #107

var.list[which(OOB.mae==min(OOB.mae)):length(OOB.mae)]

#double check with rounded error#
OOB.rmse_1 <- round(OOB.rmse, 2)
OOB.mae_1 <- round(OOB.mae, 2)
which(OOB.rmse_1==min(OOB.rmse_1)) #108
which(OOB.mae_1==min(OOB.mae_1)) #107

var.list[which(OOB.mae_1==min(OOB.mae_1)):length(OOB.mae_1)] #same

# Save output - name files in a way that lets you know what differs between different runs of the algorithm

saveRDS(VI.list, "cf_VI_at_iter_thresh90_ria_F.rds") # R object (list) with the variable importances at each iteration of the algorithm

RFE_info <- data.frame(Worst_Var=var.list, OOB_rmse=OOB.rmse, OOB_mae=OOB.mae, Train_rmse=train.rmse, Train_mae=train.mae)

write.csv(RFE_info, "cf_RFE_info_thresh90_ria_F.csv", row.names = FALSE)

#set mtry = # of predictors left in model / 3 (nump/3) - and round it up to solid number 
# Refit final model
ria_F_final <- ria_F %>% dplyr::select(Effect, "SCN-Inhibition", "Tramadol", "SLC-Inhibition")
set.seed(3)
cf.final <- cforest(Effect~., data=ria_F_final, control=cforest_unbiased(mtry=1, ntree=5000))

# OOB predictions
ria_F_final.pred.oob <- data.frame(predict(cf.final, OOB=T)) %>% rename(Predicted_OOB=Effect)
ria_F_final_wPred <- cbind(ria_F_final, ria_F_final.pred.oob)

# OOB error
mean(abs(ria_F_final$Effect-ria_F_final.pred.oob$Predicted_OOB)) #1.897099
sqrt(mean((ria_F_final$Effect-ria_F_final.pred.oob$Predicted_OOB)^2)) #2.655092


#CYP1A1_Intestine####
names(cyp_1a1_intestine)
cyp_1a1_intestine_sub <- cyp_1a1_intestine[,c(2:111)]
cyp_1a1_intestine_sub[,1] <- round(cyp_1a1_intestine_sub[,1],2)
View(cyp_1a1_intestine_sub)

p <- ncol(cyp_1a1_intestine_sub) - 1

#cforest - all variables
set.seed(3)
cf.all <- cforest(Effect~., data =cyp_1a1_intestine_sub, control = cforest_unbiased(mtry =ceiling(p/3), ntree = 5000))
cf.pred <- predict(cf.all, OOB = T)
sqrt(mean((cyp_1a1_intestine_sub$Effect-cf.pred)^2)) #13.24071
mean(abs(cyp_1a1_intestine_sub$Effect - cf.pred)) #9.190668

#plot predicted vs. observed
plot(y=cf.pred, x= cyp_1a1_intestine_sub$Effect)
abline(a=0, b=1)

# Unconditional variable importance (threshold=1)
set.seed(3) 
cf.upi <- permimp(cf.all, conditional=TRUE, threshold=1, thresholdDiagnostics = TRUE)
as.matrix(sort(cf.upi$values, decreasing=T)) # CF, unconditional

# Conditional variable importance (threshold=0.9)
set.seed(3) 
cf.cpi <- permimp(cf.all, conditional=TRUE, threshold=0.90, thresholdDiagnostics = TRUE)
as.matrix(sort(cf.cpi$values, decreasing=T)) # CF, conditional

# cf.cpi$perTree # VI values per tree
plot(cf.cpi, type="bar")

### Loop to recursively eliminate least informative variables ##

# Pull out response and predictors - columns 5, 7-53
cyp_1a1_intestine_preds <- cyp_1a1_intestine[, c(2:111)] # original predictors
cyp_1a1_intestine_preds[,1] <- round(cyp_1a1_intestine_preds[,1], 2)
cyp_1a1_intestine_preds_loop <- cyp_1a1_intestine_preds 

p <- ncol(cyp_1a1_intestine_preds)-1

# Save OOB error, variable importance
VI.list <- NULL
var.list <- rep(NA,p)

OOB.rmse <- rep(NA,p)
OOB.mae <- rep(NA,p)

train.rmse <- rep(NA,p)
train.mae <- rep(NA,p)

#run RF - conditional
set.seed(3) # Set seed for reproducibility

for(i in 1:p){
  print(i)
  
  nump <- ncol(cyp_1a1_intestine_preds_loop)-1 # Number predictors
  
  set.seed(3) 
  cf <- cforest(Effect~., data=cyp_1a1_intestine_preds_loop, control=cforest_unbiased(mtry=ceiling(nump/3), ntree=5000))
  cf.pred <- predict(cf, OOB=T)
  cf.pred.tr <- predict(cf, OOB=F)
  
  # OOB error
  OOB.rmse[i] <- sqrt(mean((cyp_1a1_intestine_preds$Effect-cf.pred)^2)) # OOB rmse
  OOB.mae[i] <- mean(abs(cyp_1a1_intestine_preds$Effect-cf.pred)) # oob mae
  
  # Training error
  train.rmse[i] <- sqrt(mean((cyp_1a1_intestine_preds$Effect-cf.pred.tr)^2)) # training rmse
  train.mae[i] <- mean(abs(cyp_1a1_intestine_preds$Effect-cf.pred.tr)) #  training mae
  
  
  # Variable importance
  set.seed(3) 
  cf.cpi <- permimp(cf, conditional=TRUE, threshold=0.9) # set threshold to 1 for unconditional
  
  VI.list[[i]] <- as.matrix(sort(cf.cpi$values, decreasing=F))
  var.list[i] <- row.names(VI.list[[i]])[1] # Least informative variable to remove
  
  
  cyp_1a1_intestine_preds_loop <- cyp_1a1_intestine_preds_loop[,!colnames(cyp_1a1_intestine_preds_loop)==var.list[i]] # remove variable
  
}
beep()

# Plot OOB error - minimum is the best subset of predictors
plot(OOB.rmse, xlab="Iteration", type='l')
points(OOB.rmse, pch=16)
plot(OOB.mae, xlab="Iteration", type='l')
points(OOB.mae, pch=16)

# training error (not to be used for optimizing number predictors)
plot(train.rmse, xlab="Iteration", type='l')
points(train.rmse, pch=16)

which(OOB.rmse==min(OOB.rmse)) #95
which(OOB.mae==min(OOB.mae)) #95

var.list[which(OOB.mae==min(OOB.mae)):length(OOB.mae)]
# [1] "PGR_Inhibition"             "Benzo[a]pyrene"             "ADORA2_ADRA_HRH-Inhibition" "Pyrene"                    
# [5] "PPARG_Activation"           "androstenedione_Inhibition" "NR1I3_Binding"              "NR1I2_Activation"          
# [9] "Anthraquinone"              "ESR2_Activation"            "NR1I3_Activation"           "RXRA_Activation"           
# [13] "AR_Inhibition" 
#double check with rounded error#
OOB.rmse_1 <- round(OOB.rmse, 2)
OOB.mae_1 <- round(OOB.mae, 2)
which(OOB.rmse_1==min(OOB.rmse_1)) #95
which(OOB.mae_1==min(OOB.mae_1)) #95

var.list[which(OOB.mae_1==min(OOB.mae_1)):length(OOB.mae_1)] #same

# Save output - name files in a way that lets you know what differs between different runs of the algorithm

saveRDS(VI.list, "cf_VI_at_iter_thresh90_cyp1a1_intest.rds") # R object (list) with the variable importances at each iteration of the algorithm

RFE_info <- data.frame(Worst_Var=var.list, OOB_rmse=OOB.rmse, OOB_mae=OOB.mae, Train_rmse=train.rmse, Train_mae=train.mae)

write.csv(RFE_info, "cf_RFE_info_thresh90_cyp1a1_intest.csv", row.names = FALSE)

#set mtry = # of predictors left in model / 3 (nump/3) - and round it up to solid number 
# Refit final model
cyp_1a1_intestine_final <- cyp_1a1_intestine %>% dplyr::select(Effect,
                                                               "PGR_Inhibition","Benzo[a]pyrene","ADORA2_ADRA_HRH-Inhibition", "Pyrene",
                                                               "PPARG_Activation", "androstenedione_Inhibition", "NR1I3_Binding", "NR1I2_Activation",
                                                               "Anthraquinone", "ESR2_Activation", "NR1I3_Activation", "RXRA_Activation", "AR_Inhibition")
set.seed(3)
cf.final <- cforest(Effect~., data=cyp_1a1_intestine_final, control=cforest_unbiased(mtry=4, ntree=5000))

# OOB predictions
cyp_1a1_intestine_final.pred.oob <- data.frame(predict(cf.final, OOB=T)) %>% rename(Predicted_OOB=Effect)
cyp_1a1_intestine_final_wPred <- cbind(cyp_1a1_intestine_final, cyp_1a1_intestine_final.pred.oob)

# OOB error
mean(abs(cyp_1a1_intestine_final$Effect-cyp_1a1_intestine_final.pred.oob$Predicted_OOB)) #9.129761
sqrt(mean((cyp_1a1_intestine_final$Effect-cyp_1a1_intestine_final.pred.oob$Predicted_OOB)^2)) #13.16639

# Predicted (OOB) vs. observed
plot(Predicted_OOB~Effect, data=cyp_1a1_intestine_final_wPred, pch=16)
abline(a=0, b=1)

#CYP1A1_Liver_2017####
names(cyp_1a1_liver_2017)
cyp_1a1_liver_2017_sub <- cyp_1a1_liver_2017[,c(2:111)]
cyp_1a1_liver_2017_sub[,1] <- round(cyp_1a1_liver_2017_sub[,1],2)
View(cyp_1a1_liver_2017_sub)

p <- ncol(cyp_1a1_liver_2017_sub) - 1

#cforest - all variables
set.seed(3)
cf.all <- cforest(Effect~., data =cyp_1a1_liver_2017_sub, control = cforest_unbiased(mtry =ceiling(p/3), ntree = 5000))
cf.pred <- predict(cf.all, OOB = T)
sqrt(mean((cyp_1a1_liver_2017_sub$Effect-cf.pred)^2)) #2.255933
mean(abs(cyp_1a1_liver_2017_sub$Effect - cf.pred)) #1.238203

#plot predicted vs. observed
plot(y=cf.pred, x= cyp_1a1_liver_2017_sub$Effect)
abline(a=0, b=1)

# Unconditional variable importance (threshold=1)
set.seed(3) 
cf.upi <- permimp(cf.all, conditional=TRUE, threshold=1, thresholdDiagnostics = TRUE)
as.matrix(sort(cf.upi$values, decreasing=T)) # CF, unconditional

# Conditional variable importance (threshold=0.9)
set.seed(3) 
cf.cpi <- permimp(cf.all, conditional=TRUE, threshold=0.90, thresholdDiagnostics = TRUE)
as.matrix(sort(cf.cpi$values, decreasing=T)) # CF, conditional

# cf.cpi$perTree # VI values per tree
plot(cf.cpi, type="bar")

### Loop to recursively eliminate least informative variables ##

# Pull out response and predictors - columns 5, 7-53
cyp_1a1_liver_2017_preds <- cyp_1a1_liver_2017[, c(2:111)] # original predictors
cyp_1a1_liver_2017_preds[,1] <- round(cyp_1a1_liver_2017_preds[,1], 2)
cyp_1a1_liver_2017_preds_loop <- cyp_1a1_liver_2017_preds 
View(cyp_1a1_liver_2017_preds_loop)

p <- ncol(cyp_1a1_liver_2017_preds)-1

# Save OOB error, variable importance
VI.list <- NULL
var.list <- rep(NA,p)

OOB.rmse <- rep(NA,p)
OOB.mae <- rep(NA,p)

train.rmse <- rep(NA,p)
train.mae <- rep(NA,p)

#run RF - conditional
set.seed(3) # Set seed for reproducibility

for(i in 1:p){
  print(i)
  
  nump <- ncol(cyp_1a1_liver_2017_preds_loop)-1 # Number predictors
  
  set.seed(3) 
  cf <- cforest(Effect~., data=cyp_1a1_liver_2017_preds_loop, control=cforest_unbiased(mtry=ceiling(nump/3), ntree=5000))
  cf.pred <- predict(cf, OOB=T)
  cf.pred.tr <- predict(cf, OOB=F)
  
  # OOB error
  OOB.rmse[i] <- sqrt(mean((cyp_1a1_liver_2017_preds$Effect-cf.pred)^2)) # OOB rmse
  OOB.mae[i] <- mean(abs(cyp_1a1_liver_2017_preds$Effect-cf.pred)) # oob mae
  
  # Training error
  train.rmse[i] <- sqrt(mean((cyp_1a1_liver_2017_preds$Effect-cf.pred.tr)^2)) # training rmse
  train.mae[i] <- mean(abs(cyp_1a1_liver_2017_preds$Effect-cf.pred.tr)) #  training mae
  
  
  # Variable importance
  set.seed(3) 
  cf.cpi <- permimp(cf, conditional=TRUE, threshold=0.9) # set threshold to 1 for unconditional
  
  VI.list[[i]] <- as.matrix(sort(cf.cpi$values, decreasing=F))
  var.list[i] <- row.names(VI.list[[i]])[1] # Least informative variable to remove
  
  
  cyp_1a1_liver_2017_preds_loop <- cyp_1a1_liver_2017_preds_loop[,!colnames(cyp_1a1_liver_2017_preds_loop)==var.list[i]] # remove variable
  
}
beep()

# Plot OOB error - minimum is the best subset of predictors
plot(OOB.rmse, xlab="Iteration", type='l')
points(OOB.rmse, pch=16)
plot(OOB.mae, xlab="Iteration", type='l')
points(OOB.mae, pch=16)

# training error (not to be used for optimizing number predictors)
plot(train.rmse, xlab="Iteration", type='l')
points(train.rmse, pch=16)

which(OOB.rmse==min(OOB.rmse)) #105
which(OOB.mae==min(OOB.mae)) #107

var.list[which(OOB.mae==min(OOB.mae)):length(OOB.mae)]

#double check with rounded error#
OOB.rmse_1 <- round(OOB.rmse, 2)
OOB.mae_1 <- round(OOB.mae, 2)
which(OOB.rmse_1==min(OOB.rmse_1)) #105
which(OOB.mae_1==min(OOB.mae_1)) #107

var.list[which(OOB.mae_1==min(OOB.mae_1)):length(OOB.mae_1)] #same

# Save output - name files in a way that lets you know what differs between different runs of the algorithm

saveRDS(VI.list, "cf_VI_at_iter_thresh90_cyp1a1_liv.rds") # R object (list) with the variable importances at each iteration of the algorithm

RFE_info <- data.frame(Worst_Var=var.list, OOB_rmse=OOB.rmse, OOB_mae=OOB.mae, Train_rmse=train.rmse, Train_mae=train.mae)

write.csv(RFE_info, "cf_RFE_info_thresh90_cyp1a1_liv.csv", row.names = FALSE)

#set mtry = # of predictors left in model / 3 (nump/3) - and round it up to solid number 
# Refit final model
cyp_1a1_liver_2017_final <- cyp_1a1_liver_2017 %>% dplyr::select(Effect, ESR2_Activation)
set.seed(3)
cf.final <- cforest(Effect~., data=cyp_1a1_liver_2017_final, control=cforest_unbiased(mtry=1, ntree=5000))

# OOB predictions
cyp_1a1_liver_2017_final.pred.oob <- data.frame(predict(cf.final, OOB=T)) %>% rename(Predicted_OOB=Effect)
cyp_1a1_liver_2017_final_wPred <- cbind(cyp_1a1_liver_2017_final, cyp_1a1_liver_2017_final.pred.oob)

# OOB error
mean(abs(cyp_1a1_liver_2017_final$Effect-cyp_1a1_liver_2017_final.pred.oob$Predicted_OOB)) #1.232205
sqrt(mean((cyp_1a1_liver_2017_final$Effect-cyp_1a1_liver_2017_final.pred.oob$Predicted_OOB)^2)) #2.268648

# Predicted (OOB) vs. observed
plot(Predicted_OOB~Effect, data=cyp_1a1_liver_2017_final_wPred, pch=16)
abline(a=0, b=1)

#CYP3A_liver####
names(cyp3a_liver)
cyp3a_liver_sub <- cyp3a_liver[,c(2:111)]
cyp3a_liver_sub[,1] <- round(cyp3a_liver_sub[,1],2)

p <- ncol(cyp3a_liver_sub) - 1

#cforest - all variables
set.seed(3)
cf.all <- cforest(Effect~., data =cyp3a_liver_sub, control = cforest_unbiased(mtry =ceiling(p/3), ntree = 5000))
cf.pred <- predict(cf.all, OOB = T)
sqrt(mean((cyp3a_liver_sub$Effect-cf.pred)^2)) #0.5759348
mean(abs(cyp3a_liver_sub$Effect - cf.pred)) #0.4592107

#plot predicted vs. observed
plot(y=cf.pred, x= cyp3a_liver_sub$Effect)
abline(a=0, b=1)

# Unconditional variable importance (threshold=1)
set.seed(3) 
cf.upi <- permimp(cf.all, conditional=TRUE, threshold=1, thresholdDiagnostics = TRUE)
as.matrix(sort(cf.upi$values, decreasing=T)) # CF, unconditional

# Conditional variable importance (threshold=0.9)
set.seed(3) 
cf.cpi <- permimp(cf.all, conditional=TRUE, threshold=0.90, thresholdDiagnostics = TRUE)
as.matrix(sort(cf.cpi$values, decreasing=T)) # CF, conditional

# cf.cpi$perTree # VI values per tree
plot(cf.cpi, type="bar")

### Loop to recursively eliminate least informative variables ##

# Pull out response and predictors - columns 5, 7-53
cyp3a_liver_preds <- cyp3a_liver[, c(2:111)] # original predictors
cyp3a_liver_preds[,1] <- round(cyp3a_liver_preds[,1], 2)
cyp3a_liver_preds_loop <- cyp3a_liver_preds 
View(cyp3a_liver_preds_loop)

p <- ncol(cyp3a_liver_preds)-1

# Save OOB error, variable importance
VI.list <- NULL
var.list <- rep(NA,p)

OOB.rmse <- rep(NA,p)
OOB.mae <- rep(NA,p)

train.rmse <- rep(NA,p)
train.mae <- rep(NA,p)

#run RF - conditional
set.seed(3) # Set seed for reproducibility

for(i in 1:p){
  print(i)
  
  nump <- ncol(cyp3a_liver_preds_loop)-1 # Number predictors
  
  set.seed(3) 
  cf <- cforest(Effect~., data=cyp3a_liver_preds_loop, control=cforest_unbiased(mtry=ceiling(nump/3), ntree=5000))
  cf.pred <- predict(cf, OOB=T)
  cf.pred.tr <- predict(cf, OOB=F)
  
  # OOB error
  OOB.rmse[i] <- sqrt(mean((cyp3a_liver_preds$Effect-cf.pred)^2)) # OOB rmse
  OOB.mae[i] <- mean(abs(cyp3a_liver_preds$Effect-cf.pred)) # oob mae
  
  # Training error
  train.rmse[i] <- sqrt(mean((cyp3a_liver_preds$Effect-cf.pred.tr)^2)) # training rmse
  train.mae[i] <- mean(abs(cyp3a_liver_preds$Effect-cf.pred.tr)) #  training mae
  
  
  # Variable importance
  set.seed(3) 
  cf.cpi <- permimp(cf, conditional=TRUE, threshold=0.9) # set threshold to 1 for unconditional
  
  VI.list[[i]] <- as.matrix(sort(cf.cpi$values, decreasing=F))
  var.list[i] <- row.names(VI.list[[i]])[1] # Least informative variable to remove
  
  
  cyp3a_liver_preds_loop <- cyp3a_liver_preds_loop[,!colnames(cyp3a_liver_preds_loop)==var.list[i]] # remove variable
  
}
beep()

# Plot OOB error - minimum is the best subset of predictors
plot(OOB.rmse, xlab="Iteration", type='l')
points(OOB.rmse, pch=16)
plot(OOB.mae, xlab="Iteration", type='l')
points(OOB.mae, pch=16)

# training error (not to be used for optimizing number predictors)
plot(train.rmse, xlab="Iteration", type='l')
points(train.rmse, pch=16)

which(OOB.rmse==min(OOB.rmse)) #109
which(OOB.mae==min(OOB.mae)) #108

var.list[which(OOB.mae==min(OOB.mae)):length(OOB.mae)]

#double check with rounded error#
OOB.rmse_1 <- round(OOB.rmse, 2)
OOB.mae_1 <- round(OOB.mae, 2)
which(OOB.rmse_1==min(OOB.rmse_1)) #109
which(OOB.mae_1==min(OOB.mae_1)) #108

var.list[which(OOB.mae_1==min(OOB.mae_1)):length(OOB.mae_1)]

# Save output - name files in a way that lets you know what differs between different runs of the algorithm

saveRDS(VI.list, "cf_VI_at_iter_thresh90_cyp3a_liv.rds") # R object (list) with the variable importances at each iteration of the algorithm

RFE_info <- data.frame(Worst_Var=var.list, OOB_rmse=OOB.rmse, OOB_mae=OOB.mae, Train_rmse=train.rmse, Train_mae=train.mae)

write.csv(RFE_info, "cf_RFE_info_thresh90_cyp3a_liv.csv", row.names = FALSE)

#set mtry = # of predictors left in model / 3 (nump/3) - and round it up to solid number 
# Refit final model
cyp3a_liver_final <- cyp3a_liver %>% dplyr::select(Effect, Indole, Carbazole)
set.seed(3)
cf.final <- cforest(Effect~., data=cyp3a_liver_final, control=cforest_unbiased(mtry=1, ntree=5000))

# OOB predictions
cyp3a_liver_final.pred.oob <- data.frame(predict(cf.final, OOB=T)) %>% rename(Predicted_OOB=Effect)
cyp3a_liver_final_final_wPred <- cbind(cyp3a_liver_final, cyp3a_liver_final.pred.oob)

# OOB error
mean(abs(cyp3a_liver_final$Effect-cyp3a_liver_final.pred.oob$Predicted_OOB)) #0.4534382
sqrt(mean((cyp3a_liver_final$Effect-cyp3a_liver_final.pred.oob$Predicted_OOB)^2)) #0.5643797

# Predicted (OOB) vs. observed
plot(Predicted_OOB~Effect, data=cyp3a_liver_final_final_wPred, pch=16)
abline(a=0, b=1)

#CYP3A_intestine####
names(CYP3A_intestine)
CYP3A_intestine_sub <- CYP3A_intestine[,c(2:111)]
CYP3A_intestine_sub[,1] <- round(CYP3A_intestine_sub[,1],2)
View(CYP3A_intestine_sub)

p <- ncol(CYP3A_intestine_sub) - 1

#cforest - all variables
set.seed(3)
cf.all <- cforest(Effect~., data =CYP3A_intestine_sub, control = cforest_unbiased(mtry =ceiling(p/3), ntree = 5000))
cf.pred <- predict(cf.all, OOB = T)
sqrt(mean((CYP3A_intestine_sub$Effect-cf.pred)^2)) #2.047725
mean(abs(CYP3A_intestine_sub$Effect - cf.pred)) #1.5011

#plot predicted vs. observed
plot(y=cf.pred, x= CYP3A_intestine_sub$Effect)
abline(a=0, b=1)

# Unconditional variable importance (threshold=1)
set.seed(3) 
cf.upi <- permimp(cf.all, conditional=TRUE, threshold=1, thresholdDiagnostics = TRUE)
as.matrix(sort(cf.upi$values, decreasing=T)) # CF, unconditional

# Conditional variable importance (threshold=0.9)
set.seed(3) 
cf.cpi <- permimp(cf.all, conditional=TRUE, threshold=0.90, thresholdDiagnostics = TRUE)
as.matrix(sort(cf.cpi$values, decreasing=T)) # CF, conditional

# cf.cpi$perTree # VI values per tree
plot(cf.cpi, type="bar")

### Loop to recursively eliminate least informative variables ##

# Pull out response and predictors - columns 5, 7-53
CYP3A_intestine_preds <- CYP3A_intestine[, c(2:111)] # original predictors
CYP3A_intestine_preds[,1] <- round(CYP3A_intestine_preds[,1], 2)
CYP3A_intestine_preds_loop <- CYP3A_intestine_preds 
View(CYP3A_intestine_preds_loop)

p <- ncol(CYP3A_intestine_preds)-1

# Save OOB error, variable importance
VI.list <- NULL
var.list <- rep(NA,p)

OOB.rmse <- rep(NA,p)
OOB.mae <- rep(NA,p)

train.rmse <- rep(NA,p)
train.mae <- rep(NA,p)

#run RF - conditional
set.seed(3) # Set seed for reproducibility

for(i in 1:p){
  print(i)
  
  nump <- ncol(CYP3A_intestine_preds_loop)-1 # Number predictors
  
  set.seed(3) 
  cf <- cforest(Effect~., data=CYP3A_intestine_preds_loop, control=cforest_unbiased(mtry=ceiling(nump/3), ntree=5000))
  cf.pred <- predict(cf, OOB=T)
  cf.pred.tr <- predict(cf, OOB=F)
  
  # OOB error
  OOB.rmse[i] <- sqrt(mean((CYP3A_intestine_preds$Effect-cf.pred)^2)) # OOB rmse
  OOB.mae[i] <- mean(abs(CYP3A_intestine_preds$Effect-cf.pred)) # oob mae
  
  # Training error
  train.rmse[i] <- sqrt(mean((CYP3A_intestine_preds$Effect-cf.pred.tr)^2)) # training rmse
  train.mae[i] <- mean(abs(CYP3A_intestine_preds$Effect-cf.pred.tr)) #  training mae
  
  
  # Variable importance
  set.seed(3) 
  cf.cpi <- permimp(cf, conditional=TRUE, threshold=0.9) # set threshold to 1 for unconditional
  
  VI.list[[i]] <- as.matrix(sort(cf.cpi$values, decreasing=F))
  var.list[i] <- row.names(VI.list[[i]])[1] # Least informative variable to remove
  
  
  CYP3A_intestine_preds_loop <- CYP3A_intestine_preds_loop[,!colnames(CYP3A_intestine_preds_loop)==var.list[i]] # remove variable
  
}
beep()

# Plot OOB error - minimum is the best subset of predictors
plot(OOB.rmse, xlab="Iteration", type='l')
points(OOB.rmse, pch=16)
plot(OOB.mae, xlab="Iteration", type='l')
points(OOB.mae, pch=16)

# training error (not to be used for optimizing number predictors)
plot(train.rmse, xlab="Iteration", type='l')
points(train.rmse, pch=16)

which(OOB.rmse==min(OOB.rmse)) #108
which(OOB.mae==min(OOB.mae)) #108

var.list[which(OOB.mae==min(OOB.mae)):length(OOB.mae)]
# "NR1I3_Activation" "AR_Inhibition" 

#double check with rounded error#
OOB.rmse_1 <- round(OOB.rmse, 2)
OOB.mae_1 <- round(OOB.mae, 2)
which(OOB.rmse_1==min(OOB.rmse_1)) #108
which(OOB.mae_1==min(OOB.mae_1)) #108


# Save output - name files in a way that lets you know what differs between different runs of the algorithm

saveRDS(VI.list, "cf_VI_at_iter_thresh90_cyp3a_int.rds") # R object (list) with the variable importances at each iteration of the algorithm

RFE_info <- data.frame(Worst_Var=var.list, OOB_rmse=OOB.rmse, OOB_mae=OOB.mae, Train_rmse=train.rmse, Train_mae=train.mae)

write.csv(RFE_info, "cf_RFE_info_thresh90_cyp3a_int.csv", row.names = FALSE)

#set mtry = # of predictors left in model / 3 (nump/3) - and round it up to solid number 
# Refit final model
CYP3A_intestine_final <- CYP3A_intestine %>% dplyr::select(Effect, "NR1I3_Activation", AR_Inhibition)
set.seed(3)
cf.final <- cforest(Effect~., data=CYP3A_intestine_final, control=cforest_unbiased(mtry=1, ntree=5000))

# OOB predictions
CYP3A_intestine_final.pred.oob <- data.frame(predict(cf.final, OOB=T)) %>% rename(Predicted_OOB=Effect)
CYP3A_intestine_final_wPred <- cbind(CYP3A_intestine_final, CYP3A_intestine_final.pred.oob)

# OOB error
mean(abs(CYP3A_intestine_final$Effect-CYP3A_intestine_final.pred.oob$Predicted_OOB)) #1.469033
sqrt(mean((CYP3A_intestine_final$Effect-CYP3A_intestine_final.pred.oob$Predicted_OOB)^2)) #2.032153

# Predicted (OOB) vs. observed
plot(Predicted_OOB~Effect, data=CYP3A_intestine_final_wPred, pch=16)
abline(a=0, b=1)

#UGT1A1_liver####
names(UGT_liver)
UGT_liver_sub <- UGT_liver[,c(2:111)]
UGT_liver_sub[,1] <- round(UGT_liver_sub[,1],2)
View(UGT_liver_sub)

p <- ncol(UGT_liver_sub) - 1

#cforest - all variables
set.seed(3)
cf.all <- cforest(Effect~., data =UGT_liver_sub, control = cforest_unbiased(mtry =ceiling(p/3), ntree = 5000))
cf.pred <- predict(cf.all, OOB = T)
sqrt(mean((UGT_liver_sub$Effect-cf.pred)^2)) #1.832809
mean(abs(UGT_liver_sub$Effect - cf.pred)) #1.242095

#plot predicted vs. observed
plot(y=cf.pred, x= UGT_liver_sub$Effect)
abline(a=0, b=1)

# Unconditional variable importance (threshold=1)
set.seed(3) 
cf.upi <- permimp(cf.all, conditional=TRUE, threshold=1, thresholdDiagnostics = TRUE)
as.matrix(sort(cf.upi$values, decreasing=T)) # CF, unconditional

# Conditional variable importance (threshold=0.9)
set.seed(3) 
cf.cpi <- permimp(cf.all, conditional=TRUE, threshold=0.90, thresholdDiagnostics = TRUE)
as.matrix(sort(cf.cpi$values, decreasing=T)) # CF, conditional

# cf.cpi$perTree # VI values per tree
plot(cf.cpi, type="bar")

### Loop to recursively eliminate least informative variables ##

# Pull out response and predictors - columns 5, 7-53
UGT_liver_preds <- UGT_liver[, c(2:111)] # original predictors
UGT_liver_preds[,1] <- round(UGT_liver_preds[,1], 2)
UGT_liver_preds_loop <- UGT_liver_preds 

p <- ncol(UGT_liver_preds)-1

# Save OOB error, variable importance
VI.list <- NULL
var.list <- rep(NA,p)

OOB.rmse <- rep(NA,p)
OOB.mae <- rep(NA,p)

train.rmse <- rep(NA,p)
train.mae <- rep(NA,p)

#run RF - conditional
set.seed(3) # Set seed for reproducibility

for(i in 1:p){
  print(i)
  
  nump <- ncol(UGT_liver_preds_loop)-1 # Number predictors
  
  set.seed(3) 
  cf <- cforest(Effect~., data=UGT_liver_preds_loop, control=cforest_unbiased(mtry=ceiling(nump/3), ntree=5000))
  cf.pred <- predict(cf, OOB=T)
  cf.pred.tr <- predict(cf, OOB=F)
  
  # OOB error
  OOB.rmse[i] <- sqrt(mean((UGT_liver_preds$Effect-cf.pred)^2)) # OOB rmse
  OOB.mae[i] <- mean(abs(UGT_liver_preds$Effect-cf.pred)) # oob mae
  
  # Training error
  train.rmse[i] <- sqrt(mean((UGT_liver_preds$Effect-cf.pred.tr)^2)) # training rmse
  train.mae[i] <- mean(abs(UGT_liver_preds$Effect-cf.pred.tr)) #  training mae
  
  
  # Variable importance
  set.seed(3) 
  cf.cpi <- permimp(cf, conditional=TRUE, threshold=0.9) # set threshold to 1 for unconditional
  
  VI.list[[i]] <- as.matrix(sort(cf.cpi$values, decreasing=F))
  var.list[i] <- row.names(VI.list[[i]])[1] # Least informative variable to remove
  
  
  UGT_liver_preds_loop <- UGT_liver_preds_loop[,!colnames(UGT_liver_preds_loop)==var.list[i]] # remove variable
  
}
beep()

# Plot OOB error - minimum is the best subset of predictors
plot(OOB.rmse, xlab="Iteration", type='l')
points(OOB.rmse, pch=16)
plot(OOB.mae, xlab="Iteration", type='l')
points(OOB.mae, pch=16)

# training error (not to be used for optimizing number predictors)
plot(train.rmse, xlab="Iteration", type='l')
points(train.rmse, pch=16)

which(OOB.rmse==min(OOB.rmse)) #86
which(OOB.mae==min(OOB.mae)) #97

var.list[which(OOB.mae==min(OOB.mae)):length(OOB.mae)]
# [1] "CYP2B6_Activation"            "Tramadol"                    
# [3] "ESRRA_Inhibition"             "Scn1a_Inhibition"            
# [5] "ESR1_Activation"              "Oprk1_Binding"               
# [7] "OP_Fire_Retardants"           "ADORA2_ADRA_HRH-Inhibition"  
# [9] "PGR_Inhibition"               "NR1I2_Binding"               
# [11] "NR1I2_Activation"             "PPARG_Activation"            
# [13] "Tris(2-butoxyethyl)phosphate"

#double check with rounded error#
OOB.rmse_1 <- round(OOB.rmse, 2)
OOB.mae_1 <- round(OOB.mae, 2)
which(OOB.rmse_1==min(OOB.rmse_1)) #109
which(OOB.mae_1==min(OOB.mae_1)) #97

var.list[which(OOB.mae_1==min(OOB.mae_1)):length(OOB.mae_1)] #same


# Save output - name files in a way that lets you know what differs between different runs of the algorithm

saveRDS(VI.list, "cf_VI_at_iter_thresh90_ugt_liv.rds") # R object (list) with the variable importances at each iteration of the algorithm

RFE_info <- data.frame(Worst_Var=var.list, OOB_rmse=OOB.rmse, OOB_mae=OOB.mae, Train_rmse=train.rmse, Train_mae=train.mae)

write.csv(RFE_info, "cf_RFE_info_thresh90_ugt_liv.csv", row.names = FALSE)

#set mtry = # of predictors left in model / 3 (nump/3) - and round it up to solid number 
# Refit final model
UGT_liver_final <- UGT_liver %>% dplyr::select(Effect, CYP2B6_Activation, Tramadol, ESRRA_Inhibition, 
                                               Scn1a_Inhibition, ESR1_Activation, Oprk1_Binding,
                                               OP_Fire_Retardants, "ADORA2_ADRA_HRH-Inhibition",
                                               PGR_Inhibition, NR1I2_Binding, NR1I2_Activation,
                                               PPARG_Activation, "Tris(2-butoxyethyl)phosphate") 
set.seed(3)
cf.final <- cforest(Effect~., data=UGT_liver_final, control=cforest_unbiased(mtry=4, ntree=5000))

# OOB predictions
UGT_liver_final.pred.oob <- data.frame(predict(cf.final, OOB=T)) %>% rename(Predicted_OOB=Effect)
UGT_liver_final_wPred <- cbind(UGT_liver_final, UGT_liver_final.pred.oob)

# OOB error
mean(abs(UGT_liver_final$Effect-UGT_liver_final.pred.oob$Predicted_OOB)) #1.233759
sqrt(mean((UGT_liver_final$Effect-UGT_liver_final.pred.oob$Predicted_OOB)^2)) #1.817499

# Predicted (OOB) vs. observed
plot(Predicted_OOB~Effect, data=UGT_liver_final_wPred, pch=16)
abline(a=0, b=1)

#UGT1A1_intestine####
names(UGT_intestine)
UGT_intestine_sub <- UGT_intestine[,c(2:111)]
UGT_intestine_sub[,1] <- round(UGT_intestine_sub[,1],2)
View(UGT_intestine_sub)

p <- ncol(UGT_intestine_sub) - 1

#cforest - all variables
set.seed(3)
cf.all <- cforest(Effect~., data =UGT_intestine_sub, control = cforest_unbiased(mtry =ceiling(p/3), ntree = 5000))
cf.pred <- predict(cf.all, OOB = T)
sqrt(mean((UGT_intestine_sub$Effect-cf.pred)^2)) #5.986977
mean(abs(UGT_intestine_sub$Effect - cf.pred)) #3.721449

#plot predicted vs. observed
plot(y=cf.pred, x= UGT_intestine_sub$Effect)
abline(a=0, b=1)

# Unconditional variable importance (threshold=1)
set.seed(3) 
cf.upi <- permimp(cf.all, conditional=TRUE, threshold=1, thresholdDiagnostics = TRUE)
as.matrix(sort(cf.upi$values, decreasing=T)) # CF, unconditional

# Conditional variable importance (threshold=0.9)
set.seed(3) 
cf.cpi <- permimp(cf.all, conditional=TRUE, threshold=0.90, thresholdDiagnostics = TRUE)
as.matrix(sort(cf.cpi$values, decreasing=T)) # CF, conditional

# cf.cpi$perTree # VI values per tree
plot(cf.cpi, type="bar")

### Loop to recursively eliminate least informative variables ##

# Pull out response and predictors - columns 5, 7-53
UGT_intestine_preds <- UGT_intestine[, c(2:111)] # original predictors
UGT_intestine_preds[,1] <- round(UGT_intestine_preds[,1], 2)
UGT_intestine_preds_loop <- UGT_intestine_preds 

p <- ncol(UGT_intestine_preds)-1

# Save OOB error, variable importance
VI.list <- NULL
var.list <- rep(NA,p)

OOB.rmse <- rep(NA,p)
OOB.mae <- rep(NA,p)

train.rmse <- rep(NA,p)
train.mae <- rep(NA,p)

#run RF - conditional
set.seed(3) # Set seed for reproducibility

for(i in 1:p){
  print(i)
  
  nump <- ncol(UGT_intestine_preds_loop)-1 # Number predictors
  
  set.seed(3) 
  cf <- cforest(Effect~., data=UGT_intestine_preds_loop, control=cforest_unbiased(mtry=ceiling(nump/3), ntree=5000))
  cf.pred <- predict(cf, OOB=T)
  cf.pred.tr <- predict(cf, OOB=F)
  
  # OOB error
  OOB.rmse[i] <- sqrt(mean((UGT_intestine_preds$Effect-cf.pred)^2)) # OOB rmse
  OOB.mae[i] <- mean(abs(UGT_intestine_preds$Effect-cf.pred)) # oob mae
  
  # Training error
  train.rmse[i] <- sqrt(mean((UGT_intestine_preds$Effect-cf.pred.tr)^2)) # training rmse
  train.mae[i] <- mean(abs(UGT_intestine_preds$Effect-cf.pred.tr)) #  training mae
  
  
  # Variable importance
  set.seed(3) 
  cf.cpi <- permimp(cf, conditional=TRUE, threshold=0.9) # set threshold to 1 for unconditional
  
  VI.list[[i]] <- as.matrix(sort(cf.cpi$values, decreasing=F))
  var.list[i] <- row.names(VI.list[[i]])[1] # Least informative variable to remove
  
  
  UGT_intestine_preds_loop <- UGT_intestine_preds_loop[,!colnames(UGT_intestine_preds_loop)==var.list[i]] # remove variable
  
}
beep()

# Plot OOB error - minimum is the best subset of predictors
plot(OOB.rmse, xlab="Iteration", type='l')
points(OOB.rmse, pch=16)
plot(OOB.mae, xlab="Iteration", type='l')
points(OOB.mae, pch=16)

# training error (not to be used for optimizing number predictors)
plot(train.rmse, xlab="Iteration", type='l')
points(train.rmse, pch=16)

which(OOB.rmse==min(OOB.rmse)) #96
which(OOB.mae==min(OOB.mae)) #22

var.list[which(OOB.mae==min(OOB.mae)):length(OOB.mae)]
# [1] "RXRB_Activation"                        "Cyp2a2_Inhibition"                     
# [3] "JUN_Activation"                         "17alpha-hydroxypregnenolone_Inhibition"
# [5] "NR1I3_Activation"                       "Meprobamate"                           
# [7] "POU2F1_Activation"                      "Carbazole"                             
# [9] "ESR2_Activation"                        "Nicotine"                              
# [11] "Bisphenol A"                            "CREB3_Activation"                      
# [13] "Tributyl phosphate"                     "NR3C1_Binding"                         
# [15] "Phenanthrene"                           "NR1I3_Binding"                         
# [17] "Anthraquinone"                          "Sterols"                               
# [19] "NR1I2_Binding"                          "androstenedione_Inhibition"            
# [21] "testosterone_Inhibition"                "estrone_Activation"                    
# [23] "AR_Activation"                          "4-tert-Octylphenol monoethoxylate"     
# [25] "SCN-Inhibition"                         "5-Methyl-1H-benzotriazole"             
# [27] "CYP2C19_Inhibition"                     "CYP1A2_Inhibition"                     
# [29] "PAHs"                                   "VDR_Activation"                        
# [31] "3,4-Dichlorophenyl isocyanate"          "Pyrene"                                
# [33] "beta-Sitosterol"                        "PGR_Binding"                           
# [35] "NR1H4_Activation"                       "TOX21_CAR_Antagonist"                  
# [37] "PPARG_Activation"                       "Industrial_PPCP_Narcotics"             
# [39] "PGR_Inhibition"                         "CYP1A2_Activation"                     
# [41] "Metolachlor"                            "ADORA2A_Binding"                       
# [43] "OPRM1_Binding"                          "Tris(2-butoxyethyl)phosphate"          
# [45] "Benzo[a]pyrene"                         "HDAC1_Inhibition"                      
# [47] "NR3C1_Activation"                       "1-Methylnaphthalene"                   
# [49] "PPCPs_Pesticides_Plasticizers"          "ESR1_Activation"                       
# [51] "HSF1_Activation"                        "Fuels_Industrial_WWIs_Pesticides"      
# [53] "Carbaryl"                               "HIF1A_Activation"                      
# [55] "Acetaminophen"                          "ESRRA_Inhibition"                      
# [57] "Oprk1_Binding"                          "Antibiotics"                           
# [59] "Methotrexate"                           "NR1I2_Activation"                      
# [61] "NFE2L2_Activation"                      "CYP24A1_Activation"                    
# [63] "TP53_Activation"                        "ESRRA_Activation"                      
# [65] "Tris(2-chloroethyl)phosphate"           "Indole"                                
# [67] "CYP2B6_Activation"                      "NR1H4_Inhibition"                      
# [69] "Cholesterol"                            "ATF6_Activation"                       
# [71] "Triamterene"                            "TSPO_Inhibition"                       
# [73] "Tramadol"                               "RXRA_Activation"                       
# [75] "Sulfamethoxazole"                       "11-deoxycortisol_Inhibition"           
# [77] "N,N-diethyl-meta-toluamide"             "PPARD_Activation"                      
# [79] "Fluoranthene"                           "Isophorone"                            
# [81] "SNRIs_WWI"                              "Carbamazepine"                         
# [83] "SLC-Inhibition"                         "Fexofenadine"                          
# [85] "Desvenlafaxine"                         "Lidocaine"                             
# [87] "Venlafaxine"                            "Metformin"  

#double check with rounded error#
OOB.rmse_1 <- round(OOB.rmse, 2)
OOB.mae_1 <- round(OOB.mae, 2)
which(OOB.rmse_1==min(OOB.rmse_1)) #96
which(OOB.mae_1==min(OOB.mae_1)) #22 - 120 as well, which was used because it reduced mae in final model

var.list[which(OOB.mae_1==min(OOB.mae_1)):length(OOB.mae_1)] 

# Save output - name files in a way that lets you know what differs between different runs of the algorithm

saveRDS(VI.list, "cf_VI_at_iter_thresh90_ugt_int.rds") # R object (list) with the variable importances at each iteration of the algorithm

RFE_info <- data.frame(Worst_Var=var.list, OOB_rmse=OOB.rmse, OOB_mae=OOB.mae, Train_rmse=train.rmse, Train_mae=train.mae)

write.csv(RFE_info, "cf_RFE_info_thresh90_ugt_int.csv", row.names = FALSE)

#set mtry = # of predictors left in model / 3 (nump/3) - and round it up to solid number 
# Refit final model
UGT_intestine_final <- UGT_intestine %>% dplyr::select(Effect,"Metformin", "Venlafaxine", "Lidocaine", "Desvenlafaxine", "Fexofenadine", "SLC-Inhibition", "Carbamazepine", "SNRIs_WWI", "Isophorone", "Fluoranthene", "PPARD_Activation", "N,N-diethyl-meta-toluamide", "11-deoxycortisol_Inhibition", "Sulfamethoxazole", "RXRA_Activation", "Tramadol", "TSPO_Inhibition", "Triamterene", "ATF6_Activation", "Cholesterol", "NR1H4_Inhibition", "CYP2B6_Activation", "Indole", "Tris(2-chloroethyl)phosphate", "ESRRA_Activation", "TP53_Activation", "CYP24A1_Activation", "NFE2L2_Activation", "NR1I2_Activation", "Methotrexate", "Antibiotics", "Oprk1_Binding", "ESRRA_Inhibition", "Acetaminophen", "HIF1A_Activation", "Carbaryl", "Fuels_Industrial_WWIs_Pesticides", "HSF1_Activation", "ESR1_Activation", "PPCPs_Pesticides_Plasticizers", "1-Methylnaphthalene", "NR3C1_Activation", "HDAC1_Inhibition", "Benzo[a]pyrene", "Tris(2-butoxyethyl)phosphate", "OPRM1_Binding", "ADORA2A_Binding", "Metolachlor", "CYP1A2_Activation", "PGR_Inhibition", "Industrial_PPCP_Narcotics", "PPARG_Activation", "TOX21_CAR_Antagonist", "NR1H4_Activation", "PGR_Binding", "beta-Sitosterol", "Pyrene", "3,4-Dichlorophenyl isocyanate", "VDR_Activation", "PAHs", "CYP1A2_Inhibition", "CYP2C19_Inhibition", "5-Methyl-1H-benzotriazole", "SCN-Inhibition", "4-tert-Octylphenol monoethoxylate", "AR_Activation", "estrone_Activation", "testosterone_Inhibition", "androstenedione_Inhibition", "NR1I2_Binding", "Sterols", "Anthraquinone", "NR1I3_Binding", "Phenanthrene", "NR3C1_Binding", "Tributyl phosphate", "CREB3_Activation", "Bisphenol A", "Nicotine", "ESR2_Activation", "Carbazole", "POU2F1_Activation", "Meprobamate", "NR1I3_Activation", "17alpha-hydroxypregnenolone_Inhibition", "JUN_Activation", "Cyp2a2_Inhibition", "RXRB_Activation"
)

set.seed(3)
cf.final <- cforest(Effect~., data=UGT_intestine_final, control=cforest_unbiased(mtry=(29), ntree=5000))

# OOB predictions
UGT_intestine_final.pred.oob <- data.frame(predict(cf.final, OOB=T)) %>% rename(Predicted_OOB=Effect)
UGT_intestine_final_wPred <- cbind(UGT_intestine_final, UGT_intestine_final.pred.oob)

# OOB error
mean(abs(UGT_intestine_final$Effect-UGT_intestine_final.pred.oob$Predicted_OOB)) #3.736375
sqrt(mean((UGT_intestine_final$Effect-UGT_intestine_final.pred.oob$Predicted_OOB)^2)) #5.998585

# Predicted (OOB) vs. observed
plot(Predicted_OOB~Effect, data=UGT_intestine_final_wPred, pch=16)
abline(a=0, b=1)


#CYP2N13_intestine####
names(CYP2N13)
CYP2N13_sub <- CYP2N13[,c(2:111)]
CYP2N13_sub[,1] <- round(CYP2N13_sub[,1],2)

p <- ncol(CYP2N13_sub) - 1

#cforest - all variables
set.seed(3)
cf.all <- cforest(Effect~., data =CYP2N13_sub, control = cforest_unbiased(mtry =ceiling(p/3), ntree = 5000))
cf.pred <- predict(cf.all, OOB = T)
sqrt(mean((CYP2N13_sub$Effect-cf.pred)^2)) #10.15187
mean(abs(CYP2N13_sub$Effect - cf.pred)) #5.725558

#plot predicted vs. observed
plot(y=cf.pred, x= CYP2N13_sub$Effect)
abline(a=0, b=1)

# Unconditional variable importance (threshold=1)
set.seed(3) 
cf.upi <- permimp(cf.all, conditional=TRUE, threshold=1, thresholdDiagnostics = TRUE)
as.matrix(sort(cf.upi$values, decreasing=T)) # CF, unconditional

# Conditional variable importance (threshold=0.9)
set.seed(3) 
cf.cpi <- permimp(cf.all, conditional=TRUE, threshold=0.90, thresholdDiagnostics = TRUE)
as.matrix(sort(cf.cpi$values, decreasing=T)) # CF, conditional

# cf.cpi$perTree # VI values per tree
plot(cf.cpi, type="bar")

### Loop to recursively eliminate least informative variables ##

# Pull out response and predictors - columns 5, 7-53
CYP2N13_preds <- CYP2N13[, c(2:111)] # original predictors
CYP2N13_preds[,1] <- round(CYP2N13_preds[,1], 2)
CYP2N13_preds_loop <- CYP2N13_preds 

p <- ncol(CYP2N13_preds)-1

# Save OOB error, variable importance
VI.list <- NULL
var.list <- rep(NA,p)

OOB.rmse <- rep(NA,p)
OOB.mae <- rep(NA,p)

train.rmse <- rep(NA,p)
train.mae <- rep(NA,p)

#run RF - conditional
set.seed(3) # Set seed for reproducibility

for(i in 1:p){
  print(i)
  
  nump <- ncol(CYP2N13_preds_loop)-1 # Number predictors
  
  set.seed(3) 
  cf <- cforest(Effect~., data=CYP2N13_preds_loop, control=cforest_unbiased(mtry=ceiling(nump/3), ntree=5000))
  cf.pred <- predict(cf, OOB=T)
  cf.pred.tr <- predict(cf, OOB=F)
  
  # OOB error
  OOB.rmse[i] <- sqrt(mean((CYP2N13_preds$Effect-cf.pred)^2)) # OOB rmse
  OOB.mae[i] <- mean(abs(CYP2N13_preds$Effect-cf.pred)) # oob mae
  
  # Training error
  train.rmse[i] <- sqrt(mean((CYP2N13_preds$Effect-cf.pred.tr)^2)) # training rmse
  train.mae[i] <- mean(abs(CYP2N13_preds$Effect-cf.pred.tr)) #  training mae
  
  
  # Variable importance
  set.seed(3) 
  cf.cpi <- permimp(cf, conditional=TRUE, threshold=0.9) # set threshold to 1 for unconditional
  
  VI.list[[i]] <- as.matrix(sort(cf.cpi$values, decreasing=F))
  var.list[i] <- row.names(VI.list[[i]])[1] # Least informative variable to remove
  
  
  CYP2N13_preds_loop <- CYP2N13_preds_loop[,!colnames(CYP2N13_preds_loop)==var.list[i]] # remove variable
  
}
beep()

# Plot OOB error - minimum is the best subset of predictors
plot(OOB.rmse, xlab="Iteration", type='l')
points(OOB.rmse, pch=16)
plot(OOB.mae, xlab="Iteration", type='l')
points(OOB.mae, pch=16)

# training error (not to be used for optimizing number predictors)
plot(train.rmse, xlab="Iteration", type='l')
points(train.rmse, pch=16)

which(OOB.rmse==min(OOB.rmse)) #99
which(OOB.mae==min(OOB.mae)) #86

var.list[which(OOB.mae==min(OOB.mae)):length(OOB.mae)]
# [1] "PGR_Inhibition"                         "SCN-Inhibition"                         "PAHs"                                  
# [4] "testosterone_Inhibition"                "NR1H4_Activation"                       "TSPO_Inhibition"                       
# [7] "NR1I2_Activation"                       "androstenedione_Inhibition"             "Fexofenadine"                          
# [10] "NR1I3_Binding"                          "NFE2L2_Activation"                      "Tramadol"                              
# [13] "17alpha-hydroxypregnenolone_Inhibition" "Metformin"                              "Carbamazepine"                         
# [16] "VDR_Activation"                         "ESRRA_Inhibition"                       "Lidocaine"                             
# [19] "RXRA_Activation"                        "ESR1_Activation"                        "Venlafaxine"                           
# [22] "Desvenlafaxine"                         "SLC-Inhibition"                         "SNRIs_WWI"
#double check with rounded error#
OOB.rmse_1 <- round(OOB.rmse, 2)
OOB.mae_1 <- round(OOB.mae, 2)
which(OOB.rmse_1==min(OOB.rmse_1)) #99
which(OOB.mae_1==min(OOB.mae_1)) #86

var.list[which(OOB.mae_1==min(OOB.mae_1)):length(OOB.mae_1)] #same

# Save output - name files in a way that lets you know what differs between different runs of the algorithm

saveRDS(VI.list, "cf_VI_at_iter_thresh90_cyp2n13.rds") # R object (list) with the variable importances at each iteration of the algorithm

RFE_info <- data.frame(Worst_Var=var.list, OOB_rmse=OOB.rmse, OOB_mae=OOB.mae, Train_rmse=train.rmse, Train_mae=train.mae)

write.csv(RFE_info, "cf_RFE_info_thresh90_cyp2n13.csv", row.names = FALSE)

#set mtry = # of predictors left in model / 3 (nump/3) - and round it up to solid number 
# Refit final model
CYP2N13_final <- CYP2N13 %>% dplyr::select(Effect, "PGR_Inhibition" , "SCN-Inhibition" , "PAHs" , "testosterone_Inhibition" , "NR1H4_Activation" , "TSPO_Inhibition" , "NR1I2_Activation" , "androstenedione_Inhibition" , "Fexofenadine" , "NR1I3_Binding" , "NFE2L2_Activation" , "Tramadol" , "17alpha-hydroxypregnenolone_Inhibition" , "Metformin" , "Carbamazepine" , "VDR_Activation" , "ESRRA_Inhibition" , "Lidocaine" , "RXRA_Activation" , "ESR1_Activation" , "Venlafaxine" , "Desvenlafaxine" , "SLC-Inhibition" , "SNRIs_WWI" 
)
set.seed(3)
cf.final <- cforest(Effect~., data=CYP2N13_final, control=cforest_unbiased(mtry=8, ntree=5000))

# OOB predictions
CYP2N13_final.pred.oob <- data.frame(predict(cf.final, OOB=T)) %>% rename(Predicted_OOB=Effect)
CYP2N13_final_wPred <- cbind(CYP2N13_final, CYP2N13_final.pred.oob)

# OOB error
mean(abs(CYP2N13_final$Effect-CYP2N13_final.pred.oob$Predicted_OOB)) #5.723969
sqrt(mean((CYP2N13_final$Effect-CYP2N13_final.pred.oob$Predicted_OOB)^2)) #10.17777

# Predicted (OOB) vs. observed
plot(Predicted_OOB~Effect, data=CYP2N13_final_wPred, pch=16)
abline(a=0, b=1)

#CYP2AD6_intestine####
names(CYP_2AD6)
CYP_2AD6_sub <- CYP_2AD6[,c(2:111)]
CYP_2AD6_sub[,1] <- round(CYP_2AD6_sub[,1],2)
View(CYP_2AD6_sub)

p <- ncol(CYP_2AD6_sub) - 1

#cforest - all variables
set.seed(3)
cf.all <- cforest(Effect~., data =CYP_2AD6_sub, control = cforest_unbiased(mtry =ceiling(p/3), ntree = 5000))
cf.pred <- predict(cf.all, OOB = T)
sqrt(mean((CYP_2AD6_sub$Effect-cf.pred)^2)) #2.557994
mean(abs(CYP_2AD6_sub$Effect - cf.pred)) #1.791248

#plot predicted vs. observed
plot(y=cf.pred, x= CYP_2AD6_sub$Effect)
abline(a=0, b=1)

# Unconditional variable importance (threshold=1)
set.seed(3) 
cf.upi <- permimp(cf.all, conditional=TRUE, threshold=1, thresholdDiagnostics = TRUE)
as.matrix(sort(cf.upi$values, decreasing=T)) # CF, unconditional

# Conditional variable importance (threshold=0.9)
set.seed(3) 
cf.cpi <- permimp(cf.all, conditional=TRUE, threshold=0.90, thresholdDiagnostics = TRUE)
as.matrix(sort(cf.cpi$values, decreasing=T)) # CF, conditional

# cf.cpi$perTree # VI values per tree
plot(cf.cpi, type="bar")

### Loop to recursively eliminate least informative variables ##

# Pull out response and predictors - columns 5, 7-53
CYP_2AD6_preds <- CYP_2AD6[, c(2:111)] # original predictors
CYP_2AD6_preds[,1] <- round(CYP_2AD6_preds[,1], 2)
CYP_2AD6_preds_loop <- CYP_2AD6_preds 
View(CYP_2AD6_preds_loop)

p <- ncol(CYP_2AD6_preds)-1

# Save OOB error, variable importance
VI.list <- NULL
var.list <- rep(NA,p)

OOB.rmse <- rep(NA,p)
OOB.mae <- rep(NA,p)

train.rmse <- rep(NA,p)
train.mae <- rep(NA,p)

#run RF - conditional
set.seed(3) # Set seed for reproducibility

for(i in 1:p){
  print(i)
  
  nump <- ncol(CYP_2AD6_preds_loop)-1 # Number predictors
  
  set.seed(3) 
  cf <- cforest(Effect~., data=CYP_2AD6_preds_loop, control=cforest_unbiased(mtry=ceiling(nump/3), ntree=5000))
  cf.pred <- predict(cf, OOB=T)
  cf.pred.tr <- predict(cf, OOB=F)
  
  # OOB error
  OOB.rmse[i] <- sqrt(mean((CYP_2AD6_preds$Effect-cf.pred)^2)) # OOB rmse
  OOB.mae[i] <- mean(abs(CYP_2AD6_preds$Effect-cf.pred)) # oob mae
  
  # Training error
  train.rmse[i] <- sqrt(mean((CYP_2AD6_preds$Effect-cf.pred.tr)^2)) # training rmse
  train.mae[i] <- mean(abs(CYP_2AD6_preds$Effect-cf.pred.tr)) #  training mae
  
  
  # Variable importance
  set.seed(3) 
  cf.cpi <- permimp(cf, conditional=TRUE, threshold=0.9) # set threshold to 1 for unconditional
  
  VI.list[[i]] <- as.matrix(sort(cf.cpi$values, decreasing=F))
  var.list[i] <- row.names(VI.list[[i]])[1] # Least informative variable to remove
  
  
  CYP_2AD6_preds_loop <- CYP_2AD6_preds_loop[,!colnames(CYP_2AD6_preds_loop)==var.list[i]] # remove variable
  
}
beep()

# Plot OOB error - minimum is the best subset of predictors
plot(OOB.rmse, xlab="Iteration", type='l')
points(OOB.rmse, pch=16)
plot(OOB.mae, xlab="Iteration", type='l')
points(OOB.mae, pch=16)

# training error (not to be used for optimizing number predictors)
plot(train.rmse, xlab="Iteration", type='l')
points(train.rmse, pch=16)

which(OOB.rmse==min(OOB.rmse)) #99
which(OOB.mae==min(OOB.mae)) #90

var.list[which(OOB.mae==min(OOB.mae)):length(OOB.mae)]

#double check with rounded error#
OOB.rmse_1 <- round(OOB.rmse, 2)
OOB.mae_1 <- round(OOB.mae, 2)
which(OOB.rmse_1==min(OOB.rmse_1)) #102
which(OOB.mae_1==min(OOB.mae_1)) #96

var.list[which(OOB.mae_1==min(OOB.mae_1)):length(OOB.mae_1)]

# Save output - name files in a way that lets you know what differs between different runs of the algorithm

saveRDS(VI.list, "cf_VI_at_iter_thresh90_cyp2ad6.rds") # R object (list) with the variable importances at each iteration of the algorithm

RFE_info <- data.frame(Worst_Var=var.list, OOB_rmse=OOB.rmse, OOB_mae=OOB.mae, Train_rmse=train.rmse, Train_mae=train.mae)

write.csv(RFE_info, "cf_RFE_info_thresh90_cyp2ad6.csv", row.names = FALSE)

#set mtry = # of predictors left in model / 3 (nump/3) - and round it up to solid number 
# Refit final model
CYP_2AD6_final <- CYP_2AD6 %>% dplyr::select(Effect, "Desvenlafaxine", "SLC-Inhibition", "SNRIs_WWI", "Carbamazepine", "SCN-Inhibition", "Tramadol", "Lidocaine", "Cholesterol", "RXRA_Activation", "Venlafaxine", "Metformin", "Fexofenadine", "CYP2B6_Activation", "N,N-diethyl-meta-toluamide"
)
set.seed(3)
cf.final <- cforest(Effect~., data=CYP_2AD6_final, control=cforest_unbiased(mtry=5, ntree=5000))

# OOB predictions
CYP_2AD6_final.pred.oob <- data.frame(predict(cf.final, OOB=T)) %>% rename(Predicted_OOB=Effect)
CYP_2AD6_final_wPred <- cbind(CYP_2AD6_final, CYP_2AD6_final.pred.oob)

# OOB error
mean(abs(CYP_2AD6_final$Effect-CYP_2AD6_final.pred.oob$Predicted_OOB)) #1.792987
sqrt(mean((CYP_2AD6_final$Effect-CYP_2AD6_final.pred.oob$Predicted_OOB)^2)) #2.552183

# Predicted (OOB) vs. observed
plot(Predicted_OOB~Effect, data=CYP_2AD6_final_wPred, pch=16)
abline(a=0, b=1)

#CYP1A1_2018####
names(cyp_1a1_2018)
cyp_1a1_2018_sub <- cyp_1a1_2018[,c(2:213)]
cyp_1a1_2018_sub[,1] <- round(cyp_1a1_2018_sub[,1],2)

p <- ncol(cyp_1a1_2018_sub) - 1

#cforest - all variables
set.seed(3)
cf.all <- cforest(Effect~., data =cyp_1a1_2018_sub, control = cforest_unbiased(mtry =ceiling(p/3), ntree = 5000))
cf.pred <- predict(cf.all, OOB = T)
sqrt(mean((cyp_1a1_2018_sub$Effect-cf.pred)^2)) #3.238356
mean(abs(cyp_1a1_2018_sub$Effect - cf.pred)) #1.777179

#plot predicted vs. observed
plot(y=cf.pred, x= cyp_1a1_2018_sub$Effect)
abline(a=0, b=1)

# Unconditional variable importance (threshold=1)
set.seed(3) 
cf.upi <- permimp(cf.all, conditional=TRUE, threshold=1, thresholdDiagnostics = TRUE)
as.matrix(sort(cf.upi$values, decreasing=T)) # CF, unconditional

# Conditional variable importance (threshold=0.9)
set.seed(3) 
cf.cpi <- permimp(cf.all, conditional=TRUE, threshold=0.90, thresholdDiagnostics = TRUE)
as.matrix(sort(cf.cpi$values, decreasing=T)) # CF, conditional

# cf.cpi$perTree # VI values per tree
plot(cf.cpi, type="bar")

### Loop to recursively eliminate least informative variables ##

# Pull out response and predictors - columns 5, 7-53
cyp_1a1_2018_preds <- cyp_1a1_2018[, c(2:213)] # original predictors
cyp_1a1_2018_preds[,1] <- round(cyp_1a1_2018_preds[,1], 2)
cyp_1a1_2018_preds_loop <- cyp_1a1_2018_preds 

p <- ncol(cyp_1a1_2018_preds)-1

# Save OOB error, variable importance
VI.list <- NULL
var.list <- rep(NA,p)

OOB.rmse <- rep(NA,p)
OOB.mae <- rep(NA,p)

train.rmse <- rep(NA,p)
train.mae <- rep(NA,p)

#run RF - conditional
set.seed(3) # Set seed for reproducibility

for(i in 1:p){
  print(i)
  
  nump <- ncol(cyp_1a1_2018_preds_loop)-1 # Number predictors
  
  set.seed(3) 
  cf <- cforest(Effect~., data=cyp_1a1_2018_preds_loop, control=cforest_unbiased(mtry=ceiling(nump/3), ntree=5000))
  cf.pred <- predict(cf, OOB=T)
  cf.pred.tr <- predict(cf, OOB=F)
  
  # OOB error
  OOB.rmse[i] <- sqrt(mean((cyp_1a1_2018_preds$Effect-cf.pred)^2)) # OOB rmse
  OOB.mae[i] <- mean(abs(cyp_1a1_2018_preds$Effect-cf.pred)) # oob mae
  
  # Training error
  train.rmse[i] <- sqrt(mean((cyp_1a1_2018_preds$Effect-cf.pred.tr)^2)) # training rmse
  train.mae[i] <- mean(abs(cyp_1a1_2018_preds$Effect-cf.pred.tr)) #  training mae
  
  
  # Variable importance
  set.seed(3) 
  cf.cpi <- permimp(cf, conditional=TRUE, threshold=0.9) # set threshold to 1 for unconditional
  
  VI.list[[i]] <- as.matrix(sort(cf.cpi$values, decreasing=F))
  var.list[i] <- row.names(VI.list[[i]])[1] # Least informative variable to remove
  
  
  cyp_1a1_2018_preds_loop <- cyp_1a1_2018_preds_loop[,!colnames(cyp_1a1_2018_preds_loop)==var.list[i]] # remove variable
  
}
beep()

# Plot OOB error - minimum is the best subset of predictors
plot(OOB.rmse, xlab="Iteration", type='l')
points(OOB.rmse, pch=16)
plot(OOB.mae, xlab="Iteration", type='l')
points(OOB.mae, pch=16)

# training error (not to be used for optimizing number predictors)
plot(train.rmse, xlab="Iteration", type='l')
points(train.rmse, pch=16)

which(OOB.rmse==min(OOB.rmse)) #143
which(OOB.mae==min(OOB.mae)) #155

var.list[which(OOB.mae==min(OOB.mae)):length(OOB.mae)]
# 
# [1] "FOS_Activation"                 "ESR2_Activation"                "Venlafaxine"                    "Bupropion"                      "NR1I2_Activation"               "Atrazine_Famotidine_Ranitidine"
# [7] "ESRRA_Inhibition"               "Fluoranthene"                   "NFE2L2_Activation"              "GLI3_Inhibition"                "Dextromethorphan"               "TP53_Activation"               
# [13] "CYP1A2_Activation"              "3-beta-Coprostanol"             "Fexofenadine"                   "PTGS2_Inhibition"               "PGR_Inhibition"                 "Ranitidine"                    
# [19] "Nicotine"                       "CYP2B6_Activation"              "RXRA_Activation"                "CYP1A2_Inhibition"              "Benzo[a]pyrene"                 "Sulfadimethoxine"              
# [25] "OPRM1_Binding"                  "Metformin"                      "Atrazine"                       "Metoprolol"                     "N,N-diethyl-meta-toluamide"     "ESR2_Inhibition"               
# [31] "beta-Sitosterol"                "RXRB_Activation"                "NR1H4_Activation"               "Metolachlor"                    "PPCPs_Pesticides_Plasticizers"  "NR1I3_Activation"              
# [37] "Antibiotics"                    "1-Methylnaphthalene"            "Anthraquinone"                  "Sterols"                        "Maoa_Inhibition"                "AR_Inhibition"                 
# [43] "JUN_Activation"                 "ESR1_Activation"                "ADORA2_ADRA_HRH-Inhibition"     "Acetaminophen"                  "FOS|JUN_Activation"             "Industrial_PPCP_Narcotics"     
# [49] "Phenanthrene"                   "Cotinine"                       "AR_Activation"                  "Cholesterol"                    "HDAC1_Inhibition"               "Scn1a_Inhibition"              
# [55] "Caffeine"                       "Tributyl phosphate"             "Methyl-1H-benzotriazole"  

#double check with rounded error#
OOB.rmse_1 <- round(OOB.rmse, 2)
OOB.mae_1 <- round(OOB.mae, 2)
which(OOB.rmse_1==min(OOB.rmse_1)) #206
which(OOB.mae_1==min(OOB.mae_1)) #181

var.list[which(OOB.mae_1==min(OOB.mae_1)):length(OOB.mae_1)]
# Cholesterol
# ADORA2A_Binding
# TRPV1-Activation
# Narcotic_PPCPs_CNS_Stimulants_Antivirals
# ADRA2A_Binding
# Tributyl phosphate
# Methyl-1H-benzotriazole


# Save output - name files in a way that lets you know what differs between different runs of the algorithm

saveRDS(VI.list, "cf_VI_at_iter_thresh90_cyp1a1_2018.rds") # R object (list) with the variable importances at each iteration of the algorithm

RFE_info <- data.frame(Worst_Var=var.list, OOB_rmse=OOB.rmse, OOB_mae=OOB.mae, Train_rmse=train.rmse, Train_mae=train.mae)

write.csv(RFE_info, "cf_RFE_info_thresh90_cyp1a1_2018.csv", row.names = FALSE)

#set mtry = # of predictors left in model / 3 (nump/3) - and round it up to solid number 
# Refit final model
cyp_1a1_2018_final <- cyp_1a1_2018 %>% dplyr::select(Effect, "Methyl-1H-benzotriazole", "Tributyl phosphate", "Caffeine", "Scn1a_Inhibition", "HDAC1_Inhibition", "Cholesterol", "AR_Activation", "Cotinine", "Phenanthrene", "Industrial_PPCP_Narcotics", "FOS|JUN_Activation", "Acetaminophen", "ADORA2_ADRA_HRH-Inhibition", "ESR1_Activation", "JUN_Activation", "AR_Inhibition", "Maoa_Inhibition", "Sterols", "Anthraquinone", "1-Methylnaphthalene", "Antibiotics", "NR1I3_Activation", "PPCPs_Pesticides_Plasticizers", "Metolachlor", "NR1H4_Activation", "RXRB_Activation", "beta-Sitosterol", "ESR2_Inhibition", "N,N-diethyl-meta-toluamide", "Metoprolol", "Atrazine"
)
set.seed(3)
cf.final <- cforest(Effect~., data=cyp_1a1_2018_final, control=cforest_unbiased(mtry=10, ntree=5000))

# OOB predictions
cyp_1a1_2018_final.pred.oob <- data.frame(predict(cf.final, OOB=T)) %>% rename(Predicted_OOB=Effect)
cyp_1a1_2018_final_wPred <- cbind(cyp_1a1_2018_final, cyp_1a1_2018_final.pred.oob)

# OOB error
mean(abs(cyp_1a1_2018_final$Effect-cyp_1a1_2018_final.pred.oob$Predicted_OOB)) #1.763817
sqrt(mean((cyp_1a1_2018_final$Effect-cyp_1a1_2018_final.pred.oob$Predicted_OOB)^2)) #3.233082

# Predicted (OOB) vs. observed
plot(Predicted_OOB~Effect, data=cyp_1a1_2018_final_wPred, pch=16)
abline(a=0, b=1)

#partial dependence plots####

#in vitro XM related#### 
ahr_final <- ahr %>% dplyr::select(Effect, Acetaminophen,
                                   AR_Activation, androstenedione_Inhibition, Anthraquinone,
                                   Cotinine, Caffeine)

set.seed(3)
ahr.final <- cforest(Effect~., data=ahr_final, control=cforest_unbiased(mtry=2, ntree=5000))

ahr_pd_1 <- partial(ahr.final, pred.var="Acetaminophen", grid.resolution = 10)
ahr_pd_4 <- partial(ahr.final, pred.var="AR_Activation", grid.resolution = 10)
ahr_pd_5 <- partial(ahr.final, pred.var="androstenedione_Inhibition", grid.resolution = 10)
ahr_pd_7 <- partial(ahr.final, pred.var="Anthraquinone", grid.resolution = 10)
ahr_pd_8 <- partial(ahr.final, pred.var="Cotinine", grid.resolution = 10)
ahr_pd_9 <- partial(ahr.final, pred.var="Caffeine", grid.resolution = 10)


ahr_pdp_1 <- autoplot(ahr_pd_1, size=1.2) + theme_classic() + xlab("Acetaminophen") + ylab("Effect") +
  theme(text=element_text(size=20)) + ggtitle("AhR")
ahr_pdp_4 <- autoplot(ahr_pd_4, size=1.2) + theme_classic() + xlab("AR_Activation") + ylab("Effect") +
  theme(text=element_text(size=20)) + ggtitle("AhR")
ahr_pdp_5 <- autoplot(ahr_pd_5, size=1.2) + theme_classic() + xlab("androstenedione_Inhibition") + ylab("Effect") +
  theme(text=element_text(size=20)) + ggtitle("AhR")
ahr_pdp_7 <- autoplot(ahr_pd_7, size=1.2) + theme_classic() + xlab("Anthraquinone") + ylab("Effect") +
  theme(text=element_text(size=20)) + ggtitle("AhR")
ahr_pdp_8 <- autoplot(ahr_pd_8, size=1.2) + theme_classic() + xlab("Cotinine") + ylab("Effect") +
  theme(text=element_text(size=20)) + ggtitle("AhR")
ahr_pdp_9 <- autoplot(ahr_pd_9, size=1.2) + theme_classic() + xlab("Caffeine") + ylab("Effect") +
  theme(text=element_text(size=20)) + ggtitle("AhR")


pxr_cis_final <- pxr_cis %>% dplyr::select(Effect,"beta-Sitosterol")
set.seed(3)
pxrcis.final <- cforest(Effect~., data=pxr_cis_final, control=cforest_unbiased(mtry=1, ntree=5000))

pxrcis_pd_1 <- partial(pxrcis.final, pred.var="beta-Sitosterol", grid.resolution = 10)

pxrcis_pdp_1 <- autoplot(pxrcis_pd_1, size=1.2) + theme_classic() + xlab("beta-Sitosterol") + ylab("Effect") +
  theme(text=element_text(size=20)) + ggtitle("PXR_cis")

pxr_trans_final <- pxr_trans %>% dplyr::select(Effect,
                                               "beta-Sitosterol")
set.seed(3)
cf.final <- cforest(Effect~., data=pxr_trans_final, control=cforest_unbiased(mtry=1, ntree=5000))

pxrtrans_pd_1 <- partial(cf.final, pred.var="beta-Sitosterol", grid.resolution = 10)
pxrtrans_pdp_1 <- autoplot(pxrtrans_pd_1, size=1.2) + theme_classic() + xlab("beta-Sitosterol") + ylab("Effect") +
  theme(text=element_text(size=20)) + ggtitle("PXR_trans")


pparg_final <- PPARg %>% dplyr::select(Effect, Anthraquinone, Cotinine)
set.seed(3)
pparg.final <- cforest(Effect~., data=pparg_final, control=cforest_unbiased(mtry=1, ntree=5000))

pparg_pd_1 <- partial(pparg.final, pred.var="Cotinine", grid.resolution = 10)
pparg_pd_9 <- partial(pparg.final, pred.var="Anthraquinone", grid.resolution = 10)

pparg_pdp_1 <- autoplot(pparg_pd_1, size=1.2) + theme_classic() + xlab("Cotinine") + ylab("Effect") +
  theme(text=element_text(size=20)) + ggtitle("PPARg")

pparg_pdp_9 <- autoplot(pparg_pd_9, size=1.2) + theme_classic() + xlab("Anthraquinone") + ylab("Effect") +
  theme(text=element_text(size=20)) + ggtitle("PPARg")


ppara_final <- PPARa %>% dplyr::select(Effect, Camphor, Anthraquinone, "Industrial_PPCP_Narcotics")
set.seed(3)
cf.final <- cforest(Effect~., data=ppara_final, control=cforest_unbiased(mtry=1, ntree=5000))

ppara_pd_1 <- partial(cf.final, pred.var="Camphor", grid.resolution = 10)
ppara_pd_2 <- partial(cf.final, pred.var="Anthraquinone", grid.resolution = 10)
ppara_pd_3 <- partial(cf.final, pred.var="Industrial_PPCP_Narcotics", grid.resolution = 10)

ppara_pdp_1 <- autoplot(ppara_pd_1, size=1.2) + theme_classic() + xlab("Camphor") + ylab("Effect") +
  theme(text=element_text(size=20)) + ggtitle("PPARa")
ppara_pdp_2 <- autoplot(ppara_pd_2, size=1.2) + theme_classic() + xlab("Anthraquinone") + ylab("Effect") +
  theme(text=element_text(size=20)) + ggtitle("PPARa")
ppara_pdp_3 <- autoplot(ppara_pd_3, size=1.2) + theme_classic() + xlab("Industrial_PPCP_Narcotics") + ylab("Effect") +
  theme(text=element_text(size=20)) + ggtitle("PPARa")

are_final <- ARE %>% dplyr::select(Effect, Caffeine)
set.seed(3)
cf.final <- cforest(Effect~., data=are_final, control=cforest_unbiased(mtry=1, ntree=5000))
are_pd_1 <- partial(cf.final, pred.var="Caffeine", grid.resolution = 10)
are_pdp_1 <- autoplot(are_pd_1, size=1.2) + theme_classic() + xlab("Caffeine") + ylab("Effect") +
  theme(text=element_text(size=20)) + ggtitle("ARE")

ppre_final <- PPRE %>% dplyr::select(Effect,AR_Activation,Cotinine)
set.seed(3)
cf.final <- cforest(Effect~., data=ppre_final, control=cforest_unbiased(mtry =1, ntree=5000))
ppre_pd_1 <- partial(cf.final, pred.var="AR_Activation", grid.resolution = 10)
ppre_pd_2 <- partial(cf.final, pred.var="Cotinine", grid.resolution = 10)

ppre_pdp_1 <- autoplot(ppre_pd_1, size=1.2) + theme_classic() + xlab("AR_Activation") + ylab("Effect") +
  theme(text=element_text(size=20)) + ggtitle("PPRE")
ppre_pdp_2 <- autoplot(ppre_pd_2, size=1.2) + theme_classic() + xlab("Cotinine") + ylab("Effect") +
  theme(text=element_text(size=20)) + ggtitle("PPRE")

gr_final <- GR %>% dplyr::select(Effect,Triamterene, Metformin,
                                 PPCPs_Pesticides_Plasticizers)
set.seed(3)
cf.final <- cforest(Effect~., data=gr_final, control=cforest_unbiased(mtry=1, ntree=5000))
gr_pd_1 <- partial(cf.final, pred.var="Triamterene", grid.resolution = 10)
gr_pd_2 <- partial(cf.final, pred.var="Metformin", grid.resolution = 10)
gr_pd_3 <- partial(cf.final, pred.var="PPCPs_Pesticides_Plasticizers", grid.resolution = 10)

gr_pdp_1 <- autoplot(gr_pd_1, size=1.2) + theme_classic() + xlab("Triamterene") + ylab("Effect") +
  theme(text=element_text(size=20)) + ggtitle("GR")
gr_pdp_2 <- autoplot(gr_pd_2, size=1.2) + theme_classic() + xlab("Metformin") + ylab("Effect") +
  theme(text=element_text(size=20)) + ggtitle("GR")

rxrb_final <- RXRb %>% dplyr::select(Effect, Acetaminophen, Anthraquinone, RXRB_Activation, Industrial_PPCP_Narcotics)
set.seed(3)
cf.final <- cforest(Effect~., data=rxrb_final, control=cforest_unbiased(mtry=1, ntree=5000))
rxrb_pd_1 <- partial(cf.final, pred.var="Acetaminophen", grid.resolution = 10)
rxrb_pd_2 <- partial(cf.final, pred.var="Anthraquinone", grid.resolution = 10)
rxrb_pd_3 <- partial(cf.final, pred.var="RXRB_Activation", grid.resolution = 10)
rxrb_pd_4 <- partial(cf.final, pred.var="Industrial_PPCP_Narcotics", grid.resolution = 10)

rxrb_pdp_1 <- autoplot(rxrb_pd_1, size=1.2) + theme_classic() + xlab("Acetaminophen") + ylab("Effect") +
  theme(text=element_text(size=20)) + ggtitle("RXRb")
rxrb_pdp_2 <- autoplot(rxrb_pd_2, size=1.2) + theme_classic() + xlab("Anthraquinone") + ylab("Effect") +
  theme(text=element_text(size=20)) + ggtitle("RXRb")
rxrb_pdp_3 <- autoplot(rxrb_pd_3, size=1.2) + theme_classic() + xlab("RXRB_Activation") + ylab("Effect") +
  theme(text=element_text(size=20)) + ggtitle("RXRb")
rxrb_pdp_4 <- autoplot(rxrb_pd_4, size=1.2) + theme_classic() + xlab("Industrial_PPCP_Narcotics") + ylab("Effect") +
  theme(text=element_text(size=20)) + ggtitle("RXRb")


invitro_XM_related_pdp <- grid.arrange(ahr_pdp_1, ahr_pdp_4, ahr_pdp_5, 
                                      ahr_pdp_7, ahr_pdp_8, ahr_pdp_9, pxrcis_pdp_1,
                                       pxrtrans_pdp_1, pparg_pdp_1, pparg_pdp_9, ppara_pdp_1, ppara_pdp_2, ppara_pdp_3, 
                                      are_pdp_1, ppre_pdp_1, ppre_pdp_2, gr_pdp_1, gr_pdp_2, 
                                      rxrb_pdp_1, rxrb_pdp_2, rxrb_pdp_3, rxrb_pdp_4,
                                      ncol = 5)

ggsave("invitro_XM_related_pdp.jpeg", invitro_XM_related_pdp, height = 30, width = 30)

#in vitro ER-related####
era_trans_final <- ERa_trans %>% dplyr::select(Effect, "Anthraquinone", Industrial_PPCP_Narcotics)
set.seed(3)
cf.final <- cforest(Effect~., data=era_trans_final, control=cforest_unbiased(mtry=1, ntree=5000))

eratrans_pd_1 <- partial(cf.final, pred.var="Industrial_PPCP_Narcotics", grid.resolution = 10)
eratrans_pd_2 <- partial(cf.final, pred.var="Anthraquinone", grid.resolution = 10)

eratrans_pdp_1 <- autoplot(eratrans_pd_1, size=1.2) + theme_classic() + xlab("Industrial_PPCP_Narcotics") + ylab("Effect") +
  theme(text=element_text(size=20)) + ggtitle("ERa")
eratrans_pdp_2 <- autoplot(eratrans_pd_2, size=1.2) + theme_classic() + xlab("Anthraquinone") + ylab("Effect") +
  theme(text=element_text(size=20)) + ggtitle("ERa")


ere_cis_final <- ERE_cis %>% dplyr::select(Effect, "Sterols","beta-Sitosterol", Industrial_PPCP_Narcotics)

set.seed(3)
cf.final <- cforest(Effect~., data=ere_cis_final, control=cforest_unbiased(mtry=1, ntree=5000))

ere_pd_1 <- partial(cf.final, pred.var="Sterols", grid.resolution = 10)
ere_pd_2 <- partial(cf.final, pred.var="Industrial_PPCP_Narcotics", grid.resolution = 10)
ere_pd_3 <- partial(cf.final, pred.var="beta-Sitosterol", grid.resolution = 10)

ere_pdp_1 <- autoplot(ere_pd_1, size=1.2) + theme_classic() + xlab("Sterols") + ylab("Effect") +
  theme(text=element_text(size=20)) + ggtitle("ERE")
ere_pdp_2 <- autoplot(ere_pd_2, size=1.2) + theme_classic() + xlab("Industrial_PPCP_Narcotics") + ylab("Effect") +
  theme(text=element_text(size=20)) + ggtitle("ERE")
ere_pdp_3 <- autoplot(ere_pd_3, size=1.2) + theme_classic() + xlab("beta-Sitosterol") + ylab("Effect") +
  theme(text=element_text(size=20)) + ggtitle("ERE")

t47_final <- t47 %>% dplyr::select(Effect, "Methyl-1H-benzotriazole", Cholesterol, 
                                   Sterols, Dextromethorphan, "beta-Sitosterol", Industrial_PPCP_Narcotics)
set.seed(3)
cf.final <- cforest(Effect~., data=t47_final, control=cforest_unbiased(mtry=2, ntree=5000))

t47_pd_1 <- partial(cf.final, pred.var="Methyl-1H-benzotriazole", grid.resolution = 10)
t47_pd_2 <- partial(cf.final, pred.var="Cholesterol", grid.resolution = 10)
t47_pd_3 <- partial(cf.final, pred.var="Sterols", grid.resolution = 10)
t47_pd_4 <- partial(cf.final, pred.var="Dextromethorphan", grid.resolution = 10)
t47_pd_5 <- partial(cf.final, pred.var="beta-Sitosterol", grid.resolution = 10)
t47_pd_6 <- partial(cf.final, pred.var="Industrial_PPCP_Narcotics", grid.resolution = 10)


t47_pdp_1 <- autoplot(t47_pd_1, size=1.2) + theme_classic() + xlab("Methyl-1H-benzotriazole") + ylab("Effect") +
  theme(text=element_text(size=20)) + ggtitle("E2-EQ")
t47_pdp_2 <- autoplot(t47_pd_2, size=1.2) + theme_classic() + xlab("Cholesterol") + ylab("Effect") +
  theme(text=element_text(size=20)) + ggtitle("E2-EQ")
t47_pdp_3 <- autoplot(t47_pd_3, size=1.2) + theme_classic() + xlab("Sterols") + ylab("Effect") +
  theme(text=element_text(size=20)) + ggtitle("E2-EQ")
t47_pdp_4 <- autoplot(t47_pd_4, size=1.2) + theme_classic() + xlab("Dextromethorphan") + ylab("Effect") +
  theme(text=element_text(size=20)) + ggtitle("E2-EQ")
t47_pdp_5 <- autoplot(t47_pd_5, size=1.2) + theme_classic() + xlab("beta-Sitosterol") + ylab("Effect") +
  theme(text=element_text(size=20)) + ggtitle("E2-EQ")
t47_pdp_6 <- autoplot(t47_pd_6, size=1.2) + theme_classic() + xlab("Industrial_PPCP_Narcotics") + ylab("Effect") +
  theme(text=element_text(size=20)) + ggtitle("E2-EQ")

invitro_ER_pdp <- grid.arrange(eratrans_pdp_1, eratrans_pdp_2,
                               ere_pdp_1, ere_pdp_2, ere_pdp_3, t47_pdp_1, t47_pdp_2, 
                               t47_pdp_3, t47_pdp_4, t47_pdp_5, t47_pdp_6, ncol = 5)

ggsave("invitro_ER_related_pdp.jpeg", invitro_ER_pdp, height = 20, width = 30)

#in vivo ER-related####
ria_final <- ria %>% dplyr::select(Effect, NR1I2_Activation,
                                   Desvenlafaxine, "SLC-Inhibition", Venlafaxine)

set.seed(3)
cf.final <- cforest(Effect~., data=ria_final, control=cforest_unbiased(mtry=ceiling(1), ntree=5000))

ria_pd_1 <- partial(cf.final, pred.var="NR1I2_Activation", grid.resolution = 10)
ria_pd_2 <- partial(cf.final, pred.var="Desvenlafaxine", grid.resolution = 10)
ria_pd_3 <- partial(cf.final, pred.var="SLC-Inhibition", grid.resolution = 10)
ria_pd_4 <- partial(cf.final, pred.var="Venlafaxine", grid.resolution = 10)

ria_pdp_1 <- autoplot(ria_pd_1, size=1.2) + theme_classic() + xlab("NR1I2_Activation") + ylab("Effect") +
  theme(text=element_text(size=20)) + ggtitle("E2-M")
ria_pdp_2 <- autoplot(ria_pd_2, size=1.2) + theme_classic() + xlab("Desvenlafaxine") + ylab("Effect") +
  theme(text=element_text(size=20)) + ggtitle("E2-M")
ria_pdp_3 <- autoplot(ria_pd_3, size=1.2) + theme_classic() + xlab("SLC-Inhibition") + ylab("Effect") +
  theme(text=element_text(size=20)) + ggtitle("E2-M")
ria_pdp_4 <- autoplot(ria_pd_4, size=1.2) + theme_classic() + xlab("Venlafaxine") + ylab("Effect") +
  theme(text=element_text(size=20)) + ggtitle("E2-M")

ria_F_final <- ria_F %>% dplyr::select(Effect, "SCN-Inhibition", "Tramadol", "SLC-Inhibition")
set.seed(3)
cf.final <- cforest(Effect~., data=ria_F_final, control=cforest_unbiased(mtry=1, ntree=5000))

riaf_pd_1 <- partial(cf.final, pred.var="SCN-Inhibition", grid.resolution = 10)
riaf_pd_2 <- partial(cf.final, pred.var="Tramadol", grid.resolution = 10)
riaf_pd_3 <- partial(cf.final, pred.var="SLC-Inhibition", grid.resolution = 10)

riaf_pdp_1 <- autoplot(riaf_pd_1, size=1.2) + theme_classic() + xlab("SCN-Inhibition") + ylab("Effect") +
  theme(text=element_text(size=20)) + ggtitle("E2-F")
riaf_pdp_2 <- autoplot(riaf_pd_1, size=1.2) + theme_classic() + xlab("Tramadol") + ylab("Effect") +
  theme(text=element_text(size=20)) + ggtitle("E2-F")
riaf_pdp_3 <- autoplot(riaf_pd_1, size=1.2) + theme_classic() + xlab("SLC-Inhibition") + ylab("Effect") +
  theme(text=element_text(size=20)) + ggtitle("E2-F")

invivo_ER_pdp <- grid.arrange(ria_pdp_1, ria_pdp_2,ria_pdp_3, ria_pdp_4, 
                              riaf_pdp_1, riaf_pdp_2, riaf_pdp_3, ncol = 5)

ggsave("invivo_ER_pdp.jpeg", invivo_ER_pdp, height = 15, width = 30)


#in vivo XM intestinal####

cyp_1a1_intestine_final <- cyp_1a1_intestine %>% dplyr::select(Effect,
                                                               "PGR_Inhibition","Benzo[a]pyrene","ADORA2_ADRA_HRH-Inhibition", "Pyrene",
                                                               "PPARG_Activation", "androstenedione_Inhibition", "NR1I3_Binding", "NR1I2_Activation",
                                                               "Anthraquinone", "ESR2_Activation", "NR1I3_Activation", "RXRA_Activation", "AR_Inhibition")
set.seed(3)
cf.final <- cforest(Effect~., data=cyp_1a1_intestine_final, control=cforest_unbiased(mtry=4, ntree=5000))

int1a1_pd_1 <- partial(cf.final, pred.var="PGR_Inhibition", grid.resolution = 10)
int1a1_pd_2 <- partial(cf.final, pred.var="Benzo[a]pyrene", grid.resolution = 10)
int1a1_pd_3 <- partial(cf.final, pred.var="ADORA2_ADRA_HRH-Inhibition", grid.resolution = 10)
int1a1_pd_4 <- partial(cf.final, pred.var="Pyrene", grid.resolution = 10)
int1a1_pd_5 <- partial(cf.final, pred.var="PPARG_Activation", grid.resolution = 10)
int1a1_pd_6 <- partial(cf.final, pred.var="androstenedione_Inhibition", grid.resolution = 10)
int1a1_pd_7 <- partial(cf.final, pred.var="NR1I3_Binding", grid.resolution = 10)
int1a1_pd_8 <- partial(cf.final, pred.var="NR1I2_Activation", grid.resolution = 10)
int1a1_pd_9 <- partial(cf.final, pred.var="Anthraquinone", grid.resolution = 10)
int1a1_pd_10 <- partial(cf.final, pred.var="ESR2_Activation", grid.resolution = 10)
int1a1_pd_11 <- partial(cf.final, pred.var="NR1I3_Activation", grid.resolution = 10)
int1a1_pd_12 <- partial(cf.final, pred.var="RXRA_Activation", grid.resolution = 10)
int1a1_pd_13 <- partial(cf.final, pred.var="AR_Inhibition", grid.resolution = 10)

int1a1_pdp_1 <- autoplot(int1a1_pd_1, size=1.2) + theme_classic() + xlab("PGR_Inhibition") + ylab("Effect") +
  theme(text=element_text(size=20)) + ggtitle("CYP 1A1")
int1a1_pdp_2 <- autoplot(int1a1_pd_2, size=1.2) + theme_classic() + xlab("Benzo[a]pyrene") + ylab("Effect") +
  theme(text=element_text(size=20)) + ggtitle("CYP 1A1")
int1a1_pdp_3 <- autoplot(int1a1_pd_3, size=1.2) + theme_classic() + xlab("ADORA2_ADRA_HRH-Inhibition") + ylab("Effect") +
  theme(text=element_text(size=20)) + ggtitle("CYP 1A1")
int1a1_pdp_4 <- autoplot(int1a1_pd_4, size=1.2) + theme_classic() + xlab("Pyrene") + ylab("Effect") +
  theme(text=element_text(size=20)) + ggtitle("CYP 1A1")
int1a1_pdp_5 <- autoplot(int1a1_pd_5, size=1.2) + theme_classic() + xlab("PPARG_Activation") + ylab("Effect") +
  theme(text=element_text(size=20)) + ggtitle("CYP 1A1")
int1a1_pdp_6 <- autoplot(int1a1_pd_6, size=1.2) + theme_classic() + xlab("androstenedione_Inhibition") + ylab("Effect") +
  theme(text=element_text(size=20)) + ggtitle("CYP 1A1")
int1a1_pdp_7 <- autoplot(int1a1_pd_7, size=1.2) + theme_classic() + xlab("NR1I3_Binding") + ylab("Effect") +
  theme(text=element_text(size=20)) + ggtitle("CYP 1A1")
int1a1_pdp_8 <- autoplot(int1a1_pd_8, size=1.2) + theme_classic() + xlab("NR1I2_Activation") + ylab("Effect") +
  theme(text=element_text(size=20)) + ggtitle("CYP 1A1")
int1a1_pdp_9 <- autoplot(int1a1_pd_9, size=1.2) + theme_classic() + xlab("Anthraquinone") + ylab("Effect") +
  theme(text=element_text(size=20)) + ggtitle("CYP 1A1")
int1a1_pdp_10 <- autoplot(int1a1_pd_10, size=1.2) + theme_classic() + xlab("ESR2_Activation") + ylab("Effect") +
  theme(text=element_text(size=20)) + ggtitle("CYP 1A1")
int1a1_pdp_11 <- autoplot(int1a1_pd_11, size=1.2) + theme_classic() + xlab("NR1I3_Activation") + ylab("Effect") +
  theme(text=element_text(size=20)) + ggtitle("CYP 1A1")
int1a1_pdp_12 <- autoplot(int1a1_pd_12, size=1.2) + theme_classic() + xlab("RXRA_Activation") + ylab("Effect") +
  theme(text=element_text(size=20)) + ggtitle("CYP 1A1")
int1a1_pdp_13 <- autoplot(int1a1_pd_13, size=1.2) + theme_classic() + xlab("AR_Inhibition") + ylab("Effect") +
  theme(text=element_text(size=20)) + ggtitle("CYP 1A1")

CYP3A_intestine_final <- CYP3A_intestine %>% dplyr::select(Effect, "NR1I3_Activation", AR_Inhibition)
set.seed(3)
cf.final <- cforest(Effect~., data=CYP3A_intestine_final, control=cforest_unbiased(mtry=1, ntree=5000))

int3A_pd_1 <- partial(cf.final, pred.var="NR1I3_Activation", grid.resolution = 10)
int3A_pd_2 <- partial(cf.final, pred.var="AR_Inhibition", grid.resolution = 10)

int3A_pdp_1 <- autoplot(int3A_pd_1, size=1.2) + theme_classic() + xlab("NR1I3_Activation") + ylab("Effect") +
  theme(text=element_text(size=20)) + ggtitle("CYP 3A")
int3A_pdp_2 <- autoplot(int3A_pd_2, size=1.2) + theme_classic() + xlab("AR_Inhibition") + ylab("Effect") +
  theme(text=element_text(size=20)) + ggtitle("CYP 3A")


# Refit final model
UGT_intestine_final <- UGT_intestine %>% dplyr::select(Effect,"Metformin", "Venlafaxine", "Lidocaine", 
                                                       "Desvenlafaxine", "Fexofenadine", "SLC-Inhibition",
                                                       "Carbamazepine", "SNRIs_WWI", "Isophorone",
                                                       "Fluoranthene", "PPARD_Activation", "N,N-diethyl-meta-toluamide",
                                                       "11-deoxycortisol_Inhibition", "Sulfamethoxazole", "RXRA_Activation",
                                                       "Tramadol", "TSPO_Inhibition", "Triamterene",
                                                       "ATF6_Activation", "Cholesterol", "NR1H4_Inhibition",
                                                       "CYP2B6_Activation", "Indole", "Tris(2-chloroethyl)phosphate",
                                                       "ESRRA_Activation", "TP53_Activation", "CYP24A1_Activation",
                                                       "NFE2L2_Activation", "NR1I2_Activation", "Methotrexate",
                                                       "Antibiotics", "Oprk1_Binding", "ESRRA_Inhibition",
                                                       "Acetaminophen", "HIF1A_Activation", "Carbaryl",
                                                       "Fuels_Industrial_WWIs_Pesticides", "HSF1_Activation", "ESR1_Activation",
                                                       "PPCPs_Pesticides_Plasticizers", "1-Methylnaphthalene", "NR3C1_Activation",
                                                       "HDAC1_Inhibition", "Benzo[a]pyrene", "Tris(2-butoxyethyl)phosphate",
                                                       "OPRM1_Binding", "ADORA2A_Binding", "Metolachlor",
                                                       "CYP1A2_Activation", "PGR_Inhibition", "Industrial_PPCP_Narcotics",
                                                       "PPARG_Activation", "TOX21_CAR_Antagonist", "NR1H4_Activation",
                                                       "PGR_Binding", "beta-Sitosterol", "Pyrene",
                                                       "3,4-Dichlorophenyl isocyanate", "VDR_Activation", "PAHs", 
                                                       "CYP1A2_Inhibition", "CYP2C19_Inhibition", "5-Methyl-1H-benzotriazole",
                                                       "SCN-Inhibition", "4-tert-Octylphenol monoethoxylate", "AR_Activation",
                                                       "estrone_Activation", "testosterone_Inhibition", "androstenedione_Inhibition",
                                                       "NR1I2_Binding", "Sterols", "Anthraquinone",
                                                       "NR1I3_Binding", "Phenanthrene", "NR3C1_Binding",
                                                       "Tributyl phosphate", "CREB3_Activation", "Bisphenol A",
                                                       "Nicotine", "ESR2_Activation", "Carbazole",
                                                       "POU2F1_Activation", "Meprobamate", "NR1I3_Activation", 
                                                       "17alpha-hydroxypregnenolone_Inhibition", "JUN_Activation", "Cyp2a2_Inhibition", "RXRB_Activation"
)

set.seed(3)
ugtint1a1.final <- cforest(Effect~., data=UGT_intestine_final, control=cforest_unbiased(mtry=(29), ntree=5000))
ugtint1a1_pd_1 <- partial(ugtint1a1.final, pred.var="Metformin", grid.resolution = 10)
ugtint1a1_pd_2 <- partial(ugtint1a1.final, pred.var="Venlafaxine", grid.resolution = 10)
ugtint1a1_pd_3 <- partial(ugtint1a1.final, pred.var="Lidocaine", grid.resolution = 10)
ugtint1a1_pd_4 <- partial(ugtint1a1.final, pred.var="Desvenlafaxine", grid.resolution = 10)
ugtint1a1_pd_5 <- partial(ugtint1a1.final, pred.var="Fexofenadine", grid.resolution = 10)
ugtint1a1_pd_6 <- partial(ugtint1a1.final, pred.var="SLC-Inhibition", grid.resolution = 10)
ugtint1a1_pd_7 <- partial(ugtint1a1.final, pred.var="Carbamazepine", grid.resolution = 10)
ugtint1a1_pd_8 <- partial(ugtint1a1.final, pred.var="SNRIs_WWI", grid.resolution = 10)
ugtint1a1_pd_9 <- partial(ugtint1a1.final, pred.var="Isophorone", grid.resolution = 10)
ugtint1a1_pd_10 <- partial(ugtint1a1.final, pred.var="Fluoranthene", grid.resolution = 10)
ugtint1a1_pd_11 <- partial(ugtint1a1.final, pred.var="PPARD_Activation", grid.resolution = 10)
ugtint1a1_pd_12 <- partial(ugtint1a1.final, pred.var="N,N-diethyl-meta-toluamide", grid.resolution = 10)
ugtint1a1_pd_13 <- partial(ugtint1a1.final, pred.var="11-deoxycortisol_Inhibition", grid.resolution = 10)
ugtint1a1_pd_14 <- partial(ugtint1a1.final, pred.var="Sulfamethoxazole", grid.resolution = 10)
ugtint1a1_pd_15 <- partial(ugtint1a1.final, pred.var="RXRA_Activation", grid.resolution = 10)
ugtint1a1_pd_16 <- partial(ugtint1a1.final, pred.var="Tramadol", grid.resolution = 10)
ugtint1a1_pd_17 <- partial(ugtint1a1.final, pred.var="TSPO_Inhibition", grid.resolution = 10)
ugtint1a1_pd_18 <- partial(ugtint1a1.final, pred.var="Triamterene", grid.resolution = 10)
ugtint1a1_pd_19 <- partial(ugtint1a1.final, pred.var="ATF6_Activation", grid.resolution = 10)
ugtint1a1_pd_20 <- partial(ugtint1a1.final, pred.var="Cholesterol", grid.resolution = 10)
ugtint1a1_pd_21 <- partial(ugtint1a1.final, pred.var="NR1H4_Inhibition", grid.resolution = 10)
ugtint1a1_pd_22 <- partial(ugtint1a1.final, pred.var="CYP2B6_Activation", grid.resolution = 10)
ugtint1a1_pd_23 <- partial(ugtint1a1.final, pred.var="Indole", grid.resolution = 10)
ugtint1a1_pd_24 <- partial(ugtint1a1.final, pred.var="Tris(2-chloroethyl)phosphate", grid.resolution = 10)
ugtint1a1_pd_25 <- partial(ugtint1a1.final, pred.var="ESRRA_Activation", grid.resolution = 10)
ugtint1a1_pd_26 <- partial(ugtint1a1.final, pred.var="TP53_Activation", grid.resolution = 10)
ugtint1a1_pd_27 <- partial(ugtint1a1.final, pred.var="CYP24A1_Activation", grid.resolution = 10)
ugtint1a1_pd_28 <- partial(ugtint1a1.final, pred.var="NFE2L2_Activation", grid.resolution = 10)
ugtint1a1_pd_29 <- partial(ugtint1a1.final, pred.var="NR1I2_Activation", grid.resolution = 10)
ugtint1a1_pd_30 <- partial(ugtint1a1.final, pred.var="Methotrexate", grid.resolution = 10)
ugtint1a1_pd_31 <- partial(ugtint1a1.final, pred.var="Antibiotics", grid.resolution = 10)
ugtint1a1_pd_32 <- partial(ugtint1a1.final, pred.var="Oprk1_Binding", grid.resolution = 10)
ugtint1a1_pd_33 <- partial(ugtint1a1.final, pred.var="ESRRA_Inhibition", grid.resolution = 10)
ugtint1a1_pd_34 <- partial(ugtint1a1.final, pred.var="Acetaminophen", grid.resolution = 10)
ugtint1a1_pd_35 <- partial(ugtint1a1.final, pred.var="HIF1A_Activation", grid.resolution = 10)
ugtint1a1_pd_36 <- partial(ugtint1a1.final, pred.var="Carbaryl", grid.resolution = 10)
ugtint1a1_pd_37 <- partial(ugtint1a1.final, pred.var="Fuels_Industrial_WWIs_Pesticides", grid.resolution = 10)
ugtint1a1_pd_38 <- partial(ugtint1a1.final, pred.var="HSF1_Activation", grid.resolution = 10)
ugtint1a1_pd_39 <- partial(ugtint1a1.final, pred.var="ESR1_Activation", grid.resolution = 10)
ugtint1a1_pd_40 <- partial(ugtint1a1.final, pred.var="PPCPs_Pesticides_Plasticizers", grid.resolution = 10)
ugtint1a1_pd_41 <- partial(ugtint1a1.final, pred.var="1-Methylnaphthalene", grid.resolution = 10)
ugtint1a1_pd_42 <- partial(ugtint1a1.final, pred.var="NR3C1_Activation", grid.resolution = 10)
ugtint1a1_pd_43 <- partial(ugtint1a1.final, pred.var="HDAC1_Inhibition", grid.resolution = 10)
ugtint1a1_pd_44 <- partial(ugtint1a1.final, pred.var="Benzo[a]pyrene", grid.resolution = 10)
ugtint1a1_pd_45 <- partial(ugtint1a1.final, pred.var="Tris(2-butoxyethyl)phosphate", grid.resolution = 10)
ugtint1a1_pd_46 <- partial(ugtint1a1.final, pred.var="OPRM1_Binding", grid.resolution = 10)
ugtint1a1_pd_47 <- partial(ugtint1a1.final, pred.var="ADORA2A_Binding", grid.resolution = 10)
ugtint1a1_pd_48 <- partial(ugtint1a1.final, pred.var="Metolachlor", grid.resolution = 10)
ugtint1a1_pd_49 <- partial(ugtint1a1.final, pred.var="CYP1A2_Activation", grid.resolution = 10)
ugtint1a1_pd_50 <- partial(ugtint1a1.final, pred.var="PGR_Inhibition", grid.resolution = 10)
ugtint1a1_pd_51 <- partial(ugtint1a1.final, pred.var="Industrial_PPCP_Narcotics", grid.resolution = 10)
ugtint1a1_pd_52 <- partial(ugtint1a1.final, pred.var="PPARG_Activation", grid.resolution = 10)
ugtint1a1_pd_53 <- partial(ugtint1a1.final, pred.var="TOX21_CAR_Antagonist", grid.resolution = 10)
ugtint1a1_pd_54 <- partial(ugtint1a1.final, pred.var="NR1H4_Activation", grid.resolution = 10)
ugtint1a1_pd_55 <- partial(ugtint1a1.final, pred.var="PGR_Binding", grid.resolution = 10)
ugtint1a1_pd_56 <- partial(ugtint1a1.final, pred.var="beta-Sitosterol", grid.resolution = 10)
ugtint1a1_pd_57 <- partial(ugtint1a1.final, pred.var="Pyrene", grid.resolution = 10)
ugtint1a1_pd_58 <- partial(ugtint1a1.final, pred.var="3,4-Dichlorophenyl isocyanate", grid.resolution = 10)
ugtint1a1_pd_59 <- partial(ugtint1a1.final, pred.var="VDR_Activation", grid.resolution = 10)
ugtint1a1_pd_60 <- partial(ugtint1a1.final, pred.var="PAHs", grid.resolution = 10)
ugtint1a1_pd_61 <- partial(ugtint1a1.final, pred.var="CYP1A2_Activation", grid.resolution = 10)
ugtint1a1_pd_62 <- partial(ugtint1a1.final, pred.var="CYP2C19_Inhibition", grid.resolution = 10)
ugtint1a1_pd_63 <- partial(ugtint1a1.final, pred.var="5-Methyl-1H-benzotriazole", grid.resolution = 10)
ugtint1a1_pd_64 <- partial(ugtint1a1.final, pred.var="SCN-Inhibition", grid.resolution = 10)
ugtint1a1_pd_65 <- partial(ugtint1a1.final, pred.var="4-tert-Octylphenol monoethoxylate", grid.resolution = 10)
ugtint1a1_pd_66 <- partial(ugtint1a1.final, pred.var="AR_Activation", grid.resolution = 10)
ugtint1a1_pd_67 <- partial(ugtint1a1.final, pred.var="estrone_Activation", grid.resolution = 10)
ugtint1a1_pd_68 <- partial(ugtint1a1.final, pred.var="testosterone_Inhibition", grid.resolution = 10)
ugtint1a1_pd_69 <- partial(ugtint1a1.final, pred.var="androstenedione_Inhibition", grid.resolution = 10)
ugtint1a1_pd_70 <- partial(ugtint1a1.final, pred.var="NR1I2_Binding", grid.resolution = 10)
ugtint1a1_pd_71 <- partial(ugtint1a1.final, pred.var="Sterols", grid.resolution = 10)
ugtint1a1_pd_72 <- partial(ugtint1a1.final, pred.var="Anthraquinone", grid.resolution = 10)
ugtint1a1_pd_73 <- partial(ugtint1a1.final, pred.var="NR1I3_Binding", grid.resolution = 10)
ugtint1a1_pd_74 <- partial(ugtint1a1.final, pred.var="Phenanthrene", grid.resolution = 10)
ugtint1a1_pd_75 <- partial(ugtint1a1.final, pred.var="NR3C1_Binding", grid.resolution = 10)
ugtint1a1_pd_76 <- partial(ugtint1a1.final, pred.var="Tributyl phosphate", grid.resolution = 10)
ugtint1a1_pd_77 <- partial(ugtint1a1.final, pred.var="CREB3_Activation", grid.resolution = 10)
ugtint1a1_pd_78 <- partial(ugtint1a1.final, pred.var="Bisphenol A", grid.resolution = 10)
ugtint1a1_pd_79 <- partial(ugtint1a1.final, pred.var="Nicotine", grid.resolution = 10)
ugtint1a1_pd_80 <- partial(ugtint1a1.final, pred.var="ESR2_Activation", grid.resolution = 10)
ugtint1a1_pd_81 <- partial(ugtint1a1.final, pred.var="Carbazole", grid.resolution = 10)
ugtint1a1_pd_82 <- partial(ugtint1a1.final, pred.var="POU2F1_Activation", grid.resolution = 10)
ugtint1a1_pd_83 <- partial(ugtint1a1.final, pred.var="Meprobamate", grid.resolution = 10)
ugtint1a1_pd_84 <- partial(ugtint1a1.final, pred.var="NR1I3_Activation", grid.resolution = 10)
ugtint1a1_pd_85 <- partial(ugtint1a1.final, pred.var="17alpha-hydroxypregnenolone_Inhibition", grid.resolution = 10)
ugtint1a1_pd_86 <- partial(ugtint1a1.final, pred.var="JUN_Activation", grid.resolution = 10)
ugtint1a1_pd_87 <- partial(ugtint1a1.final, pred.var="Cyp2a2_Inhibition", grid.resolution = 10)
ugtint1a1_pd_88 <- partial(ugtint1a1.final, pred.var="RXRB_Activation", grid.resolution = 10)

ugtint1a1_pdp_1 <- autoplot(ugtint1a1_pd_1, size=1.2) + theme_classic() + xlab("Metformin") + ylab("Effect") +
  theme(text=element_text(size=20)) + ggtitle("UGT 1A1")
ugtint1a1_pdp_2 <- autoplot(ugtint1a1_pd_2, size=1.2) + theme_classic() + xlab("Venlafaxine") + ylab("Effect") +
  theme(text=element_text(size=20)) + ggtitle("UGT 1A1")
ugtint1a1_pdp_3 <- autoplot(ugtint1a1_pd_3, size=1.2) + theme_classic() + xlab("Lidocaine") + ylab("Effect") +
  theme(text=element_text(size=20)) + ggtitle("UGT 1A1")
ugtint1a1_pdp_4 <- autoplot(ugtint1a1_pd_4, size=1.2) + theme_classic() + xlab("Desvenlafaxine") + ylab("Effect") +
  theme(text=element_text(size=20)) + ggtitle("UGT 1A1")
ugtint1a1_pdp_5 <- autoplot(ugtint1a1_pd_5, size=1.2) + theme_classic() + xlab("Fexofenadine") + ylab("Effect") +
  theme(text=element_text(size=20)) + ggtitle("UGT 1A1")
ugtint1a1_pdp_6 <- autoplot(ugtint1a1_pd_6, size=1.2) + theme_classic() + xlab("SLC-Inhibition") + ylab("Effect") +
  theme(text=element_text(size=20)) + ggtitle("UGT 1A1")
ugtint1a1_pdp_7 <- autoplot(ugtint1a1_pd_7, size=1.2) + theme_classic() + xlab("Carbamazepine") + ylab("Effect") +
  theme(text=element_text(size=20)) + ggtitle("UGT 1A1")
ugtint1a1_pdp_8 <- autoplot(ugtint1a1_pd_8, size=1.2) + theme_classic() + xlab("SNRIs_WWI") + ylab("Effect") +
  theme(text=element_text(size=20)) + ggtitle("UGT 1A1")
ugtint1a1_pdp_9 <- autoplot(ugtint1a1_pd_9, size=1.2) + theme_classic() + xlab("Isophorone") + ylab("Effect") +
  theme(text=element_text(size=20)) + ggtitle("UGT 1A1")
ugtint1a1_pdp_10 <- autoplot(ugtint1a1_pd_10, size=1.2) + theme_classic() + xlab("Fluoranthene") + ylab("Effect") +
  theme(text=element_text(size=20)) + ggtitle("UGT 1A1")
ugtint1a1_pdp_11 <- autoplot(ugtint1a1_pd_11, size=1.2) + theme_classic() + xlab("PPARD_Activation") + ylab("Effect") +
  theme(text=element_text(size=20)) + ggtitle("UGT 1A1")
ugtint1a1_pdp_12 <- autoplot(ugtint1a1_pd_12, size=1.2) + theme_classic() + xlab("N,N-diethyl-meta-toluamide") + ylab("Effect") +
  theme(text=element_text(size=20)) + ggtitle("UGT 1A1")
ugtint1a1_pdp_13 <- autoplot(ugtint1a1_pd_13, size=1.2) + theme_classic() + xlab("11-deoxycortisol_Inhibition") + ylab("Effect") +
  theme(text=element_text(size=20)) + ggtitle("UGT 1A1")
ugtint1a1_pdp_14 <- autoplot(ugtint1a1_pd_14, size=1.2) + theme_classic() + xlab("Sulfamethoxazole") + ylab("Effect") +
  theme(text=element_text(size=20)) + ggtitle("UGT 1A1")
ugtint1a1_pdp_15 <- autoplot(ugtint1a1_pd_15, size=1.2) + theme_classic() + xlab("RXRA_Activation") + ylab("Effect") +
  theme(text=element_text(size=20)) + ggtitle("UGT 1A1")
ugtint1a1_pdp_16 <- autoplot(ugtint1a1_pd_16, size=1.2) + theme_classic() + xlab("Tramadol") + ylab("Effect") +
  theme(text=element_text(size=20)) + ggtitle("UGT 1A1")
ugtint1a1_pdp_17 <- autoplot(ugtint1a1_pd_17, size=1.2) + theme_classic() + xlab("TSPO_Inhibition") + ylab("Effect") +
  theme(text=element_text(size=20)) + ggtitle("UGT 1A1")
ugtint1a1_pdp_18 <- autoplot(ugtint1a1_pd_18, size=1.2) + theme_classic() + xlab("Triamterene") + ylab("Effect") +
  theme(text=element_text(size=20)) + ggtitle("UGT 1A1")
ugtint1a1_pdp_19 <- autoplot(ugtint1a1_pd_19, size=1.2) + theme_classic() + xlab("ATF6_Activation") + ylab("Effect") +
  theme(text=element_text(size=20)) + ggtitle("UGT 1A1")
ugtint1a1_pdp_20 <- autoplot(ugtint1a1_pd_20, size=1.2) + theme_classic() + xlab("Cholesterol") + ylab("Effect") +
  theme(text=element_text(size=20)) + ggtitle("UGT 1A1")
ugtint1a1_pdp_21 <- autoplot(ugtint1a1_pd_21, size=1.2) + theme_classic() + xlab("NR1H4_Inhibition") + ylab("Effect") +
  theme(text=element_text(size=20)) + ggtitle("UGT 1A1")
ugtint1a1_pdp_22 <- autoplot(ugtint1a1_pd_22, size=1.2) + theme_classic() + xlab("CYP2B6_Activation") + ylab("Effect") +
  theme(text=element_text(size=20)) + ggtitle("UGT 1A1")
ugtint1a1_pdp_23 <- autoplot(ugtint1a1_pd_23, size=1.2) + theme_classic() + xlab("Indole") + ylab("Effect") +
  theme(text=element_text(size=20)) + ggtitle("UGT 1A1")
ugtint1a1_pdp_24 <- autoplot(ugtint1a1_pd_24, size=1.2) + theme_classic() + xlab("Tris(2-chloroethyl)phosphate") + ylab("Effect") +
  theme(text=element_text(size=20)) + ggtitle("UGT 1A1")
ugtint1a1_pdp_25 <- autoplot(ugtint1a1_pd_25, size=1.2) + theme_classic() + xlab("ESRRA_Activation") + ylab("Effect") +
  theme(text=element_text(size=20)) + ggtitle("UGT 1A1")
ugtint1a1_pdp_26 <- autoplot(ugtint1a1_pd_26, size=1.2) + theme_classic() + xlab("TP53_Activation") + ylab("Effect") +
  theme(text=element_text(size=20)) + ggtitle("UGT 1A1")
ugtint1a1_pdp_27 <- autoplot(ugtint1a1_pd_27, size=1.2) + theme_classic() + xlab("CYP24A1_Activation") + ylab("Effect") +
  theme(text=element_text(size=20)) + ggtitle("UGT 1A1")
ugtint1a1_pdp_28 <- autoplot(ugtint1a1_pd_28, size=1.2) + theme_classic() + xlab("NFE2L2_Activation") + ylab("Effect") +
  theme(text=element_text(size=20)) + ggtitle("UGT 1A1")
ugtint1a1_pdp_29 <- autoplot(ugtint1a1_pd_29, size=1.2) + theme_classic() + xlab("NR1I2_Activation") + ylab("Effect") +
  theme(text=element_text(size=20)) + ggtitle("UGT 1A1")
ugtint1a1_pdp_30 <- autoplot(ugtint1a1_pd_30, size=1.2) + theme_classic() + xlab("Methotrexate") + ylab("Effect") +
  theme(text=element_text(size=20)) + ggtitle("UGT 1A1")
ugtint1a1_pdp_31 <- autoplot(ugtint1a1_pd_31, size=1.2) + theme_classic() + xlab("Antibiotics") + ylab("Effect") +
  theme(text=element_text(size=20)) + ggtitle("UGT 1A1")
ugtint1a1_pdp_33 <- autoplot(ugtint1a1_pd_32, size=1.2) + theme_classic() + xlab("Oprk1_Binding") + ylab("Effect") +
  theme(text=element_text(size=20)) + ggtitle("UGT 1A1")
ugtint1a1_pdp_23 <- autoplot(ugtint1a1_pd_33, size=1.2) + theme_classic() + xlab("ESRRA_Inhibition") + ylab("Effect") +
  theme(text=element_text(size=20)) + ggtitle("UGT 1A1")
ugtint1a1_pdp_34 <- autoplot(ugtint1a1_pd_34, size=1.2) + theme_classic() + xlab("Acetaminophen") + ylab("Effect") +
  theme(text=element_text(size=20)) + ggtitle("UGT 1A1")
ugtint1a1_pdp_35 <- autoplot(ugtint1a1_pd_35, size=1.2) + theme_classic() + xlab("HIF1A_Activation") + ylab("Effect") +
  theme(text=element_text(size=20)) + ggtitle("UGT 1A1")
ugtint1a1_pdp_36 <- autoplot(ugtint1a1_pd_36, size=1.2) + theme_classic() + xlab("Carbaryl") + ylab("Effect") +
  theme(text=element_text(size=20)) + ggtitle("UGT 1A1")
ugtint1a1_pdp_37 <- autoplot(ugtint1a1_pd_37, size=1.2) + theme_classic() + xlab("Fuels_Industrial_WWIs_Pesticides") + ylab("Effect") +
  theme(text=element_text(size=20)) + ggtitle("UGT 1A1")
ugtint1a1_pdp_38 <- autoplot(ugtint1a1_pd_38, size=1.2) + theme_classic() + xlab("HSF1_Activation") + ylab("Effect") +
  theme(text=element_text(size=20)) + ggtitle("UGT 1A1")
ugtint1a1_pdp_39 <- autoplot(ugtint1a1_pd_39, size=1.2) + theme_classic() + xlab("ESR1_Activation") + ylab("Effect") +
  theme(text=element_text(size=20)) + ggtitle("UGT 1A1")
ugtint1a1_pdp_40 <- autoplot(ugtint1a1_pd_40, size=1.2) + theme_classic() + xlab("PPCPs_Pesticides_Plasticizers") + ylab("Effect") +
  theme(text=element_text(size=20)) + ggtitle("UGT 1A1")
ugtint1a1_pdp_41 <- autoplot(ugtint1a1_pd_41, size=1.2) + theme_classic() + xlab("1-Methylnaphthalene") + ylab("Effect") +
  theme(text=element_text(size=20)) + ggtitle("UGT 1A1")
ugtint1a1_pdp_42 <- autoplot(ugtint1a1_pd_42, size=1.2) + theme_classic() + xlab("NR3C1_Activation") + ylab("Effect") +
  theme(text=element_text(size=20)) + ggtitle("UGT 1A1")
ugtint1a1_pdp_43 <- autoplot(ugtint1a1_pd_43, size=1.2) + theme_classic() + xlab("HDAC1_Inhibition") + ylab("Effect") +
  theme(text=element_text(size=20)) + ggtitle("UGT 1A1")
ugtint1a1_pdp_44 <- autoplot(ugtint1a1_pd_44, size=1.2) + theme_classic() + xlab("Benzo[a]pyrene") + ylab("Effect") +
  theme(text=element_text(size=20)) + ggtitle("UGT 1A1")
ugtint1a1_pdp_45 <- autoplot(ugtint1a1_pd_45, size=1.2) + theme_classic() + xlab("Tris(2-butoxyethyl)phosphate") + ylab("Effect") +
  theme(text=element_text(size=20)) + ggtitle("UGT 1A1")
ugtint1a1_pdp_46 <- autoplot(ugtint1a1_pd_46, size=1.2) + theme_classic() + xlab("OPRM1_Binding") + ylab("Effect") +
  theme(text=element_text(size=20)) + ggtitle("UGT 1A1")
ugtint1a1_pdp_47 <- autoplot(ugtint1a1_pd_47, size=1.2) + theme_classic() + xlab("ADORA2A_Binding") + ylab("Effect") +
  theme(text=element_text(size=20)) + ggtitle("UGT 1A1")
ugtint1a1_pdp_48 <- autoplot(ugtint1a1_pd_48, size=1.2) + theme_classic() + xlab("Metolachlor") + ylab("Effect") +
  theme(text=element_text(size=20)) + ggtitle("UGT 1A1")
ugtint1a1_pdp_49 <- autoplot(ugtint1a1_pd_49, size=1.2) + theme_classic() + xlab("CYP1A2_Activation") + ylab("Effect") +
  theme(text=element_text(size=20)) + ggtitle("UGT 1A1")
ugtint1a1_pdp_50 <- autoplot(ugtint1a1_pd_50, size=1.2) + theme_classic() + xlab("PGR_Inhibition") + ylab("Effect") +
  theme(text=element_text(size=20)) + ggtitle("UGT 1A1")
ugtint1a1_pdp_51 <- autoplot(ugtint1a1_pd_51, size=1.2) + theme_classic() + xlab("Industrial_PPCP_Narcotics") + ylab("Effect") +
  theme(text=element_text(size=20)) + ggtitle("UGT 1A1")
ugtint1a1_pdp_52 <- autoplot(ugtint1a1_pd_52, size=1.2) + theme_classic() + xlab("PPARG_Activation") + ylab("Effect") +
  theme(text=element_text(size=20)) + ggtitle("UGT 1A1")
ugtint1a1_pdp_53 <- autoplot(ugtint1a1_pd_53, size=1.2) + theme_classic() + xlab("TOX21_CAR_Antagonist") + ylab("Effect") +
  theme(text=element_text(size=20)) + ggtitle("UGT 1A1")
ugtint1a1_pdp_54 <- autoplot(ugtint1a1_pd_54, size=1.2) + theme_classic() + xlab("NR1H4_Activation") + ylab("Effect") +
  theme(text=element_text(size=20)) + ggtitle("UGT 1A1")
ugtint1a1_pdp_55 <- autoplot(ugtint1a1_pd_55, size=1.2) + theme_classic() + xlab("PGR_Binding") + ylab("Effect") +
  theme(text=element_text(size=20)) + ggtitle("UGT 1A1")
ugtint1a1_pdp_56 <- autoplot(ugtint1a1_pd_56, size=1.2) + theme_classic() + xlab("beta-Sitosterol") + ylab("Effect") +
  theme(text=element_text(size=20)) + ggtitle("UGT 1A1")
ugtint1a1_pdp_57 <- autoplot(ugtint1a1_pd_57, size=1.2) + theme_classic() + xlab("Pyrene") + ylab("Effect") +
  theme(text=element_text(size=20)) + ggtitle("UGT 1A1")
ugtint1a1_pdp_58 <- autoplot(ugtint1a1_pd_58, size=1.2) + theme_classic() + xlab("3,4-Dichlorophenyl isocyanate") + ylab("Effect") +
  theme(text=element_text(size=20)) + ggtitle("UGT 1A1")
ugtint1a1_pdp_59 <- autoplot(ugtint1a1_pd_59, size=1.2) + theme_classic() + xlab("VDR_Activation") + ylab("Effect") +
  theme(text=element_text(size=20)) + ggtitle("UGT 1A1")
ugtint1a1_pdp_60 <- autoplot(ugtint1a1_pd_60, size=1.2) + theme_classic() + xlab("PAHs") + ylab("Effect") +
  theme(text=element_text(size=20)) + ggtitle("UGT 1A1")
ugtint1a1_pdp_61 <- autoplot(ugtint1a1_pd_61, size=1.2) + theme_classic() + xlab("CYP1A2_Inhibition") + ylab("Effect") +
  theme(text=element_text(size=20)) + ggtitle("UGT 1A1")
ugtint1a1_pdp_62 <- autoplot(ugtint1a1_pd_62, size=1.2) + theme_classic() + xlab("CYP2C19_Inhibition") + ylab("Effect") +
  theme(text=element_text(size=20)) + ggtitle("UGT 1A1")
ugtint1a1_pdp_63 <- autoplot(ugtint1a1_pd_63, size=1.2) + theme_classic() + xlab("5-Methyl-1H-benzotriazole") + ylab("Effect") +
  theme(text=element_text(size=20)) + ggtitle("UGT 1A1")
ugtint1a1_pdp_64 <- autoplot(ugtint1a1_pd_64, size=1.2) + theme_classic() + xlab("SCN-Inhibition") + ylab("Effect") +
  theme(text=element_text(size=20)) + ggtitle("UGT 1A1")
ugtint1a1_pdp_65 <- autoplot(ugtint1a1_pd_65, size=1.2) + theme_classic() + xlab("4-tert-Octylphenol monoethoxylate") + ylab("Effect") +
  theme(text=element_text(size=20)) + ggtitle("UGT 1A1")
ugtint1a1_pdp_66 <- autoplot(ugtint1a1_pd_66, size=1.2) + theme_classic() + xlab("AR_Activation") + ylab("Effect") +
  theme(text=element_text(size=20)) + ggtitle("UGT 1A1")
ugtint1a1_pdp_67 <- autoplot(ugtint1a1_pd_67, size=1.2) + theme_classic() + xlab("estrone_Activation") + ylab("Effect") +
  theme(text=element_text(size=20)) + ggtitle("UGT 1A1")
ugtint1a1_pdp_68 <- autoplot(ugtint1a1_pd_68, size=1.2) + theme_classic() + xlab("testosterone_Inhibition") + ylab("Effect") +
  theme(text=element_text(size=20)) + ggtitle("UGT 1A1")
ugtint1a1_pdp_69 <- autoplot(ugtint1a1_pd_69, size=1.2) + theme_classic() + xlab("androstenedione_Inhibition") + ylab("Effect") +
  theme(text=element_text(size=20)) + ggtitle("UGT 1A1")
ugtint1a1_pdp_70 <- autoplot(ugtint1a1_pd_70, size=1.2) + theme_classic() + xlab("NR1I2_Binding") + ylab("Effect") +
  theme(text=element_text(size=20)) + ggtitle("UGT 1A1")
ugtint1a1_pdp_71 <- autoplot(ugtint1a1_pd_71, size=1.2) + theme_classic() + xlab("Sterols") + ylab("Effect") +
  theme(text=element_text(size=20)) + ggtitle("UGT 1A1")
ugtint1a1_pdp_72 <- autoplot(ugtint1a1_pd_72, size=1.2) + theme_classic() + xlab("Anthraquinone") + ylab("Effect") +
  theme(text=element_text(size=20)) + ggtitle("UGT 1A1")
ugtint1a1_pdp_73 <- autoplot(ugtint1a1_pd_73, size=1.2) + theme_classic() + xlab("NR1I3_Binding") + ylab("Effect") +
  theme(text=element_text(size=20)) + ggtitle("UGT 1A1")
ugtint1a1_pdp_74 <- autoplot(ugtint1a1_pd_74, size=1.2) + theme_classic() + xlab("Phenanthrene") + ylab("Effect") +
  theme(text=element_text(size=20)) + ggtitle("UGT 1A1")
ugtint1a1_pdp_75 <- autoplot(ugtint1a1_pd_75, size=1.2) + theme_classic() + xlab("NR3C1_Binding") + ylab("Effect") +
  theme(text=element_text(size=20)) + ggtitle("UGT 1A1")
ugtint1a1_pdp_76 <- autoplot(ugtint1a1_pd_76, size=1.2) + theme_classic() + xlab("Tributyl phosphate") + ylab("Effect") +
  theme(text=element_text(size=20)) + ggtitle("UGT 1A1")
ugtint1a1_pdp_77 <- autoplot(ugtint1a1_pd_77, size=1.2) + theme_classic() + xlab("CREB3_Activation") + ylab("Effect") +
  theme(text=element_text(size=20)) + ggtitle("UGT 1A1")
ugtint1a1_pdp_78 <- autoplot(ugtint1a1_pd_78, size=1.2) + theme_classic() + xlab("Bisphenol A") + ylab("Effect") +
  theme(text=element_text(size=20)) + ggtitle("UGT 1A1")
ugtint1a1_pdp_79 <- autoplot(ugtint1a1_pd_79, size=1.2) + theme_classic() + xlab("Nicotine") + ylab("Effect") +
  theme(text=element_text(size=20)) + ggtitle("UGT 1A1")
ugtint1a1_pdp_80 <- autoplot(ugtint1a1_pd_80, size=1.2) + theme_classic() + xlab("ESR2_Activation") + ylab("Effect") +
  theme(text=element_text(size=20)) + ggtitle("UGT 1A1")
ugtint1a1_pdp_81 <- autoplot(ugtint1a1_pd_81, size=1.2) + theme_classic() + xlab("Carbazole") + ylab("Effect") +
  theme(text=element_text(size=20)) + ggtitle("UGT 1A1")
ugtint1a1_pdp_82 <- autoplot(ugtint1a1_pd_82, size=1.2) + theme_classic() + xlab("POU2F1_Activation") + ylab("Effect") +
  theme(text=element_text(size=20)) + ggtitle("UGT 1A1")
ugtint1a1_pdp_83 <- autoplot(ugtint1a1_pd_83, size=1.2) + theme_classic() + xlab("Meprobamate") + ylab("Effect") +
  theme(text=element_text(size=20)) + ggtitle("UGT 1A1")
ugtint1a1_pdp_84 <- autoplot(ugtint1a1_pd_84, size=1.2) + theme_classic() + xlab("NR1I3_Activation") + ylab("Effect") +
  theme(text=element_text(size=20)) + ggtitle("UGT 1A1")
ugtint1a1_pdp_85 <- autoplot(ugtint1a1_pd_85, size=1.2) + theme_classic() + xlab("17alpha-hydroxypregnenolone_Inhibition") + ylab("Effect") +
  theme(text=element_text(size=20)) + ggtitle("UGT 1A1")
ugtint1a1_pdp_86 <- autoplot(ugtint1a1_pd_86, size=1.2) + theme_classic() + xlab("JUN_Activation") + ylab("Effect") +
  theme(text=element_text(size=20)) + ggtitle("UGT 1A1")
ugtint1a1_pdp_87 <- autoplot(ugtint1a1_pd_87, size=1.2) + theme_classic() + xlab("Cyp2a2_Inhibition") + ylab("Effect") +
  theme(text=element_text(size=20)) + ggtitle("UGT 1A1")
ugtint1a1_pdp_88 <- autoplot(ugtint1a1_pd_88, size=1.2) + theme_classic() + xlab("RXRB_Activation") + ylab("Effect") +
  theme(text=element_text(size=20)) + ggtitle("UGT 1A1")


CYP_2AD6_final <- CYP_2AD6 %>% dplyr::select(Effect, "Desvenlafaxine", "SLC-Inhibition", "SNRIs_WWI",
                                             "Carbamazepine", "SCN-Inhibition", "Tramadol",
                                             "Lidocaine", "Cholesterol", "RXRA_Activation",
                                             "Venlafaxine", "Metformin", "Fexofenadine",
                                             "CYP2B6_Activation", "N,N-diethyl-meta-toluamide"
)
set.seed(3)
int2ad6.final <- cforest(Effect~., data=CYP_2AD6_final, control=cforest_unbiased(mtry=5, ntree=5000))

int2ad6_pd_1 <- partial(int2ad6.final, pred.var="Desvenlafaxine", grid.resolution = 10)
int2ad6_pd_2 <- partial(int2ad6.final, pred.var="SLC-Inhibition", grid.resolution = 10)
int2ad6_pd_3 <- partial(int2ad6.final, pred.var="SNRIs_WWI", grid.resolution = 10)
int2ad6_pd_4 <- partial(int2ad6.final, pred.var="Carbamazepine", grid.resolution = 10)
int2ad6_pd_5 <- partial(int2ad6.final, pred.var="Tramadol", grid.resolution = 10)
int2ad6_pd_6 <- partial(int2ad6.final, pred.var="Lidocaine", grid.resolution = 10)
int2ad6_pd_7 <- partial(int2ad6.final, pred.var="Cholesterol", grid.resolution = 10)
int2ad6_pd_8 <- partial(int2ad6.final, pred.var="RXRA_Activation", grid.resolution = 10)
int2ad6_pd_9 <- partial(int2ad6.final, pred.var="Venlafaxine", grid.resolution = 10)
int2ad6_pd_10 <- partial(int2ad6.final, pred.var="Metformin", grid.resolution = 10)
int2ad6_pd_11 <- partial(int2ad6.final, pred.var="Fexofenadine", grid.resolution = 10)
int2ad6_pd_12 <- partial(int2ad6.final, pred.var="CYP2B6_Activation", grid.resolution = 10)
int2ad6_pd_13 <- partial(int2ad6.final, pred.var="N,N-diethyl-meta-toluamide", grid.resolution = 10)

int2ad6_pdp_1 <- autoplot(int2ad6_pd_1, size=1.2) + theme_classic() + xlab("Desvenlafaxine") + ylab("Effect") +
  theme(text=element_text(size=20)) + ggtitle("CYP 2AD6")
int2ad6_pdp_2 <- autoplot(int2ad6_pd_2, size=1.2) + theme_classic() + xlab("SLC-Inhibition") + ylab("Effect") +
  theme(text=element_text(size=20)) + ggtitle("CYP 2AD6")
int2ad6_pdp_3 <- autoplot(int2ad6_pd_3, size=1.2) + theme_classic() + xlab("SNRIs_WWI") + ylab("Effect") +
  theme(text=element_text(size=20)) + ggtitle("CYP 2AD6")
int2ad6_pdp_4 <- autoplot(int2ad6_pd_4, size=1.2) + theme_classic() + xlab("Carbamazepine") + ylab("Effect") +
  theme(text=element_text(size=20)) + ggtitle("CYP 2AD6")
int2ad6_pdp_5 <- autoplot(int2ad6_pd_5, size=1.2) + theme_classic() + xlab("Tramadol") + ylab("Effect") +
  theme(text=element_text(size=20)) + ggtitle("CYP 2AD6")
int2ad6_pdp_6 <- autoplot(int2ad6_pd_6, size=1.2) + theme_classic() + xlab("Lidocaine") + ylab("Effect") +
  theme(text=element_text(size=20)) + ggtitle("CYP 2AD6")
int2ad6_pdp_7 <- autoplot(int2ad6_pd_7, size=1.2) + theme_classic() + xlab("Cholesterol") + ylab("Effect") +
  theme(text=element_text(size=20)) + ggtitle("CYP 2AD6")
int2ad6_pdp_8 <- autoplot(int2ad6_pd_8, size=1.2) + theme_classic() + xlab("RXRA_Activation") + ylab("Effect") +
  theme(text=element_text(size=20)) + ggtitle("CYP 2AD6")
int2ad6_pdp_9 <- autoplot(int2ad6_pd_9, size=1.2) + theme_classic() + xlab("Venlafaxine") + ylab("Effect") +
  theme(text=element_text(size=20)) + ggtitle("CYP 2AD6")
int2ad6_pdp_10 <- autoplot(int2ad6_pd_10, size=1.2) + theme_classic() + xlab("Metformin") + ylab("Effect") +
  theme(text=element_text(size=20)) + ggtitle("CYP 2AD6")
int2ad6_pdp_11 <- autoplot(int2ad6_pd_11, size=1.2) + theme_classic() + xlab("Fexofenadine") + ylab("Effect") +
  theme(text=element_text(size=20)) + ggtitle("CYP 2AD6")
int2ad6_pdp_12 <- autoplot(int2ad6_pd_12, size=1.2) + theme_classic() + xlab("CYP2B6_Activation") + ylab("Effect") +
  theme(text=element_text(size=20)) + ggtitle("CYP 2AD6")
int2ad6_pdp_13 <- autoplot(int2ad6_pd_13, size=1.2) + theme_classic() + xlab("N,N-diethyl-meta-toluamide") + ylab("Effect") +
  theme(text=element_text(size=20)) + ggtitle("CYP 2AD6")



CYP2N13_final <- CYP2N13 %>% dplyr::select(Effect, "PGR_Inhibition" , "SCN-Inhibition" , "PAHs" ,
                                           "testosterone_Inhibition" , "NR1H4_Activation" , "TSPO_Inhibition" ,
                                           "NR1I2_Activation" , "androstenedione_Inhibition" , "Fexofenadine" ,
                                           "NR1I3_Binding" , "NFE2L2_Activation" , "Tramadol" ,
                                           "17alpha-hydroxypregnenolone_Inhibition" , "Metformin" , "Carbamazepine" ,
                                           "VDR_Activation" , "ESRRA_Inhibition" , "Lidocaine" ,
                                           "RXRA_Activation" , "ESR1_Activation" , "Venlafaxine" ,
                                           "Desvenlafaxine" , "SLC-Inhibition" , "SNRIs_WWI" 
)
set.seed(3)
int2n13.final <- cforest(Effect~., data=CYP2N13_final, control=cforest_unbiased(mtry=8, ntree=5000))

int2n13_pd_1 <- partial(int2n13.final, pred.var="PGR_Inhibition", grid.resolution = 10)
int2n13_pd_2 <- partial(int2n13.final, pred.var="SCN-Inhibition", grid.resolution = 10)
int2n13_pd_3 <- partial(int2n13.final, pred.var="PAHs", grid.resolution = 10)
int2n13_pd_4 <- partial(int2n13.final, pred.var="testosterone_Inhibition", grid.resolution = 10)
int2n13_pd_5 <- partial(int2n13.final, pred.var="NR1H4_Activation", grid.resolution = 10)
int2n13_pd_6 <- partial(int2n13.final, pred.var="TSPO_Inhibition", grid.resolution = 10)
int2n13_pd_7 <- partial(int2n13.final, pred.var="NR1I3_Binding", grid.resolution = 10)
int2n13_pd_8 <- partial(int2n13.final, pred.var="NFE2L2_Activation", grid.resolution = 10)
int2n13_pd_9 <- partial(int2n13.final, pred.var="Tramadol", grid.resolution = 10)
int2n13_pd_10 <- partial(int2n13.final, pred.var="17alpha-hydroxypregnenolone_Inhibition", grid.resolution = 10)
int2n13_pd_11 <- partial(int2n13.final, pred.var="Metformin", grid.resolution = 10)
int2n13_pd_12 <- partial(int2n13.final, pred.var="Carbamazepine", grid.resolution = 10)
int2n13_pd_13 <- partial(int2n13.final, pred.var="VDR_Activation", grid.resolution = 10)
int2n13_pd_14 <- partial(int2n13.final, pred.var="ESRRA_Inhibition", grid.resolution = 10)
int2n13_pd_15 <- partial(int2n13.final, pred.var="Lidocaine", grid.resolution = 10)
int2n13_pd_16 <- partial(int2n13.final, pred.var="RXRA_Activation", grid.resolution = 10)
int2n13_pd_17 <- partial(int2n13.final, pred.var="ESR1_Activation", grid.resolution = 10)
int2n13_pd_18 <- partial(int2n13.final, pred.var="Venlafaxine", grid.resolution = 10)
int2n13_pd_19 <- partial(int2n13.final, pred.var="Desvenlafaxine", grid.resolution = 10)
int2n13_pd_20 <- partial(int2n13.final, pred.var="SLC-Inhibition", grid.resolution = 10)
int2n13_pd_21 <- partial(int2n13.final, pred.var="SNRIs_WWI", grid.resolution = 10)

int2n13_pdp_1 <- autoplot(int2n13_pd_1, size=1.2) + theme_classic() + xlab("PGR_Inhibition") + ylab("Effect") +
  theme(text=element_text(size=20)) + ggtitle("CYP 2N13")
int2n13_pdp_2 <- autoplot(int2n13_pd_2, size=1.2) + theme_classic() + xlab("SCN-Inhibition") + ylab("Effect") +
  theme(text=element_text(size=20)) + ggtitle("CYP 2N13")
int2n13_pdp_3 <- autoplot(int2n13_pd_3, size=1.2) + theme_classic() + xlab("PAHs") + ylab("Effect") +
  theme(text=element_text(size=20)) + ggtitle("CYP 2N13")
int2n13_pdp_4 <- autoplot(int2n13_pd_4, size=1.2) + theme_classic() + xlab("testosterone_Inhibition") + ylab("Effect") +
  theme(text=element_text(size=20)) + ggtitle("CYP 2N13")
int2n13_pdp_5 <- autoplot(int2n13_pd_5, size=1.2) + theme_classic() + xlab("NR1H4_Activation") + ylab("Effect") +
  theme(text=element_text(size=20)) + ggtitle("CYP 2N13")
int2n13_pdp_6 <- autoplot(int2n13_pd_6, size=1.2) + theme_classic() + xlab("TSPO_Inhibition") + ylab("Effect") +
  theme(text=element_text(size=20)) + ggtitle("CYP 2N13")
int2n13_pdp_7 <- autoplot(int2n13_pd_7, size=1.2) + theme_classic() + xlab("NR1I3_Binding") + ylab("Effect") +
  theme(text=element_text(size=20)) + ggtitle("CYP 2N13")
int2n13_pdp_8 <- autoplot(int2n13_pd_8, size=1.2) + theme_classic() + xlab("NFE2L2_Activation") + ylab("Effect") +
  theme(text=element_text(size=20)) + ggtitle("CYP 2N13")
int2n13_pdp_9 <- autoplot(int2n13_pd_9, size=1.2) + theme_classic() + xlab("Tramadol") + ylab("Effect") +
  theme(text=element_text(size=20)) + ggtitle("CYP 2N13")
int2n13_pdp_10 <- autoplot(int2n13_pd_10, size=1.2) + theme_classic() + xlab("17alpha-hydroxypregnenolone_Inhibition") + ylab("Effect") +
  theme(text=element_text(size=20)) + ggtitle("CYP 2N13")
int2n13_pdp_11 <- autoplot(int2n13_pd_11, size=1.2) + theme_classic() + xlab("Metformin") + ylab("Effect") +
  theme(text=element_text(size=20)) + ggtitle("CYP 2N13")
int2n13_pdp_12 <- autoplot(int2n13_pd_12, size=1.2) + theme_classic() + xlab("Carbamazepine") + ylab("Effect") +
  theme(text=element_text(size=20)) + ggtitle("CYP 2N13")
int2n13_pdp_13 <- autoplot(int2n13_pd_13, size=1.2) + theme_classic() + xlab("VDR_Activation") + ylab("Effect") +
  theme(text=element_text(size=20)) + ggtitle("CYP 2N13")
int2n13_pdp_14 <- autoplot(int2n13_pd_14, size=1.2) + theme_classic() + xlab("ESRRA_Inhibition") + ylab("Effect") +
  theme(text=element_text(size=20)) + ggtitle("CYP 2N13")
int2n13_pdp_15 <- autoplot(int2n13_pd_15, size=1.2) + theme_classic() + xlab("Lidocaine") + ylab("Effect") +
  theme(text=element_text(size=20)) + ggtitle("CYP 2N13")
int2n13_pdp_16 <- autoplot(int2n13_pd_16, size=1.2) + theme_classic() + xlab("RXRA_Activation") + ylab("Effect") +
  theme(text=element_text(size=20)) + ggtitle("CYP 2N13")
int2n13_pdp_17 <- autoplot(int2n13_pd_17, size=1.2) + theme_classic() + xlab("ESR1_Activation") + ylab("Effect") +
  theme(text=element_text(size=20)) + ggtitle("CYP 2N13")
int2n13_pdp_18 <- autoplot(int2n13_pd_18, size=1.2) + theme_classic() + xlab("Venlafaxine") + ylab("Effect") +
  theme(text=element_text(size=20)) + ggtitle("CYP 2N13")
int2n13_pdp_19 <- autoplot(int2n13_pd_19, size=1.2) + theme_classic() + xlab("Desvenlafaxine") + ylab("Effect") +
  theme(text=element_text(size=20)) + ggtitle("CYP 2N13")
int2n13_pdp_20 <- autoplot(int2n13_pd_20, size=1.2) + theme_classic() + xlab("SLC-Inhibition") + ylab("Effect") +
  theme(text=element_text(size=20)) + ggtitle("CYP 2N13")
int2n13_pdp_21 <- autoplot(int2n13_pd_20, size=1.2) + theme_classic() + xlab("SNRIs_WWI") + ylab("Effect") +
  theme(text=element_text(size=20)) + ggtitle("CYP 2N13")

#whole plot
invivo_XM_intest_1 <- grid.arrange(int1a1_pdp_1, int1a1_pdp_2, int1a1_pdp_3, int1a1_pdp_4, int1a1_pdp_5, int1a1_pdp_6, 
                                   int1a1_pdp_7, int1a1_pdp_8, int1a1_pdp_9, int1a1_pdp_10,int1a1_pdp_11,
                                   int1a1_pdp_12, int1a1_pdp_13, 
                                   int3A_pdp_1, int3A_pdp_2,int2ad6_pdp_1, int2ad6_pdp_2, 
                                   int2ad6_pdp_3, int2ad6_pdp_4, int2ad6_pdp_5, int2ad6_pdp_6, int2ad6_pdp_7, int2ad6_pdp_8, 
                                   int2ad6_pdp_9, int2ad6_pdp_10, int2ad6_pdp_11, int2ad6_pdp_12, int2ad6_pdp_13,
                                   int2n13_pdp_1, int2n13_pdp_2, int2n13_pdp_3, int2n13_pdp_4, int2n13_pdp_5, 
                                   int2n13_pdp_6, int2n13_pdp_7, int2n13_pdp_8, int2n13_pdp_9, int2n13_pdp_10, int2n13_pdp_11, 
                                   int2n13_pdp_12, int2n13_pdp_13, int2n13_pdp_14,
                                   int2n13_pdp_15, int2n13_pdp_16, int2n13_pdp_17, int2n13_pdp_18,
                                   int2n13_pdp_19, int2n13_pdp_20, int2n13_pdp_21, ncol = 6)

invivo_XM_intest_2 <- grid.arrange(ugtint1a1_pdp_1, ugtint1a1_pdp_2, ugtint1a1_pdp_3, ugtint1a1_pdp_4, ugtint1a1_pdp_5,
                                   ugtint1a1_pdp_6, ugtint1a1_pdp_7, ugtint1a1_pdp_8, ugtint1a1_pdp_9, ugtint1a1_pdp_10, 
                                   ugtint1a1_pdp_11, ugtint1a1_pdp_12, ugtint1a1_pdp_13, ugtint1a1_pdp_14, ugtint1a1_pdp_15,
                                   ugtint1a1_pdp_16, ugtint1a1_pdp_17, ugtint1a1_pdp_18, ugtint1a1_pdp_19, ugtint1a1_pdp_20, 
                                   ugtint1a1_pdp_21, ugtint1a1_pdp_22, ugtint1a1_pdp_23, ugtint1a1_pdp_24,
                                   ugtint1a1_pdp_25, ugtint1a1_pdp_26, ugtint1a1_pdp_27, ugtint1a1_pdp_28,
                                   ugtint1a1_pdp_29, ugtint1a1_pdp_30, ugtint1a1_pdp_31, 
                                   ugtint1a1_pdp_33, ugtint1a1_pdp_34, ugtint1a1_pdp_35, ugtint1a1_pdp_36,
                                   ugtint1a1_pdp_37, ugtint1a1_pdp_38, ugtint1a1_pdp_39, ugtint1a1_pdp_40,
                                   ugtint1a1_pdp_41, ugtint1a1_pdp_42, ugtint1a1_pdp_43, ugtint1a1_pdp_44, 
                                   ugtint1a1_pdp_45, ugtint1a1_pdp_46, ugtint1a1_pdp_47, ugtint1a1_pdp_48, 
                                   ugtint1a1_pdp_49, ugtint1a1_pdp_50, ugtint1a1_pdp_51, ugtint1a1_pdp_52, 
                                   ugtint1a1_pdp_53, ugtint1a1_pdp_54, ugtint1a1_pdp_55, ugtint1a1_pdp_56, 
                                   ugtint1a1_pdp_57, ugtint1a1_pdp_58, ugtint1a1_pdp_59, ugtint1a1_pdp_60, 
                                   ugtint1a1_pdp_61, ugtint1a1_pdp_62, ugtint1a1_pdp_63, ugtint1a1_pdp_64,
                                   ugtint1a1_pdp_65, ugtint1a1_pdp_66, ugtint1a1_pdp_67, ugtint1a1_pdp_68,
                                   ugtint1a1_pdp_69, ugtint1a1_pdp_70, ugtint1a1_pdp_71, ugtint1a1_pdp_72, 
                                   ugtint1a1_pdp_73, ugtint1a1_pdp_74, ugtint1a1_pdp_75, ugtint1a1_pdp_76, 
                                   ugtint1a1_pdp_77, ugtint1a1_pdp_78, ugtint1a1_pdp_79, ugtint1a1_pdp_80, 
                                   ugtint1a1_pdp_81, ugtint1a1_pdp_82, ugtint1a1_pdp_83, ugtint1a1_pdp_84, 
                                   ugtint1a1_pdp_85, ugtint1a1_pdp_86, ugtint1a1_pdp_87, ugtint1a1_pdp_88,
                                   ncol = 6)


ggsave("invivo_XM_intestine_pdp_1.jpeg", invivo_XM_intest_1, height = 40, width = 30)
ggsave("invivo_XM_intestine_pdp_2.jpeg", invivo_XM_intest_2, height = 40, width = 30)
beep()

#in vivo XM hepatic####
cyp_1a1_liver_2017_final <- cyp_1a1_liver_2017 %>% dplyr::select(Effect, ESR2_Activation)
set.seed(3)
hep1a1_2017.final <- cforest(Effect~., data=cyp_1a1_liver_2017_final, control=cforest_unbiased(mtry=1, ntree=5000))

hep1a1_2017_pd_1 <- partial(hep1a1_2017.final, pred.var="ESR2_Activation", grid.resolution = 10)
hep1a1_2017_pdp_1 <- autoplot(hep1a1_2017_pd_1, size=1.2) + theme_classic() + xlab("ESR2_Activation") + ylab("Effect") +
  theme(text=element_text(size=20)) + ggtitle("CYP 1A1 (2017)")


cyp_1a1_2018_final <- cyp_1a1_2018 %>% dplyr::select(Effect, "Methyl-1H-benzotriazole", "Tributyl phosphate", "Caffeine",
                                                     "Scn1a_Inhibition", "HDAC1_Inhibition", "Cholesterol",
                                                     "AR_Activation", "Cotinine", "Phenanthrene",
                                                     "Industrial_PPCP_Narcotics", "FOS|JUN_Activation", "Acetaminophen",
                                                     "ADORA2_ADRA_HRH-Inhibition", "ESR1_Activation", "JUN_Activation",
                                                     "AR_Inhibition", "Maoa_Inhibition", "Sterols",
                                                     "Anthraquinone", "1-Methylnaphthalene", "Antibiotics",
                                                     "NR1I3_Activation", "PPCPs_Pesticides_Plasticizers", "Metolachlor",
                                                     "NR1H4_Activation", "RXRB_Activation", "beta-Sitosterol",
                                                     "ESR2_Inhibition", "N,N-diethyl-meta-toluamide", "Metoprolol", "Atrazine"
)
set.seed(3)
hep1a1_2018.final <- cforest(Effect~., data=cyp_1a1_2018_final, control=cforest_unbiased(mtry=11, ntree=5000))

hep1a1_2018_pd_1 <- partial(hep1a1_2018.final, pred.var="Methyl-1H-benzotriazole", grid.resolution = 10)
hep1a1_2018_pdp_1 <- autoplot(hep1a1_2018_pd_1, size=1.2) + theme_classic() + xlab("Methyl-1H-benzotriazole") + ylab("Effect") +
  theme(text=element_text(size=20)) + ggtitle("CYP 1A1 (2018)")
hep1a1_2018_pd_2 <- partial(hep1a1_2018.final, pred.var="Tributyl phosphate", grid.resolution = 10)
hep1a1_2018_pdp_2 <- autoplot(hep1a1_2018_pd_2, size=1.2) + theme_classic() + xlab("Tributyl phosphate") + ylab("Effect") +
  theme(text=element_text(size=20)) + ggtitle("CYP 1A1 (2018)")
hep1a1_2018_pd_3 <- partial(hep1a1_2018.final, pred.var="Caffeine", grid.resolution = 10)
hep1a1_2018_pdp_3 <- autoplot(hep1a1_2018_pd_3, size=1.2) + theme_classic() + xlab("Caffeine") + ylab("Effect") +
  theme(text=element_text(size=20)) + ggtitle("CYP 1A1 (2018)")
hep1a1_2018_pd_4 <- partial(hep1a1_2018.final, pred.var="Scn1a_Inhibition", grid.resolution = 10)
hep1a1_2018_pdp_4 <- autoplot(hep1a1_2018_pd_4, size=1.2) + theme_classic() + xlab("Scn1a_Inhibition") + ylab("Effect") +
  theme(text=element_text(size=20)) + ggtitle("CYP 1A1 (2018)")
hep1a1_2018_pd_5 <- partial(hep1a1_2018.final, pred.var="HDAC1_Inhibition", grid.resolution = 10)
hep1a1_2018_pdp_5 <- autoplot(hep1a1_2018_pd_5, size=1.2) + theme_classic() + xlab("HDAC1_Inhibition") + ylab("Effect") +
  theme(text=element_text(size=20)) + ggtitle("CYP 1A1 (2018)")
hep1a1_2018_pd_6 <- partial(hep1a1_2018.final, pred.var="Cholesterol", grid.resolution = 10)
hep1a1_2018_pdp_6 <- autoplot(hep1a1_2018_pd_6, size=1.2) + theme_classic() + xlab("Cholesterol") + ylab("Effect") +
  theme(text=element_text(size=20)) + ggtitle("CYP 1A1 (2018)")
hep1a1_2018_pd_7 <- partial(hep1a1_2018.final, pred.var="AR_Activation", grid.resolution = 10)
hep1a1_2018_pdp_7 <- autoplot(hep1a1_2018_pd_7, size=1.2) + theme_classic() + xlab("AR_Activation") + ylab("Effect") +
  theme(text=element_text(size=20)) + ggtitle("CYP 1A1 (2018)")
hep1a1_2018_pd_8 <- partial(hep1a1_2018.final, pred.var="Cotinine", grid.resolution = 10)
hep1a1_2018_pdp_8 <- autoplot(hep1a1_2018_pd_8, size=1.2) + theme_classic() + xlab("Cotinine") + ylab("Effect") +
  theme(text=element_text(size=20)) + ggtitle("CYP 1A1 (2018)")
hep1a1_2018_pd_9 <- partial(hep1a1_2018.final, pred.var="Phenanthrene", grid.resolution = 10)
hep1a1_2018_pdp_9 <- autoplot(hep1a1_2018_pd_9, size=1.2) + theme_classic() + xlab("Phenanthrene") + ylab("Effect") +
  theme(text=element_text(size=20)) + ggtitle("CYP 1A1 (2018)")
hep1a1_2018_pd_10 <- partial(hep1a1_2018.final, pred.var="Industrial_PPCP_Narcotics", grid.resolution = 10)
hep1a1_2018_pdp_10 <- autoplot(hep1a1_2018_pd_10, size=1.2) + theme_classic() + xlab("Industrial_PPCP_Narcotics") + ylab("Effect") +
  theme(text=element_text(size=20)) + ggtitle("CYP 1A1 (2018)")
hep1a1_2018_pd_11 <- partial(hep1a1_2018.final, pred.var="FOS|JUN_Activation", grid.resolution = 10)
hep1a1_2018_pdp_11 <- autoplot(hep1a1_2018_pd_11, size=1.2) + theme_classic() + xlab("FOS|JUN_Activation") + ylab("Effect") +
  theme(text=element_text(size=20)) + ggtitle("CYP 1A1 (2018)")
hep1a1_2018_pd_12 <- partial(hep1a1_2018.final, pred.var="Acetaminophen", grid.resolution = 10)
hep1a1_2018_pdp_12 <- autoplot(hep1a1_2018_pd_12, size=1.2) + theme_classic() + xlab("Acetaminophen") + ylab("Effect") +
  theme(text=element_text(size=20)) + ggtitle("CYP 1A1 (2018)")
hep1a1_2018_pd_13 <- partial(hep1a1_2018.final, pred.var="ADORA2_ADRA_HRH-Inhibition", grid.resolution = 10)
hep1a1_2018_pdp_13 <- autoplot(hep1a1_2018_pd_13, size=1.2) + theme_classic() + xlab("ADORA2_ADRA_HRH-Inhibition") + ylab("Effect") +
  theme(text=element_text(size=20)) + ggtitle("CYP 1A1 (2018)")
hep1a1_2018_pd_14 <- partial(hep1a1_2018.final, pred.var="ESR1_Activation", grid.resolution = 10)
hep1a1_2018_pdp_14 <- autoplot(hep1a1_2018_pd_14, size=1.2) + theme_classic() + xlab("ESR1_Activation") + ylab("Effect") +
  theme(text=element_text(size=20)) + ggtitle("CYP 1A1 (2018)")
hep1a1_2018_pd_15 <- partial(hep1a1_2018.final, pred.var="JUN_Activation", grid.resolution = 10)
hep1a1_2018_pdp_15 <- autoplot(hep1a1_2018_pd_15, size=1.2) + theme_classic() + xlab("JUN_Activation") + ylab("Effect") +
  theme(text=element_text(size=20)) + ggtitle("CYP 1A1 (2018)")
hep1a1_2018_pd_16 <- partial(hep1a1_2018.final, pred.var="AR_Inhibition", grid.resolution = 10)
hep1a1_2018_pdp_16 <- autoplot(hep1a1_2018_pd_16, size=1.2) + theme_classic() + xlab("AR_Inhibition") + ylab("Effect") +
  theme(text=element_text(size=20)) + ggtitle("CYP 1A1 (2018)")
hep1a1_2018_pd_17 <- partial(hep1a1_2018.final, pred.var="Maoa_Inhibition", grid.resolution = 10)
hep1a1_2018_pdp_17 <- autoplot(hep1a1_2018_pd_17, size=1.2) + theme_classic() + xlab("Maoa_Inhibition") + ylab("Effect") +
  theme(text=element_text(size=20)) + ggtitle("CYP 1A1 (2018)")
hep1a1_2018_pd_18 <- partial(hep1a1_2018.final, pred.var="Sterols", grid.resolution = 10)
hep1a1_2018_pdp_18 <- autoplot(hep1a1_2018_pd_18, size=1.2) + theme_classic() + xlab("Sterols") + ylab("Effect") +
  theme(text=element_text(size=20)) + ggtitle("CYP 1A1 (2018)")
hep1a1_2018_pd_19 <- partial(hep1a1_2018.final, pred.var="Anthraquinone", grid.resolution = 10)
hep1a1_2018_pdp_19 <- autoplot(hep1a1_2018_pd_19, size=1.2) + theme_classic() + xlab("Anthraquinone") + ylab("Effect") +
  theme(text=element_text(size=20)) + ggtitle("CYP 1A1 (2018)")
hep1a1_2018_pd_20 <- partial(hep1a1_2018.final, pred.var="1-Methylnaphthalene", grid.resolution = 10)
hep1a1_2018_pdp_20 <- autoplot(hep1a1_2018_pd_20, size=1.2) + theme_classic() + xlab("1-Methylnaphthalene") + ylab("Effect") +
  theme(text=element_text(size=20)) + ggtitle("CYP 1A1 (2018)")
hep1a1_2018_pd_21 <- partial(hep1a1_2018.final, pred.var="Antibiotics", grid.resolution = 10)
hep1a1_2018_pdp_21 <- autoplot(hep1a1_2018_pd_21, size=1.2) + theme_classic() + xlab("Antibiotics") + ylab("Effect") +
  theme(text=element_text(size=20)) + ggtitle("CYP 1A1 (2018)")
hep1a1_2018_pd_22 <- partial(hep1a1_2018.final, pred.var="NR1I3_Activation", grid.resolution = 10)
hep1a1_2018_pdp_22 <- autoplot(hep1a1_2018_pd_22, size=1.2) + theme_classic() + xlab("NR1I3_Activation") + ylab("Effect") +
  theme(text=element_text(size=20)) + ggtitle("CYP 1A1 (2018)")
hep1a1_2018_pd_23 <- partial(hep1a1_2018.final, pred.var="PPCPs_Pesticides_Plasticizers", grid.resolution = 10)
hep1a1_2018_pdp_23 <- autoplot(hep1a1_2018_pd_23, size=1.2) + theme_classic() + xlab("PPCPs_Pesticides_Plasticizers") + ylab("Effect") +
  theme(text=element_text(size=20)) + ggtitle("CYP 1A1 (2018)")
hep1a1_2018_pd_24 <- partial(hep1a1_2018.final, pred.var="Metolachlor", grid.resolution = 10)
hep1a1_2018_pdp_24 <- autoplot(hep1a1_2018_pd_24, size=1.2) + theme_classic() + xlab("Metolachlor") + ylab("Effect") +
  theme(text=element_text(size=20)) + ggtitle("CYP 1A1 (2018)")
hep1a1_2018_pd_25 <- partial(hep1a1_2018.final, pred.var="NR1H4_Activation", grid.resolution = 10)
hep1a1_2018_pdp_25 <- autoplot(hep1a1_2018_pd_25, size=1.2) + theme_classic() + xlab("NR1H4_Activation") + ylab("Effect") +
  theme(text=element_text(size=20)) + ggtitle("CYP 1A1 (2018)")
hep1a1_2018_pd_26 <- partial(hep1a1_2018.final, pred.var="RXRB_Activation", grid.resolution = 10)
hep1a1_2018_pdp_26 <- autoplot(hep1a1_2018_pd_26, size=1.2) + theme_classic() + xlab("RXRB_Activation") + ylab("Effect") +
  theme(text=element_text(size=20)) + ggtitle("CYP 1A1 (2018)")
hep1a1_2018_pd_27 <- partial(hep1a1_2018.final, pred.var="beta-Sitosterol", grid.resolution = 10)
hep1a1_2018_pdp_27 <- autoplot(hep1a1_2018_pd_27, size=1.2) + theme_classic() + xlab("beta-Sitosterol") + ylab("Effect") +
  theme(text=element_text(size=20)) + ggtitle("CYP 1A1 (2018)")
hep1a1_2018_pd_28 <- partial(hep1a1_2018.final, pred.var="ESR2_Inhibition", grid.resolution = 10)
hep1a1_2018_pdp_28 <- autoplot(hep1a1_2018_pd_28, size=1.2) + theme_classic() + xlab("ESR2_Inhibition") + ylab("Effect") +
  theme(text=element_text(size=20)) + ggtitle("CYP 1A1 (2018)")
hep1a1_2018_pd_29 <- partial(hep1a1_2018.final, pred.var="N,N-diethyl-meta-toluamide", grid.resolution = 10)
hep1a1_2018_pdp_29 <- autoplot(hep1a1_2018_pd_29, size=1.2) + theme_classic() + xlab("N,N-diethyl-meta-toluamide") + ylab("Effect") +
  theme(text=element_text(size=20)) + ggtitle("CYP 1A1 (2018)")
hep1a1_2018_pd_30 <- partial(hep1a1_2018.final, pred.var="Metoprolol", grid.resolution = 10)
hep1a1_2018_pdp_30 <- autoplot(hep1a1_2018_pd_30, size=1.2) + theme_classic() + xlab("Metoprolol") + ylab("Effect") +
  theme(text=element_text(size=20)) + ggtitle("CYP 1A1 (2018)")
hep1a1_2018_pd_31 <- partial(hep1a1_2018.final, pred.var="Atrazine", grid.resolution = 10)
hep1a1_2018_pdp_31 <- autoplot(hep1a1_2018_pd_31, size=1.2) + theme_classic() + xlab("Atrazine") + ylab("Effect") +
  theme(text=element_text(size=20)) + ggtitle("CYP 1A1 (2018)")


cyp3a_liver_final <- cyp3a_liver %>% dplyr::select(Effect, Indole, Carbazole)
set.seed(3)
cf.final <- cforest(Effect~., data=cyp3a_liver_final, control=cforest_unbiased(mtry=1, ntree=5000))

cyp3ahep_pd_1 <- partial(cf.final, pred.var="Indole", grid.resolution = 10)
cyp3ahep_pdp_1 <- autoplot(cyp3ahep_pd_1, size=1.2) + theme_classic() + xlab("Indole") + ylab("Effect") +
  theme(text=element_text(size=20)) + ggtitle("CYP 3A")
cyp3ahep_pd_2 <- partial(cf.final, pred.var="Carbazole", grid.resolution = 10)
cyp3ahep_pdp_2 <- autoplot(cyp3ahep_pd_2, size=1.2) + theme_classic() + xlab("Indole") + ylab("Effect") +
  theme(text=element_text(size=20)) + ggtitle("CYP 3A")


UGT_liver_final <- UGT_liver %>% dplyr::select(Effect, CYP2B6_Activation, Tramadol, ESRRA_Inhibition, 
                                               Scn1a_Inhibition, ESR1_Activation, Oprk1_Binding,
                                               OP_Fire_Retardants, "ADORA2_ADRA_HRH-Inhibition",
                                               PGR_Inhibition, NR1I2_Binding, NR1I2_Activation,
                                               PPARG_Activation, "Tris(2-butoxyethyl)phosphate") 
set.seed(3)
ugthep.final <- cforest(Effect~., data=UGT_liver_final, control=cforest_unbiased(mtry=4, ntree=5000))

ugthep_pd_1 <- partial(ugthep.final, pred.var="CYP2B6_Activation", grid.resolution = 10)
ugthep_pdp_1 <- autoplot(ugthep_pd_1, size=1.2) + theme_classic() + xlab("CYP2B6_Activation") + ylab("Effect") +
  theme(text=element_text(size=20)) + ggtitle("UGT 1A1")
ugthep_pd_2 <- partial(ugthep.final, pred.var="Tramadol", grid.resolution = 10)
ugthep_pdp_2 <- autoplot(ugthep_pd_2, size=1.2) + theme_classic() + xlab("Tramadol") + ylab("Effect") +
  theme(text=element_text(size=20)) + ggtitle("UGT 1A1")
ugthep_pd_3 <- partial(ugthep.final, pred.var="ESRRA_Inhibition", grid.resolution = 10)
ugthep_pdp_3 <- autoplot(ugthep_pd_3, size=1.2) + theme_classic() + xlab("ESRRA_Inhibition") + ylab("Effect") +
  theme(text=element_text(size=20)) + ggtitle("UGT 1A1")
ugthep_pd_4 <- partial(ugthep.final, pred.var="Scn1a_Inhibition", grid.resolution = 10)
ugthep_pdp_4 <- autoplot(ugthep_pd_4, size=1.2) + theme_classic() + xlab("ESRRA_Inhibition") + ylab("Effect") +
  theme(text=element_text(size=20)) + ggtitle("UGT 1A1")
ugthep_pd_5 <- partial(ugthep.final, pred.var="ESR1_Activation", grid.resolution = 10)
ugthep_pdp_5 <- autoplot(ugthep_pd_5, size=1.2) + theme_classic() + xlab("ESRRA_Inhibition") + ylab("Effect") +
  theme(text=element_text(size=20)) + ggtitle("UGT 1A1")
ugthep_pd_6 <- partial(ugthep.final, pred.var="Oprk1_Binding", grid.resolution = 10)
ugthep_pdp_6 <- autoplot(ugthep_pd_6, size=1.2) + theme_classic() + xlab("Oprk1_Binding") + ylab("Effect") +
  theme(text=element_text(size=20)) + ggtitle("UGT 1A1")
ugthep_pd_7 <- partial(ugthep.final, pred.var="OP_Fire_Retardants", grid.resolution = 10)
ugthep_pdp_7 <- autoplot(ugthep_pd_7, size=1.2) + theme_classic() + xlab("OP_Fire_Retardants") + ylab("Effect") +
  theme(text=element_text(size=20)) + ggtitle("UGT 1A1")
ugthep_pd_8 <- partial(ugthep.final, pred.var="ADORA2_ADRA_HRH-Inhibition", grid.resolution = 10)
ugthep_pdp_8 <- autoplot(ugthep_pd_8, size=1.2) + theme_classic() + xlab("ADORA2_ADRA_HRH-Inhibition") + ylab("Effect") +
  theme(text=element_text(size=20)) + ggtitle("UGT 1A1")
ugthep_pd_9 <- partial(ugthep.final, pred.var="PGR_Inhibition", grid.resolution = 10)
ugthep_pdp_9 <- autoplot(ugthep_pd_9, size=1.2) + theme_classic() + xlab("PGR_Inhibition") + ylab("Effect") +
  theme(text=element_text(size=20)) + ggtitle("UGT 1A1")
ugthep_pd_10 <- partial(ugthep.final, pred.var="NR1I2_Binding", grid.resolution = 10)
ugthep_pdp_10 <- autoplot(ugthep_pd_10, size=1.2) + theme_classic() + xlab("NR1I2_Binding") + ylab("Effect") +
  theme(text=element_text(size=20)) + ggtitle("UGT 1A1")
ugthep_pd_11 <- partial(ugthep.final, pred.var="NR1I2_Activation", grid.resolution = 10)
ugthep_pdp_11 <- autoplot(ugthep_pd_11, size=1.2) + theme_classic() + xlab("NR1I2_Activation") + ylab("Effect") +
  theme(text=element_text(size=20)) + ggtitle("UGT 1A1")
ugthep_pd_11 <- partial(ugthep.final, pred.var="PPARG_Activation", grid.resolution = 10)
ugthep_pdp_11 <- autoplot(ugthep_pd_11, size=1.2) + theme_classic() + xlab("PPARG_Activation") + ylab("Effect") +
  theme(text=element_text(size=20)) + ggtitle("UGT 1A1")
ugthep_pd_12 <- partial(ugthep.final, pred.var="Tris(2-butoxyethyl)phosphate", grid.resolution = 10)
ugthep_pdp_12 <- autoplot(ugthep_pd_12, size=1.2) + theme_classic() + xlab("Tris(2-butoxyethyl)phosphate") + ylab("Effect") +
  theme(text=element_text(size=20)) + ggtitle("UGT 1A1")



invivo_hep_ppdp <- grid.arrange(hep1a1_2017_pdp_1, hep1a1_2018_pdp_1, hep1a1_2018_pdp_2, hep1a1_2018_pdp_3, 
                                hep1a1_2018_pdp_4, hep1a1_2018_pdp_5, hep1a1_2018_pdp_6, hep1a1_2018_pdp_7,
                                hep1a1_2018_pdp_8, hep1a1_2018_pdp_9, hep1a1_2018_pdp_10, hep1a1_2018_pdp_11,
                                hep1a1_2018_pdp_12, hep1a1_2018_pdp_13, hep1a1_2018_pdp_14, hep1a1_2018_pdp_15,
                                hep1a1_2018_pdp_16, hep1a1_2018_pdp_17, hep1a1_2018_pdp_18, hep1a1_2018_pdp_19,
                                hep1a1_2018_pdp_20, hep1a1_2018_pdp_21, hep1a1_2018_pdp_22, hep1a1_2018_pdp_23,
                                hep1a1_2018_pdp_24, hep1a1_2018_pdp_25, hep1a1_2018_pdp_26, hep1a1_2018_pdp_27,
                                hep1a1_2018_pdp_28, hep1a1_2018_pdp_29, hep1a1_2018_pdp_30, hep1a1_2018_pdp_31,
                                cyp3ahep_pdp_1, cyp3ahep_pdp_2, 
                                ugthep_pdp_1, ugthep_pdp_2, ugthep_pdp_3, ugthep_pdp_4, ugthep_pdp_5, ugthep_pdp_6,
                                ugthep_pdp_7, ugthep_pdp_8, ugthep_pdp_9, ugthep_pdp_10, ugthep_pdp_11, ugthep_pdp_12,
                                ncol = 6) 
ggsave("invivo_hep_ppdp.jpeg", invivo_hep_ppdp, height = 30, width = 20)
beep()
