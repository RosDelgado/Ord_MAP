
###########################################
###########################################
#
#  EXPERIMENTAL PHASE   (Section 5)
# 
################################################################################
#
# PROCEDURE 4): tuning Random Forest 
#
#########################################
############   train function, from caret library
############
############ This function sets up a grid of tuning parameters for a number of 
############ classification and regression routines, fits each model and 
############ calculates a resampling based performance measure.
############ Uses "trainControl" argument from caret
###########
########################################

library(caret)

source("mat_square.R")

####   FUNCTIONS

# Distance between two classes: |i-j| # used for MAE computation
distance <- function(M) {  # M squared confusion matrix
  r <- nrow(M)
  d <- as.data.frame(which(M !='NA', arr.ind = T))
  d <- abs(d$row-d$col)
  d <- matrix(d, nrow = r)
  d
}


## MAE: Mean Absolute Error
# Given a confusion matrix C, the function returns the Mean Absolute Error of C
MAE <- function(C) {
  # Check parameters
  stopifnot("The argument of the function should be a matrix" = class(C)[1] %in% c('matrix','table') ) # == 'matrix'
  stopifnot("The argument of the function should be a square matrix" = nrow(C) == ncol(C))
  stopifnot("The argument of the function is not a confusion matrix (not all its elements are integer numbers)" = all(C == floor(C)))
  
  # Variables
  r <- nrow(C)
  N <- sum(C)
  
  # Cost matrix m.ij
  m <- distance(C)/N
  
  # MAE
  MAE <- sum(C*m)
  MAE
}



#########.  DATASETS


# DATASET (a) World Values Survey (WVS) Dataset (carData package)
#
# Dataset from the World Values Surveys for Australia, 
# Norway, Sweden, and the United States from "carData" package in R.
#
# Poverty is the multi-class ordered dependent variable with categories
# ‘Too Little’, ‘About Right’ and ‘Too Much’. 
# We have the following five independent variables
#
# Religion: member of a religion -no or yes
# Degree: held a university degree -no or yes
# Country: Australia, Norway, Sweden or the USA
# Age: age (years)
# Gender: male or female
#
library(carData)

data(WVS) 
str(WVS)    # 5381 obs. 6 variables.
head(WVS)  # output to predict: "poverty", with 3 categories
table(WVS$poverty) # Too Little, About Right, Too Much
#
levels.poverty.encod<-sort(unique(as.numeric(WVS$poverty)))
levels(WVS$poverty)<-levels.poverty.encod
table(WVS$poverty) # 1, 2, 3
#
features<-c(2:6)
################################################################################
####### preparing for k-fold cross-validation with k=10
#######

L=dim(WVS)[1]

set.seed(12345)
fold<-sample(c(1:10),L,replace=TRUE)
table(fold)

training<-list()
test<-list()

for (i in 1:10)
{test[[i]]<-WVS[which(fold==i),]
training[[i]]<-WVS[-which(fold==i),]}


categories<-names(table(WVS$poverty))
r<-length(categories)

################################################################################
####### Error functions for summaryFunction argument, trainControl function
#######

standard.MAE.ord<-function(data,lev=1:length(categories),model=NULL)
  # data=dataframe with columns "obs" and "pred" of character/factor type
  # lev=character string with outcome factor levels
{Conf.mat<- mat.square(table(as.numeric(data$pred),as.numeric(data$obs)),lev)
value<-MAE(Conf.mat)
c(MAE.ord=value)
}



################################################################################
################################################################################
########## caret::train. Resampling method: cross-validation
################################################################################

mtry<-sqrt(ncol(WVS)-4)
ntree<-3

fitControl.Accuracy <- trainControl(
  method = "cv",
  number = 3,
  search="random")



fitControl.MAE <- trainControl(
  method = "cv",
  number = 3,
  ## Evaluate performance using 
  ## the following function
  summaryFunction = standard.MAE.ord)




tuned.rf.caret.Accuracy.WVS<-list()
tuned.rf.caret.MAE.WVS<-list()


for (i in 1:10)
{ set.seed(12345)
  tuned.rf.caret.Accuracy.WVS[[i]] <- caret::train(training[[i]][ ,features], 
                                                         training[[i]][ ,1], 
                                                         method="rf", 
                                                         metric="Accuracy",
                                                         tuneLength=10,
                                                         trControl=fitControl.Accuracy)
  
  print(i)}




##############

for (i in 1:10)
{ set.seed(12345)
  tuned.rf.caret.MAE.WVS[[i]] <- caret::train(training[[i]][ ,features], 
                                                    training[[i]][ ,1], 
                                                    method="rf",
                                                    metric="MAE.ord",
                                                    maximize=FALSE,
                                                    tuneLength=10,
                                                    trControl=fitControl.MAE)
  
  print(i)}


################################################################################
#
# Predictions of any of the test set
# 
#

pred.test.rf.caret.Accuracy.WVS<-list()
for (i in 1:10)
{pred.test.rf.caret.Accuracy.WVS[[i]] <- predict(tuned.rf.caret.Accuracy.WVS[[i]],test[[i]][,features],type="prob")
}


pred.test.rf.caret.MAE.WVS<-list()
for (i in 1:10)
{pred.test.rf.caret.MAE.WVS[[i]] <- predict(tuned.rf.caret.MAE.WVS[[i]],test[[i]][,features],type="prob")
}

###

################################################################################
### predictions with standard MAP criterion:
###

pred.MAP.rf.caret.Accuracy.WVS<-list()
for (i in 1:10)
{pred.MAP.rf.caret.Accuracy.WVS[[i]]<-vector()
for (j in 1:dim(pred.test.rf.caret.Accuracy.WVS[[i]])[1])
{pred.MAP.rf.caret.Accuracy.WVS[[i]][j]<-categories[which.max(pred.test.rf.caret.Accuracy.WVS[[i]][j,])]
}
}

pred.MAP.rf.caret.MAE.WVS<-list()
for (i in 1:10)
{pred.MAP.rf.caret.MAE.WVS[[i]]<-vector()
for (j in 1:dim(pred.test.rf.caret.MAE.WVS[[i]])[1])
{pred.MAP.rf.caret.MAE.WVS[[i]][j]<-categories[which.max(pred.test.rf.caret.MAE.WVS[[i]][j,])]
}
}


# confusion matrices


conf.mat.MAP.rf.caret.Accuracy.WVS<-list()
for(i in 1:10)
{
  conf.mat.MAP.rf.caret.Accuracy.WVS[[i]]<-mat.square(table(pred.MAP.rf.caret.Accuracy.WVS[[i]],test[[i]][[1]]),categories)
}


conf.mat.MAP.rf.caret.MAE.WVS<-list()
for(i in 1:10)
{
  conf.mat.MAP.rf.caret.MAE.WVS[[i]]<-mat.square(table(pred.MAP.rf.caret.MAE.WVS[[i]],test[[i]][[1]]),categories)
}



# MAE 

MAE.MAP.rf.caret.Accuracy.WVS<-vector()
for (i in 1:10)
{
  MAE.MAP.rf.caret.Accuracy.WVS[i]<-MAE(conf.mat.MAP.rf.caret.Accuracy.WVS[[i]])
}


MAE.MAP.rf.caret.MAE.WVS<-vector()
for (i in 1:10)
{
  MAE.MAP.rf.caret.MAE.WVS[i]<-MAE(conf.mat.MAP.rf.caret.MAE.WVS[[i]])
}


###

################################################################################
### predictions with Ord.MAP criterion:
###

pred.Ord.MAP.rf.caret.Accuracy.WVS<-list()
for (i in 1:10)
{pred.Ord.MAP.rf.caret.Accuracy.WVS[[i]]<-vector()
for (j in 1:dim(pred.test.rf.caret.Accuracy.WVS[[i]])[1])
{h<-min(which(cumsum(as.vector(pred.test.rf.caret.Accuracy.WVS[[i]][j,]))>=0.5))
pred.Ord.MAP.rf.caret.Accuracy.WVS[[i]][j]<-categories[h]
}
}


pred.Ord.MAP.rf.caret.MAE.WVS<-list()
for (i in 1:10)
{pred.Ord.MAP.rf.caret.MAE.WVS[[i]]<-vector()
for (j in 1:dim(pred.test.rf.caret.MAE.WVS[[i]])[1])
{h<-min(which(cumsum(as.vector(pred.test.rf.caret.MAE.WVS[[i]][j,]))>=0.5))
pred.Ord.MAP.rf.caret.MAE.WVS[[i]][j]<-categories[h]
}
}


###
# confusion matrices


conf.mat.Ord.MAP.rf.caret.Accuracy.WVS<-list()
for(i in 1:10)
{
  conf.mat.Ord.MAP.rf.caret.Accuracy.WVS[[i]]<-mat.square(table(pred.Ord.MAP.rf.caret.Accuracy.WVS[[i]],test[[i]][[1]]),categories)
}



conf.mat.Ord.MAP.rf.caret.MAE.WVS<-list()
for(i in 1:10)
{
  conf.mat.Ord.MAP.rf.caret.MAE.WVS[[i]]<-mat.square(table(pred.Ord.MAP.rf.caret.MAE.WVS[[i]],test[[i]][[1]]),categories)
}


# MAE 


MAE.Ord.MAP.rf.caret.Accuracy.WVS<-vector()
for (i in 1:10)
{
  MAE.Ord.MAP.rf.caret.Accuracy.WVS[i]<-MAE(conf.mat.Ord.MAP.rf.caret.Accuracy.WVS[[i]])
}


MAE.Ord.MAP.rf.caret.MAE.WVS<-vector()
for (i in 1:10)
{
  MAE.Ord.MAP.rf.caret.MAE.WVS[i]<-MAE(conf.mat.Ord.MAP.rf.caret.MAE.WVS[[i]])
}



###
################################################################################
###. Comparison MAP vs Ord.MAP
###
###.  Tests MAE

p.val.test.MAE.MAP.greater.MAE.ORD.MAP.rf.caret.Accuracy.WVS<-wilcox.test(MAE.MAP.rf.caret.Accuracy.WVS,
                                                                          MAE.Ord.MAP.rf.caret.Accuracy.WVS,paired = TRUE,alternative="greater")$p.value
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.rf.caret.Accuracy.WVS.less<-wilcox.test(MAE.MAP.rf.caret.Accuracy.WVS,
                                                                          MAE.Ord.MAP.rf.caret.Accuracy.WVS,paired = TRUE,alternative="less")$p.value


p.val.test.MAE.MAP.greater.MAE.ORD.MAP.rf.caret.MAE.WVS<-wilcox.test(MAE.MAP.rf.caret.MAE.WVS,
                                                                     MAE.Ord.MAP.rf.caret.MAE.WVS,paired = TRUE,alternative="greater")$p.value
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.rf.caret.MAE.WVS.less<-wilcox.test(MAE.MAP.rf.caret.MAE.WVS,
                                                                     MAE.Ord.MAP.rf.caret.MAE.WVS,paired = TRUE,alternative="less")$p.value




p.val.test.MAE.MAP.greater.MAE.ORD.MAP.rf.caret.Accuracy.WVS
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.rf.caret.Accuracy.WVS.less

p.val.test.MAE.MAP.greater.MAE.ORD.MAP.rf.caret.MAE.WVS
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.rf.caret.MAE.WVS.less


###############################
###############################
###   Boxplots   MAE

MAE.rf.caret.WVS.results<-as.data.frame(cbind(MAE.MAP.rf.caret.Accuracy.WVS-MAE.Ord.MAP.rf.caret.Accuracy.WVS,
                                                    MAE.MAP.rf.caret.MAE.WVS-MAE.Ord.MAP.rf.caret.MAE.WVS))

boxplot(MAE.rf.caret.WVS.results,
        main=paste("(a) World Values Surveys (WVS) dataset"), 
        xlab="Summary Function", ylab="MAE increment (MAP - Ord.MAP)",names=c("Accuracy", "MAE"),
        col="white")

abline(h = 0, col = 2, lwd = 2, lty=2)
# Dots
stripchart(MAE.rf.caret.WVS.results,
           #method = "jitter",
           pch = 19,
           col = 1:2,
           vertical = TRUE,
           add = TRUE)




#
##
###
##
#


# DATASET (b) Wine (package "ordinal")
#
# The wine data set is adopted from Randall(1989) and from a factorial experiment 
# on factors determining the bitterness of wine. 
# Two treatment factors (temperature and contact) each have two
# levels. Temperature and contact between juice and skins can be controlled when cruching grapes
# during wine production. Nine judges each assessed wine from two bottles from each of the four
# treatment conditions, hence there are 72 observations in all.
#

library(ordinal)

data(wine) 
str(wine)    # 72 obs. 6 variables.
head(wine)  # output to predict: "rating", with 5 categories
table(wine$rating)   # 1:5

#
#


################################################################################
####### preparing for k-fold cross-validation with k=10
#######

L=dim(wine)[1]

set.seed(12345)
fold<-sample(c(1:10),L,replace=TRUE)
table(fold)

training<-list()
test<-list()

for (i in 1:10)
{test[[i]]<-wine[which(fold==i),]
training[[i]]<-wine[-which(fold==i),]}


categories<-names(table(wine$rating))
r<-length(categories)

features<-c(3,4)


################################################################################
################################################################################
########## caret::train. Resampling method: cross-validation
################################################################################

mtry<-sqrt(ncol(wine)-4)
ntree<-3

fitControl.Accuracy <- trainControl(
  method = "cv",
  number = 3,
  search="random")



fitControl.MAE <- trainControl(
  method = "cv",
  number = 3,
  ## Evaluate performance using 
  ## the following function
  summaryFunction = standard.MAE.ord)




tuned.rf.caret.Accuracy.wine<-list()
tuned.rf.caret.MAE.wine<-list()


for (i in 1:10)
{ set.seed(12345)
  tuned.rf.caret.Accuracy.wine[[i]] <- caret::train(training[[i]][ ,features], 
                                                   training[[i]][ ,2], 
                                                   method="rf", 
                                                   metric="Accuracy",
                                                   tuneLength=10,
                                                   trControl=fitControl.Accuracy)
  
  print(i)}




##############

for (i in 1:10)
{ set.seed(12345)
  tuned.rf.caret.MAE.wine[[i]] <- caret::train(training[[i]][ ,features], 
                                              training[[i]][ ,2], 
                                              method="rf",
                                              metric="MAE.ord",
                                              maximize=FALSE,
                                              tuneLength=10,
                                              trControl=fitControl.MAE)
  
  print(i)}


################################################################################
#
# Predictions of any of the test set
# 
#

pred.test.rf.caret.Accuracy.wine<-list()
for (i in 1:10)
{pred.test.rf.caret.Accuracy.wine[[i]] <- predict(tuned.rf.caret.Accuracy.wine[[i]],test[[i]][,features],type="prob")
}


pred.test.rf.caret.MAE.wine<-list()
for (i in 1:10)
{pred.test.rf.caret.MAE.wine[[i]] <- predict(tuned.rf.caret.MAE.wine[[i]],test[[i]][,features],type="prob")
}

###

################################################################################
### predictions with standard MAP criterion:
###

pred.MAP.rf.caret.Accuracy.wine<-list()
for (i in 1:10)
{pred.MAP.rf.caret.Accuracy.wine[[i]]<-vector()
for (j in 1:dim(pred.test.rf.caret.Accuracy.wine[[i]])[1])
{pred.MAP.rf.caret.Accuracy.wine[[i]][j]<-categories[which.max(pred.test.rf.caret.Accuracy.wine[[i]][j,])]
}
}

pred.MAP.rf.caret.MAE.wine<-list()
for (i in 1:10)
{pred.MAP.rf.caret.MAE.wine[[i]]<-vector()
for (j in 1:dim(pred.test.rf.caret.MAE.wine[[i]])[1])
{pred.MAP.rf.caret.MAE.wine[[i]][j]<-categories[which.max(pred.test.rf.caret.MAE.wine[[i]][j,])]
}
}


# confusion matrices


conf.mat.MAP.rf.caret.Accuracy.wine<-list()
for(i in 1:10)
{
  conf.mat.MAP.rf.caret.Accuracy.wine[[i]]<-mat.square(table(pred.MAP.rf.caret.Accuracy.wine[[i]],test[[i]][[2]]),categories)
}


conf.mat.MAP.rf.caret.MAE.wine<-list()
for(i in 1:10)
{
  conf.mat.MAP.rf.caret.MAE.wine[[i]]<-mat.square(table(pred.MAP.rf.caret.MAE.wine[[i]],test[[i]][[2]]),categories)
}



# MAE 

MAE.MAP.rf.caret.Accuracy.wine<-vector()
for (i in 1:10)
{
  MAE.MAP.rf.caret.Accuracy.wine[i]<-MAE(conf.mat.MAP.rf.caret.Accuracy.wine[[i]])
}


MAE.MAP.rf.caret.MAE.wine<-vector()
for (i in 1:10)
{
  MAE.MAP.rf.caret.MAE.wine[i]<-MAE(conf.mat.MAP.rf.caret.MAE.wine[[i]])
}


###

################################################################################
### predictions with Ord.MAP criterion:
###

pred.Ord.MAP.rf.caret.Accuracy.wine<-list()
for (i in 1:10)
{pred.Ord.MAP.rf.caret.Accuracy.wine[[i]]<-vector()
for (j in 1:dim(pred.test.rf.caret.Accuracy.wine[[i]])[1])
{h<-min(which(cumsum(as.vector(pred.test.rf.caret.Accuracy.wine[[i]][j,]))>=0.5))
pred.Ord.MAP.rf.caret.Accuracy.wine[[i]][j]<-categories[h]
}
}


pred.Ord.MAP.rf.caret.MAE.wine<-list()
for (i in 1:10)
{pred.Ord.MAP.rf.caret.MAE.wine[[i]]<-vector()
for (j in 1:dim(pred.test.rf.caret.MAE.wine[[i]])[1])
{h<-min(which(cumsum(as.vector(pred.test.rf.caret.MAE.wine[[i]][j,]))>=0.5))
pred.Ord.MAP.rf.caret.MAE.wine[[i]][j]<-categories[h]
}
}


###
# confusion matrices


conf.mat.Ord.MAP.rf.caret.Accuracy.wine<-list()
for(i in 1:10)
{
  conf.mat.Ord.MAP.rf.caret.Accuracy.wine[[i]]<-mat.square(table(pred.Ord.MAP.rf.caret.Accuracy.wine[[i]],test[[i]][[2]]),categories)
}



conf.mat.Ord.MAP.rf.caret.MAE.wine<-list()
for(i in 1:10)
{
  conf.mat.Ord.MAP.rf.caret.MAE.wine[[i]]<-mat.square(table(pred.Ord.MAP.rf.caret.MAE.wine[[i]],test[[i]][[2]]),categories)
}


# MAE 


MAE.Ord.MAP.rf.caret.Accuracy.wine<-vector()
for (i in 1:10)
{
  MAE.Ord.MAP.rf.caret.Accuracy.wine[i]<-MAE(conf.mat.Ord.MAP.rf.caret.Accuracy.wine[[i]])
}


MAE.Ord.MAP.rf.caret.MAE.wine<-vector()
for (i in 1:10)
{
  MAE.Ord.MAP.rf.caret.MAE.wine[i]<-MAE(conf.mat.Ord.MAP.rf.caret.MAE.wine[[i]])
}



###
################################################################################
###. Comparison MAP vs Ord.MAP
###
###.  Tests MAE

p.val.test.MAE.MAP.greater.MAE.ORD.MAP.rf.caret.Accuracy.wine<-wilcox.test(MAE.MAP.rf.caret.Accuracy.wine,
                                                                          MAE.Ord.MAP.rf.caret.Accuracy.wine,paired = TRUE,alternative="greater")$p.value
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.rf.caret.Accuracy.wine.less<-wilcox.test(MAE.MAP.rf.caret.Accuracy.wine,
                                                                               MAE.Ord.MAP.rf.caret.Accuracy.wine,paired = TRUE,alternative="less")$p.value


p.val.test.MAE.MAP.greater.MAE.ORD.MAP.rf.caret.MAE.wine<-wilcox.test(MAE.MAP.rf.caret.MAE.wine,
                                                                     MAE.Ord.MAP.rf.caret.MAE.wine,paired = TRUE,alternative="greater")$p.value
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.rf.caret.MAE.wine.less<-wilcox.test(MAE.MAP.rf.caret.MAE.wine,
                                                                          MAE.Ord.MAP.rf.caret.MAE.wine,paired = TRUE,alternative="less")$p.value




p.val.test.MAE.MAP.greater.MAE.ORD.MAP.rf.caret.Accuracy.wine
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.rf.caret.Accuracy.wine.less

p.val.test.MAE.MAP.greater.MAE.ORD.MAP.rf.caret.MAE.wine
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.rf.caret.MAE.wine.less


###############################
###############################
###   Boxplots   MAE

MAE.rf.caret.wine.results<-as.data.frame(cbind(MAE.MAP.rf.caret.Accuracy.wine-MAE.Ord.MAP.rf.caret.Accuracy.wine,
                                              MAE.MAP.rf.caret.MAE.wine-MAE.Ord.MAP.rf.caret.MAE.wine))

boxplot(MAE.rf.caret.wine.results,
        main=paste("(b) Wine dataset"), 
        xlab="Summary Function", ylab="MAE increment (MAP - Ord.MAP)",names=c("Accuracy", "MAE"),
        col="white")

abline(h = 0, col = 2, lwd = 2, lty=2)
# Dots
stripchart(MAE.rf.caret.WVS.results,
           #method = "jitter",
           pch = 19,
           col = 1:2,
           vertical = TRUE,
           add = TRUE)




#
##
###
##
#


# DATASET (c) Hearth (ordinalForest package)
#
# Dataset heart included in the package ordinalForest.
# This data includes 294 patients undergoing angiography at the Hungarian Institute of Cardiology in
# Budapest between 1983 and 1987.
# Is a  data frame with 294 observations, ten covariates and one ordinal target variable
# age. numeric. Age in years
# sex. factor. Sex (1 = male; 0 = female)
# chest_pain. factor. Chest pain type (1 = typical angina; 2 = atypical angina; 3 = non-anginal
#                                       pain; 4 = asymptomatic)
# trestbps. numeric. Resting blood pressure (in mm Hg on admission to the hospital)
# chol. numeric. Serum cholestoral in mg/dl
# fbs. factor. Fasting blood sugar > 120 mg/dl (1 = true; 0 = false)
# restecg. factor. Resting electrocardiographic results (1 = having ST-T wave abnormality (T
#                    wave inversions and/or ST elevation or depression of > 0.05 mV); 0 = normal)
# thalach. numeric. Maximum heart rate achieved
# exang. factor. Exercise induced angina (1 = yes; 0 = no)
# oldpeak. numeric. ST depression induced by exercise relative to rest
# Class. factor. Ordinal target variable - severity of coronary artery disease (determined using
# angiograms) (1 = no disease; 2 = degree 1; 3 = degree 2; 4 = degree 3; 5 = degree 4)


library(ordinalForest)

data(hearth)
str(hearth)

# Inspect the data:
table(hearth$Class)
dim(hearth)

head(hearth) 


################################################################################
####### preparing for k-fold cross-validation with k=10
#######

L=dim(hearth)[1]

set.seed(12345)
fold<-sample(c(1:10),L,replace=TRUE)
table(fold)

training<-list()
test<-list()

for (i in 1:10)
{test[[i]]<-hearth[which(fold==i),]
training[[i]]<-hearth[-which(fold==i),]}


categories<-names(table(hearth$Class))
r<-length(categories)


features<-c(1:10)

################################################################################
################################################################################
########## caret::train. Resampling method: cross-validation
################################################################################

mtry<-sqrt(ncol(hearth)-4)
ntree<-3

fitControl.Accuracy <- trainControl(
  method = "cv",
  number = 3,
  search="random")



fitControl.MAE <- trainControl(
  method = "cv",
  number = 3,
  ## Evaluate performance using 
  ## the following function
  summaryFunction = standard.MAE.ord)




tuned.rf.caret.Accuracy<-list()
tuned.rf.caret.MAE<-list()


for (i in 1:10)
{ set.seed(12345)
  tuned.rf.caret.Accuracy[[i]] <- caret::train(training[[i]][ ,features], 
                                                    training[[i]][ ,11], 
                                                    method="rf", 
                                                    metric="Accuracy",
                                                    tuneLength=10,
                                                    trControl=fitControl.Accuracy)
  
  print(i)}




##############

for (i in 1:10)
{ set.seed(12345)
  tuned.rf.caret.MAE[[i]] <- caret::train(training[[i]][ ,features], 
                                               training[[i]][ ,11], 
                                               method="rf",
                                               metric="MAE.ord",
                                               maximize=FALSE,
                                               tuneLength=10,
                                               trControl=fitControl.MAE)
  
  print(i)}


################################################################################
#
# Predictions of any of the test set
# 
#

pred.test.rf.caret.Accuracy<-list()
for (i in 1:10)
{pred.test.rf.caret.Accuracy[[i]] <- predict(tuned.rf.caret.Accuracy[[i]],test[[i]][,features],type="prob")
}


pred.test.rf.caret.MAE<-list()
for (i in 1:10)
{pred.test.rf.caret.MAE[[i]] <- predict(tuned.rf.caret.MAE[[i]],test[[i]][,features],type="prob")
}

###

################################################################################
### predictions with standard MAP criterion:
###

pred.MAP.rf.caret.Accuracy<-list()
for (i in 1:10)
{pred.MAP.rf.caret.Accuracy[[i]]<-vector()
for (j in 1:dim(pred.test.rf.caret.Accuracy[[i]])[1])
{pred.MAP.rf.caret.Accuracy[[i]][j]<-categories[which.max(pred.test.rf.caret.Accuracy[[i]][j,])]
}
}

pred.MAP.rf.caret.MAE<-list()
for (i in 1:10)
{pred.MAP.rf.caret.MAE[[i]]<-vector()
for (j in 1:dim(pred.test.rf.caret.MAE[[i]])[1])
{pred.MAP.rf.caret.MAE[[i]][j]<-categories[which.max(pred.test.rf.caret.MAE[[i]][j,])]
}
}


# confusion matrices


conf.mat.MAP.rf.caret.Accuracy<-list()
for(i in 1:10)
{
  conf.mat.MAP.rf.caret.Accuracy[[i]]<-mat.square(table(pred.MAP.rf.caret.Accuracy[[i]],test[[i]][[11]]),categories)
}


conf.mat.MAP.rf.caret.MAE<-list()
for(i in 1:10)
{
  conf.mat.MAP.rf.caret.MAE[[i]]<-mat.square(table(pred.MAP.rf.caret.MAE[[i]],test[[i]][[11]]),categories)
}



# MAE 

MAE.MAP.rf.caret.Accuracy<-vector()
for (i in 1:10)
{
  MAE.MAP.rf.caret.Accuracy[i]<-MAE(conf.mat.MAP.rf.caret.Accuracy[[i]])
}


MAE.MAP.rf.caret.MAE<-vector()
for (i in 1:10)
{
  MAE.MAP.rf.caret.MAE[i]<-MAE(conf.mat.MAP.rf.caret.MAE[[i]])
}


###

################################################################################
### predictions with Ord.MAP criterion:
###

pred.Ord.MAP.rf.caret.Accuracy<-list()
for (i in 1:10)
{pred.Ord.MAP.rf.caret.Accuracy[[i]]<-vector()
for (j in 1:dim(pred.test.rf.caret.Accuracy[[i]])[1])
{h<-min(which(cumsum(as.vector(pred.test.rf.caret.Accuracy[[i]][j,]))>=0.5))
pred.Ord.MAP.rf.caret.Accuracy[[i]][j]<-categories[h]
}
}


pred.Ord.MAP.rf.caret.MAE<-list()
for (i in 1:10)
{pred.Ord.MAP.rf.caret.MAE[[i]]<-vector()
for (j in 1:dim(pred.test.rf.caret.MAE[[i]])[1])
{h<-min(which(cumsum(as.vector(pred.test.rf.caret.MAE[[i]][j,]))>=0.5))
pred.Ord.MAP.rf.caret.MAE[[i]][j]<-categories[h]
}
}


###
# confusion matrices


conf.mat.Ord.MAP.rf.caret.Accuracy<-list()
for(i in 1:10)
{
  conf.mat.Ord.MAP.rf.caret.Accuracy[[i]]<-mat.square(table(pred.Ord.MAP.rf.caret.Accuracy[[i]],test[[i]][[11]]),categories)
}



conf.mat.Ord.MAP.rf.caret.MAE<-list()
for(i in 1:10)
{
  conf.mat.Ord.MAP.rf.caret.MAE[[i]]<-mat.square(table(pred.Ord.MAP.rf.caret.MAE[[i]],test[[i]][[11]]),categories)
}


# MAE 


MAE.Ord.MAP.rf.caret.Accuracy<-vector()
for (i in 1:10)
{
  MAE.Ord.MAP.rf.caret.Accuracy[i]<-MAE(conf.mat.Ord.MAP.rf.caret.Accuracy[[i]])
}


MAE.Ord.MAP.rf.caret.MAE<-vector()
for (i in 1:10)
{
  MAE.Ord.MAP.rf.caret.MAE[i]<-MAE(conf.mat.Ord.MAP.rf.caret.MAE[[i]])
}



###
################################################################################
###. Comparison MAP vs Ord.MAP
###
###.  Tests MAE

p.val.test.MAE.MAP.greater.MAE.ORD.MAP.rf.caret.Accuracy<-wilcox.test(MAE.MAP.rf.caret.Accuracy,
                                                                           MAE.Ord.MAP.rf.caret.Accuracy,paired = TRUE,alternative="greater")$p.value
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.rf.caret.Accuracy.less<-wilcox.test(MAE.MAP.rf.caret.Accuracy,
                                                                                MAE.Ord.MAP.rf.caret.Accuracy,paired = TRUE,alternative="less")$p.value


p.val.test.MAE.MAP.greater.MAE.ORD.MAP.rf.caret.MAE<-wilcox.test(MAE.MAP.rf.caret.MAE,
                                                                      MAE.Ord.MAP.rf.caret.MAE,paired = TRUE,alternative="greater")$p.value
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.rf.caret.MAE.less<-wilcox.test(MAE.MAP.rf.caret.MAE,
                                                                           MAE.Ord.MAP.rf.caret.MAE,paired = TRUE,alternative="less")$p.value




p.val.test.MAE.MAP.greater.MAE.ORD.MAP.rf.caret.Accuracy
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.rf.caret.Accuracy.less

p.val.test.MAE.MAP.greater.MAE.ORD.MAP.rf.caret.MAE
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.rf.caret.MAE.less


###############################
###############################
###   Boxplots   MAE

MAE.rf.caret.results<-as.data.frame(cbind(MAE.MAP.rf.caret.Accuracy-MAE.Ord.MAP.rf.caret.Accuracy,
                                               MAE.MAP.rf.caret.MAE-MAE.Ord.MAP.rf.caret.MAE))

boxplot(MAE.rf.caret.results,
        main=paste("(c) Hearth dataset"), 
        xlab="Summary Function", ylab="MAE increment (MAP - Ord.MAP)",names=c("Accuracy", "MAE"),
        col="white")

abline(h = 0, col = 2, lwd = 2, lty=2)
# Dots
stripchart(MAE.rf.caret.WVS.results,
           #method = "jitter",
           pch = 19,
           col = 1:2,
           vertical = TRUE,
           add = TRUE)


#
##
###
##
#



# DATASET (d) Parkinson (Output: V5 and V6, both with 4, 5 and 6 categories) 
#
parkinsons <- read.csv("parkinsons_updrs.data", header=FALSE)
#View(parkinsons)
str(parkinsons)

parkinson<-parkinsons[-1,] # first row are column names

### objective: predict "V5: motor_UPDRS" and "V6: total_UPDRS" scores from the 16 voice measures + sex + age + test_time
### The unified (total)Parkinson’s disease rating scale (UPDRS) reflects the presence and severity of symptoms 
### (but does not quantify their underlying causes). For untreated patients, it spans the range
### 0–176, with 0 representing healthy state and 176 representing total disabilities, and consists of three sections: 
### 1) mentation,behavior, and mood; 
### 2) activities of daily living; and 
### 3) motor.
### The motor UPDRS ranges from 0 to 108, with 0 denoting symptom free and 108 denoting severe motor impairment, 
### and encompasses tasks such as speech, facial expression, tremor, and rigidity. 
### Speech has two explicit headings and ranges between 0 and 8 with 8 being unintelligible.

### From: "Accurate Telemonitoring of Parkinson’s Disease Progression by Noninvasive Speech Tests" 
### Athanasios Tsanas, Max A. Little, Member, IEEE, Patrick E. McSharry, Senior Member, IEEE, and Lorraine O. Ramig
### IEEE TRANSACTIONS ON BIOMEDICAL ENGINEERING, VOL. 57, NO. 4, APRIL 2010


# V1 = patient identification

parkinson<-as.data.frame(lapply(parkinson,as.numeric))

features<-c(7:22)

#### FIRST: OUTPUT VARIABLE TO PREDICT: V5: motor_UPDRS (ranges from 0 to 108)

########## V5 with 5 categories:

summary(parkinson$V5)
cut.points<-c(0,13,18,24,29,1000)

V5.bin <- arules::discretize(parkinson$V5, method = "fixed", breaks=cut.points, infinity=TRUE)
table(V5.bin)
parkinson<-as.data.frame(cbind(parkinson,V5.bin))


levels.V5.int<-names(table(parkinson$V5.bin))
levels(parkinson$V5.bin)<-c("<13","[13,18)","[18,24)","[24,29)",">=29")
levels.V5.ordinal.encod<-sort(unique(as.numeric(parkinson$V5.bin)))

parkinson$V5.bin.num<-as.numeric(parkinson$V5.bin)
table(parkinson$V5.bin)
table(parkinson$V5.bin.num)

parkinson$V5.bin.num.factor<-as.factor(parkinson$V5.bin.num)   # classes in number but factor type 

parkinson$V5.bin.num.factor<-ordered(parkinson$V5.bin.num.factor, levels = c("1", "2", "3", "4","5"))

categories<-names(table(parkinson$V5.bin.num.factor))
r<-length(categories)

########## V5 with 4 categories:

summary(parkinson$V5)
cut.points<-c(0,15,20,30,1000)

V5.bin <- arules::discretize(parkinson$V5, method = "fixed", breaks=cut.points, infinity=TRUE)
table(V5.bin)
parkinson<-as.data.frame(cbind(parkinson,V5.bin))


levels.V5.int<-names(table(parkinson$V5.bin))
levels(parkinson$V5.bin)<-c("<15","[15,20)","[20,30)",">=30")
levels.V5.ordinal.encod<-sort(unique(as.numeric(parkinson$V5.bin)))


parkinson$V5.bin.num<-as.numeric(parkinson$V5.bin)
table(parkinson$V5.bin)
table(parkinson$V5.bin.num)

parkinson$V5.bin.num.factor<-as.factor(parkinson$V5.bin.num)   # classes in number but factor type 

parkinson$V5.bin.num.factor<-ordered(parkinson$V5.bin.num.factor, levels = c("1", "2", "3", "4"))

categories<-names(table(parkinson$V5.bin.num.factor))
r<-length(categories)


########## V5 with 6 categories:

summary(parkinson$V5)
cut.points<-c(0,10,15,20,25,30,1000)

V5.bin <- arules::discretize(parkinson$V5, method = "fixed", breaks=cut.points, infinity=TRUE)
table(V5.bin)
parkinson<-as.data.frame(cbind(parkinson,V5.bin))


levels.V5.int<-names(table(parkinson$V5.bin))
levels(parkinson$V5.bin)<-c("<10","[10,15)","[15,20)","[20,25)","[25,30)",">=30")
levels.V5.ordinal.encod<-sort(unique(as.numeric(parkinson$V5.bin)))

parkinson$V5.bin.num<-as.numeric(parkinson$V5.bin)
table(parkinson$V5.bin)
table(parkinson$V5.bin.num)

parkinson$V5.bin.num.factor<-as.factor(parkinson$V5.bin.num)   # classes in number but factor type 

parkinson$V5.bin.num.factor<-ordered(parkinson$V5.bin.num.factor, levels = c("1", "2", "3", "4", "5", "6"))

categories<-names(table(parkinson$V5.bin.num.factor))
r<-length(categories)


#### SECOND: OUTPUT VARIABLE TO PREDICT: V6: total_UPDRS (ranges from 0 to 176)

########## V6 with 5 categories:

cut.points.v6<-c(0,20,25,30,40,1000)

V6.bin <- arules::discretize(parkinson$V6, method = "fixed", breaks=cut.points.v6, infinity=TRUE)
table(V6.bin)
parkinson<-as.data.frame(cbind(parkinson,V6.bin))


levels.V6.int<-names(table(parkinson$V6.bin))
levels(parkinson$V6.bin)<-c("<20","[20,25)","[25,30)","[30,40)",">=40")
levels.V6.ordinal.encod<-sort(unique(as.numeric(parkinson$V6.bin)))


parkinson$V6.bin.num<-as.numeric(parkinson$V6.bin)
table(parkinson$V6.bin)
table(parkinson$V6.bin.num)

parkinson$V6.bin.num.factor<-as.factor(parkinson$V6.bin.num)   # classes in number but factor type 

parkinson$V6.bin.num.factor<-ordered(parkinson$V6.bin.num.factor, levels = c("1", "2", "3", "4","5"))


categories<-names(table(parkinson$V6.bin.num.factor))
r<-length(categories)


########## V6 with 4 categories:

cut.points.v6<-c(0,22,30,40,1000)

V6.bin <- arules::discretize(parkinson$V6, method = "fixed", breaks=cut.points.v6, infinity=TRUE)
table(V6.bin)
parkinson<-as.data.frame(cbind(parkinson,V6.bin))


levels.V6.int<-names(table(parkinson$V6.bin))
levels(parkinson$V6.bin)<-c("<22","[22,30)","[30,40)",">=40")
levels.V6.ordinal.encod<-sort(unique(as.numeric(parkinson$V6.bin)))

parkinson$V6.bin.num<-as.numeric(parkinson$V6.bin)
table(parkinson$V6.bin)
table(parkinson$V6.bin.num)

parkinson$V6.bin.num.factor<-as.factor(parkinson$V6.bin.num)   # classes in number but factor type 

parkinson$V6.bin.num.factor<-ordered(parkinson$V6.bin.num.factor, levels = c("1", "2", "3", "4"))


categories<-names(table(parkinson$V6.bin.num.factor))
r<-length(categories)

########## V6 with 6 categories:

cut.points.v6<-c(0,20,25,30,35,45,1000)

V6.bin <- arules::discretize(parkinson$V6, method = "fixed", breaks=cut.points.v6, infinity=TRUE)
table(V6.bin)
parkinson<-as.data.frame(cbind(parkinson,V6.bin))


levels.V6.int<-names(table(parkinson$V6.bin))
levels(parkinson$V6.bin)<-c("<20","[20,25)","[25,30)","[30,35)","[35,45)",">=45")
levels.V6.ordinal.encod<-sort(unique(as.numeric(parkinson$V6.bin)))

parkinson$V6.bin.num<-as.numeric(parkinson$V6.bin)
table(parkinson$V6.bin)
table(parkinson$V6.bin.num)

parkinson$V6.bin.num.factor<-as.factor(parkinson$V6.bin.num)   # classes in number but factor type 

parkinson$V6.bin.num.factor<-ordered(parkinson$V6.bin.num.factor, levels = c("1", "2", "3", "4","5","6"))


categories<-names(table(parkinson$V6.bin.num.factor))
r<-length(categories)



################################################################################
####### preparing for k-fold cross-validation with k=10
#######

N=dim(parkinson)[1]
n=round(N/10)

set.seed(12345)
fold<-sample(c(1:10),N,replace=TRUE)
table(fold)

training<-list()
test<-list()
sub.train<-list()

for (i in 1:10)
{test[[i]]<-parkinson[which(fold==i),]
training[[i]]<-parkinson[-which(fold==i),]}

###############################################################################
################################################################################
########## caret::train. Resampling method: cross-validation
################################################################################

mtry<-sqrt(ncol(parkinson)-4)
ntree<-3

fitControl.Accuracy <- trainControl(
  method = "cv",
  number = 3,
  search="random")



fitControl.MAE <- trainControl(
  method = "cv",
  number = 3,
  ## Evaluate performance using 
  ## the following function
  summaryFunction = standard.MAE.ord)




tuned.rf.caret.Accuracy<-list()
tuned.rf.caret.MAE<-list()


for (i in 1:10)
{ set.seed(12345)
  tuned.rf.caret.Accuracy[[i]] <- caret::train(training[[i]][ ,features], 
                                               training[[i]][ ,25], 
                                               method="rf", 
                                               metric="Accuracy",
                                               tuneLength=10,
                                               trControl=fitControl.Accuracy)
  
  print(i)}




##############

for (i in 1:10)
{ set.seed(12345)
  tuned.rf.caret.MAE[[i]] <- caret::train(training[[i]][ ,features], 
                                          training[[i]][ ,25], 
                                          method="rf",
                                          metric="MAE.ord",
                                          maximize=FALSE,
                                          tuneLength=10,
                                          trControl=fitControl.MAE)
  
  print(i)}


################################################################################
#
# Predictions of any of the test set
# 
#

pred.test.rf.caret.Accuracy<-list()
for (i in 1:10)
{pred.test.rf.caret.Accuracy[[i]] <- predict(tuned.rf.caret.Accuracy[[i]],test[[i]][,features],type="prob")
}


pred.test.rf.caret.MAE<-list()
for (i in 1:10)
{pred.test.rf.caret.MAE[[i]] <- predict(tuned.rf.caret.MAE[[i]],test[[i]][,features],type="prob")
}

###

################################################################################
### predictions with standard MAP criterion:
###

pred.MAP.rf.caret.Accuracy<-list()
for (i in 1:10)
{pred.MAP.rf.caret.Accuracy[[i]]<-vector()
for (j in 1:dim(pred.test.rf.caret.Accuracy[[i]])[1])
{pred.MAP.rf.caret.Accuracy[[i]][j]<-categories[which.max(pred.test.rf.caret.Accuracy[[i]][j,])]
}
}

pred.MAP.rf.caret.MAE<-list()
for (i in 1:10)
{pred.MAP.rf.caret.MAE[[i]]<-vector()
for (j in 1:dim(pred.test.rf.caret.MAE[[i]])[1])
{pred.MAP.rf.caret.MAE[[i]][j]<-categories[which.max(pred.test.rf.caret.MAE[[i]][j,])]
}
}


# confusion matrices


conf.mat.MAP.rf.caret.Accuracy<-list()
for(i in 1:10)
{
  conf.mat.MAP.rf.caret.Accuracy[[i]]<-mat.square(table(pred.MAP.rf.caret.Accuracy[[i]],test[[i]][[25]]),categories)
}


conf.mat.MAP.rf.caret.MAE<-list()
for(i in 1:10)
{
  conf.mat.MAP.rf.caret.MAE[[i]]<-mat.square(table(pred.MAP.rf.caret.MAE[[i]],test[[i]][[25]]),categories)
}



# MAE 

MAE.MAP.rf.caret.Accuracy<-vector()
for (i in 1:10)
{
  MAE.MAP.rf.caret.Accuracy[i]<-MAE(conf.mat.MAP.rf.caret.Accuracy[[i]])
}


MAE.MAP.rf.caret.MAE<-vector()
for (i in 1:10)
{
  MAE.MAP.rf.caret.MAE[i]<-MAE(conf.mat.MAP.rf.caret.MAE[[i]])
}


###

################################################################################
### predictions with Ord.MAP criterion:
###

pred.Ord.MAP.rf.caret.Accuracy<-list()
for (i in 1:10)
{pred.Ord.MAP.rf.caret.Accuracy[[i]]<-vector()
for (j in 1:dim(pred.test.rf.caret.Accuracy[[i]])[1])
{h<-min(which(cumsum(as.vector(pred.test.rf.caret.Accuracy[[i]][j,]))>=0.5))
pred.Ord.MAP.rf.caret.Accuracy[[i]][j]<-categories[h]
}
}


pred.Ord.MAP.rf.caret.MAE<-list()
for (i in 1:10)
{pred.Ord.MAP.rf.caret.MAE[[i]]<-vector()
for (j in 1:dim(pred.test.rf.caret.MAE[[i]])[1])
{h<-min(which(cumsum(as.vector(pred.test.rf.caret.MAE[[i]][j,]))>=0.5))
pred.Ord.MAP.rf.caret.MAE[[i]][j]<-categories[h]
}
}


###
# confusion matrices


conf.mat.Ord.MAP.rf.caret.Accuracy<-list()
for(i in 1:10)
{
  conf.mat.Ord.MAP.rf.caret.Accuracy[[i]]<-mat.square(table(pred.Ord.MAP.rf.caret.Accuracy[[i]],test[[i]][[25]]),categories)
}



conf.mat.Ord.MAP.rf.caret.MAE<-list()
for(i in 1:10)
{
  conf.mat.Ord.MAP.rf.caret.MAE[[i]]<-mat.square(table(pred.Ord.MAP.rf.caret.MAE[[i]],test[[i]][[25]]),categories)
}


# MAE 


MAE.Ord.MAP.rf.caret.Accuracy<-vector()
for (i in 1:10)
{
  MAE.Ord.MAP.rf.caret.Accuracy[i]<-MAE(conf.mat.Ord.MAP.rf.caret.Accuracy[[i]])
}


MAE.Ord.MAP.rf.caret.MAE<-vector()
for (i in 1:10)
{
  MAE.Ord.MAP.rf.caret.MAE[i]<-MAE(conf.mat.Ord.MAP.rf.caret.MAE[[i]])
}



###
################################################################################
###. Comparison MAP vs Ord.MAP
###
###.  Tests MAE

p.val.test.MAE.MAP.greater.MAE.ORD.MAP.rf.caret.Accuracy<-wilcox.test(MAE.MAP.rf.caret.Accuracy,
                                                                      MAE.Ord.MAP.rf.caret.Accuracy,paired = TRUE,alternative="greater")$p.value
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.rf.caret.Accuracy.less<-wilcox.test(MAE.MAP.rf.caret.Accuracy,
                                                                           MAE.Ord.MAP.rf.caret.Accuracy,paired = TRUE,alternative="less")$p.value


p.val.test.MAE.MAP.greater.MAE.ORD.MAP.rf.caret.MAE<-wilcox.test(MAE.MAP.rf.caret.MAE,
                                                                 MAE.Ord.MAP.rf.caret.MAE,paired = TRUE,alternative="greater")$p.value
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.rf.caret.MAE.less<-wilcox.test(MAE.MAP.rf.caret.MAE,
                                                                      MAE.Ord.MAP.rf.caret.MAE,paired = TRUE,alternative="less")$p.value




p.val.test.MAE.MAP.greater.MAE.ORD.MAP.rf.caret.Accuracy
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.rf.caret.Accuracy.less

p.val.test.MAE.MAP.greater.MAE.ORD.MAP.rf.caret.MAE
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.rf.caret.MAE.less


###############################
###############################
###   Boxplots   MAE

MAE.rf.caret.results<-as.data.frame(cbind(MAE.MAP.rf.caret.Accuracy-MAE.Ord.MAP.rf.caret.Accuracy,
                                          MAE.MAP.rf.caret.MAE-MAE.Ord.MAP.rf.caret.MAE))

boxplot(MAE.rf.caret.results,
        main=paste("(d) Parkinson dataset: output V6 (6 classes)"), 
        xlab="Summary Function", ylab="MAE increment (MAP - Ord.MAP)",names=c("Accuracy", "MAE"),
        col="white")

abline(h = 0, col = 2, lwd = 2, lty=2)
# Dots
stripchart(MAE.rf.caret.WVS.results,
           #method = "jitter",
           pch = 19,
           col = 1:2,
           vertical = TRUE,
           add = TRUE)


#
##
###
##
#

# DATASET (e) CES11 (package carData) (Output: importance (of the religion), with 4 categories) 
#

library(carData)

data(CES11)
str(CES11)

df<-CES11[,-1]
str(df)

table(df$importance)

df$Class<-factor(df$importance,order=TRUE)


#
levels.Class.encod<-sort(unique(as.numeric(df$Class)))
levels(df$Class)<-levels.Class.encod
table(df$Class) # 1, 2, 3, 4
#
features<-c(1,3:5,7,8)

categories<-names(table(df$Class))
r<-length(categories)



################################################################################
####### preparing for k-fold cross-validation with k=10
#######

N=dim(df)[1]
n=round(N/10)

set.seed(12345)
fold<-sample(c(1:10),N,replace=TRUE)
table(fold)

training<-list()
test<-list()
sub.train<-list()

for (i in 1:10)
{test[[i]]<-df[which(fold==i),]
training[[i]]<-df[-which(fold==i),]}


################################################################################
################################################################################
########## caret::train. Resampling method: cross-validation
################################################################################

mtry<-sqrt(ncol(df)-4)
ntree<-3

fitControl.Accuracy <- trainControl(
  method = "cv",
  number = 3,
  search="random")



fitControl.MAE <- trainControl(
  method = "cv",
  number = 3,
  ## Evaluate performance using 
  ## the following function
  summaryFunction = standard.MAE.ord)




tuned.rf.caret.Accuracy<-list()
tuned.rf.caret.MAE<-list()


for (i in 1:10)
{ set.seed(12345)
  tuned.rf.caret.Accuracy[[i]] <- caret::train(training[[i]][ ,features], 
                                               training[[i]][ ,9], 
                                               method="rf", 
                                               metric="Accuracy",
                                               tuneLength=10,
                                               trControl=fitControl.Accuracy)
  
  print(i)}




##############

for (i in 1:10)
{ set.seed(12345)
  tuned.rf.caret.MAE[[i]] <- caret::train(training[[i]][ ,features], 
                                          training[[i]][ ,9], 
                                          method="rf",
                                          metric="MAE.ord",
                                          maximize=FALSE,
                                          tuneLength=10,
                                          trControl=fitControl.MAE)
  
  print(i)}


################################################################################
#
# Predictions of any of the test set
# 
#

pred.test.rf.caret.Accuracy<-list()
for (i in 1:10)
{pred.test.rf.caret.Accuracy[[i]] <- predict(tuned.rf.caret.Accuracy[[i]],test[[i]][,features],type="prob")
}


pred.test.rf.caret.MAE<-list()
for (i in 1:10)
{pred.test.rf.caret.MAE[[i]] <- predict(tuned.rf.caret.MAE[[i]],test[[i]][,features],type="prob")
}

###

################################################################################
### predictions with standard MAP criterion:
###

pred.MAP.rf.caret.Accuracy<-list()
for (i in 1:10)
{pred.MAP.rf.caret.Accuracy[[i]]<-vector()
for (j in 1:dim(pred.test.rf.caret.Accuracy[[i]])[1])
{pred.MAP.rf.caret.Accuracy[[i]][j]<-categories[which.max(pred.test.rf.caret.Accuracy[[i]][j,])]
}
}

pred.MAP.rf.caret.MAE<-list()
for (i in 1:10)
{pred.MAP.rf.caret.MAE[[i]]<-vector()
for (j in 1:dim(pred.test.rf.caret.MAE[[i]])[1])
{pred.MAP.rf.caret.MAE[[i]][j]<-categories[which.max(pred.test.rf.caret.MAE[[i]][j,])]
}
}


# confusion matrices


conf.mat.MAP.rf.caret.Accuracy<-list()
for(i in 1:10)
{
  conf.mat.MAP.rf.caret.Accuracy[[i]]<-mat.square(table(pred.MAP.rf.caret.Accuracy[[i]],test[[i]][[9]]),categories)
}


conf.mat.MAP.rf.caret.MAE<-list()
for(i in 1:10)
{
  conf.mat.MAP.rf.caret.MAE[[i]]<-mat.square(table(pred.MAP.rf.caret.MAE[[i]],test[[i]][[9]]),categories)
}



# MAE 

MAE.MAP.rf.caret.Accuracy<-vector()
for (i in 1:10)
{
  MAE.MAP.rf.caret.Accuracy[i]<-MAE(conf.mat.MAP.rf.caret.Accuracy[[i]])
}


MAE.MAP.rf.caret.MAE<-vector()
for (i in 1:10)
{
  MAE.MAP.rf.caret.MAE[i]<-MAE(conf.mat.MAP.rf.caret.MAE[[i]])
}


###

################################################################################
### predictions with Ord.MAP criterion:
###

pred.Ord.MAP.rf.caret.Accuracy<-list()
for (i in 1:10)
{pred.Ord.MAP.rf.caret.Accuracy[[i]]<-vector()
for (j in 1:dim(pred.test.rf.caret.Accuracy[[i]])[1])
{h<-min(which(cumsum(as.vector(pred.test.rf.caret.Accuracy[[i]][j,]))>=0.5))
pred.Ord.MAP.rf.caret.Accuracy[[i]][j]<-categories[h]
}
}


pred.Ord.MAP.rf.caret.MAE<-list()
for (i in 1:10)
{pred.Ord.MAP.rf.caret.MAE[[i]]<-vector()
for (j in 1:dim(pred.test.rf.caret.MAE[[i]])[1])
{h<-min(which(cumsum(as.vector(pred.test.rf.caret.MAE[[i]][j,]))>=0.5))
pred.Ord.MAP.rf.caret.MAE[[i]][j]<-categories[h]
}
}


###
# confusion matrices


conf.mat.Ord.MAP.rf.caret.Accuracy<-list()
for(i in 1:10)
{
  conf.mat.Ord.MAP.rf.caret.Accuracy[[i]]<-mat.square(table(pred.Ord.MAP.rf.caret.Accuracy[[i]],test[[i]][[9]]),categories)
}



conf.mat.Ord.MAP.rf.caret.MAE<-list()
for(i in 1:10)
{
  conf.mat.Ord.MAP.rf.caret.MAE[[i]]<-mat.square(table(pred.Ord.MAP.rf.caret.MAE[[i]],test[[i]][[9]]),categories)
}


# MAE 


MAE.Ord.MAP.rf.caret.Accuracy<-vector()
for (i in 1:10)
{
  MAE.Ord.MAP.rf.caret.Accuracy[i]<-MAE(conf.mat.Ord.MAP.rf.caret.Accuracy[[i]])
}


MAE.Ord.MAP.rf.caret.MAE<-vector()
for (i in 1:10)
{
  MAE.Ord.MAP.rf.caret.MAE[i]<-MAE(conf.mat.Ord.MAP.rf.caret.MAE[[i]])
}



###
################################################################################
###. Comparison MAP vs Ord.MAP
###
###.  Tests MAE

p.val.test.MAE.MAP.greater.MAE.ORD.MAP.rf.caret.Accuracy<-wilcox.test(MAE.MAP.rf.caret.Accuracy,
                                                                      MAE.Ord.MAP.rf.caret.Accuracy,paired = TRUE,alternative="greater")$p.value
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.rf.caret.Accuracy.less<-wilcox.test(MAE.MAP.rf.caret.Accuracy,
                                                                           MAE.Ord.MAP.rf.caret.Accuracy,paired = TRUE,alternative="less")$p.value


p.val.test.MAE.MAP.greater.MAE.ORD.MAP.rf.caret.MAE<-wilcox.test(MAE.MAP.rf.caret.MAE,
                                                                 MAE.Ord.MAP.rf.caret.MAE,paired = TRUE,alternative="greater")$p.value
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.rf.caret.MAE.less<-wilcox.test(MAE.MAP.rf.caret.MAE,
                                                                      MAE.Ord.MAP.rf.caret.MAE,paired = TRUE,alternative="less")$p.value




p.val.test.MAE.MAP.greater.MAE.ORD.MAP.rf.caret.Accuracy
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.rf.caret.Accuracy.less

p.val.test.MAE.MAP.greater.MAE.ORD.MAP.rf.caret.MAE
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.rf.caret.MAE.less


###############################
###############################
###   Boxplots   MAE

MAE.rf.caret.results<-as.data.frame(cbind(MAE.MAP.rf.caret.Accuracy-MAE.Ord.MAP.rf.caret.Accuracy,
                                          MAE.MAP.rf.caret.MAE-MAE.Ord.MAP.rf.caret.MAE))

boxplot(MAE.rf.caret.results,
        main=paste("(e) CES11 dataset"), 
        xlab="Summary Function", ylab="MAE increment (MAP - Ord.MAP)",names=c("Accuracy", "MAE"),
        col="white")

abline(h = 0, col = 2, lwd = 2, lty=2)
# Dots
stripchart(MAE.rf.caret.WVS.results,
           #method = "jitter",
           pch = 19,
           col = 1:2,
           vertical = TRUE,
           add = TRUE)


#
##
###
##
#


