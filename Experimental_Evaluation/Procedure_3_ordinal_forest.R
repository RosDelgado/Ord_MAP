
###########################################
###########################################
#
#  EXPERIMENTAL PHASE  (Section 5) 
# 
################################################################################
#
# PROCEDURE 3): Ordinal Forest (OF) method allows ordinal regression with high-dimensional 
# and low-dimensional data. After having constructed an OF prediction rule using 
# a training dataset, it can be used to predict the values of the ordinal target 
# variable for new observations. Moreover, by means of the (permutation-based) 
# variable importance measure of OF, it is also possible to rank the covariates 
# with respect to their importances in the prediction of the values of the ordinal target variable.
# OF is presented in Hornung (2020).

# ordinal random forest using ordfor function from the ordinalForest package. 
#


library(ordinalForest)

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
table(WVS$poverty)   # Too Little, About Right, Too Much
#
#
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
#
# Construct OF prediction rule using the training dataset (default 
# perffunction = "probability" corresponding to the 
# (negative) ranked probability score as performance function):

v.a<-1000   # values for the hyperparameter nsets (default=1000)
v.b<-seq(from=8, to=12, by=1)    # values for the hyperparameter nbest (default=10)



model.ord.forest<-list()
for (k in 1:length(v.b))
{ model.ord.forest[[k]]<-list()
for (i in 1:10)
{model.ord.forest[[k]][[i]] <- ordfor(depvar="poverty", data=training[[i]], nsets=v.a, nbest=v.b[k],ntreeperdiv=100, 
                                      ntreefinal=1000, npermtrial=100,perffunction = "probability")
print("i")
print(i)
}
print("k")
print(k)
}


################################################################################
#
# Predictions of any of the test set
# 
#

pred.test.ord.forest<-list()
for (k in 1:length(v.b))
{ pred.test.ord.forest[[k]]<-list()
for (i in 1:10)
{pred.test.ord.forest[[k]][[i]] <- predict(model.ord.forest[[k]][[i]],test[[i]][,-1])
print("i")
print(i)
}
print("k")
print(k)
}

#pred.test.ord.forest[[k]][[i]]$ypred
#pred.test.ord.forest[[k]][[i]]$classprobs


################################################################################
### predictions with standard MAP criterion:
###

pred.MAP.ord.forest<-list()
for(k in 1:length(v.b))
{ pred.MAP.ord.forest[[k]]<-list()
for (i in 1:10)
{pred.MAP.ord.forest[[k]][[i]]<-vector()
for (j in 1:dim(pred.test.ord.forest[[k]][[i]]$classprobs)[1])
{pred.MAP.ord.forest[[k]][[i]][j]<-categories[which.max(pred.test.ord.forest[[k]][[i]]$classprobs[j,])]
}
print("i")
print(i)
}
print("k")
print(k)
}


#################

# confusion matrices

conf.mat.MAP.ord.forest<-list()
for (k in 1:length(v.b))
{ conf.mat.MAP.ord.forest[[k]]<-list()
for(i in 1:10)
{conf.mat.MAP.ord.forest[[k]][[i]]<-mat.square(table(pred.MAP.ord.forest[[k]][[i]],test[[i]]$poverty),categories)
print("i")
print(i)
}
print("k")
print(k)
}


# MAE 

MAE.MAP.ord.forest<-list()
for(k in 1:length(v.b))
{ MAE.MAP.ord.forest[[k]]<-vector()
for (i in 1:10)
{MAE.MAP.ord.forest[[k]][i]<-MAE(conf.mat.MAP.ord.forest[[k]][[i]])
}
}


################################################################################
### predictions with Ord.MAP criterion:
###

pred.Ord.MAP.ord.forest<-list()
for (k in 1:length(v.b))
{ pred.Ord.MAP.ord.forest[[k]]<-list()
for (i in 1:10)
{ pred.Ord.MAP.ord.forest[[k]][[i]]<-vector()
for (j in 1:dim(pred.test.ord.forest[[k]][[i]]$classprobs)[1])
{h<-min(which(cumsum(pred.test.ord.forest[[k]][[i]]$classprobs[j,])>=0.5))
pred.Ord.MAP.ord.forest[[k]][[i]][j]<-categories[h]
}
print("i")
print(i)
}
print("k")
print(k)
}



# confusion matrices

conf.mat.Ord.MAP.ord.forest<-list()
for (k in 1:length(v.b))
{ conf.mat.Ord.MAP.ord.forest[[k]]<-list()
for(i in 1:10)
{conf.mat.Ord.MAP.ord.forest[[k]][[i]]<-mat.square(table(pred.Ord.MAP.ord.forest[[k]][[i]],test[[i]]$poverty),categories)
}
}


# MAE 

MAE.Ord.MAP.ord.forest<-list()
for (k in 1:length(v.b))
{ MAE.Ord.MAP.ord.forest[[k]]<-vector()
for (i in 1:10)
{MAE.Ord.MAP.ord.forest[[k]][i]<-MAE(conf.mat.Ord.MAP.ord.forest[[k]][[i]])
}
}



################################################################################
###. Comparison MAP vs Ord.MAP
###
###.  Tests MAE

p.val.test.MAE.MAP.greater.MAE.ORD.MAP.ord.forest<-vector()
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.ord.forest.less<-vector()

for (k in 1:length(v.b))
{p.val.test.MAE.MAP.greater.MAE.ORD.MAP.ord.forest[k]<-wilcox.test(MAE.MAP.ord.forest[[k]],MAE.Ord.MAP.ord.forest[[k]],paired = TRUE,alternative="greater")$p.value
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.ord.forest.less[k]<-wilcox.test(MAE.MAP.ord.forest[[k]],MAE.Ord.MAP.ord.forest[[k]],paired = TRUE,alternative="less")$p.value
}


p.val.test.MAE.MAP.greater.MAE.ORD.MAP.ord.forest
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.ord.forest.less


###############################
###############################
###   Boxplots   MAE


MAE.results.ord.forest<-as.data.frame(cbind(MAE.MAP.ord.forest[[1]]-MAE.Ord.MAP.ord.forest[[1]],
                                            MAE.MAP.ord.forest[[2]]-MAE.Ord.MAP.ord.forest[[2]],
                                            MAE.MAP.ord.forest[[3]]-MAE.Ord.MAP.ord.forest[[3]],
                                            MAE.MAP.ord.forest[[4]]-MAE.Ord.MAP.ord.forest[[4]],
                                            MAE.MAP.ord.forest[[5]]-MAE.Ord.MAP.ord.forest[[5]]
))

boxplot(MAE.results.ord.forest,
        main=paste("(a) World Values Surveys (WVS) dataset"), 
        xlab="nbest hyperparmeter", ylab="MAE increment (MAP - Ord.MAP)",names=c("8", "9", "10", "11", "12"),
        col="white")
abline(h = 0, col = 2, lwd = 2, lty=2)
# Dots
stripchart(MAE.results.ord.forest,
           #method = "jitter",
           pch = 19,
           col = 1:5,
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

################################################################################
#
# Construct OF prediction rule using the training dataset (default 
# perffunction = "probability" corresponding to the 
# (negative) ranked probability score as performance function):

v.a<-1000   # values for the hyperparameter nsets (default=1000)
v.b<-seq(from=8, to=12, by=1)    # values for the hyperparameter nbest (default=10)


model.ord.forest<-list()
for (k in 1:length(v.b))
{ model.ord.forest[[k]]<-list()
for (i in 1:10)
{model.ord.forest[[k]][[i]] <- ordfor(depvar="rating", data=training[[i]], nsets=v.a, nbest=v.b[k],ntreeperdiv=100, 
                                      ntreefinal=1000, npermtrial=100,perffunction = "probability")
print("i")
print(i)
}
print("k")
print(k)
}


################################################################################
#
# Predictions of any of the test set
# 
#

pred.test.ord.forest<-list()
for (k in 1:length(v.b))
{ pred.test.ord.forest[[k]]<-list()
for (i in 1:10)
{pred.test.ord.forest[[k]][[i]] <- predict(model.ord.forest[[k]][[i]],test[[i]][,-2])
print("i")
print(i)
}
print("k")
print(k)
}

#pred.test.ord.forest[[k]][[i]]$ypred
#pred.test.ord.forest[[k]][[i]]$classprobs


################################################################################
### predictions with standard MAP criterion:
###

pred.MAP.ord.forest<-list()
for(k in 1:length(v.b))
{ pred.MAP.ord.forest[[k]]<-list()
for (i in 1:10)
{pred.MAP.ord.forest[[k]][[i]]<-vector()
for (j in 1:dim(pred.test.ord.forest[[k]][[i]]$classprobs)[1])
{pred.MAP.ord.forest[[k]][[i]][j]<-categories[which.max(pred.test.ord.forest[[k]][[i]]$classprobs[j,])]
}
print("i")
print(i)
}
print("k")
print(k)
}


#################

# confusion matrices

conf.mat.MAP.ord.forest<-list()
for (k in 1:length(v.b))
{ conf.mat.MAP.ord.forest[[k]]<-list()
for(i in 1:10)
{conf.mat.MAP.ord.forest[[k]][[i]]<-mat.square(table(pred.MAP.ord.forest[[k]][[i]],test[[i]]$rating),categories)
print("i")
print(i)
}
print("k")
print(k)
}


# MAE 

MAE.MAP.ord.forest<-list()
for(k in 1:length(v.b))
{ MAE.MAP.ord.forest[[k]]<-vector()
for (i in 1:10)
{MAE.MAP.ord.forest[[k]][i]<-MAE(conf.mat.MAP.ord.forest[[k]][[i]])
}
}


################################################################################
### predictions with Ord.MAP criterion:
###

pred.Ord.MAP.ord.forest<-list()
for (k in 1:length(v.b))
{ pred.Ord.MAP.ord.forest[[k]]<-list()
for (i in 1:10)
{ pred.Ord.MAP.ord.forest[[k]][[i]]<-vector()
for (j in 1:dim(pred.test.ord.forest[[k]][[i]]$classprobs)[1])
{h<-min(which(cumsum(pred.test.ord.forest[[k]][[i]]$classprobs[j,])>=0.5))
pred.Ord.MAP.ord.forest[[k]][[i]][j]<-categories[h]
}
print("i")
print(i)
}
print("k")
print(k)
}



# confusion matrices

conf.mat.Ord.MAP.ord.forest<-list()
for (k in 1:length(v.b))
{ conf.mat.Ord.MAP.ord.forest[[k]]<-list()
for(i in 1:10)
{conf.mat.Ord.MAP.ord.forest[[k]][[i]]<-mat.square(table(pred.Ord.MAP.ord.forest[[k]][[i]],test[[i]]$rating),categories)
}
}


# MAE 

MAE.Ord.MAP.ord.forest<-list()
for (k in 1:length(v.b))
{ MAE.Ord.MAP.ord.forest[[k]]<-vector()
for (i in 1:10)
{MAE.Ord.MAP.ord.forest[[k]][i]<-MAE(conf.mat.Ord.MAP.ord.forest[[k]][[i]])
}
}



################################################################################
###. Comparison MAP vs Ord.MAP
###
###.  Tests MAE

p.val.test.MAE.MAP.greater.MAE.ORD.MAP.ord.forest<-vector()
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.ord.forest.less<-vector()

for (k in 1:length(v.b))
{p.val.test.MAE.MAP.greater.MAE.ORD.MAP.ord.forest[k]<-wilcox.test(MAE.MAP.ord.forest[[k]],MAE.Ord.MAP.ord.forest[[k]],paired = TRUE,alternative="greater")$p.value
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.ord.forest.less[k]<-wilcox.test(MAE.MAP.ord.forest[[k]],MAE.Ord.MAP.ord.forest[[k]],paired = TRUE,alternative="less")$p.value
}


p.val.test.MAE.MAP.greater.MAE.ORD.MAP.ord.forest
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.ord.forest.less



###############################
###############################
###   Boxplots   MAE


MAE.results.ord.forest<-as.data.frame(cbind(MAE.MAP.ord.forest[[1]]-MAE.Ord.MAP.ord.forest[[1]],
                                            MAE.MAP.ord.forest[[2]]-MAE.Ord.MAP.ord.forest[[2]],
                                            MAE.MAP.ord.forest[[3]]-MAE.Ord.MAP.ord.forest[[3]],
                                            MAE.MAP.ord.forest[[4]]-MAE.Ord.MAP.ord.forest[[4]],
                                            MAE.MAP.ord.forest[[5]]-MAE.Ord.MAP.ord.forest[[5]]
))

boxplot(MAE.results.ord.forest,
        main=paste("(b) Wine dataset"), 
        xlab="nbest hyperparmeter", ylab="MAE increment (MAP - Ord.MAP)",names=c("8", "9", "10", "11", "12"),
        col="white")
abline(h = 0, col = 2, lwd = 2, lty=2)
# Dots
stripchart(MAE.results.ord.forest,
           #method = "jitter",
           pch = 19,
           col = 1:5,
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

################################################################################
#
# Construct OF prediction rule using the training dataset (default 
# perffunction = "probability" corresponding to the 
# (negative) ranked probability score as performance function):

v.a<-1000   # values for the hyperparameter nsets (default=1000)
v.b<-seq(from=8, to=12, by=1)    # values for the hyperparameter nbest (default=10)


model.ord.forest<-list()
for (k in 1:length(v.b))
{ model.ord.forest[[k]]<-list()
for (i in 1:10)
{model.ord.forest[[k]][[i]] <- ordfor(depvar="Class", data=training[[i]], nsets=v.a, nbest=v.b[k],ntreeperdiv=100, 
                                      ntreefinal=1000, npermtrial=100,perffunction = "probability")
print("i")
print(i)
}
print("k")
print(k)
}


################################################################################
#
# Predictions of any of the test set
# 
#

pred.test.ord.forest<-list()
for (k in 1:length(v.b))
{ pred.test.ord.forest[[k]]<-list()
for (i in 1:10)
{pred.test.ord.forest[[k]][[i]] <- predict(model.ord.forest[[k]][[i]],test[[i]][,-11])
print("i")
print(i)
}
print("k")
print(k)
}

#pred.test.ord.forest[[k]][[i]]$ypred
#pred.test.ord.forest[[k]][[i]]$classprobs


################################################################################
### predictions with standard MAP criterion:
###

pred.MAP.ord.forest<-list()
for(k in 1:length(v.b))
{ pred.MAP.ord.forest[[k]]<-list()
for (i in 1:10)
{pred.MAP.ord.forest[[k]][[i]]<-vector()
for (j in 1:dim(pred.test.ord.forest[[k]][[i]]$classprobs)[1])
{pred.MAP.ord.forest[[k]][[i]][j]<-categories[which.max(pred.test.ord.forest[[k]][[i]]$classprobs[j,])]
}
print("i")
print(i)
}
print("k")
print(k)
}


#################

# confusion matrices

conf.mat.MAP.ord.forest<-list()
for (k in 1:length(v.b))
{ conf.mat.MAP.ord.forest[[k]]<-list()
for(i in 1:10)
{conf.mat.MAP.ord.forest[[k]][[i]]<-mat.square(table(pred.MAP.ord.forest[[k]][[i]],test[[i]]$Class),categories)
print("i")
print(i)
}
print("k")
print(k)
}


# MAE 

MAE.MAP.ord.forest<-list()
for(k in 1:length(v.b))
{ MAE.MAP.ord.forest[[k]]<-vector()
for (i in 1:10)
{MAE.MAP.ord.forest[[k]][i]<-MAE(conf.mat.MAP.ord.forest[[k]][[i]])
}
}


################################################################################
### predictions with Ord.MAP criterion:
###

pred.Ord.MAP.ord.forest<-list()
for (k in 1:length(v.b))
{ pred.Ord.MAP.ord.forest[[k]]<-list()
for (i in 1:10)
{ pred.Ord.MAP.ord.forest[[k]][[i]]<-vector()
for (j in 1:dim(pred.test.ord.forest[[k]][[i]]$classprobs)[1])
{h<-min(which(cumsum(pred.test.ord.forest[[k]][[i]]$classprobs[j,])>=0.5))
pred.Ord.MAP.ord.forest[[k]][[i]][j]<-categories[h]
}
print("i")
print(i)
}
print("k")
print(k)
}



# confusion matrices

conf.mat.Ord.MAP.ord.forest<-list()
for (k in 1:length(v.b))
{ conf.mat.Ord.MAP.ord.forest[[k]]<-list()
for(i in 1:10)
{conf.mat.Ord.MAP.ord.forest[[k]][[i]]<-mat.square(table(pred.Ord.MAP.ord.forest[[k]][[i]],test[[i]]$Class),categories)
}
}


# MAE 

MAE.Ord.MAP.ord.forest<-list()
for (k in 1:length(v.b))
{ MAE.Ord.MAP.ord.forest[[k]]<-vector()
for (i in 1:10)
{MAE.Ord.MAP.ord.forest[[k]][i]<-MAE(conf.mat.Ord.MAP.ord.forest[[k]][[i]])
}
}



################################################################################
###. Comparison MAP vs Ord.MAP
###
###.  Tests MAE

p.val.test.MAE.MAP.greater.MAE.ORD.MAP.ord.forest<-vector()
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.ord.forest.less<-vector()

for (k in 1:length(v.b))
{p.val.test.MAE.MAP.greater.MAE.ORD.MAP.ord.forest[k]<-wilcox.test(MAE.MAP.ord.forest[[k]],MAE.Ord.MAP.ord.forest[[k]],paired = TRUE,alternative="greater")$p.value
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.ord.forest.less[k]<-wilcox.test(MAE.MAP.ord.forest[[k]],MAE.Ord.MAP.ord.forest[[k]],paired = TRUE,alternative="less")$p.value
}


p.val.test.MAE.MAP.greater.MAE.ORD.MAP.ord.forest
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.ord.forest.less



###############################
###############################
###   Boxplots   MAE


MAE.results.ord.forest<-as.data.frame(cbind(MAE.MAP.ord.forest[[1]]-MAE.Ord.MAP.ord.forest[[1]],
                                            MAE.MAP.ord.forest[[2]]-MAE.Ord.MAP.ord.forest[[2]],
                                            MAE.MAP.ord.forest[[3]]-MAE.Ord.MAP.ord.forest[[3]],
                                            MAE.MAP.ord.forest[[4]]-MAE.Ord.MAP.ord.forest[[4]],
                                            MAE.MAP.ord.forest[[5]]-MAE.Ord.MAP.ord.forest[[5]]
))

boxplot(MAE.results.ord.forest,
        main=paste("(c) Hearth dataset"), 
        xlab="nbest hyperparmeter", ylab="MAE increment (MAP - Ord.MAP)",names=c("8", "9", "10", "11", "12"),
        col="white")
abline(h = 0, col = 2, lwd = 2, lty=2)
# Dots
stripchart(MAE.results.ord.forest,
           #method = "jitter",
           pch = 19,
           col = 1:5,
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
# View(parkinsons)
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

install.packages("arules")
library(arules)

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


################################################################################
#
# Construct OF prediction rule using the training dataset (default 
# perffunction = "probability" corresponding to the 
# (negative) ranked probability score as performance function):

v.a<-1000   # values for the hyperparameter nsets (default=1000)
v.b<-seq(from=8, to=12, by=1)    # values for the hyperparameter nbest (default=10)


#
#### FIRST: OUTPUT VARIABLE TO PREDICT: V5.bin.num.factor

model.ord.forest<-list()
for (k in 1:length(v.b))
{ model.ord.forest[[k]]<-list()
for (i in 1:10)
{model.ord.forest[[k]][[i]] <- ordfor(depvar="V5.bin.num.factor", data=training[[i]], nsets=v.a, nbest=v.b[k],ntreeperdiv=5, 
                                      ntreefinal=50, npermtrial=5, perffunction = "probability")
print("i")
print(i)
}
print("k")
print(k)
}


################################################################################
#
# Predictions of any of the test set
# 
#

pred.test.ord.forest<-list()
for (k in 1:length(v.b))
{ pred.test.ord.forest[[k]]<-list()
for (i in 1:10)
{pred.test.ord.forest[[k]][[i]] <- predict(model.ord.forest[[k]][[i]],test[[i]][,])
print("i")
print(i)
}
print("k")
print(k)
}

#pred.test.ord.forest[[k]][[i]]$ypred
#pred.test.ord.forest[[k]][[i]]$classprobs


################################################################################
### predictions with standard MAP criterion:
###

pred.MAP.ord.forest<-list()
for(k in 1:length(v.b))
{ pred.MAP.ord.forest[[k]]<-list()
for (i in 1:10)
{pred.MAP.ord.forest[[k]][[i]]<-vector()
for (j in 1:dim(pred.test.ord.forest[[k]][[i]]$classprobs)[1])
{pred.MAP.ord.forest[[k]][[i]][j]<-categories[which.max(pred.test.ord.forest[[k]][[i]]$classprobs[j,])]
}
print("i")
print(i)
}
print("k")
print(k)
}


#################

# confusion matrices

conf.mat.MAP.ord.forest<-list()
for (k in 1:length(v.b))
{ conf.mat.MAP.ord.forest[[k]]<-list()
for(i in 1:10)
{conf.mat.MAP.ord.forest[[k]][[i]]<-mat.square(table(pred.MAP.ord.forest[[k]][[i]],test[[i]]$V5.bin.num.factor),categories)
print("i")
print(i)
}
print("k")
print(k)
}


# MAE 

MAE.MAP.ord.forest<-list()
for(k in 1:length(v.b))
{ MAE.MAP.ord.forest[[k]]<-vector()
for (i in 1:10)
{MAE.MAP.ord.forest[[k]][i]<-MAE(conf.mat.MAP.ord.forest[[k]][[i]])
}
}


################################################################################
### predictions with Ord.MAP criterion:
###

pred.Ord.MAP.ord.forest<-list()
for (k in 1:length(v.b))
{ pred.Ord.MAP.ord.forest[[k]]<-list()
for (i in 1:10)
{ pred.Ord.MAP.ord.forest[[k]][[i]]<-vector()
for (j in 1:dim(pred.test.ord.forest[[k]][[i]]$classprobs)[1])
{h<-min(which(cumsum(pred.test.ord.forest[[k]][[i]]$classprobs[j,])>=0.5))
pred.Ord.MAP.ord.forest[[k]][[i]][j]<-categories[h]
}
print("i")
print(i)
}
print("k")
print(k)
}



# confusion matrices

conf.mat.Ord.MAP.ord.forest<-list()
for (k in 1:length(v.b))
{ conf.mat.Ord.MAP.ord.forest[[k]]<-list()
for(i in 1:10)
{conf.mat.Ord.MAP.ord.forest[[k]][[i]]<-mat.square(table(pred.Ord.MAP.ord.forest[[k]][[i]],test[[i]]$V5.bin.num.factor),categories)
}
}


# MAE 

MAE.Ord.MAP.ord.forest<-list()
for (k in 1:length(v.b))
{ MAE.Ord.MAP.ord.forest[[k]]<-vector()
for (i in 1:10)
{MAE.Ord.MAP.ord.forest[[k]][i]<-MAE(conf.mat.Ord.MAP.ord.forest[[k]][[i]])
}
}



################################################################################
###. Comparison MAP vs Ord.MAP
###
###.  Tests MAE

p.val.test.MAE.MAP.greater.MAE.ORD.MAP.ord.forest<-vector()
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.ord.forest.less<-vector()

for (k in 1:length(v.b))
{p.val.test.MAE.MAP.greater.MAE.ORD.MAP.ord.forest[k]<-wilcox.test(MAE.MAP.ord.forest[[k]],MAE.Ord.MAP.ord.forest[[k]],paired = TRUE,alternative="greater")$p.value
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.ord.forest.less[k]<-wilcox.test(MAE.MAP.ord.forest[[k]],MAE.Ord.MAP.ord.forest[[k]],paired = TRUE,alternative="less")$p.value
}


p.val.test.MAE.MAP.greater.MAE.ORD.MAP.ord.forest
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.ord.forest.less



###############################
###############################
###   Boxplots   MAE


MAE.results.ord.forest<-as.data.frame(cbind(MAE.MAP.ord.forest[[1]]-MAE.Ord.MAP.ord.forest[[1]],
                                            MAE.MAP.ord.forest[[2]]-MAE.Ord.MAP.ord.forest[[2]],
                                            MAE.MAP.ord.forest[[3]]-MAE.Ord.MAP.ord.forest[[3]],
                                            MAE.MAP.ord.forest[[4]]-MAE.Ord.MAP.ord.forest[[4]],
                                            MAE.MAP.ord.forest[[5]]-MAE.Ord.MAP.ord.forest[[5]]
))

boxplot(MAE.results.ord.forest,
        main=paste("(d) Parkinson dataset: output V5 (5 classes)"), 
        xlab="nbest hyperparmeter", ylab="MAE increment (MAP - Ord.MAP)",names=c("8", "9", "10", "11", "12"),
        col="white")
abline(h = 0, col = 2, lwd = 2, lty=2)
# Dots
stripchart(MAE.results.ord.forest,
           #method = "jitter",
           pch = 19,
           col = 1:5,
           vertical = TRUE,
           add = TRUE)





#
#### SECOND: OUTPUT VARIABLE TO PREDICT: V6.bin.num.factor

model.ord.forest<-list()
for (k in 1:length(v.b))
{ model.ord.forest[[k]]<-list()
for (i in 1:10)
{model.ord.forest[[k]][[i]] <- ordfor(depvar="V6.bin.num.factor", data=training[[i]], nsets=v.a, nbest=v.b[k],ntreeperdiv=5, 
                                      ntreefinal=50, npermtrial=5, perffunction = "probability")
print("i")
print(i)
}
print("k")
print(k)
}


################################################################################
#
# Predictions of any of the test set
# 
#

pred.test.ord.forest<-list()
for (k in 1:length(v.b))
{ pred.test.ord.forest[[k]]<-list()
for (i in 1:10)
{pred.test.ord.forest[[k]][[i]] <- predict(model.ord.forest[[k]][[i]],test[[i]][,])
print("i")
print(i)
}
print("k")
print(k)
}

#pred.test.ord.forest[[k]][[i]]$ypred
#pred.test.ord.forest[[k]][[i]]$classprobs


################################################################################
### predictions with standard MAP criterion:
###

pred.MAP.ord.forest<-list()
for(k in 1:length(v.b))
{ pred.MAP.ord.forest[[k]]<-list()
for (i in 1:10)
{pred.MAP.ord.forest[[k]][[i]]<-vector()
for (j in 1:dim(pred.test.ord.forest[[k]][[i]]$classprobs)[1])
{pred.MAP.ord.forest[[k]][[i]][j]<-categories[which.max(pred.test.ord.forest[[k]][[i]]$classprobs[j,])]
}
print("i")
print(i)
}
print("k")
print(k)
}


#################

# confusion matrices

conf.mat.MAP.ord.forest<-list()
for (k in 1:length(v.b))
{ conf.mat.MAP.ord.forest[[k]]<-list()
for(i in 1:10)
{conf.mat.MAP.ord.forest[[k]][[i]]<-mat.square(table(pred.MAP.ord.forest[[k]][[i]],test[[i]]$V6.bin.num.factor),categories)
print("i")
print(i)
}
print("k")
print(k)
}


# MAE 

MAE.MAP.ord.forest<-list()
for(k in 1:length(v.b))
{ MAE.MAP.ord.forest[[k]]<-vector()
for (i in 1:10)
{MAE.MAP.ord.forest[[k]][i]<-MAE(conf.mat.MAP.ord.forest[[k]][[i]])
}
}


################################################################################
### predictions with Ord.MAP criterion:
###

pred.Ord.MAP.ord.forest<-list()
for (k in 1:length(v.b))
{ pred.Ord.MAP.ord.forest[[k]]<-list()
for (i in 1:10)
{ pred.Ord.MAP.ord.forest[[k]][[i]]<-vector()
for (j in 1:dim(pred.test.ord.forest[[k]][[i]]$classprobs)[1])
{h<-min(which(cumsum(pred.test.ord.forest[[k]][[i]]$classprobs[j,])>=0.5))
pred.Ord.MAP.ord.forest[[k]][[i]][j]<-categories[h]
}
print("i")
print(i)
}
print("k")
print(k)
}



# confusion matrices

conf.mat.Ord.MAP.ord.forest<-list()
for (k in 1:length(v.b))
{ conf.mat.Ord.MAP.ord.forest[[k]]<-list()
for(i in 1:10)
{conf.mat.Ord.MAP.ord.forest[[k]][[i]]<-mat.square(table(pred.Ord.MAP.ord.forest[[k]][[i]],test[[i]]$V6.bin.num.factor),categories)
}
}


# MAE 

MAE.Ord.MAP.ord.forest<-list()
for (k in 1:length(v.b))
{ MAE.Ord.MAP.ord.forest[[k]]<-vector()
for (i in 1:10)
{MAE.Ord.MAP.ord.forest[[k]][i]<-MAE(conf.mat.Ord.MAP.ord.forest[[k]][[i]])
}
}



################################################################################
###. Comparison MAP vs Ord.MAP
###
###.  Tests MAE

p.val.test.MAE.MAP.greater.MAE.ORD.MAP.ord.forest<-vector()
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.ord.forest.less<-vector()

for (k in 1:length(v.b))
{p.val.test.MAE.MAP.greater.MAE.ORD.MAP.ord.forest[k]<-wilcox.test(MAE.MAP.ord.forest[[k]],MAE.Ord.MAP.ord.forest[[k]],paired = TRUE,alternative="greater")$p.value
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.ord.forest.less[k]<-wilcox.test(MAE.MAP.ord.forest[[k]],MAE.Ord.MAP.ord.forest[[k]],paired = TRUE,alternative="less")$p.value
}


p.val.test.MAE.MAP.greater.MAE.ORD.MAP.ord.forest
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.ord.forest.less



###############################
###############################
###   Boxplots   MAE


MAE.results.ord.forest<-as.data.frame(cbind(MAE.MAP.ord.forest[[1]]-MAE.Ord.MAP.ord.forest[[1]],
                                            MAE.MAP.ord.forest[[2]]-MAE.Ord.MAP.ord.forest[[2]],
                                            MAE.MAP.ord.forest[[3]]-MAE.Ord.MAP.ord.forest[[3]],
                                            MAE.MAP.ord.forest[[4]]-MAE.Ord.MAP.ord.forest[[4]],
                                            MAE.MAP.ord.forest[[5]]-MAE.Ord.MAP.ord.forest[[5]]
))

boxplot(MAE.results.ord.forest,
        main=paste("(d) Parkinson dataset: output V6 (5 classes)"), 
        xlab="nbest hyperparmeter", ylab="MAE increment (MAP - Ord.MAP)",names=c("8", "9", "10", "11", "12"),
        col="white")
abline(h = 0, col = 2, lwd = 2, lty=2)
# Dots
stripchart(MAE.results.ord.forest,
           #method = "jitter",
           pch = 19,
           col = 1:5,
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
#
# Construct OF prediction rule using the training dataset (default 
# perffunction = "probability" corresponding to the 
# (negative) ranked probability score as performance function):

v.a<-1000   # values for the hyperparameter nsets (default=1000)
v.b<-seq(from=8, to=12, by=1)    # values for the hyperparameter nbest (default=10)

model.ord.forest<-list()
for (k in 1:length(v.b))
{ model.ord.forest[[k]]<-list()
for (i in 1:10)
{model.ord.forest[[k]][[i]] <- ordfor(depvar="Class", data=training[[i]], nsets=v.a, nbest=v.b[k],ntreeperdiv=100, 
                                      ntreefinal=1000, npermtrial=100,perffunction = "probability")
print("i")
print(i)
}
print("k")
print(k)
}


################################################################################
#
# Predictions of any of the test set
# 
#

pred.test.ord.forest<-list()
for (k in 1:length(v.b))
{ pred.test.ord.forest[[k]]<-list()
for (i in 1:10)
{pred.test.ord.forest[[k]][[i]] <- predict(model.ord.forest[[k]][[i]],test[[i]][, ])
print("i")
print(i)
}
print("k")
print(k)
}

#pred.test.ord.forest[[k]][[i]]$ypred
#pred.test.ord.forest[[k]][[i]]$classprobs


################################################################################
### predictions with standard MAP criterion:
###

pred.MAP.ord.forest<-list()
for(k in 1:length(v.b))
{ pred.MAP.ord.forest[[k]]<-list()
for (i in 1:10)
{pred.MAP.ord.forest[[k]][[i]]<-vector()
for (j in 1:dim(pred.test.ord.forest[[k]][[i]]$classprobs)[1])
{pred.MAP.ord.forest[[k]][[i]][j]<-categories[which.max(pred.test.ord.forest[[k]][[i]]$classprobs[j,])]
}
print("i")
print(i)
}
print("k")
print(k)
}


#################

# confusion matrices

conf.mat.MAP.ord.forest<-list()
for (k in 1:length(v.b))
{ conf.mat.MAP.ord.forest[[k]]<-list()
for(i in 1:10)
{conf.mat.MAP.ord.forest[[k]][[i]]<-mat.square(table(pred.MAP.ord.forest[[k]][[i]],test[[i]]$Class),categories)
print("i")
print(i)
}
print("k")
print(k)
}


# MAE 

MAE.MAP.ord.forest<-list()
for(k in 1:length(v.b))
{ MAE.MAP.ord.forest[[k]]<-vector()
for (i in 1:10)
{MAE.MAP.ord.forest[[k]][i]<-MAE(conf.mat.MAP.ord.forest[[k]][[i]])
}
}


################################################################################
### predictions with Ord.MAP criterion:
###

pred.Ord.MAP.ord.forest<-list()
for (k in 1:length(v.b))
{ pred.Ord.MAP.ord.forest[[k]]<-list()
for (i in 1:10)
{ pred.Ord.MAP.ord.forest[[k]][[i]]<-vector()
for (j in 1:dim(pred.test.ord.forest[[k]][[i]]$classprobs)[1])
{h<-min(which(cumsum(pred.test.ord.forest[[k]][[i]]$classprobs[j,])>=0.5))
pred.Ord.MAP.ord.forest[[k]][[i]][j]<-categories[h]
}
print("i")
print(i)
}
print("k")
print(k)
}



# confusion matrices

conf.mat.Ord.MAP.ord.forest<-list()
for (k in 1:length(v.b))
{ conf.mat.Ord.MAP.ord.forest[[k]]<-list()
for(i in 1:10)
{conf.mat.Ord.MAP.ord.forest[[k]][[i]]<-mat.square(table(pred.Ord.MAP.ord.forest[[k]][[i]],test[[i]]$Class),categories)
}
}


# MAE 

MAE.Ord.MAP.ord.forest<-list()
for (k in 1:length(v.b))
{ MAE.Ord.MAP.ord.forest[[k]]<-vector()
for (i in 1:10)
{MAE.Ord.MAP.ord.forest[[k]][i]<-MAE(conf.mat.Ord.MAP.ord.forest[[k]][[i]])
}
}



################################################################################
###. Comparison MAP vs Ord.MAP
###
###.  Tests MAE

p.val.test.MAE.MAP.greater.MAE.ORD.MAP.ord.forest<-vector()
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.ord.forest.less<-vector()

for (k in 1:length(v.b))
{p.val.test.MAE.MAP.greater.MAE.ORD.MAP.ord.forest[k]<-wilcox.test(MAE.MAP.ord.forest[[k]],MAE.Ord.MAP.ord.forest[[k]],paired = TRUE,alternative="greater")$p.value
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.ord.forest.less[k]<-wilcox.test(MAE.MAP.ord.forest[[k]],MAE.Ord.MAP.ord.forest[[k]],paired = TRUE,alternative="less")$p.value
}


p.val.test.MAE.MAP.greater.MAE.ORD.MAP.ord.forest
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.ord.forest.less



###############################
###############################
###   Boxplots   MAE


MAE.results.ord.forest<-as.data.frame(cbind(MAE.MAP.ord.forest[[1]]-MAE.Ord.MAP.ord.forest[[1]],
                                            MAE.MAP.ord.forest[[2]]-MAE.Ord.MAP.ord.forest[[2]],
                                            MAE.MAP.ord.forest[[3]]-MAE.Ord.MAP.ord.forest[[3]],
                                            MAE.MAP.ord.forest[[4]]-MAE.Ord.MAP.ord.forest[[4]],
                                            MAE.MAP.ord.forest[[5]]-MAE.Ord.MAP.ord.forest[[5]]
))

boxplot(MAE.results.ord.forest,
        main=paste("(e) CES11 dataset"), 
        xlab="nbest hyperparmeter", ylab="MAE increment (MAP - Ord.MAP)",names=c("8", "9", "10", "11", "12"),
        col="white")
abline(h = 0, col = 2, lwd = 2, lty=2)
# Dots
stripchart(MAE.results.ord.forest,
           #method = "jitter",
           pch = 19,
           col = 1:5,
           vertical = TRUE,
           add = TRUE)



#
##
###
##
#
