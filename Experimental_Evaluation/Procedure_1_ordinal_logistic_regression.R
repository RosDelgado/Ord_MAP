
###########################################
###########################################
#
#  EXPERIMENTAL PHASE  (Section 5)
# 
############################################
#
# PROCEDURE 1): Ordinal Logistic Regression (OLR)
# using polr function from the MASS package.
#
# 

library(MASS)

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
# We’ll now fit the Proportional Odds Logistic Regression model 
# using polr function from the MASS package.
#

model.logistic<-list()
for (i in 1:10)
{model.logistic[[i]] <- polr(poverty~religion+degree+country+age+gender, data = training[[i]], Hess = TRUE)
}

model.probit<-list()
for (i in 1:10)
{model.probit[[i]] <- polr(poverty~religion+degree+country+age+gender, data = training[[i]], Hess = TRUE,method="probit")
}

model.loglog<-list()
for (i in 1:10)
{model.loglog[[i]] <- polr(poverty~religion+degree+country+age+gender, data = training[[i]], Hess = TRUE,method="loglog")
}

model.cloglog<-list()
for (i in 1:10)
{model.cloglog[[i]] <- polr(poverty~religion+degree+country+age+gender, data = training[[i]], Hess = TRUE,method="cloglog")
}

model.cauchit<-list()
for (i in 1:10)
{model.cauchit[[i]] <- polr(poverty~religion+degree+country+age+gender, data = training[[i]], Hess = TRUE,method="cauchit")
}

# polr method = "logistic" by default. Other: "probit", "loglog", "cloglog", "cauchit".

# Some link functions in ordinal regression models may fail to converge or produce non-finite results 
# for certain datasets. This is often due to the inherent asymmetry of the underlying distributions 
# that these link functions assume. Specifically, the probit link assumes a normal distribution, 
# which can be problematic when the data is heavily skewed or has extreme outliers. 
# Similarly, the cauchit link, which is based on the inverse of the tangent of the logistic function, 
# is more sensitive to extreme values and can lead to numerical instability if the data does not conform 
# well to this distribution. The loglog and cloglog links, which are based on extreme value distributions, 
# can also be problematic in datasets with low variability or when the outcome categories are imbalanced, 
# as they assume a more severe skew in the response variable. These issues typically arise when there 
# is multicollinearity among predictors or when predictors have highly different scales.


# summary(model_fit[[2]])


################################################################################
#
# Predictions of any of the test set
# 
#

pred.test.logistic<-list()
for (i in 1:10)
{pred.test.logistic[[i]] <- predict(model.logistic[[i]],test[[i]][,-1],type = "p")
}


pred.test.probit<-list()
for (i in 1:10)
{pred.test.probit[[i]] <- predict(model.probit[[i]],test[[i]][,-1],type = "p")
}

pred.test.loglog<-list()
for (i in 1:10)
{pred.test.loglog[[i]] <- predict(model.loglog[[i]],test[[i]][,-1],type = "p")
}

pred.test.cloglog<-list()
for (i in 1:10)
{pred.test.cloglog[[i]] <- predict(model.cloglog[[i]],test[[i]][,-1],type = "p")
}

pred.test.cauchit<-list()
for (i in 1:10)
{pred.test.cauchit[[i]] <- predict(model.cauchit[[i]],test[[i]][,-1],type = "p")
}



################################################################################
### predictions with standard MAP criterion:
###

pred.MAP.logistic<-list()
for (i in 1:10)
{pred.MAP.logistic[[i]]<-vector()
 for (j in 1:dim(pred.test.logistic[[i]])[1])
  {pred.MAP.logistic[[i]][j]<-categories[which.max(pred.test.logistic[[i]][j,])]
  }
}


pred.MAP.probit<-list()
for (i in 1:10)
{pred.MAP.probit[[i]]<-vector()
for (j in 1:dim(pred.test.probit[[i]])[1])
{pred.MAP.probit[[i]][j]<-categories[which.max(pred.test.probit[[i]][j,])]
}
}

pred.MAP.loglog<-list()
for (i in 1:10)
{pred.MAP.loglog[[i]]<-vector()
for (j in 1:dim(pred.test.loglog[[i]])[1])
{pred.MAP.loglog[[i]][j]<-categories[which.max(pred.test.loglog[[i]][j,])]
}
}

pred.MAP.cloglog<-list()
for (i in 1:10)
{pred.MAP.cloglog[[i]]<-vector()
for (j in 1:dim(pred.test.cloglog[[i]])[1])
{pred.MAP.cloglog[[i]][j]<-categories[which.max(pred.test.cloglog[[i]][j,])]
}
}

pred.MAP.cauchit<-list()
for (i in 1:10)
{pred.MAP.cauchit[[i]]<-vector()
for (j in 1:dim(pred.test.cauchit[[i]])[1])
{pred.MAP.cauchit[[i]][j]<-categories[which.max(pred.test.cauchit[[i]][j,])]
}
}


# confusion matrices


conf.mat.MAP.logistic<-list()
for(i in 1:10)
{
conf.mat.MAP.logistic[[i]]<-mat.square(table(pred.MAP.logistic[[i]],test[[i]]$poverty),categories)
}

conf.mat.MAP.probit<-list()
for(i in 1:10)
{
  conf.mat.MAP.probit[[i]]<-mat.square(table(pred.MAP.probit[[i]],test[[i]]$poverty),categories)
}

conf.mat.MAP.loglog<-list()
for(i in 1:10)
{
  conf.mat.MAP.loglog[[i]]<-mat.square(table(pred.MAP.loglog[[i]],test[[i]]$poverty),categories)
}

conf.mat.MAP.cloglog<-list()
for(i in 1:10)
{
  conf.mat.MAP.cloglog[[i]]<-mat.square(table(pred.MAP.cloglog[[i]],test[[i]]$poverty),categories)
}

conf.mat.MAP.cauchit<-list()
for(i in 1:10)
{
  conf.mat.MAP.cauchit[[i]]<-mat.square(table(pred.MAP.cauchit[[i]],test[[i]]$poverty),categories)
}


# MAE 

MAE.MAP.logistic<-vector()
for (i in 1:10)
{
  MAE.MAP.logistic[i]<-MAE(conf.mat.MAP.logistic[[i]])
}

MAE.MAP.probit<-vector()
for (i in 1:10)
{
  MAE.MAP.probit[i]<-MAE(conf.mat.MAP.probit[[i]])
}

MAE.MAP.loglog<-vector()
for (i in 1:10)
{
  MAE.MAP.loglog[i]<-MAE(conf.mat.MAP.loglog[[i]])
}

MAE.MAP.cloglog<-vector()
for (i in 1:10)
{
  MAE.MAP.cloglog[i]<-MAE(conf.mat.MAP.cloglog[[i]])
}

MAE.MAP.cauchit<-vector()
for (i in 1:10)
{
  MAE.MAP.cauchit[i]<-MAE(conf.mat.MAP.cauchit[[i]])
}


################################################################################
### predictions with Ord.MAP criterion:
###

pred.Ord.MAP.logistic<-list()
for (i in 1:10)
{ pred.Ord.MAP.logistic[[i]]<-vector()
 for (j in 1:dim(pred.test.logistic[[i]])[1])
 {h<-min(which(cumsum(pred.test.logistic[[i]][j,])>=0.5))
  pred.Ord.MAP.logistic[[i]][j]<-categories[h]
}
}


pred.Ord.MAP.probit<-list()
for (i in 1:10)
{ pred.Ord.MAP.probit[[i]]<-vector()
for (j in 1:dim(pred.test.probit[[i]])[1])
{h<-min(which(cumsum(pred.test.probit[[i]][j,])>=0.5))
pred.Ord.MAP.probit[[i]][j]<-categories[h]
}
}

pred.Ord.MAP.loglog<-list()
for (i in 1:10)
{ pred.Ord.MAP.loglog[[i]]<-vector()
for (j in 1:dim(pred.test.loglog[[i]])[1])
{h<-min(which(cumsum(pred.test.loglog[[i]][j,])>=0.5))
pred.Ord.MAP.loglog[[i]][j]<-categories[h]
}
}

pred.Ord.MAP.cloglog<-list()
for (i in 1:10)
{ pred.Ord.MAP.cloglog[[i]]<-vector()
for (j in 1:dim(pred.test.cloglog[[i]])[1])
{h<-min(which(cumsum(pred.test.cloglog[[i]][j,])>=0.5))
pred.Ord.MAP.cloglog[[i]][j]<-categories[h]
}
}

pred.Ord.MAP.cauchit<-list()
for (i in 1:10)
{ pred.Ord.MAP.cauchit[[i]]<-vector()
for (j in 1:dim(pred.test.cauchit[[i]])[1])
{h<-min(which(cumsum(pred.test.cauchit[[i]][j,])>=0.5))
pred.Ord.MAP.cauchit[[i]][j]<-categories[h]
}
}



# confusion matrices

conf.mat.Ord.MAP.logistic<-list()
for(i in 1:10)
{
  conf.mat.Ord.MAP.logistic[[i]]<-mat.square(table(pred.Ord.MAP.logistic[[i]],test[[i]]$poverty),categories)
}

conf.mat.Ord.MAP.probit<-list()
for(i in 1:10)
{
  conf.mat.Ord.MAP.probit[[i]]<-mat.square(table(pred.Ord.MAP.probit[[i]],test[[i]]$poverty),categories)
}

conf.mat.Ord.MAP.loglog<-list()
for(i in 1:10)
{
  conf.mat.Ord.MAP.loglog[[i]]<-mat.square(table(pred.Ord.MAP.loglog[[i]],test[[i]]$poverty),categories)
}

conf.mat.Ord.MAP.cloglog<-list()
for(i in 1:10)
{
  conf.mat.Ord.MAP.cloglog[[i]]<-mat.square(table(pred.Ord.MAP.cloglog[[i]],test[[i]]$poverty),categories)
}

conf.mat.Ord.MAP.cauchit<-list()
for(i in 1:10)
{
  conf.mat.Ord.MAP.cauchit[[i]]<-mat.square(table(pred.Ord.MAP.cauchit[[i]],test[[i]]$poverty),categories)
}


# MAE 

MAE.Ord.MAP.logistic<-vector()
for (i in 1:10)
{
  MAE.Ord.MAP.logistic[i]<-MAE(conf.mat.Ord.MAP.logistic[[i]])
}

MAE.Ord.MAP.probit<-vector()
for (i in 1:10)
{
  MAE.Ord.MAP.probit[i]<-MAE(conf.mat.Ord.MAP.probit[[i]])
}

MAE.Ord.MAP.loglog<-vector()
for (i in 1:10)
{
  MAE.Ord.MAP.loglog[i]<-MAE(conf.mat.Ord.MAP.loglog[[i]])
}

MAE.Ord.MAP.cloglog<-vector()
for (i in 1:10)
{
  MAE.Ord.MAP.cloglog[i]<-MAE(conf.mat.Ord.MAP.cloglog[[i]])
}

MAE.Ord.MAP.cauchit<-vector()
for (i in 1:10)
{
  MAE.Ord.MAP.cauchit[i]<-MAE(conf.mat.Ord.MAP.cauchit[[i]])
}


################################################################################
###. Comparison MAP vs Ord.MAP
###
###.  Tests MAE

p.val.test.MAE.MAP.greater.MAE.ORD.MAP.logistic<-wilcox.test(MAE.MAP.logistic,MAE.Ord.MAP.logistic,paired = TRUE,alternative="greater")$p.value
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.probit<-wilcox.test(MAE.MAP.probit,MAE.Ord.MAP.probit,paired = TRUE,alternative="greater")$p.value
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.loglog<-wilcox.test(MAE.MAP.loglog,MAE.Ord.MAP.loglog,paired = TRUE,alternative="greater")$p.value
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.cloglog<-wilcox.test(MAE.MAP.cloglog,MAE.Ord.MAP.cloglog,paired = TRUE,alternative="greater")$p.value
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.cauchit<-wilcox.test(MAE.MAP.cauchit,MAE.Ord.MAP.cauchit,paired = TRUE,alternative="greater")$p.value

p.val.test.MAE.MAP.greater.MAE.ORD.MAP.logistic.less<-wilcox.test(MAE.MAP.logistic,MAE.Ord.MAP.logistic,paired = TRUE,alternative="less")$p.value
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.probit.less<-wilcox.test(MAE.MAP.probit,MAE.Ord.MAP.probit,paired = TRUE,alternative="less")$p.value
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.loglog.less<-wilcox.test(MAE.MAP.loglog,MAE.Ord.MAP.loglog,paired = TRUE,alternative="less")$p.value
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.cloglog.less<-wilcox.test(MAE.MAP.cloglog,MAE.Ord.MAP.cloglog,paired = TRUE,alternative="less")$p.value
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.cauchit.less<-wilcox.test(MAE.MAP.cauchit,MAE.Ord.MAP.cauchit,paired = TRUE,alternative="less")$p.value



p.val.test.MAE.MAP.greater.MAE.ORD.MAP.logistic
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.probit
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.loglog
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.cloglog
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.cauchit

p.val.test.MAE.MAP.greater.MAE.ORD.MAP.logistic.less
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.probit.less
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.loglog.less
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.cloglog.less
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.cauchit.less



###############################
###############################
###   Boxplots   MAE

MAE.results<-as.data.frame(cbind(MAE.MAP.logistic-MAE.Ord.MAP.logistic,
                                 MAE.MAP.probit-MAE.Ord.MAP.probit,
                                 MAE.MAP.loglog-MAE.Ord.MAP.loglog,
                                 MAE.MAP.cloglog-MAE.Ord.MAP.cloglog,
                                 MAE.MAP.cauchit-MAE.Ord.MAP.cauchit))

boxplot(MAE.results,
        main=paste("(a) World Values Surveys (WVS) dataset"), 
        xlab=" ", ylab="MAE increment (MAP - Ord.MAP)",names=c("logistic", "probit","loglog","cloglog","cauchit"),
        col="white")
abline(h = 0, col = 2, lwd = 2, lty=2)
# Dots
stripchart(MAE.results,
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
# We’ll now fit the Proportional Odds Logistic Regression model 
# using polr function from the MASS package.
#

model.logistic<-list()
for (i in 1:10)
{model.logistic[[i]] <- polr(rating ~ temp + contact, data = training[[i]], Hess = TRUE)
}

model.probit<-list()
for (i in 1:10)
{model.probit[[i]] <- polr(rating ~ temp + contact, data = training[[i]], Hess = TRUE,method="probit")
}

model.loglog<-list()
for (i in 1:10)
{model.loglog[[i]] <- polr(rating ~ temp + contact, data = training[[i]], Hess = TRUE,method="loglog")
}

model.cloglog<-list()
for (i in 1:10)
{model.cloglog[[i]] <- polr(rating ~ temp + contact, data = training[[i]], Hess = TRUE,method="cloglog")
}

model.cauchit<-list()
for (i in 1:10)
{model.cauchit[[i]] <- polr(rating ~ temp + contact, data = training[[i]], Hess = TRUE,method="cauchit")
}

# polr method = "logistic" by default. Other: "probit", "loglog", "cloglog", "cauchit"

# summary(model_fit[[2]])


################################################################################
#
# Predictions of any of the test set
# 
#

pred.test.logistic<-list()
for (i in 1:10)
{pred.test.logistic[[i]] <- predict(model.logistic[[i]],test[[i]][,-2],type = "p")
}


pred.test.probit<-list()
for (i in 1:10)
{pred.test.probit[[i]] <- predict(model.probit[[i]],test[[i]][,-2],type = "p")
}

pred.test.loglog<-list()
for (i in 1:10)
{pred.test.loglog[[i]] <- predict(model.loglog[[i]],test[[i]][,-2],type = "p")
}

pred.test.cloglog<-list()
for (i in 1:10)
{pred.test.cloglog[[i]] <- predict(model.cloglog[[i]],test[[i]][,-2],type = "p")
}

pred.test.cauchit<-list()
for (i in 1:10)
{pred.test.cauchit[[i]] <- predict(model.cauchit[[i]],test[[i]][,-2],type = "p")
}



################################################################################
### predictions with standard MAP criterion:
###

pred.MAP.logistic<-list()
for (i in 1:10)
{pred.MAP.logistic[[i]]<-vector()
for (j in 1:dim(pred.test.logistic[[i]])[1])
{pred.MAP.logistic[[i]][j]<-categories[which.max(pred.test.logistic[[i]][j,])]
}
}


pred.MAP.probit<-list()
for (i in 1:10)
{pred.MAP.probit[[i]]<-vector()
for (j in 1:dim(pred.test.probit[[i]])[1])
{pred.MAP.probit[[i]][j]<-categories[which.max(pred.test.probit[[i]][j,])]
}
}

pred.MAP.loglog<-list()
for (i in 1:10)
{pred.MAP.loglog[[i]]<-vector()
for (j in 1:dim(pred.test.loglog[[i]])[1])
{pred.MAP.loglog[[i]][j]<-categories[which.max(pred.test.loglog[[i]][j,])]
}
}

pred.MAP.cloglog<-list()
for (i in 1:10)
{pred.MAP.cloglog[[i]]<-vector()
for (j in 1:dim(pred.test.cloglog[[i]])[1])
{pred.MAP.cloglog[[i]][j]<-categories[which.max(pred.test.cloglog[[i]][j,])]
}
}

pred.MAP.cauchit<-list()
for (i in 1:10)
{pred.MAP.cauchit[[i]]<-vector()
for (j in 1:dim(pred.test.cauchit[[i]])[1])
{pred.MAP.cauchit[[i]][j]<-categories[which.max(pred.test.cauchit[[i]][j,])]
}
}


# confusion matrices


conf.mat.MAP.logistic<-list()
for(i in 1:10)
{
  conf.mat.MAP.logistic[[i]]<-mat.square(table(pred.MAP.logistic[[i]],test[[i]]$rating),categories)
}

conf.mat.MAP.probit<-list()
for(i in 1:10)
{
  conf.mat.MAP.probit[[i]]<-mat.square(table(pred.MAP.probit[[i]],test[[i]]$rating),categories)
}

conf.mat.MAP.loglog<-list()
for(i in 1:10)
{
  conf.mat.MAP.loglog[[i]]<-mat.square(table(pred.MAP.loglog[[i]],test[[i]]$rating),categories)
}

conf.mat.MAP.cloglog<-list()
for(i in 1:10)
{
  conf.mat.MAP.cloglog[[i]]<-mat.square(table(pred.MAP.cloglog[[i]],test[[i]]$rating),categories)
}

conf.mat.MAP.cauchit<-list()
for(i in 1:10)
{
  conf.mat.MAP.cauchit[[i]]<-mat.square(table(pred.MAP.cauchit[[i]],test[[i]]$rating),categories)
}


# MAE 

MAE.MAP.logistic<-vector()
for (i in 1:10)
{
  MAE.MAP.logistic[i]<-MAE(conf.mat.MAP.logistic[[i]])
}

MAE.MAP.probit<-vector()
for (i in 1:10)
{
  MAE.MAP.probit[i]<-MAE(conf.mat.MAP.probit[[i]])
}

MAE.MAP.loglog<-vector()
for (i in 1:10)
{
  MAE.MAP.loglog[i]<-MAE(conf.mat.MAP.loglog[[i]])
}

MAE.MAP.cloglog<-vector()
for (i in 1:10)
{
  MAE.MAP.cloglog[i]<-MAE(conf.mat.MAP.cloglog[[i]])
}

MAE.MAP.cauchit<-vector()
for (i in 1:10)
{
  MAE.MAP.cauchit[i]<-MAE(conf.mat.MAP.cauchit[[i]])
}




################################################################################
### predictions with Ord.MAP criterion:
###

pred.Ord.MAP.logistic<-list()
for (i in 1:10)
{ pred.Ord.MAP.logistic[[i]]<-vector()
for (j in 1:dim(pred.test.logistic[[i]])[1])
{h<-min(which(cumsum(pred.test.logistic[[i]][j,])>=0.5))
pred.Ord.MAP.logistic[[i]][j]<-categories[h]
}
}


pred.Ord.MAP.probit<-list()
for (i in 1:10)
{ pred.Ord.MAP.probit[[i]]<-vector()
for (j in 1:dim(pred.test.probit[[i]])[1])
{h<-min(which(cumsum(pred.test.probit[[i]][j,])>=0.5))
pred.Ord.MAP.probit[[i]][j]<-categories[h]
}
}

pred.Ord.MAP.loglog<-list()
for (i in 1:10)
{ pred.Ord.MAP.loglog[[i]]<-vector()
for (j in 1:dim(pred.test.loglog[[i]])[1])
{h<-min(which(cumsum(pred.test.loglog[[i]][j,])>=0.5))
pred.Ord.MAP.loglog[[i]][j]<-categories[h]
}
}

pred.Ord.MAP.cloglog<-list()
for (i in 1:10)
{ pred.Ord.MAP.cloglog[[i]]<-vector()
for (j in 1:dim(pred.test.cloglog[[i]])[1])
{h<-min(which(cumsum(pred.test.cloglog[[i]][j,])>=0.5))
pred.Ord.MAP.cloglog[[i]][j]<-categories[h]
}
}

pred.Ord.MAP.cauchit<-list()
for (i in 1:10)
{ pred.Ord.MAP.cauchit[[i]]<-vector()
for (j in 1:dim(pred.test.cauchit[[i]])[1])
{h<-min(which(cumsum(pred.test.cauchit[[i]][j,])>=0.5))
pred.Ord.MAP.cauchit[[i]][j]<-categories[h]
}
}



# confusion matrices

conf.mat.Ord.MAP.logistic<-list()
for(i in 1:10)
{
  conf.mat.Ord.MAP.logistic[[i]]<-mat.square(table(pred.Ord.MAP.logistic[[i]],test[[i]]$rating),categories)
}

conf.mat.Ord.MAP.probit<-list()
for(i in 1:10)
{
  conf.mat.Ord.MAP.probit[[i]]<-mat.square(table(pred.Ord.MAP.probit[[i]],test[[i]]$rating),categories)
}

conf.mat.Ord.MAP.loglog<-list()
for(i in 1:10)
{
  conf.mat.Ord.MAP.loglog[[i]]<-mat.square(table(pred.Ord.MAP.loglog[[i]],test[[i]]$rating),categories)
}

conf.mat.Ord.MAP.cloglog<-list()
for(i in 1:10)
{
  conf.mat.Ord.MAP.cloglog[[i]]<-mat.square(table(pred.Ord.MAP.cloglog[[i]],test[[i]]$rating),categories)
}

conf.mat.Ord.MAP.cauchit<-list()
for(i in 1:10)
{
  conf.mat.Ord.MAP.cauchit[[i]]<-mat.square(table(pred.Ord.MAP.cauchit[[i]],test[[i]]$rating),categories)
}


# MAE 

MAE.Ord.MAP.logistic<-vector()
for (i in 1:10)
{
  MAE.Ord.MAP.logistic[i]<-MAE(conf.mat.Ord.MAP.logistic[[i]])
}

MAE.Ord.MAP.probit<-vector()
for (i in 1:10)
{
  MAE.Ord.MAP.probit[i]<-MAE(conf.mat.Ord.MAP.probit[[i]])
}

MAE.Ord.MAP.loglog<-vector()
for (i in 1:10)
{
  MAE.Ord.MAP.loglog[i]<-MAE(conf.mat.Ord.MAP.loglog[[i]])
}

MAE.Ord.MAP.cloglog<-vector()
for (i in 1:10)
{
  MAE.Ord.MAP.cloglog[i]<-MAE(conf.mat.Ord.MAP.cloglog[[i]])
}

MAE.Ord.MAP.cauchit<-vector()
for (i in 1:10)
{
  MAE.Ord.MAP.cauchit[i]<-MAE(conf.mat.Ord.MAP.cauchit[[i]])
}


################################################################################
###. Comparison MAP vs Ord.MAP
###
###.  Tests MAE


p.val.test.MAE.MAP.greater.MAE.ORD.MAP.logistic<-wilcox.test(MAE.MAP.logistic,MAE.Ord.MAP.logistic,paired = TRUE,alternative="greater")$p.value
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.probit<-wilcox.test(MAE.MAP.probit,MAE.Ord.MAP.probit,paired = TRUE,alternative="greater")$p.value
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.loglog<-wilcox.test(MAE.MAP.loglog,MAE.Ord.MAP.loglog,paired = TRUE,alternative="greater")$p.value
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.cloglog<-wilcox.test(MAE.MAP.cloglog,MAE.Ord.MAP.cloglog,paired = TRUE,alternative="greater")$p.value
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.cauchit<-wilcox.test(MAE.MAP.cauchit,MAE.Ord.MAP.cauchit,paired = TRUE,alternative="greater")$p.value

p.val.test.MAE.MAP.greater.MAE.ORD.MAP.logistic.less<-wilcox.test(MAE.MAP.logistic,MAE.Ord.MAP.logistic,paired = TRUE,alternative="less")$p.value
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.probit.less<-wilcox.test(MAE.MAP.probit,MAE.Ord.MAP.probit,paired = TRUE,alternative="less")$p.value
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.loglog.less<-wilcox.test(MAE.MAP.loglog,MAE.Ord.MAP.loglog,paired = TRUE,alternative="less")$p.value
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.cloglog.less<-wilcox.test(MAE.MAP.cloglog,MAE.Ord.MAP.cloglog,paired = TRUE,alternative="less")$p.value
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.cauchit.less<-wilcox.test(MAE.MAP.cauchit,MAE.Ord.MAP.cauchit,paired = TRUE,alternative="less")$p.value



p.val.test.MAE.MAP.greater.MAE.ORD.MAP.logistic
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.probit
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.loglog
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.cloglog
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.cauchit

p.val.test.MAE.MAP.greater.MAE.ORD.MAP.logistic.less
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.probit.less
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.loglog.less
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.cloglog.less
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.cauchit.less


###############################
###############################
###   Boxplots   MAE

MAE.results<-as.data.frame(cbind(MAE.MAP.logistic-MAE.Ord.MAP.logistic,
                                 MAE.MAP.probit-MAE.Ord.MAP.probit,
                                 MAE.MAP.loglog-MAE.Ord.MAP.loglog,
                                 MAE.MAP.cloglog-MAE.Ord.MAP.cloglog,
                                 MAE.MAP.cauchit-MAE.Ord.MAP.cauchit))

boxplot(MAE.results,
        main=paste("(b) Wine dataset"), 
        xlab=" ", ylab="MAE increment (MAP - Ord.MAP)",names=c("logistic", "probit","loglog","cloglog","cauchit"),
        col="white")
abline(h = 0, col = 2, lwd = 2, lty=2)
# Dots
stripchart(MAE.results,
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
# We’ll now fit the Proportional Odds Logistic Regression model 
# using polr function from the MASS package.
#

model.logistic<-list()
for (i in 1:10)
{model.logistic[[i]] <- polr(Class ~ age + sex + chest_pain + trestbps + chol +
                               fbs + restecg + thalach + exang + oldpeak, data = training[[i]], Hess = TRUE)
}

model.probit<-list()
for (i in 1:10)
{model.probit[[i]] <- polr(Class ~ age + sex + chest_pain + trestbps + chol +
                             fbs + restecg + thalach + exang + oldpeak, data = training[[i]], Hess = TRUE,method="probit")
}

model.loglog<-list()
for (i in 1:10)
{model.loglog[[i]] <- polr(Class ~ age + sex + chest_pain + trestbps + chol +
                             fbs + restecg + thalach + exang + oldpeak, data = training[[i]], Hess = TRUE,method="loglog")
}

model.cloglog<-list()
for (i in 1:10)
{model.cloglog[[i]] <- polr(Class ~ age + sex + chest_pain + trestbps + chol +
                              fbs + restecg + thalach + exang + oldpeak, data = training[[i]], Hess = TRUE,method="cloglog")
}

model.cauchit<-list()
for (i in 1:10)
{model.cauchit[[i]] <- polr(Class ~ age + sex + chest_pain + trestbps + chol +
                              fbs + restecg + thalach + exang + oldpeak, data = training[[i]], Hess = TRUE,method="cauchit")
}

# polr method = "logistic" by default. Other: "probit", "loglog", "cloglog", "cauchit"

# summary(model_fit[[2]])


################################################################################
#
# Predictions of any of the test set
# 
#

pred.test.logistic<-list()
for (i in 1:10)
{pred.test.logistic[[i]] <- predict(model.logistic[[i]],test[[i]][,-11],type = "p")
}


pred.test.probit<-list()
for (i in 1:10)
{pred.test.probit[[i]] <- predict(model.probit[[i]],test[[i]][,-11],type = "p")
}

pred.test.loglog<-list()
for (i in 1:10)
{pred.test.loglog[[i]] <- predict(model.loglog[[i]],test[[i]][,-11],type = "p")
}

pred.test.cloglog<-list()
for (i in 1:10)
{pred.test.cloglog[[i]] <- predict(model.cloglog[[i]],test[[i]][,-11],type = "p")
}

pred.test.cauchit<-list()
for (i in 1:10)
{pred.test.cauchit[[i]] <- predict(model.cauchit[[i]],test[[i]][,-11],type = "p")
}



################################################################################
### predictions with standard MAP criterion:
###

pred.MAP.logistic<-list()
for (i in 1:10)
{pred.MAP.logistic[[i]]<-vector()
for (j in 1:dim(pred.test.logistic[[i]])[1])
{pred.MAP.logistic[[i]][j]<-categories[which.max(pred.test.logistic[[i]][j,])]
}
}


pred.MAP.probit<-list()
for (i in 1:10)
{pred.MAP.probit[[i]]<-vector()
for (j in 1:dim(pred.test.probit[[i]])[1])
{pred.MAP.probit[[i]][j]<-categories[which.max(pred.test.probit[[i]][j,])]
}
}

pred.MAP.loglog<-list()
for (i in 1:10)
{pred.MAP.loglog[[i]]<-vector()
for (j in 1:dim(pred.test.loglog[[i]])[1])
{pred.MAP.loglog[[i]][j]<-categories[which.max(pred.test.loglog[[i]][j,])]
}
}

pred.MAP.cloglog<-list()
for (i in 1:10)
{pred.MAP.cloglog[[i]]<-vector()
for (j in 1:dim(pred.test.cloglog[[i]])[1])
{pred.MAP.cloglog[[i]][j]<-categories[which.max(pred.test.cloglog[[i]][j,])]
}
}

pred.MAP.cauchit<-list()
for (i in 1:10)
{pred.MAP.cauchit[[i]]<-vector()
for (j in 1:dim(pred.test.cauchit[[i]])[1])
{pred.MAP.cauchit[[i]][j]<-categories[which.max(pred.test.cauchit[[i]][j,])]
}
}


# confusion matrices


conf.mat.MAP.logistic<-list()
for(i in 1:10)
{
  conf.mat.MAP.logistic[[i]]<-mat.square(table(pred.MAP.logistic[[i]],test[[i]]$Class),categories)
}

conf.mat.MAP.probit<-list()
for(i in 1:10)
{
  conf.mat.MAP.probit[[i]]<-mat.square(table(pred.MAP.probit[[i]],test[[i]]$Class),categories)
}

conf.mat.MAP.loglog<-list()
for(i in 1:10)
{
  conf.mat.MAP.loglog[[i]]<-mat.square(table(pred.MAP.loglog[[i]],test[[i]]$Class),categories)
}

conf.mat.MAP.cloglog<-list()
for(i in 1:10)
{
  conf.mat.MAP.cloglog[[i]]<-mat.square(table(pred.MAP.cloglog[[i]],test[[i]]$Class),categories)
}

conf.mat.MAP.cauchit<-list()
for(i in 1:10)
{
  conf.mat.MAP.cauchit[[i]]<-mat.square(table(pred.MAP.cauchit[[i]],test[[i]]$Class),categories)
}


# MAE 

MAE.MAP.logistic<-vector()
for (i in 1:10)
{
  MAE.MAP.logistic[i]<-MAE(conf.mat.MAP.logistic[[i]])
}

MAE.MAP.probit<-vector()
for (i in 1:10)
{
  MAE.MAP.probit[i]<-MAE(conf.mat.MAP.probit[[i]])
}

MAE.MAP.loglog<-vector()
for (i in 1:10)
{
  MAE.MAP.loglog[i]<-MAE(conf.mat.MAP.loglog[[i]])
}

MAE.MAP.cloglog<-vector()
for (i in 1:10)
{
  MAE.MAP.cloglog[i]<-MAE(conf.mat.MAP.cloglog[[i]])
}

MAE.MAP.cauchit<-vector()
for (i in 1:10)
{
  MAE.MAP.cauchit[i]<-MAE(conf.mat.MAP.cauchit[[i]])
}




################################################################################
### predictions with Ord.MAP criterion:
###

pred.Ord.MAP.logistic<-list()
for (i in 1:10)
{ pred.Ord.MAP.logistic[[i]]<-vector()
for (j in 1:dim(pred.test.logistic[[i]])[1])
{h<-min(which(cumsum(pred.test.logistic[[i]][j,])>=0.5))
pred.Ord.MAP.logistic[[i]][j]<-categories[h]
}
}


pred.Ord.MAP.probit<-list()
for (i in 1:10)
{ pred.Ord.MAP.probit[[i]]<-vector()
for (j in 1:dim(pred.test.probit[[i]])[1])
{h<-min(which(cumsum(pred.test.probit[[i]][j,])>=0.5))
pred.Ord.MAP.probit[[i]][j]<-categories[h]
}
}

pred.Ord.MAP.loglog<-list()
for (i in 1:10)
{ pred.Ord.MAP.loglog[[i]]<-vector()
for (j in 1:dim(pred.test.loglog[[i]])[1])
{h<-min(which(cumsum(pred.test.loglog[[i]][j,])>=0.5))
pred.Ord.MAP.loglog[[i]][j]<-categories[h]
}
}

pred.Ord.MAP.cloglog<-list()
for (i in 1:10)
{ pred.Ord.MAP.cloglog[[i]]<-vector()
for (j in 1:dim(pred.test.cloglog[[i]])[1])
{h<-min(which(cumsum(pred.test.cloglog[[i]][j,])>=0.5))
pred.Ord.MAP.cloglog[[i]][j]<-categories[h]
}
}

pred.Ord.MAP.cauchit<-list()
for (i in 1:10)
{ pred.Ord.MAP.cauchit[[i]]<-vector()
for (j in 1:dim(pred.test.cauchit[[i]])[1])
{h<-min(which(cumsum(pred.test.cauchit[[i]][j,])>=0.5))
pred.Ord.MAP.cauchit[[i]][j]<-categories[h]
}
}



# confusion matrices

conf.mat.Ord.MAP.logistic<-list()
for(i in 1:10)
{
  conf.mat.Ord.MAP.logistic[[i]]<-mat.square(table(pred.Ord.MAP.logistic[[i]],test[[i]]$Class),categories)
}

conf.mat.Ord.MAP.probit<-list()
for(i in 1:10)
{
  conf.mat.Ord.MAP.probit[[i]]<-mat.square(table(pred.Ord.MAP.probit[[i]],test[[i]]$Class),categories)
}

conf.mat.Ord.MAP.loglog<-list()
for(i in 1:10)
{
  conf.mat.Ord.MAP.loglog[[i]]<-mat.square(table(pred.Ord.MAP.loglog[[i]],test[[i]]$Class),categories)
}

conf.mat.Ord.MAP.cloglog<-list()
for(i in 1:10)
{
  conf.mat.Ord.MAP.cloglog[[i]]<-mat.square(table(pred.Ord.MAP.cloglog[[i]],test[[i]]$Class),categories)
}

conf.mat.Ord.MAP.cauchit<-list()
for(i in 1:10)
{
  conf.mat.Ord.MAP.cauchit[[i]]<-mat.square(table(pred.Ord.MAP.cauchit[[i]],test[[i]]$Class),categories)
}


# MAE 

MAE.Ord.MAP.logistic<-vector()
for (i in 1:10)
{
  MAE.Ord.MAP.logistic[i]<-MAE(conf.mat.Ord.MAP.logistic[[i]])
}

MAE.Ord.MAP.probit<-vector()
for (i in 1:10)
{
  MAE.Ord.MAP.probit[i]<-MAE(conf.mat.Ord.MAP.probit[[i]])
}

MAE.Ord.MAP.loglog<-vector()
for (i in 1:10)
{
  MAE.Ord.MAP.loglog[i]<-MAE(conf.mat.Ord.MAP.loglog[[i]])
}

MAE.Ord.MAP.cloglog<-vector()
for (i in 1:10)
{
  MAE.Ord.MAP.cloglog[i]<-MAE(conf.mat.Ord.MAP.cloglog[[i]])
}

MAE.Ord.MAP.cauchit<-vector()
for (i in 1:10)
{
  MAE.Ord.MAP.cauchit[i]<-MAE(conf.mat.Ord.MAP.cauchit[[i]])
}


################################################################################
###. Comparison MAP vs Ord.MAP
###
###.  Tests MAE


p.val.test.MAE.MAP.greater.MAE.ORD.MAP.logistic<-wilcox.test(MAE.MAP.logistic,MAE.Ord.MAP.logistic,paired = TRUE,alternative="greater")$p.value
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.probit<-wilcox.test(MAE.MAP.probit,MAE.Ord.MAP.probit,paired = TRUE,alternative="greater")$p.value
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.loglog<-wilcox.test(MAE.MAP.loglog,MAE.Ord.MAP.loglog,paired = TRUE,alternative="greater")$p.value
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.cloglog<-wilcox.test(MAE.MAP.cloglog,MAE.Ord.MAP.cloglog,paired = TRUE,alternative="greater")$p.value
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.cauchit<-wilcox.test(MAE.MAP.cauchit,MAE.Ord.MAP.cauchit,paired = TRUE,alternative="greater")$p.value

p.val.test.MAE.MAP.greater.MAE.ORD.MAP.logistic.less<-wilcox.test(MAE.MAP.logistic,MAE.Ord.MAP.logistic,paired = TRUE,alternative="less")$p.value
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.probit.less<-wilcox.test(MAE.MAP.probit,MAE.Ord.MAP.probit,paired = TRUE,alternative="less")$p.value
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.loglog.less<-wilcox.test(MAE.MAP.loglog,MAE.Ord.MAP.loglog,paired = TRUE,alternative="less")$p.value
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.cloglog.less<-wilcox.test(MAE.MAP.cloglog,MAE.Ord.MAP.cloglog,paired = TRUE,alternative="less")$p.value
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.cauchit.less<-wilcox.test(MAE.MAP.cauchit,MAE.Ord.MAP.cauchit,paired = TRUE,alternative="less")$p.value



p.val.test.MAE.MAP.greater.MAE.ORD.MAP.logistic
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.probit
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.loglog
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.cloglog
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.cauchit

p.val.test.MAE.MAP.greater.MAE.ORD.MAP.logistic.less
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.probit.less
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.loglog.less
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.cloglog.less
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.cauchit.less


###############################
###############################
###   Boxplots   MAE

MAE.results<-as.data.frame(cbind(MAE.MAP.logistic-MAE.Ord.MAP.logistic,
                                 MAE.MAP.probit-MAE.Ord.MAP.probit,
                                 MAE.MAP.loglog-MAE.Ord.MAP.loglog,
                                 MAE.MAP.cloglog-MAE.Ord.MAP.cloglog,
                                 MAE.MAP.cauchit-MAE.Ord.MAP.cauchit))

boxplot(MAE.results,
        main=paste("(c) Hearth dataset"), 
        xlab=" ", ylab="MAE increment (MAP - Ord.MAP)",names=c("logistic", "probit","loglog","cloglog","cauchit"),
        col="white")
abline(h = 0, col = 2, lwd = 2, lty=2)
# Dots
stripchart(MAE.results,
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
View(parkinsons)
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
# We’ll now fit the Proportional Odds Logistic Regression model 
# using polr function from the MASS package.
#
#### FIRST: OUTPUT VARIABLE TO PREDICT: V5.bin.num.factor

model.logistic<-list()
for (i in 1:10)
{model.logistic[[i]] <- polr(V5.bin.num.factor ~ V2 + V3 + V4 + V7 + V8 + V9 + V10 +
                               V11 + V12 + V13 + V14 + V15 + V16 + V17 + V18 + V19 + 
                               V20 + V21 +V22, data = training[[i]], Hess = TRUE)
}

model.probit<-list()
for (i in 1:10)
{model.probit[[i]] <- polr(V5.bin.num.factor ~ V2 + V3 + V4 + V7 + V8 + V9 + V10 +
                             V11 + V12 + V13 + V14 + V15 + V16 + V17 + V18 + V19 + 
                             V20 + V21 +V22, data = training[[i]], Hess = TRUE,method="probit")
}

model.loglog<-list()  
for (i in 1:10)
{
  model.loglog[[i]] <- polr(V5.bin.num.factor ~ V2 + V3 + V4 + V7 + V8 + V9 + V10 +
                             V11 + V12 + V13 + V14 + V15 + V16 + V17 + V18 + V19 + 
                             V20 + V21 +V22, data = training[[i]], Hess = TRUE,method="loglog")
   }


model.cloglog<-list() 
for (i in 1:10)
{model.cloglog[[i]] <- polr(V5.bin.num.factor ~ V2 + V3 + V4 + V7 + V8 + V9 + V10 +
                              V11 + V12 + V13 + V14 + V15 + V16 + V17 + V18 + V19 + 
                              V20 + V21 +V22, data = training[[i]], Hess = TRUE,method="cloglog")
}

model.cauchit<-list()
for (i in 1:10)
{model.cauchit[[i]] <- polr(V5.bin.num.factor ~ V2 + V3 + V4 + V7 + V8 + V9 + V10 +
                              V11 + V12 + V13 + V14 + V15 + V16 + V17 + V18 + V19 + 
                              V20 + V21 +V22, data = training[[i]], Hess = TRUE,method="cauchit")
}

# polr method = "logistic" by default. Other: "probit", "loglog", "cloglog", "cauchit"
# some of the link functions do not work!!


################################################################################
#
# Predictions of any of the test set
# 
#

vector.no<-c(1,5,6,23:28)

pred.test.logistic<-list()
for (i in 1:10)
{pred.test.logistic[[i]] <- predict(model.logistic[[i]],test[[i]][,-vector.no],type = "p")
}


pred.test.probit<-list()
for (i in 1:10)
{pred.test.probit[[i]] <- predict(model.probit[[i]],test[[i]][,-vector.no],type = "p")
}

pred.test.loglog<-list()
for (i in 1:10)
{pred.test.loglog[[i]] <- predict(model.loglog[[i]],test[[i]][,-vector.no],type = "p")
}

pred.test.cloglog<-list()
for (i in 1:10)
{pred.test.cloglog[[i]] <- predict(model.cloglog[[i]],test[[i]][,-vector.no],type = "p")
}

pred.test.cauchit<-list()
for (i in 1:10)
{pred.test.cauchit[[i]] <- predict(model.cauchit[[i]],test[[i]][,-vector.no],type = "p")
}



################################################################################
### predictions with standard MAP criterion:
###

pred.MAP.logistic<-list()
for (i in 1:10)
{pred.MAP.logistic[[i]]<-vector()
for (j in 1:dim(pred.test.logistic[[i]])[1])
{pred.MAP.logistic[[i]][j]<-categories[which.max(pred.test.logistic[[i]][j,])]
}
}


pred.MAP.probit<-list()
for (i in 1:10)
{pred.MAP.probit[[i]]<-vector()
for (j in 1:dim(pred.test.probit[[i]])[1])
{pred.MAP.probit[[i]][j]<-categories[which.max(pred.test.probit[[i]][j,])]
}
}

pred.MAP.loglog<-list()
for (i in 1:10)
{pred.MAP.loglog[[i]]<-vector()
for (j in 1:dim(pred.test.loglog[[i]])[1])
{pred.MAP.loglog[[i]][j]<-categories[which.max(pred.test.loglog[[i]][j,])]
}
}

pred.MAP.cloglog<-list()
for (i in 1:10)
{pred.MAP.cloglog[[i]]<-vector()
for (j in 1:dim(pred.test.cloglog[[i]])[1])
{pred.MAP.cloglog[[i]][j]<-categories[which.max(pred.test.cloglog[[i]][j,])]
}
}

pred.MAP.cauchit<-list()
for (i in 1:10)
{pred.MAP.cauchit[[i]]<-vector()
for (j in 1:dim(pred.test.cauchit[[i]])[1])
{pred.MAP.cauchit[[i]][j]<-categories[which.max(pred.test.cauchit[[i]][j,])]
}
}


# confusion matrices


conf.mat.MAP.logistic<-list()
for(i in 1:10)
{
  conf.mat.MAP.logistic[[i]]<-mat.square(table(pred.MAP.logistic[[i]],test[[i]]$V5.bin.num.factor),categories)
}

conf.mat.MAP.probit<-list()
for(i in 1:10)
{
  conf.mat.MAP.probit[[i]]<-mat.square(table(pred.MAP.probit[[i]],test[[i]]$V5.bin.num.factor),categories)
}

conf.mat.MAP.loglog<-list()
for(i in 1:10)
{
  conf.mat.MAP.loglog[[i]]<-mat.square(table(pred.MAP.loglog[[i]],test[[i]]$V5.bin.num.factor),categories)
}

conf.mat.MAP.cloglog<-list()
for(i in 1:10)
{
  conf.mat.MAP.cloglog[[i]]<-mat.square(table(pred.MAP.cloglog[[i]],test[[i]]$V5.bin.num.factor),categories)
}

conf.mat.MAP.cauchit<-list()
for(i in 1:10)
{
  conf.mat.MAP.cauchit[[i]]<-mat.square(table(pred.MAP.cauchit[[i]],test[[i]]$V5.bin.num.factor),categories)
}


# MAE 

MAE.MAP.logistic<-vector()
for (i in 1:10)
{
  MAE.MAP.logistic[i]<-MAE(conf.mat.MAP.logistic[[i]])
}

MAE.MAP.probit<-vector()
for (i in 1:10)
{
  MAE.MAP.probit[i]<-MAE(conf.mat.MAP.probit[[i]])
}

MAE.MAP.loglog<-vector()
for (i in 1:10)
{
  MAE.MAP.loglog[i]<-MAE(conf.mat.MAP.loglog[[i]])
}

MAE.MAP.cloglog<-vector()
for (i in 1:10)
{
  MAE.MAP.cloglog[i]<-MAE(conf.mat.MAP.cloglog[[i]])
}

MAE.MAP.cauchit<-vector()
for (i in 1:10)
{
  MAE.MAP.cauchit[i]<-MAE(conf.mat.MAP.cauchit[[i]])
}




################################################################################
### predictions with Ord.MAP criterion:
###

pred.Ord.MAP.logistic<-list()
for (i in 1:10)
{ pred.Ord.MAP.logistic[[i]]<-vector()
for (j in 1:dim(pred.test.logistic[[i]])[1])
{h<-min(which(cumsum(pred.test.logistic[[i]][j,])>=0.5))
pred.Ord.MAP.logistic[[i]][j]<-categories[h]
}
}


pred.Ord.MAP.probit<-list()
for (i in 1:10)
{ pred.Ord.MAP.probit[[i]]<-vector()
for (j in 1:dim(pred.test.probit[[i]])[1])
{h<-min(which(cumsum(pred.test.probit[[i]][j,])>=0.5))
pred.Ord.MAP.probit[[i]][j]<-categories[h]
}
}

pred.Ord.MAP.loglog<-list()
for (i in 1:10)
{ pred.Ord.MAP.loglog[[i]]<-vector()
for (j in 1:dim(pred.test.loglog[[i]])[1])
{h<-min(which(cumsum(pred.test.loglog[[i]][j,])>=0.5))
pred.Ord.MAP.loglog[[i]][j]<-categories[h]
}
}

pred.Ord.MAP.cloglog<-list()
for (i in 1:10)
{ pred.Ord.MAP.cloglog[[i]]<-vector()
for (j in 1:dim(pred.test.cloglog[[i]])[1])
{h<-min(which(cumsum(pred.test.cloglog[[i]][j,])>=0.5))
pred.Ord.MAP.cloglog[[i]][j]<-categories[h]
}
}

pred.Ord.MAP.cauchit<-list()
for (i in 1:10)
{ pred.Ord.MAP.cauchit[[i]]<-vector()
for (j in 1:dim(pred.test.cauchit[[i]])[1])
{h<-min(which(cumsum(pred.test.cauchit[[i]][j,])>=0.5))
pred.Ord.MAP.cauchit[[i]][j]<-categories[h]
}
}



# confusion matrices

conf.mat.Ord.MAP.logistic<-list()
for(i in 1:10)
{
  conf.mat.Ord.MAP.logistic[[i]]<-mat.square(table(pred.Ord.MAP.logistic[[i]],test[[i]]$V5.bin.num.factor),categories)
}

conf.mat.Ord.MAP.probit<-list()
for(i in 1:10)
{
  conf.mat.Ord.MAP.probit[[i]]<-mat.square(table(pred.Ord.MAP.probit[[i]],test[[i]]$V5.bin.num.factor),categories)
}

conf.mat.Ord.MAP.loglog<-list()
for(i in 1:10)
{
  conf.mat.Ord.MAP.loglog[[i]]<-mat.square(table(pred.Ord.MAP.loglog[[i]],test[[i]]$V5.bin.num.factor),categories)
}

conf.mat.Ord.MAP.cloglog<-list()
for(i in 1:10)
{
  conf.mat.Ord.MAP.cloglog[[i]]<-mat.square(table(pred.Ord.MAP.cloglog[[i]],test[[i]]$V5.bin.num.factor),categories)
}

conf.mat.Ord.MAP.cauchit<-list()
for(i in 1:10)
{
  conf.mat.Ord.MAP.cauchit[[i]]<-mat.square(table(pred.Ord.MAP.cauchit[[i]],test[[i]]$V5.bin.num.factor),categories)
}


# MAE 

MAE.Ord.MAP.logistic<-vector()
for (i in 1:10)
{
  MAE.Ord.MAP.logistic[i]<-MAE(conf.mat.Ord.MAP.logistic[[i]])
}

MAE.Ord.MAP.probit<-vector()
for (i in 1:10)
{
  MAE.Ord.MAP.probit[i]<-MAE(conf.mat.Ord.MAP.probit[[i]])
}

MAE.Ord.MAP.loglog<-vector()
for (i in 1:10)
{
  MAE.Ord.MAP.loglog[i]<-MAE(conf.mat.Ord.MAP.loglog[[i]])
}

MAE.Ord.MAP.cloglog<-vector()
for (i in 1:10)
{
  MAE.Ord.MAP.cloglog[i]<-MAE(conf.mat.Ord.MAP.cloglog[[i]])
}

MAE.Ord.MAP.cauchit<-vector()
for (i in 1:10)
{
  MAE.Ord.MAP.cauchit[i]<-MAE(conf.mat.Ord.MAP.cauchit[[i]])
}


################################################################################
###. Comparison MAP vs Ord.MAP
###
###.  Tests MAE


p.val.test.MAE.MAP.greater.MAE.ORD.MAP.logistic<-wilcox.test(MAE.MAP.logistic,MAE.Ord.MAP.logistic,paired = TRUE,alternative="greater")$p.value
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.probit<-wilcox.test(MAE.MAP.probit,MAE.Ord.MAP.probit,paired = TRUE,alternative="greater")$p.value
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.loglog<-wilcox.test(MAE.MAP.loglog,MAE.Ord.MAP.loglog,paired = TRUE,alternative="greater")$p.value
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.cloglog<-wilcox.test(MAE.MAP.cloglog,MAE.Ord.MAP.cloglog,paired = TRUE,alternative="greater")$p.value
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.cauchit<-wilcox.test(MAE.MAP.cauchit,MAE.Ord.MAP.cauchit,paired = TRUE,alternative="greater")$p.value

p.val.test.MAE.MAP.greater.MAE.ORD.MAP.logistic.less<-wilcox.test(MAE.MAP.logistic,MAE.Ord.MAP.logistic,paired = TRUE,alternative="less")$p.value
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.probit.less<-wilcox.test(MAE.MAP.probit,MAE.Ord.MAP.probit,paired = TRUE,alternative="less")$p.value
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.loglog.less<-wilcox.test(MAE.MAP.loglog,MAE.Ord.MAP.loglog,paired = TRUE,alternative="less")$p.value
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.cloglog.less<-wilcox.test(MAE.MAP.cloglog,MAE.Ord.MAP.cloglog,paired = TRUE,alternative="less")$p.value
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.cauchit.less<-wilcox.test(MAE.MAP.cauchit,MAE.Ord.MAP.cauchit,paired = TRUE,alternative="less")$p.value



p.val.test.MAE.MAP.greater.MAE.ORD.MAP.logistic
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.probit
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.loglog
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.cloglog
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.cauchit

p.val.test.MAE.MAP.greater.MAE.ORD.MAP.logistic.less
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.probit.less
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.loglog.less
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.cloglog.less
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.cauchit.less


###############################
###############################
###   Boxplots   MAE

MAE.results<-as.data.frame(cbind(MAE.MAP.logistic-MAE.Ord.MAP.logistic,
                                 #MAE.MAP.probit-MAE.Ord.MAP.probit,
                                 #MAE.MAP.loglog-MAE.Ord.MAP.loglog,
                                 #MAE.MAP.cloglog-MAE.Ord.MAP.cloglog,
                                 MAE.MAP.cauchit-MAE.Ord.MAP.cauchit))

boxplot(MAE.results,
        main=paste("(d) Parkinson dataset: output V5 (6 classes)"), 
        xlab=" ", ylab="MAE increment (MAP - Ord.MAP)",names=c("logistic", 
                                                               #"probit",
                                                               #"loglog","cloglog",
                                                               "cauchit"),
        col="white")
abline(h = 0, col = 2, lwd = 2, lty=2)
# Dots
stripchart(MAE.results,
           #method = "jitter",
           pch = 19,
           # col = 1:5,
           col = 1:2,
           vertical = TRUE,
           add = TRUE)


#
##
###
##
#
#### SECOND: OUTPUT VARIABLE TO PREDICT: V6.bin.num.factor

model.logistic<-list()
for (i in 1:10)
{model.logistic[[i]] <- polr(V6.bin.num.factor ~ V2 + V3 + V4 + V7 + V8 + V9 + V10 +
                               V11 + V12 + V13 + V14 + V15 + V16 + V17 + V18 + V19 + 
                               V20 + V21 +V22, data = training[[i]], Hess = TRUE)
}

model.probit<-list()
for (i in 1:10)
{model.probit[[i]] <- polr(V6.bin.num.factor ~ V2 + V3 + V4 + V7 + V8 + V9 + V10 +
                             V11 + V12 + V13 + V14 + V15 + V16 + V17 + V18 + V19 + 
                             V20 + V21 +V22, data = training[[i]], Hess = TRUE,method="probit")
}

model.loglog<-list()  
for (i in 1:10)
{
  model.loglog[[i]] <- polr(V6.bin.num.factor ~ V2 + V3 + V4 + V7 + V8 + V9 + V10 +
                              V11 + V12 + V13 + V14 + V15 + V16 + V17 + V18 + V19 + 
                              V20 + V21 +V22, data = training[[i]], Hess = TRUE,method="loglog")
}


model.cloglog<-list()  
for (i in 1:10)
{model.cloglog[[i]] <- polr(V6.bin.num.factor ~ V2 + V3 + V4 + V7 + V8 + V9 + V10 +
                              V11 + V12 + V13 + V14 + V15 + V16 + V17 + V18 + V19 + 
                              V20 + V21 +V22, data = training[[i]], Hess = TRUE,method="cloglog")
}

model.cauchit<-list()
for (i in 1:10)
{model.cauchit[[i]] <- polr(V6.bin.num.factor ~ V2 + V3 + V4 + V7 + V8 + V9 + V10 +
                              V11 + V12 + V13 + V14 + V15 + V16 + V17 + V18 + V19 + 
                              V20 + V21 +V22, data = training[[i]], Hess = TRUE,method="cauchit")
}

# polr method = "logistic" by default. Other: "probit", "loglog", "cloglog", "cauchit"
# some of the link functions do not work!!


################################################################################
#
# Predictions of any of the test set
# 
#

vector.no<-c(1,5,6,23:28)

pred.test.logistic<-list()
for (i in 1:10)
{pred.test.logistic[[i]] <- predict(model.logistic[[i]],test[[i]][,-vector.no],type = "p")
}


pred.test.probit<-list()
for (i in 1:10)
{pred.test.probit[[i]] <- predict(model.probit[[i]],test[[i]][,-vector.no],type = "p")
}

pred.test.loglog<-list()
for (i in 1:10)
{pred.test.loglog[[i]] <- predict(model.loglog[[i]],test[[i]][,-vector.no],type = "p")
}

pred.test.cloglog<-list()
for (i in 1:10)
{pred.test.cloglog[[i]] <- predict(model.cloglog[[i]],test[[i]][,-vector.no],type = "p")
}

pred.test.cauchit<-list()
for (i in 1:10)
{pred.test.cauchit[[i]] <- predict(model.cauchit[[i]],test[[i]][,-vector.no],type = "p")
}



################################################################################
### predictions with standard MAP criterion:
###

pred.MAP.logistic<-list()
for (i in 1:10)
{pred.MAP.logistic[[i]]<-vector()
for (j in 1:dim(pred.test.logistic[[i]])[1])
{pred.MAP.logistic[[i]][j]<-categories[which.max(pred.test.logistic[[i]][j,])]
}
}


pred.MAP.probit<-list()
for (i in 1:10)
{pred.MAP.probit[[i]]<-vector()
for (j in 1:dim(pred.test.probit[[i]])[1])
{pred.MAP.probit[[i]][j]<-categories[which.max(pred.test.probit[[i]][j,])]
}
}

pred.MAP.loglog<-list()
for (i in 1:10)
{pred.MAP.loglog[[i]]<-vector()
for (j in 1:dim(pred.test.loglog[[i]])[1])
{pred.MAP.loglog[[i]][j]<-categories[which.max(pred.test.loglog[[i]][j,])]
}
}

pred.MAP.cloglog<-list()
for (i in 1:10)
{pred.MAP.cloglog[[i]]<-vector()
for (j in 1:dim(pred.test.cloglog[[i]])[1])
{pred.MAP.cloglog[[i]][j]<-categories[which.max(pred.test.cloglog[[i]][j,])]
}
}

pred.MAP.cauchit<-list()
for (i in 1:10)
{pred.MAP.cauchit[[i]]<-vector()
for (j in 1:dim(pred.test.cauchit[[i]])[1])
{pred.MAP.cauchit[[i]][j]<-categories[which.max(pred.test.cauchit[[i]][j,])]
}
}


# confusion matrices


conf.mat.MAP.logistic<-list()
for(i in 1:10)
{
  conf.mat.MAP.logistic[[i]]<-mat.square(table(pred.MAP.logistic[[i]],test[[i]]$V6.bin.num.factor),categories)
}

conf.mat.MAP.probit<-list()
for(i in 1:10)
{
  conf.mat.MAP.probit[[i]]<-mat.square(table(pred.MAP.probit[[i]],test[[i]]$V6.bin.num.factor),categories)
}

conf.mat.MAP.loglog<-list()
for(i in 1:10)
{
  conf.mat.MAP.loglog[[i]]<-mat.square(table(pred.MAP.loglog[[i]],test[[i]]$V6.bin.num.factor),categories)
}

conf.mat.MAP.cloglog<-list()
for(i in 1:10)
{
  conf.mat.MAP.cloglog[[i]]<-mat.square(table(pred.MAP.cloglog[[i]],test[[i]]$V6.bin.num.factor),categories)
}

conf.mat.MAP.cauchit<-list()
for(i in 1:10)
{
  conf.mat.MAP.cauchit[[i]]<-mat.square(table(pred.MAP.cauchit[[i]],test[[i]]$V6.bin.num.factor),categories)
}


# MAE 

MAE.MAP.logistic<-vector()
for (i in 1:10)
{
  MAE.MAP.logistic[i]<-MAE(conf.mat.MAP.logistic[[i]])
}

MAE.MAP.probit<-vector()
for (i in 1:10)
{
  MAE.MAP.probit[i]<-MAE(conf.mat.MAP.probit[[i]])
}

MAE.MAP.loglog<-vector()
for (i in 1:10)
{
  MAE.MAP.loglog[i]<-MAE(conf.mat.MAP.loglog[[i]])
}

MAE.MAP.cloglog<-vector()
for (i in 1:10)
{
  MAE.MAP.cloglog[i]<-MAE(conf.mat.MAP.cloglog[[i]])
}

MAE.MAP.cauchit<-vector()
for (i in 1:10)
{
  MAE.MAP.cauchit[i]<-MAE(conf.mat.MAP.cauchit[[i]])
}




################################################################################
### predictions with Ord.MAP criterion:
###

pred.Ord.MAP.logistic<-list()
for (i in 1:10)
{ pred.Ord.MAP.logistic[[i]]<-vector()
for (j in 1:dim(pred.test.logistic[[i]])[1])
{h<-min(which(cumsum(pred.test.logistic[[i]][j,])>=0.5))
pred.Ord.MAP.logistic[[i]][j]<-categories[h]
}
}


pred.Ord.MAP.probit<-list()
for (i in 1:10)
{ pred.Ord.MAP.probit[[i]]<-vector()
for (j in 1:dim(pred.test.probit[[i]])[1])
{h<-min(which(cumsum(pred.test.probit[[i]][j,])>=0.5))
pred.Ord.MAP.probit[[i]][j]<-categories[h]
}
}

pred.Ord.MAP.loglog<-list()
for (i in 1:10)
{ pred.Ord.MAP.loglog[[i]]<-vector()
for (j in 1:dim(pred.test.loglog[[i]])[1])
{h<-min(which(cumsum(pred.test.loglog[[i]][j,])>=0.5))
pred.Ord.MAP.loglog[[i]][j]<-categories[h]
}
}

pred.Ord.MAP.cloglog<-list()
for (i in 1:10)
{ pred.Ord.MAP.cloglog[[i]]<-vector()
for (j in 1:dim(pred.test.cloglog[[i]])[1])
{h<-min(which(cumsum(pred.test.cloglog[[i]][j,])>=0.5))
pred.Ord.MAP.cloglog[[i]][j]<-categories[h]
}
}

pred.Ord.MAP.cauchit<-list()
for (i in 1:10)
{ pred.Ord.MAP.cauchit[[i]]<-vector()
for (j in 1:dim(pred.test.cauchit[[i]])[1])
{h<-min(which(cumsum(pred.test.cauchit[[i]][j,])>=0.5))
pred.Ord.MAP.cauchit[[i]][j]<-categories[h]
}
}



# confusion matrices

conf.mat.Ord.MAP.logistic<-list()
for(i in 1:10)
{
  conf.mat.Ord.MAP.logistic[[i]]<-mat.square(table(pred.Ord.MAP.logistic[[i]],test[[i]]$V6.bin.num.factor),categories)
}

conf.mat.Ord.MAP.probit<-list()
for(i in 1:10)
{
  conf.mat.Ord.MAP.probit[[i]]<-mat.square(table(pred.Ord.MAP.probit[[i]],test[[i]]$V6.bin.num.factor),categories)
}

conf.mat.Ord.MAP.loglog<-list()
for(i in 1:10)
{
  conf.mat.Ord.MAP.loglog[[i]]<-mat.square(table(pred.Ord.MAP.loglog[[i]],test[[i]]$V6.bin.num.factor),categories)
}

conf.mat.Ord.MAP.cloglog<-list()
for(i in 1:10)
{
  conf.mat.Ord.MAP.cloglog[[i]]<-mat.square(table(pred.Ord.MAP.cloglog[[i]],test[[i]]$V6.bin.num.factor),categories)
}

conf.mat.Ord.MAP.cauchit<-list()
for(i in 1:10)
{
  conf.mat.Ord.MAP.cauchit[[i]]<-mat.square(table(pred.Ord.MAP.cauchit[[i]],test[[i]]$V6.bin.num.factor),categories)
}


# MAE 

MAE.Ord.MAP.logistic<-vector()
for (i in 1:10)
{
  MAE.Ord.MAP.logistic[i]<-MAE(conf.mat.Ord.MAP.logistic[[i]])
}

MAE.Ord.MAP.probit<-vector()
for (i in 1:10)
{
  MAE.Ord.MAP.probit[i]<-MAE(conf.mat.Ord.MAP.probit[[i]])
}

MAE.Ord.MAP.loglog<-vector()
for (i in 1:10)
{
  MAE.Ord.MAP.loglog[i]<-MAE(conf.mat.Ord.MAP.loglog[[i]])
}

MAE.Ord.MAP.cloglog<-vector()
for (i in 1:10)
{
  MAE.Ord.MAP.cloglog[i]<-MAE(conf.mat.Ord.MAP.cloglog[[i]])
}

MAE.Ord.MAP.cauchit<-vector()
for (i in 1:10)
{
  MAE.Ord.MAP.cauchit[i]<-MAE(conf.mat.Ord.MAP.cauchit[[i]])
}


################################################################################
###. Comparison MAP vs Ord.MAP
###
###.  Tests MAE


p.val.test.MAE.MAP.greater.MAE.ORD.MAP.logistic<-wilcox.test(MAE.MAP.logistic,MAE.Ord.MAP.logistic,paired = TRUE,alternative="greater")$p.value
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.probit<-wilcox.test(MAE.MAP.probit,MAE.Ord.MAP.probit,paired = TRUE,alternative="greater")$p.value
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.loglog<-wilcox.test(MAE.MAP.loglog,MAE.Ord.MAP.loglog,paired = TRUE,alternative="greater")$p.value
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.cloglog<-wilcox.test(MAE.MAP.cloglog,MAE.Ord.MAP.cloglog,paired = TRUE,alternative="greater")$p.value
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.cauchit<-wilcox.test(MAE.MAP.cauchit,MAE.Ord.MAP.cauchit,paired = TRUE,alternative="greater")$p.value

p.val.test.MAE.MAP.greater.MAE.ORD.MAP.logistic.less<-wilcox.test(MAE.MAP.logistic,MAE.Ord.MAP.logistic,paired = TRUE,alternative="less")$p.value
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.probit.less<-wilcox.test(MAE.MAP.probit,MAE.Ord.MAP.probit,paired = TRUE,alternative="less")$p.value
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.loglog.less<-wilcox.test(MAE.MAP.loglog,MAE.Ord.MAP.loglog,paired = TRUE,alternative="less")$p.value
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.cloglog.less<-wilcox.test(MAE.MAP.cloglog,MAE.Ord.MAP.cloglog,paired = TRUE,alternative="less")$p.value
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.cauchit.less<-wilcox.test(MAE.MAP.cauchit,MAE.Ord.MAP.cauchit,paired = TRUE,alternative="less")$p.value



p.val.test.MAE.MAP.greater.MAE.ORD.MAP.logistic
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.probit
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.loglog
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.cloglog
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.cauchit

p.val.test.MAE.MAP.greater.MAE.ORD.MAP.logistic.less
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.probit.less
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.loglog.less
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.cloglog.less
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.cauchit.less


###############################
###############################
###   Boxplots   MAE

MAE.results<-as.data.frame(cbind(MAE.MAP.logistic-MAE.Ord.MAP.logistic,
                                 MAE.MAP.probit-MAE.Ord.MAP.probit,
                                 MAE.MAP.loglog-MAE.Ord.MAP.loglog,
                                 #MAE.MAP.cloglog-MAE.Ord.MAP.cloglog,
                                 MAE.MAP.cauchit-MAE.Ord.MAP.cauchit))

boxplot(MAE.results,
        main=paste("(d) Parkinson dataset: output V6 (6 classes)"), 
        xlab=" ", ylab="MAE increment (MAP - Ord.MAP)",names=c("logistic", "probit",
                                                               "loglog",
                                                               #"cloglog",
                                                               "cauchit"),
        col="white")
abline(h = 0, col = 2, lwd = 2, lty=2)
# Dots
stripchart(MAE.results,
           #method = "jitter",
           pch = 19,
           # col = 1:5,
           col = 1:4,
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
# We’ll now fit the Proportional Odds Logistic Regression model 
# using polr function from the MASS package.
#

model.logistic<-list()
for (i in 1:10)
{model.logistic[[i]] <- polr(Class ~ . - importance - population, data = training[[i]], Hess = TRUE)
}

model.probit<-list()
for (i in 1:10)
{model.probit[[i]] <- polr(Class ~ . - importance - population, data = training[[i]], Hess = TRUE,method="probit")
}


model.loglog<-list()   
for (i in 1:10)
{
  model.loglog[[i]] <- polr(Class ~ . - importance - population, data = training[[i]], Hess = TRUE,method="loglog")
}


model.cloglog<-list()  
for (i in 1:10)
{model.cloglog[[i]] <- polr(Class ~ . - importance - population, data = training[[i]], Hess = TRUE,method="cloglog")
}

model.cauchit<-list()
for (i in 1:10)
{model.cauchit[[i]] <- polr(Class ~ . - importance - population, data = training[[i]], Hess = TRUE,method="cauchit")
}

# polr method = "logistic" by default. Other: "probit", "loglog", "cloglog", "cauchit"

# summary(model_fit[[2]])


################################################################################
#
# Predictions of any of the test set
# 
#


pred.test.logistic<-list()
for (i in 1:10)
{pred.test.logistic[[i]] <- predict(model.logistic[[i]],test[[i]][,],type = "p")
}


pred.test.probit<-list()
for (i in 1:10)
{pred.test.probit[[i]] <- predict(model.probit[[i]],test[[i]][,],type = "p")
}

pred.test.loglog<-list()
for (i in 1:10)
{pred.test.loglog[[i]] <- predict(model.loglog[[i]],test[[i]][,],type = "p")
}

pred.test.cloglog<-list()
for (i in 1:10)
{pred.test.cloglog[[i]] <- predict(model.cloglog[[i]],test[[i]][,],type = "p")
}

pred.test.cauchit<-list()
for (i in 1:10)
{pred.test.cauchit[[i]] <- predict(model.cauchit[[i]],test[[i]][,],type = "p")
}



################################################################################
### predictions with standard MAP criterion:
###

pred.MAP.logistic<-list()
for (i in 1:10)
{pred.MAP.logistic[[i]]<-vector()
for (j in 1:dim(pred.test.logistic[[i]])[1])
{pred.MAP.logistic[[i]][j]<-categories[which.max(pred.test.logistic[[i]][j,])]
}
}


pred.MAP.probit<-list()
for (i in 1:10)
{pred.MAP.probit[[i]]<-vector()
for (j in 1:dim(pred.test.probit[[i]])[1])
{pred.MAP.probit[[i]][j]<-categories[which.max(pred.test.probit[[i]][j,])]
}
}

pred.MAP.loglog<-list()
for (i in 1:10)
{pred.MAP.loglog[[i]]<-vector()
for (j in 1:dim(pred.test.loglog[[i]])[1])
{pred.MAP.loglog[[i]][j]<-categories[which.max(pred.test.loglog[[i]][j,])]
}
}

pred.MAP.cloglog<-list()
for (i in 1:10)
{pred.MAP.cloglog[[i]]<-vector()
for (j in 1:dim(pred.test.cloglog[[i]])[1])
{pred.MAP.cloglog[[i]][j]<-categories[which.max(pred.test.cloglog[[i]][j,])]
}
}

pred.MAP.cauchit<-list()
for (i in 1:10)
{pred.MAP.cauchit[[i]]<-vector()
for (j in 1:dim(pred.test.cauchit[[i]])[1])
{pred.MAP.cauchit[[i]][j]<-categories[which.max(pred.test.cauchit[[i]][j,])]
}
}


# confusion matrices


conf.mat.MAP.logistic<-list()
for(i in 1:10)
{
  conf.mat.MAP.logistic[[i]]<-mat.square(table(pred.MAP.logistic[[i]],test[[i]]$Class),categories)
}

conf.mat.MAP.probit<-list()
for(i in 1:10)
{
  conf.mat.MAP.probit[[i]]<-mat.square(table(pred.MAP.probit[[i]],test[[i]]$Class),categories)
}

conf.mat.MAP.loglog<-list()
for(i in 1:10)
{
  conf.mat.MAP.loglog[[i]]<-mat.square(table(pred.MAP.loglog[[i]],test[[i]]$Class),categories)
}

conf.mat.MAP.cloglog<-list()
for(i in 1:10)
{
  conf.mat.MAP.cloglog[[i]]<-mat.square(table(pred.MAP.cloglog[[i]],test[[i]]$Class),categories)
}

conf.mat.MAP.cauchit<-list()
for(i in 1:10)
{
  conf.mat.MAP.cauchit[[i]]<-mat.square(table(pred.MAP.cauchit[[i]],test[[i]]$Class),categories)
}


# MAE 

MAE.MAP.logistic<-vector()
for (i in 1:10)
{
  MAE.MAP.logistic[i]<-MAE(conf.mat.MAP.logistic[[i]])
}

MAE.MAP.probit<-vector()
for (i in 1:10)
{
  MAE.MAP.probit[i]<-MAE(conf.mat.MAP.probit[[i]])
}

MAE.MAP.loglog<-vector()
for (i in 1:10)
{
  MAE.MAP.loglog[i]<-MAE(conf.mat.MAP.loglog[[i]])
}

MAE.MAP.cloglog<-vector()
for (i in 1:10)
{
  MAE.MAP.cloglog[i]<-MAE(conf.mat.MAP.cloglog[[i]])
}

MAE.MAP.cauchit<-vector()
for (i in 1:10)
{
  MAE.MAP.cauchit[i]<-MAE(conf.mat.MAP.cauchit[[i]])
}




################################################################################
### predictions with Ord.MAP criterion:
###

pred.Ord.MAP.logistic<-list()
for (i in 1:10)
{ pred.Ord.MAP.logistic[[i]]<-vector()
for (j in 1:dim(pred.test.logistic[[i]])[1])
{h<-min(which(cumsum(pred.test.logistic[[i]][j,])>=0.5))
pred.Ord.MAP.logistic[[i]][j]<-categories[h]
}
}


pred.Ord.MAP.probit<-list()
for (i in 1:10)
{ pred.Ord.MAP.probit[[i]]<-vector()
for (j in 1:dim(pred.test.probit[[i]])[1])
{h<-min(which(cumsum(pred.test.probit[[i]][j,])>=0.5))
pred.Ord.MAP.probit[[i]][j]<-categories[h]
}
}

pred.Ord.MAP.loglog<-list()
for (i in 1:10)
{ pred.Ord.MAP.loglog[[i]]<-vector()
for (j in 1:dim(pred.test.loglog[[i]])[1])
{h<-min(which(cumsum(pred.test.loglog[[i]][j,])>=0.5))
pred.Ord.MAP.loglog[[i]][j]<-categories[h]
}
}

pred.Ord.MAP.cloglog<-list()
for (i in 1:10)
{ pred.Ord.MAP.cloglog[[i]]<-vector()
for (j in 1:dim(pred.test.cloglog[[i]])[1])
{h<-min(which(cumsum(pred.test.cloglog[[i]][j,])>=0.5))
pred.Ord.MAP.cloglog[[i]][j]<-categories[h]
}
}

pred.Ord.MAP.cauchit<-list()
for (i in 1:10)
{ pred.Ord.MAP.cauchit[[i]]<-vector()
for (j in 1:dim(pred.test.cauchit[[i]])[1])
{h<-min(which(cumsum(pred.test.cauchit[[i]][j,])>=0.5))
pred.Ord.MAP.cauchit[[i]][j]<-categories[h]
}
}



# confusion matrices

conf.mat.Ord.MAP.logistic<-list()
for(i in 1:10)
{
  conf.mat.Ord.MAP.logistic[[i]]<-mat.square(table(pred.Ord.MAP.logistic[[i]],test[[i]]$Class),categories)
}

conf.mat.Ord.MAP.probit<-list()
for(i in 1:10)
{
  conf.mat.Ord.MAP.probit[[i]]<-mat.square(table(pred.Ord.MAP.probit[[i]],test[[i]]$Class),categories)
}

conf.mat.Ord.MAP.loglog<-list()
for(i in 1:10)
{
  conf.mat.Ord.MAP.loglog[[i]]<-mat.square(table(pred.Ord.MAP.loglog[[i]],test[[i]]$Class),categories)
}

conf.mat.Ord.MAP.cloglog<-list()
for(i in 1:10)
{
  conf.mat.Ord.MAP.cloglog[[i]]<-mat.square(table(pred.Ord.MAP.cloglog[[i]],test[[i]]$Class),categories)
}

conf.mat.Ord.MAP.cauchit<-list()
for(i in 1:10)
{
  conf.mat.Ord.MAP.cauchit[[i]]<-mat.square(table(pred.Ord.MAP.cauchit[[i]],test[[i]]$Class),categories)
}


# MAE 

MAE.Ord.MAP.logistic<-vector()
for (i in 1:10)
{
  MAE.Ord.MAP.logistic[i]<-MAE(conf.mat.Ord.MAP.logistic[[i]])
}

MAE.Ord.MAP.probit<-vector()
for (i in 1:10)
{
  MAE.Ord.MAP.probit[i]<-MAE(conf.mat.Ord.MAP.probit[[i]])
}

MAE.Ord.MAP.loglog<-vector()
for (i in 1:10)
{
  MAE.Ord.MAP.loglog[i]<-MAE(conf.mat.Ord.MAP.loglog[[i]])
}

MAE.Ord.MAP.cloglog<-vector()
for (i in 1:10)
{
  MAE.Ord.MAP.cloglog[i]<-MAE(conf.mat.Ord.MAP.cloglog[[i]])
}

MAE.Ord.MAP.cauchit<-vector()
for (i in 1:10)
{
  MAE.Ord.MAP.cauchit[i]<-MAE(conf.mat.Ord.MAP.cauchit[[i]])
}


################################################################################
###. Comparison MAP vs Ord.MAP
###
###.  Tests MAE


p.val.test.MAE.MAP.greater.MAE.ORD.MAP.logistic<-wilcox.test(MAE.MAP.logistic,MAE.Ord.MAP.logistic,paired = TRUE,alternative="greater")$p.value
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.probit<-wilcox.test(MAE.MAP.probit,MAE.Ord.MAP.probit,paired = TRUE,alternative="greater")$p.value
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.loglog<-wilcox.test(MAE.MAP.loglog,MAE.Ord.MAP.loglog,paired = TRUE,alternative="greater")$p.value
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.cloglog<-wilcox.test(MAE.MAP.cloglog,MAE.Ord.MAP.cloglog,paired = TRUE,alternative="greater")$p.value
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.cauchit<-wilcox.test(MAE.MAP.cauchit,MAE.Ord.MAP.cauchit,paired = TRUE,alternative="greater")$p.value

p.val.test.MAE.MAP.greater.MAE.ORD.MAP.logistic.less<-wilcox.test(MAE.MAP.logistic,MAE.Ord.MAP.logistic,paired = TRUE,alternative="less")$p.value
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.probit.less<-wilcox.test(MAE.MAP.probit,MAE.Ord.MAP.probit,paired = TRUE,alternative="less")$p.value
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.loglog.less<-wilcox.test(MAE.MAP.loglog,MAE.Ord.MAP.loglog,paired = TRUE,alternative="less")$p.value
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.cloglog.less<-wilcox.test(MAE.MAP.cloglog,MAE.Ord.MAP.cloglog,paired = TRUE,alternative="less")$p.value
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.cauchit.less<-wilcox.test(MAE.MAP.cauchit,MAE.Ord.MAP.cauchit,paired = TRUE,alternative="less")$p.value



p.val.test.MAE.MAP.greater.MAE.ORD.MAP.logistic
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.probit
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.loglog
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.cloglog
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.cauchit

p.val.test.MAE.MAP.greater.MAE.ORD.MAP.logistic.less
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.probit.less
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.loglog.less
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.cloglog.less
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.cauchit.less


###############################
###############################
###   Boxplots   MAE

MAE.results<-as.data.frame(cbind(MAE.MAP.logistic-MAE.Ord.MAP.logistic,
                                 MAE.MAP.probit-MAE.Ord.MAP.probit,
                                 MAE.MAP.loglog-MAE.Ord.MAP.loglog,
                                 MAE.MAP.cloglog-MAE.Ord.MAP.cloglog,
                                 MAE.MAP.cauchit-MAE.Ord.MAP.cauchit))

boxplot(MAE.results,
        main=paste("(e) CES11 dataset"), 
        xlab=" ", ylab="MAE increment (MAP - Ord.MAP)",names=c("logistic", 
                                                               "probit",
                                                               "loglog",
                                                               "cloglog",
                                                               "cauchit"),
        col="white")
abline(h = 0, col = 2, lwd = 2, lty=2)
# Dots
stripchart(MAE.results,
           #method = "jitter",
           pch = 19,
           col = 1:5,
           vertical = TRUE,
           add = TRUE)


