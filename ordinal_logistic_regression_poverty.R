
###########################################
###########################################
#
#  EXPERIMENTAL PHASE 
#
# Ordinal Logistic Regression
# using polr function from the MASS package.
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


#
#

library(carData)
library(MASS)
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
# Fitting the Proportional Odds Logistic Regression model 
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

p.normality.MAE.logistic<-shapiro.test(MAE.MAP.logistic-MAE.Ord.MAP.logistic)$p.value

if (p.normality.MAE.logistic >= 0.05)
{p.val.test.MAE.MAP.greater.MAE.ORD.MAP.logistic<-t.test(MAE.MAP.logistic,MAE.Ord.MAP.logistic,paired = TRUE,alternative="greater")$p.value
} else {
  p.val.test.MAE.MAP.greater.MAE.ORD.MAP.logistic<-wilcox.test(MAE.MAP.logistic,MAE.Ord.MAP.logistic,paired = TRUE,alternative="greater")$p.value
}


p.normality.MAE.probit<-shapiro.test(MAE.MAP.probit-MAE.Ord.MAP.probit)$p.value

if (p.normality.MAE.probit >= 0.05)
{p.val.test.MAE.MAP.greater.MAE.ORD.MAP.probit<-t.test(MAE.MAP.probit,MAE.Ord.MAP.probit,paired = TRUE,alternative="greater")$p.value
} else {
  p.val.test.MAE.MAP.greater.MAE.ORD.MAP.probit<-wilcox.test(MAE.MAP.probit,MAE.Ord.MAP.probit,paired = TRUE,alternative="greater")$p.value
}


p.normality.MAE.loglog<-shapiro.test(MAE.MAP.loglog-MAE.Ord.MAP.loglog)$p.value

if (p.normality.MAE.loglog >= 0.05)
{p.val.test.MAE.MAP.greater.MAE.ORD.MAP.loglog<-t.test(MAE.MAP.loglog,MAE.Ord.MAP.loglog,paired = TRUE,alternative="greater")$p.value
} else {
  p.val.test.MAE.MAP.greater.MAE.ORD.MAP.loglog<-wilcox.test(MAE.MAP.loglog,MAE.Ord.MAP.loglog,paired = TRUE,alternative="greater")$p.value
}


p.normality.MAE.cloglog<-shapiro.test(MAE.MAP.cloglog-MAE.Ord.MAP.cloglog)$p.value

if (p.normality.MAE.cloglog >= 0.05)
{p.val.test.MAE.MAP.greater.MAE.ORD.MAP.cloglog<-t.test(MAE.MAP.cloglog,MAE.Ord.MAP.cloglog,paired = TRUE,alternative="greater")$p.value
} else {
  p.val.test.MAE.MAP.greater.MAE.ORD.MAP.cloglog<-wilcox.test(MAE.MAP.cloglog,MAE.Ord.MAP.cloglog,paired = TRUE,alternative="greater")$p.value
}


p.normality.MAE.cauchit<-shapiro.test(MAE.MAP.cauchit-MAE.Ord.MAP.cauchit)$p.value

if (p.normality.MAE.cauchit >= 0.05)
{p.val.test.MAE.MAP.greater.MAE.ORD.MAP.cauchit<-t.test(MAE.MAP.cauchit,MAE.Ord.MAP.cauchit,paired = TRUE,alternative="greater")$p.value
} else {
  p.val.test.MAE.MAP.greater.MAE.ORD.MAP.cauchit<-wilcox.test(MAE.MAP.cauchit,MAE.Ord.MAP.cauchit,paired = TRUE,alternative="greater")$p.value
}




p.val.test.MAE.MAP.greater.MAE.ORD.MAP.logistic
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.probit
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.loglog
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.cloglog
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.cauchit


###############################
###############################
###   Boxplots   MAE

MAE.results<-as.data.frame(cbind(MAE.MAP.logistic-MAE.Ord.MAP.logistic,
                           MAE.MAP.probit-MAE.Ord.MAP.probit,
                           MAE.MAP.loglog-MAE.Ord.MAP.loglog,
                           MAE.MAP.cloglog-MAE.Ord.MAP.cloglog,
                           MAE.MAP.cauchit-MAE.Ord.MAP.cauchit))

boxplot(MAE.results,
        main=paste(" "), 
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


#########.  END