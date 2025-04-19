
###########################################
###########################################
#
#  EXPERIMENTAL PHASE (Section 5) 
# 
################################################################################
#
# PROCEDURE 2): cumulative link (mixed) models (CLM)
# also known as ordered regression models, proportional odds models, proportional
# hazards models for grouped survival times and ordered logit/probit/...
# models, with the "ordinal" package of R.
#
# Estimation is via maximum likelihood and mixed models are fitted
# with the Laplace approximation and adaptive Gauss-Hermite quadrature.
#



library(ordinal)

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

ata(WVS) 
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
# clm, from "ordinal" package, rfits cumulative link models (CLMs) such as the propotional odds model. 
# The model allows for various link functions and structured thresholds that restricts 
# the thresholds or cut-points to be e.g.,
# equidistant or symmetrically arranged around the central threshold(s). 
#
# clm(formula, scale, nominal, data, weights, start, subset, doFit = TRUE,
#     na.action, contrasts, model = TRUE, control=list(),
#     link = c("logit", "probit", "cloglog", "loglog", "cauchit",
#              "Aranda-Ordaz", "log-gamma"),
#     threshold = c("flexible", "symmetric", "symmetric2", "equidistant"), ...)

model.clm.logit.flexible<-list()
for (i in 1:10)
{model.clm.logit.flexible[[i]] <- clm(poverty~religion+degree+country+age+gender, data = training[[i]],
                       link="logit", threshold="flexible")}


model.clm.logit.symmetric<-list()
for (i in 1:10)
{model.clm.logit.symmetric[[i]] <- clm(poverty~religion+degree+country+age+gender, data = training[[i]],
                                      link="logit", threshold="symmetric")}


model.clm.logit.symmetric2<-list()
for (i in 1:10)
{model.clm.logit.symmetric2[[i]] <- clm(poverty~religion+degree+country+age+gender, data = training[[i]],
                                       link="logit", threshold="symmetric2")}

model.clm.logit.equidistant<-list()
for (i in 1:10)
{model.clm.logit.equidistant[[i]] <- clm(poverty~religion+degree+country+age+gender, data = training[[i]],
                                       link="logit", threshold="equidistant")}


###

model.clm.probit.flexible<-list()
for (i in 1:10)
{model.clm.probit.flexible[[i]] <- clm(poverty~religion+degree+country+age+gender, data = training[[i]],
                                      link="probit", threshold="flexible")}


model.clm.probit.symmetric<-list()
for (i in 1:10)
{model.clm.probit.symmetric[[i]] <- clm(poverty~religion+degree+country+age+gender, data = training[[i]],
                                       link="probit", threshold="symmetric")}


model.clm.probit.symmetric2<-list()
for (i in 1:10)
{model.clm.probit.symmetric2[[i]] <- clm(poverty~religion+degree+country+age+gender, data = training[[i]],
                                        link="probit", threshold="symmetric2")}

model.clm.probit.equidistant<-list()
for (i in 1:10)
{model.clm.probit.equidistant[[i]] <- clm(poverty~religion+degree+country+age+gender, data = training[[i]],
                                         link="probit", threshold="equidistant")}


###

model.clm.loglog.flexible<-list()
for (i in 1:10)
{model.clm.loglog.flexible[[i]] <- clm(poverty~religion+degree+country+age+gender, data = training[[i]],
                                      link="loglog", threshold="flexible")}


model.clm.loglog.symmetric<-list()
for (i in 1:10)
{model.clm.loglog.symmetric[[i]] <- clm(poverty~religion+degree+country+age+gender, data = training[[i]],
                                       link="loglog", threshold="symmetric")}


model.clm.loglog.symmetric2<-list()
for (i in 1:10)
{model.clm.loglog.symmetric2[[i]] <- clm(poverty~religion+degree+country+age+gender, data = training[[i]],
                                        link="loglog", threshold="symmetric2")}

model.clm.loglog.equidistant<-list()
for (i in 1:10)
{model.clm.loglog.equidistant[[i]] <- clm(poverty~religion+degree+country+age+gender, data = training[[i]],
                                         link="loglog", threshold="equidistant")}


###

model.clm.cloglog.flexible<-list()
for (i in 1:10)
{model.clm.cloglog.flexible[[i]] <- clm(poverty~religion+degree+country+age+gender, data = training[[i]],
                                       link="cloglog", threshold="flexible")}


model.clm.cloglog.symmetric<-list()
for (i in 1:10)
{model.clm.cloglog.symmetric[[i]] <- clm(poverty~religion+degree+country+age+gender, data = training[[i]],
                                        link="cloglog", threshold="symmetric")}


model.clm.cloglog.symmetric2<-list()
for (i in 1:10)
{model.clm.cloglog.symmetric2[[i]] <- clm(poverty~religion+degree+country+age+gender, data = training[[i]],
                                         link="cloglog", threshold="symmetric2")}

model.clm.cloglog.equidistant<-list()
for (i in 1:10)
{model.clm.cloglog.equidistant[[i]] <- clm(poverty~religion+degree+country+age+gender, data = training[[i]],
                                          link="cloglog", threshold="equidistant")}


###

model.clm.cauchit.flexible<-list()
for (i in 1:10)
{model.clm.cauchit.flexible[[i]] <- clm(poverty~religion+degree+country+age+gender, data = training[[i]],
                                        link="cauchit", threshold="flexible")}


model.clm.cauchit.symmetric<-list()
for (i in 1:10)
{model.clm.cauchit.symmetric[[i]] <- clm(poverty~religion+degree+country+age+gender, data = training[[i]],
                                         link="cauchit", threshold="symmetric")}


model.clm.cauchit.symmetric2<-list()
for (i in 1:10)
{model.clm.cauchit.symmetric2[[i]] <- clm(poverty~religion+degree+country+age+gender, data = training[[i]],
                                          link="cauchit", threshold="symmetric2")}

model.clm.cauchit.equidistant<-list()
for (i in 1:10)
{model.clm.cauchit.equidistant[[i]] <- clm(poverty~religion+degree+country+age+gender, data = training[[i]],
                                           link="cauchit", threshold="equidistant")}



################################################################################
#
# Predictions of any of the test set
# 
#

pred.test.clm.logit.flexible<-list()
for (i in 1:10)
{pred.test.clm.logit.flexible[[i]] <- predict(model.clm.logit.flexible[[i]],test[[i]][,-1],type = "p")$fit
}

pred.test.clm.logit.symmetric<-list()
for (i in 1:10)
{pred.test.clm.logit.symmetric[[i]] <- predict(model.clm.logit.symmetric[[i]],test[[i]][,-1],type = "p")$fit
}

pred.test.clm.logit.symmetric2<-list()
for (i in 1:10)
{pred.test.clm.logit.symmetric2[[i]] <- predict(model.clm.logit.symmetric2[[i]],test[[i]][,-1],type = "p")$fit
}

pred.test.clm.logit.equidistant<-list()
for (i in 1:10)
{pred.test.clm.logit.equidistant[[i]] <- predict(model.clm.logit.equidistant[[i]],test[[i]][,-1],type = "p")$fit
}

###

pred.test.clm.probit.flexible<-list()
for (i in 1:10)
{pred.test.clm.probit.flexible[[i]] <- predict(model.clm.probit.flexible[[i]],test[[i]][,-1],type = "p")$fit
}

pred.test.clm.probit.symmetric<-list()
for (i in 1:10)
{pred.test.clm.probit.symmetric[[i]] <- predict(model.clm.probit.symmetric[[i]],test[[i]][,-1],type = "p")$fit
}

pred.test.clm.probit.symmetric2<-list()
for (i in 1:10)
{pred.test.clm.probit.symmetric2[[i]] <- predict(model.clm.probit.symmetric2[[i]],test[[i]][,-1],type = "p")$fit
}

pred.test.clm.probit.equidistant<-list()
for (i in 1:10)
{pred.test.clm.probit.equidistant[[i]] <- predict(model.clm.probit.equidistant[[i]],test[[i]][,-1],type = "p")$fit
}

###

pred.test.clm.loglog.flexible<-list()
for (i in 1:10)
{pred.test.clm.loglog.flexible[[i]] <- predict(model.clm.loglog.flexible[[i]],test[[i]][,-1],type = "p")$fit
}

pred.test.clm.loglog.symmetric<-list()
for (i in 1:10)
{pred.test.clm.loglog.symmetric[[i]] <- predict(model.clm.loglog.symmetric[[i]],test[[i]][,-1],type = "p")$fit
}

pred.test.clm.loglog.symmetric2<-list()
for (i in 1:10)
{pred.test.clm.loglog.symmetric2[[i]] <- predict(model.clm.loglog.symmetric2[[i]],test[[i]][,-1],type = "p")$fit
}

pred.test.clm.loglog.equidistant<-list()
for (i in 1:10)
{pred.test.clm.loglog.equidistant[[i]] <- predict(model.clm.loglog.equidistant[[i]],test[[i]][,-1],type = "p")$fit
}

###


pred.test.clm.cloglog.flexible<-list()
for (i in 1:10)
{pred.test.clm.cloglog.flexible[[i]] <- predict(model.clm.cloglog.flexible[[i]],test[[i]][,-1],type = "p")$fit
}

pred.test.clm.cloglog.symmetric<-list()
for (i in 1:10)
{pred.test.clm.cloglog.symmetric[[i]] <- predict(model.clm.cloglog.symmetric[[i]],test[[i]][,-1],type = "p")$fit
}

pred.test.clm.cloglog.symmetric2<-list()
for (i in 1:10)
{pred.test.clm.cloglog.symmetric2[[i]] <- predict(model.clm.cloglog.symmetric2[[i]],test[[i]][,-1],type = "p")$fit
}

pred.test.clm.cloglog.equidistant<-list()
for (i in 1:10)
{pred.test.clm.cloglog.equidistant[[i]] <- predict(model.clm.cloglog.equidistant[[i]],test[[i]][,-1],type = "p")$fit
}

###


pred.test.clm.cauchit.flexible<-list()
for (i in 1:10)
{pred.test.clm.cauchit.flexible[[i]] <- predict(model.clm.cauchit.flexible[[i]],test[[i]][,-1],type = "p")$fit
}

pred.test.clm.cauchit.symmetric<-list()
for (i in 1:10)
{pred.test.clm.cauchit.symmetric[[i]] <- predict(model.clm.cauchit.symmetric[[i]],test[[i]][,-1],type = "p")$fit
}

pred.test.clm.cauchit.symmetric2<-list()
for (i in 1:10)
{pred.test.clm.cauchit.symmetric2[[i]] <- predict(model.clm.cauchit.symmetric2[[i]],test[[i]][,-1],type = "p")$fit
}

pred.test.clm.cauchit.equidistant<-list()
for (i in 1:10)
{pred.test.clm.cauchit.equidistant[[i]] <- predict(model.clm.cauchit.equidistant[[i]],test[[i]][,-1],type = "p")$fit
}

###

################################################################################
### predictions with standard MAP criterion:
###

pred.MAP.clm.logit.flexible<-list()
for (i in 1:10)
{pred.MAP.clm.logit.flexible[[i]]<-vector()
 for (j in 1:dim(pred.test.clm.logit.flexible[[i]])[1])
  {pred.MAP.clm.logit.flexible[[i]][j]<-categories[which.max(pred.test.clm.logit.flexible[[i]][j,])]
  }
}

pred.MAP.clm.logit.symmetric<-list()
for (i in 1:10)
{pred.MAP.clm.logit.symmetric[[i]]<-vector()
for (j in 1:dim(pred.test.clm.logit.symmetric[[i]])[1])
{pred.MAP.clm.logit.symmetric[[i]][j]<-categories[which.max(pred.test.clm.logit.symmetric[[i]][j,])]
}
}

pred.MAP.clm.logit.symmetric2<-list()
for (i in 1:10)
{pred.MAP.clm.logit.symmetric2[[i]]<-vector()
for (j in 1:dim(pred.test.clm.logit.symmetric2[[i]])[1])
{pred.MAP.clm.logit.symmetric2[[i]][j]<-categories[which.max(pred.test.clm.logit.symmetric2[[i]][j,])]
}
}

pred.MAP.clm.logit.equidistant<-list()
for (i in 1:10)
{pred.MAP.clm.logit.equidistant[[i]]<-vector()
for (j in 1:dim(pred.test.clm.logit.equidistant[[i]])[1])
{pred.MAP.clm.logit.equidistant[[i]][j]<-categories[which.max(pred.test.clm.logit.equidistant[[i]][j,])]
}
}

###


pred.MAP.clm.probit.flexible<-list()
for (i in 1:10)
{pred.MAP.clm.probit.flexible[[i]]<-vector()
for (j in 1:dim(pred.test.clm.probit.flexible[[i]])[1])
{pred.MAP.clm.probit.flexible[[i]][j]<-categories[which.max(pred.test.clm.probit.flexible[[i]][j,])]
}
}

pred.MAP.clm.probit.symmetric<-list()
for (i in 1:10)
{pred.MAP.clm.probit.symmetric[[i]]<-vector()
for (j in 1:dim(pred.test.clm.probit.symmetric[[i]])[1])
{pred.MAP.clm.probit.symmetric[[i]][j]<-categories[which.max(pred.test.clm.probit.symmetric[[i]][j,])]
}
}

pred.MAP.clm.probit.symmetric2<-list()
for (i in 1:10)
{pred.MAP.clm.probit.symmetric2[[i]]<-vector()
for (j in 1:dim(pred.test.clm.probit.symmetric2[[i]])[1])
{pred.MAP.clm.probit.symmetric2[[i]][j]<-categories[which.max(pred.test.clm.probit.symmetric2[[i]][j,])]
}
}

pred.MAP.clm.probit.equidistant<-list()
for (i in 1:10)
{pred.MAP.clm.probit.equidistant[[i]]<-vector()
for (j in 1:dim(pred.test.clm.probit.equidistant[[i]])[1])
{pred.MAP.clm.probit.equidistant[[i]][j]<-categories[which.max(pred.test.clm.probit.equidistant[[i]][j,])]
}
}

###


pred.MAP.clm.loglog.flexible<-list()
for (i in 1:10)
{pred.MAP.clm.loglog.flexible[[i]]<-vector()
for (j in 1:dim(pred.test.clm.loglog.flexible[[i]])[1])
{pred.MAP.clm.loglog.flexible[[i]][j]<-categories[which.max(pred.test.clm.loglog.flexible[[i]][j,])]
}
}

pred.MAP.clm.loglog.symmetric<-list()
for (i in 1:10)
{pred.MAP.clm.loglog.symmetric[[i]]<-vector()
for (j in 1:dim(pred.test.clm.loglog.symmetric[[i]])[1])
{pred.MAP.clm.loglog.symmetric[[i]][j]<-categories[which.max(pred.test.clm.loglog.symmetric[[i]][j,])]
}
}

pred.MAP.clm.loglog.symmetric2<-list()
for (i in 1:10)
{pred.MAP.clm.loglog.symmetric2[[i]]<-vector()
for (j in 1:dim(pred.test.clm.loglog.symmetric2[[i]])[1])
{pred.MAP.clm.loglog.symmetric2[[i]][j]<-categories[which.max(pred.test.clm.loglog.symmetric2[[i]][j,])]
}
}

pred.MAP.clm.loglog.equidistant<-list()
for (i in 1:10)
{pred.MAP.clm.loglog.equidistant[[i]]<-vector()
for (j in 1:dim(pred.test.clm.loglog.equidistant[[i]])[1])
{pred.MAP.clm.loglog.equidistant[[i]][j]<-categories[which.max(pred.test.clm.loglog.equidistant[[i]][j,])]
}
}

###


pred.MAP.clm.cloglog.flexible<-list()
for (i in 1:10)
{pred.MAP.clm.cloglog.flexible[[i]]<-vector()
for (j in 1:dim(pred.test.clm.cloglog.flexible[[i]])[1])
{pred.MAP.clm.cloglog.flexible[[i]][j]<-categories[which.max(pred.test.clm.cloglog.flexible[[i]][j,])]
}
}

pred.MAP.clm.cloglog.symmetric<-list()
for (i in 1:10)
{pred.MAP.clm.cloglog.symmetric[[i]]<-vector()
for (j in 1:dim(pred.test.clm.cloglog.symmetric[[i]])[1])
{pred.MAP.clm.cloglog.symmetric[[i]][j]<-categories[which.max(pred.test.clm.cloglog.symmetric[[i]][j,])]
}
}

pred.MAP.clm.cloglog.symmetric2<-list()
for (i in 1:10)
{pred.MAP.clm.cloglog.symmetric2[[i]]<-vector()
for (j in 1:dim(pred.test.clm.cloglog.symmetric2[[i]])[1])
{pred.MAP.clm.cloglog.symmetric2[[i]][j]<-categories[which.max(pred.test.clm.cloglog.symmetric2[[i]][j,])]
}
}

pred.MAP.clm.cloglog.equidistant<-list()
for (i in 1:10)
{pred.MAP.clm.cloglog.equidistant[[i]]<-vector()
for (j in 1:dim(pred.test.clm.cloglog.equidistant[[i]])[1])
{pred.MAP.clm.cloglog.equidistant[[i]][j]<-categories[which.max(pred.test.clm.cloglog.equidistant[[i]][j,])]
}
}

###


pred.MAP.clm.cauchit.flexible<-list()
for (i in 1:10)
{pred.MAP.clm.cauchit.flexible[[i]]<-vector()
for (j in 1:dim(pred.test.clm.cauchit.flexible[[i]])[1])
{pred.MAP.clm.cauchit.flexible[[i]][j]<-categories[which.max(pred.test.clm.cauchit.flexible[[i]][j,])]
}
}

pred.MAP.clm.cauchit.symmetric<-list()
for (i in 1:10)
{pred.MAP.clm.cauchit.symmetric[[i]]<-vector()
for (j in 1:dim(pred.test.clm.cauchit.symmetric[[i]])[1])
{pred.MAP.clm.cauchit.symmetric[[i]][j]<-categories[which.max(pred.test.clm.cauchit.symmetric[[i]][j,])]
}
}

pred.MAP.clm.cauchit.symmetric2<-list()
for (i in 1:10)
{pred.MAP.clm.cauchit.symmetric2[[i]]<-vector()
for (j in 1:dim(pred.test.clm.cauchit.symmetric2[[i]])[1])
{pred.MAP.clm.cauchit.symmetric2[[i]][j]<-categories[which.max(pred.test.clm.cauchit.symmetric2[[i]][j,])]
}
}

pred.MAP.clm.cauchit.equidistant<-list()
for (i in 1:10)
{pred.MAP.clm.cauchit.equidistant[[i]]<-vector()
for (j in 1:dim(pred.test.clm.cauchit.equidistant[[i]])[1])
{pred.MAP.clm.cauchit.equidistant[[i]][j]<-categories[which.max(pred.test.clm.cauchit.equidistant[[i]][j,])]
}
}

# confusion matrices


conf.mat.MAP.clm.logit.flexible<-list()
for(i in 1:10)
{
conf.mat.MAP.clm.logit.flexible[[i]]<-mat.square(table(pred.MAP.clm.logit.flexible[[i]],test[[i]]$poverty),categories)
}


conf.mat.MAP.clm.logit.symmetric<-list()
for(i in 1:10)
{
  conf.mat.MAP.clm.logit.symmetric[[i]]<-mat.square(table(pred.MAP.clm.logit.symmetric[[i]],test[[i]]$poverty),categories)
}


conf.mat.MAP.clm.logit.symmetric2<-list()
for(i in 1:10)
{
  conf.mat.MAP.clm.logit.symmetric2[[i]]<-mat.square(table(pred.MAP.clm.logit.symmetric2[[i]],test[[i]]$poverty),categories)
}


conf.mat.MAP.clm.logit.equidistant<-list()
for(i in 1:10)
{
  conf.mat.MAP.clm.logit.equidistant[[i]]<-mat.square(table(pred.MAP.clm.logit.equidistant[[i]],test[[i]]$poverty),categories)
}

###

conf.mat.MAP.clm.probit.flexible<-list()
for(i in 1:10)
{
  conf.mat.MAP.clm.probit.flexible[[i]]<-mat.square(table(pred.MAP.clm.probit.flexible[[i]],test[[i]]$poverty),categories)
}


conf.mat.MAP.clm.probit.symmetric<-list()
for(i in 1:10)
{
  conf.mat.MAP.clm.probit.symmetric[[i]]<-mat.square(table(pred.MAP.clm.probit.symmetric[[i]],test[[i]]$poverty),categories)
}


conf.mat.MAP.clm.probit.symmetric2<-list()
for(i in 1:10)
{
  conf.mat.MAP.clm.probit.symmetric2[[i]]<-mat.square(table(pred.MAP.clm.probit.symmetric2[[i]],test[[i]]$poverty),categories)
}


conf.mat.MAP.clm.probit.equidistant<-list()
for(i in 1:10)
{
  conf.mat.MAP.clm.probit.equidistant[[i]]<-mat.square(table(pred.MAP.clm.probit.equidistant[[i]],test[[i]]$poverty),categories)
}

###

conf.mat.MAP.clm.loglog.flexible<-list()
for(i in 1:10)
{
  conf.mat.MAP.clm.loglog.flexible[[i]]<-mat.square(table(pred.MAP.clm.loglog.flexible[[i]],test[[i]]$poverty),categories)
}


conf.mat.MAP.clm.loglog.symmetric<-list()
for(i in 1:10)
{
  conf.mat.MAP.clm.loglog.symmetric[[i]]<-mat.square(table(pred.MAP.clm.loglog.symmetric[[i]],test[[i]]$poverty),categories)
}


conf.mat.MAP.clm.loglog.symmetric2<-list()
for(i in 1:10)
{
  conf.mat.MAP.clm.loglog.symmetric2[[i]]<-mat.square(table(pred.MAP.clm.loglog.symmetric2[[i]],test[[i]]$poverty),categories)
}


conf.mat.MAP.clm.loglog.equidistant<-list()
for(i in 1:10)
{
  conf.mat.MAP.clm.loglog.equidistant[[i]]<-mat.square(table(pred.MAP.clm.loglog.equidistant[[i]],test[[i]]$poverty),categories)
}

###

conf.mat.MAP.clm.cloglog.flexible<-list()
for(i in 1:10)
{
  conf.mat.MAP.clm.cloglog.flexible[[i]]<-mat.square(table(pred.MAP.clm.cloglog.flexible[[i]],test[[i]]$poverty),categories)
}


conf.mat.MAP.clm.cloglog.symmetric<-list()
for(i in 1:10)
{
  conf.mat.MAP.clm.cloglog.symmetric[[i]]<-mat.square(table(pred.MAP.clm.cloglog.symmetric[[i]],test[[i]]$poverty),categories)
}


conf.mat.MAP.clm.cloglog.symmetric2<-list()
for(i in 1:10)
{
  conf.mat.MAP.clm.cloglog.symmetric2[[i]]<-mat.square(table(pred.MAP.clm.cloglog.symmetric2[[i]],test[[i]]$poverty),categories)
}


conf.mat.MAP.clm.cloglog.equidistant<-list()
for(i in 1:10)
{
  conf.mat.MAP.clm.cloglog.equidistant[[i]]<-mat.square(table(pred.MAP.clm.cloglog.equidistant[[i]],test[[i]]$poverty),categories)
}

###

conf.mat.MAP.clm.cauchit.flexible<-list()
for(i in 1:10)
{
  conf.mat.MAP.clm.cauchit.flexible[[i]]<-mat.square(table(pred.MAP.clm.cauchit.flexible[[i]],test[[i]]$poverty),categories)
}


conf.mat.MAP.clm.cauchit.symmetric<-list()
for(i in 1:10)
{
  conf.mat.MAP.clm.cauchit.symmetric[[i]]<-mat.square(table(pred.MAP.clm.cauchit.symmetric[[i]],test[[i]]$poverty),categories)
}


conf.mat.MAP.clm.cauchit.symmetric2<-list()
for(i in 1:10)
{
  conf.mat.MAP.clm.cauchit.symmetric2[[i]]<-mat.square(table(pred.MAP.clm.cauchit.symmetric2[[i]],test[[i]]$poverty),categories)
}


conf.mat.MAP.clm.cauchit.equidistant<-list()
for(i in 1:10)
{
  conf.mat.MAP.clm.cauchit.equidistant[[i]]<-mat.square(table(pred.MAP.clm.cauchit.equidistant[[i]],test[[i]]$poverty),categories)
}


# MAE 

MAE.MAP.clm.logit.flexible<-vector()
for (i in 1:10)
{
  MAE.MAP.clm.logit.flexible[i]<-MAE(conf.mat.MAP.clm.logit.flexible[[i]])
}

MAE.MAP.clm.probit.flexible<-vector()
for (i in 1:10)
{
  MAE.MAP.clm.probit.flexible[i]<-MAE(conf.mat.MAP.clm.probit.flexible[[i]])
}

MAE.MAP.clm.loglog.flexible<-vector()
for (i in 1:10)
{
  MAE.MAP.clm.loglog.flexible[i]<-MAE(conf.mat.MAP.clm.loglog.flexible[[i]])
}

MAE.MAP.clm.cloglog.flexible<-vector()
for (i in 1:10)
{
  MAE.MAP.clm.cloglog.flexible[i]<-MAE(conf.mat.MAP.clm.cloglog.flexible[[i]])
}

MAE.MAP.clm.cauchit.flexible<-vector()
for (i in 1:10)
{
  MAE.MAP.clm.cauchit.flexible[i]<-MAE(conf.mat.MAP.clm.cauchit.flexible[[i]])
}

###

MAE.MAP.clm.logit.symmetric<-vector()
for (i in 1:10)
{
  MAE.MAP.clm.logit.symmetric[i]<-MAE(conf.mat.MAP.clm.logit.symmetric[[i]])
}

MAE.MAP.clm.probit.symmetric<-vector()
for (i in 1:10)
{
  MAE.MAP.clm.probit.symmetric[i]<-MAE(conf.mat.MAP.clm.probit.symmetric[[i]])
}

MAE.MAP.clm.loglog.symmetric<-vector()
for (i in 1:10)
{
  MAE.MAP.clm.loglog.symmetric[i]<-MAE(conf.mat.MAP.clm.loglog.symmetric[[i]])
}

MAE.MAP.clm.cloglog.symmetric<-vector()
for (i in 1:10)
{
  MAE.MAP.clm.cloglog.symmetric[i]<-MAE(conf.mat.MAP.clm.cloglog.symmetric[[i]])
}

MAE.MAP.clm.cauchit.symmetric<-vector()
for (i in 1:10)
{
  MAE.MAP.clm.cauchit.symmetric[i]<-MAE(conf.mat.MAP.clm.cauchit.symmetric[[i]])
}

###

MAE.MAP.clm.logit.symmetric2<-vector()
for (i in 1:10)
{
  MAE.MAP.clm.logit.symmetric2[i]<-MAE(conf.mat.MAP.clm.logit.symmetric2[[i]])
}

MAE.MAP.clm.probit.symmetric2<-vector()
for (i in 1:10)
{
  MAE.MAP.clm.probit.symmetric2[i]<-MAE(conf.mat.MAP.clm.probit.symmetric2[[i]])
}

MAE.MAP.clm.loglog.symmetric2<-vector()
for (i in 1:10)
{
  MAE.MAP.clm.loglog.symmetric2[i]<-MAE(conf.mat.MAP.clm.loglog.symmetric2[[i]])
}

MAE.MAP.clm.cloglog.symmetric2<-vector()
for (i in 1:10)
{
  MAE.MAP.clm.cloglog.symmetric2[i]<-MAE(conf.mat.MAP.clm.cloglog.symmetric2[[i]])
}

MAE.MAP.clm.cauchit.symmetric2<-vector()
for (i in 1:10)
{
  MAE.MAP.clm.cauchit.symmetric2[i]<-MAE(conf.mat.MAP.clm.cauchit.symmetric2[[i]])
}

###

MAE.MAP.clm.logit.equidistant<-vector()
for (i in 1:10)
{
  MAE.MAP.clm.logit.equidistant[i]<-MAE(conf.mat.MAP.clm.logit.equidistant[[i]])
}

MAE.MAP.clm.probit.equidistant<-vector()
for (i in 1:10)
{
  MAE.MAP.clm.probit.equidistant[i]<-MAE(conf.mat.MAP.clm.probit.equidistant[[i]])
}

MAE.MAP.clm.loglog.equidistant<-vector()
for (i in 1:10)
{
  MAE.MAP.clm.loglog.equidistant[i]<-MAE(conf.mat.MAP.clm.loglog.equidistant[[i]])
}

MAE.MAP.clm.cloglog.equidistant<-vector()
for (i in 1:10)
{
  MAE.MAP.clm.cloglog.equidistant[i]<-MAE(conf.mat.MAP.clm.cloglog.equidistant[[i]])
}

MAE.MAP.clm.cauchit.equidistant<-vector()
for (i in 1:10)
{
  MAE.MAP.clm.cauchit.equidistant[i]<-MAE(conf.mat.MAP.clm.cauchit.equidistant[[i]])
}



################################################################################
### predictions with Ord.MAP criterion:
###

pred.Ord.MAP.clm.logit.flexible<-list()
for (i in 1:10)
{ pred.Ord.MAP.clm.logit.flexible[[i]]<-vector()
 for (j in 1:dim(pred.test.clm.logit.flexible[[i]])[1])
 {h<-min(which(cumsum(pred.test.clm.logit.flexible[[i]][j,])>=0.5))
  pred.Ord.MAP.clm.logit.flexible[[i]][j]<-categories[h]
}
}

pred.Ord.MAP.clm.logit.symmetric<-list()
for (i in 1:10)
{ pred.Ord.MAP.clm.logit.symmetric[[i]]<-vector()
for (j in 1:dim(pred.test.clm.logit.symmetric[[i]])[1])
{h<-min(which(cumsum(pred.test.clm.logit.symmetric[[i]][j,])>=0.5))
pred.Ord.MAP.clm.logit.symmetric[[i]][j]<-categories[h]
}
}

pred.Ord.MAP.clm.logit.symmetric2<-list()
for (i in 1:10)
{ pred.Ord.MAP.clm.logit.symmetric2[[i]]<-vector()
for (j in 1:dim(pred.test.clm.logit.symmetric2[[i]])[1])
{h<-min(which(cumsum(pred.test.clm.logit.symmetric2[[i]][j,])>=0.5))
pred.Ord.MAP.clm.logit.symmetric2[[i]][j]<-categories[h]
}
}

pred.Ord.MAP.clm.logit.equidistant<-list()
for (i in 1:10)
{ pred.Ord.MAP.clm.logit.equidistant[[i]]<-vector()
for (j in 1:dim(pred.test.clm.logit.equidistant[[i]])[1])
{h<-min(which(cumsum(pred.test.clm.logit.equidistant[[i]][j,])>=0.5))
pred.Ord.MAP.clm.logit.equidistant[[i]][j]<-categories[h]
}
}

###

pred.Ord.MAP.clm.probit.flexible<-list()
for (i in 1:10)
{ pred.Ord.MAP.clm.probit.flexible[[i]]<-vector()
for (j in 1:dim(pred.test.clm.probit.flexible[[i]])[1])
{h<-min(which(cumsum(pred.test.clm.probit.flexible[[i]][j,])>=0.5))
pred.Ord.MAP.clm.probit.flexible[[i]][j]<-categories[h]
}
}

pred.Ord.MAP.clm.probit.symmetric<-list()
for (i in 1:10)
{ pred.Ord.MAP.clm.probit.symmetric[[i]]<-vector()
for (j in 1:dim(pred.test.clm.probit.symmetric[[i]])[1])
{h<-min(which(cumsum(pred.test.clm.probit.symmetric[[i]][j,])>=0.5))
pred.Ord.MAP.clm.probit.symmetric[[i]][j]<-categories[h]
}
}

pred.Ord.MAP.clm.probit.symmetric2<-list()
for (i in 1:10)
{ pred.Ord.MAP.clm.probit.symmetric2[[i]]<-vector()
for (j in 1:dim(pred.test.clm.probit.symmetric2[[i]])[1])
{h<-min(which(cumsum(pred.test.clm.probit.symmetric2[[i]][j,])>=0.5))
pred.Ord.MAP.clm.probit.symmetric2[[i]][j]<-categories[h]
}
}

pred.Ord.MAP.clm.probit.equidistant<-list()
for (i in 1:10)
{ pred.Ord.MAP.clm.probit.equidistant[[i]]<-vector()
for (j in 1:dim(pred.test.clm.probit.equidistant[[i]])[1])
{h<-min(which(cumsum(pred.test.clm.probit.equidistant[[i]][j,])>=0.5))
pred.Ord.MAP.clm.probit.equidistant[[i]][j]<-categories[h]
}
}

###

pred.Ord.MAP.clm.loglog.flexible<-list()
for (i in 1:10)
{ pred.Ord.MAP.clm.loglog.flexible[[i]]<-vector()
for (j in 1:dim(pred.test.clm.loglog.flexible[[i]])[1])
{h<-min(which(cumsum(pred.test.clm.loglog.flexible[[i]][j,])>=0.5))
pred.Ord.MAP.clm.loglog.flexible[[i]][j]<-categories[h]
}
}

pred.Ord.MAP.clm.loglog.symmetric<-list()
for (i in 1:10)
{ pred.Ord.MAP.clm.loglog.symmetric[[i]]<-vector()
for (j in 1:dim(pred.test.clm.loglog.symmetric[[i]])[1])
{h<-min(which(cumsum(pred.test.clm.loglog.symmetric[[i]][j,])>=0.5))
pred.Ord.MAP.clm.loglog.symmetric[[i]][j]<-categories[h]
}
}

pred.Ord.MAP.clm.loglog.symmetric2<-list()
for (i in 1:10)
{ pred.Ord.MAP.clm.loglog.symmetric2[[i]]<-vector()
for (j in 1:dim(pred.test.clm.loglog.symmetric2[[i]])[1])
{h<-min(which(cumsum(pred.test.clm.loglog.symmetric2[[i]][j,])>=0.5))
pred.Ord.MAP.clm.loglog.symmetric2[[i]][j]<-categories[h]
}
}

pred.Ord.MAP.clm.loglog.equidistant<-list()
for (i in 1:10)
{ pred.Ord.MAP.clm.loglog.equidistant[[i]]<-vector()
for (j in 1:dim(pred.test.clm.loglog.equidistant[[i]])[1])
{h<-min(which(cumsum(pred.test.clm.loglog.equidistant[[i]][j,])>=0.5))
pred.Ord.MAP.clm.loglog.equidistant[[i]][j]<-categories[h]
}
}

###

pred.Ord.MAP.clm.cloglog.flexible<-list()
for (i in 1:10)
{ pred.Ord.MAP.clm.cloglog.flexible[[i]]<-vector()
for (j in 1:dim(pred.test.clm.cloglog.flexible[[i]])[1])
{h<-min(which(cumsum(pred.test.clm.cloglog.flexible[[i]][j,])>=0.5))
pred.Ord.MAP.clm.cloglog.flexible[[i]][j]<-categories[h]
}
}

pred.Ord.MAP.clm.cloglog.symmetric<-list()
for (i in 1:10)
{ pred.Ord.MAP.clm.cloglog.symmetric[[i]]<-vector()
for (j in 1:dim(pred.test.clm.cloglog.symmetric[[i]])[1])
{h<-min(which(cumsum(pred.test.clm.cloglog.symmetric[[i]][j,])>=0.5))
pred.Ord.MAP.clm.cloglog.symmetric[[i]][j]<-categories[h]
}
}

pred.Ord.MAP.clm.cloglog.symmetric2<-list()
for (i in 1:10)
{ pred.Ord.MAP.clm.cloglog.symmetric2[[i]]<-vector()
for (j in 1:dim(pred.test.clm.cloglog.symmetric2[[i]])[1])
{h<-min(which(cumsum(pred.test.clm.cloglog.symmetric2[[i]][j,])>=0.5))
pred.Ord.MAP.clm.cloglog.symmetric2[[i]][j]<-categories[h]
}
}

pred.Ord.MAP.clm.cloglog.equidistant<-list()
for (i in 1:10)
{ pred.Ord.MAP.clm.cloglog.equidistant[[i]]<-vector()
for (j in 1:dim(pred.test.clm.cloglog.equidistant[[i]])[1])
{h<-min(which(cumsum(pred.test.clm.cloglog.equidistant[[i]][j,])>=0.5))
pred.Ord.MAP.clm.cloglog.equidistant[[i]][j]<-categories[h]
}
}

###

pred.Ord.MAP.clm.cauchit.flexible<-list()
for (i in 1:10)
{ pred.Ord.MAP.clm.cauchit.flexible[[i]]<-vector()
for (j in 1:dim(pred.test.clm.cauchit.flexible[[i]])[1])
{h<-min(which(cumsum(pred.test.clm.cauchit.flexible[[i]][j,])>=0.5))
pred.Ord.MAP.clm.cauchit.flexible[[i]][j]<-categories[h]
}
}

pred.Ord.MAP.clm.cauchit.symmetric<-list()
for (i in 1:10)
{ pred.Ord.MAP.clm.cauchit.symmetric[[i]]<-vector()
for (j in 1:dim(pred.test.clm.cauchit.symmetric[[i]])[1])
{h<-min(which(cumsum(pred.test.clm.cauchit.symmetric[[i]][j,])>=0.5))
pred.Ord.MAP.clm.cauchit.symmetric[[i]][j]<-categories[h]
}
}

pred.Ord.MAP.clm.cauchit.symmetric2<-list()
for (i in 1:10)
{ pred.Ord.MAP.clm.cauchit.symmetric2[[i]]<-vector()
for (j in 1:dim(pred.test.clm.cauchit.symmetric2[[i]])[1])
{h<-min(which(cumsum(pred.test.clm.cauchit.symmetric2[[i]][j,])>=0.5))
pred.Ord.MAP.clm.cauchit.symmetric2[[i]][j]<-categories[h]
}
}

pred.Ord.MAP.clm.cauchit.equidistant<-list()
for (i in 1:10)
{ pred.Ord.MAP.clm.cauchit.equidistant[[i]]<-vector()
for (j in 1:dim(pred.test.clm.cauchit.equidistant[[i]])[1])
{h<-min(which(cumsum(pred.test.clm.cauchit.equidistant[[i]][j,])>=0.5))
pred.Ord.MAP.clm.cauchit.equidistant[[i]][j]<-categories[h]
}
}

###
# confusion matrices

conf.mat.Ord.MAP.clm.logit.flexible<-list()
for(i in 1:10)
{
  conf.mat.Ord.MAP.clm.logit.flexible[[i]]<-mat.square(table(pred.Ord.MAP.clm.logit.flexible[[i]],test[[i]]$poverty),categories)
}

conf.mat.Ord.MAP.clm.logit.symmetric<-list()
for(i in 1:10)
{
  conf.mat.Ord.MAP.clm.logit.symmetric[[i]]<-mat.square(table(pred.Ord.MAP.clm.logit.symmetric[[i]],test[[i]]$poverty),categories)
}

conf.mat.Ord.MAP.clm.logit.symmetric2<-list()
for(i in 1:10)
{
  conf.mat.Ord.MAP.clm.logit.symmetric2[[i]]<-mat.square(table(pred.Ord.MAP.clm.logit.symmetric2[[i]],test[[i]]$poverty),categories)
}

conf.mat.Ord.MAP.clm.logit.equidistant<-list()
for(i in 1:10)
{
  conf.mat.Ord.MAP.clm.logit.equidistant[[i]]<-mat.square(table(pred.Ord.MAP.clm.logit.equidistant[[i]],test[[i]]$poverty),categories)
}

###

conf.mat.Ord.MAP.clm.probit.flexible<-list()
for(i in 1:10)
{
  conf.mat.Ord.MAP.clm.probit.flexible[[i]]<-mat.square(table(pred.Ord.MAP.clm.probit.flexible[[i]],test[[i]]$poverty),categories)
}

conf.mat.Ord.MAP.clm.probit.symmetric<-list()
for(i in 1:10)
{
  conf.mat.Ord.MAP.clm.probit.symmetric[[i]]<-mat.square(table(pred.Ord.MAP.clm.probit.symmetric[[i]],test[[i]]$poverty),categories)
}

conf.mat.Ord.MAP.clm.probit.symmetric2<-list()
for(i in 1:10)
{
  conf.mat.Ord.MAP.clm.probit.symmetric2[[i]]<-mat.square(table(pred.Ord.MAP.clm.probit.symmetric2[[i]],test[[i]]$poverty),categories)
}

conf.mat.Ord.MAP.clm.probit.equidistant<-list()
for(i in 1:10)
{
  conf.mat.Ord.MAP.clm.probit.equidistant[[i]]<-mat.square(table(pred.Ord.MAP.clm.probit.equidistant[[i]],test[[i]]$poverty),categories)
}

###

conf.mat.Ord.MAP.clm.loglog.flexible<-list()
for(i in 1:10)
{
  conf.mat.Ord.MAP.clm.loglog.flexible[[i]]<-mat.square(table(pred.Ord.MAP.clm.loglog.flexible[[i]],test[[i]]$poverty),categories)
}

conf.mat.Ord.MAP.clm.loglog.symmetric<-list()
for(i in 1:10)
{
  conf.mat.Ord.MAP.clm.loglog.symmetric[[i]]<-mat.square(table(pred.Ord.MAP.clm.loglog.symmetric[[i]],test[[i]]$poverty),categories)
}

conf.mat.Ord.MAP.clm.loglog.symmetric2<-list()
for(i in 1:10)
{
  conf.mat.Ord.MAP.clm.loglog.symmetric2[[i]]<-mat.square(table(pred.Ord.MAP.clm.loglog.symmetric2[[i]],test[[i]]$poverty),categories)
}

conf.mat.Ord.MAP.clm.loglog.equidistant<-list()
for(i in 1:10)
{
  conf.mat.Ord.MAP.clm.loglog.equidistant[[i]]<-mat.square(table(pred.Ord.MAP.clm.loglog.equidistant[[i]],test[[i]]$poverty),categories)
}

###


conf.mat.Ord.MAP.clm.cloglog.flexible<-list()
for(i in 1:10)
{
  conf.mat.Ord.MAP.clm.cloglog.flexible[[i]]<-mat.square(table(pred.Ord.MAP.clm.cloglog.flexible[[i]],test[[i]]$poverty),categories)
}

conf.mat.Ord.MAP.clm.cloglog.symmetric<-list()
for(i in 1:10)
{
  conf.mat.Ord.MAP.clm.cloglog.symmetric[[i]]<-mat.square(table(pred.Ord.MAP.clm.cloglog.symmetric[[i]],test[[i]]$poverty),categories)
}

conf.mat.Ord.MAP.clm.cloglog.symmetric2<-list()
for(i in 1:10)
{
  conf.mat.Ord.MAP.clm.cloglog.symmetric2[[i]]<-mat.square(table(pred.Ord.MAP.clm.cloglog.symmetric2[[i]],test[[i]]$poverty),categories)
}

conf.mat.Ord.MAP.clm.cloglog.equidistant<-list()
for(i in 1:10)
{
  conf.mat.Ord.MAP.clm.cloglog.equidistant[[i]]<-mat.square(table(pred.Ord.MAP.clm.cloglog.equidistant[[i]],test[[i]]$poverty),categories)
}

###

conf.mat.Ord.MAP.clm.cauchit.flexible<-list()
for(i in 1:10)
{
  conf.mat.Ord.MAP.clm.cauchit.flexible[[i]]<-mat.square(table(pred.Ord.MAP.clm.cauchit.flexible[[i]],test[[i]]$poverty),categories)
}

conf.mat.Ord.MAP.clm.cauchit.symmetric<-list()
for(i in 1:10)
{
  conf.mat.Ord.MAP.clm.cauchit.symmetric[[i]]<-mat.square(table(pred.Ord.MAP.clm.cauchit.symmetric[[i]],test[[i]]$poverty),categories)
}

conf.mat.Ord.MAP.clm.cauchit.symmetric2<-list()
for(i in 1:10)
{
  conf.mat.Ord.MAP.clm.cauchit.symmetric2[[i]]<-mat.square(table(pred.Ord.MAP.clm.cauchit.symmetric2[[i]],test[[i]]$poverty),categories)
}

conf.mat.Ord.MAP.clm.cauchit.equidistant<-list()
for(i in 1:10)
{
  conf.mat.Ord.MAP.clm.cauchit.equidistant[[i]]<-mat.square(table(pred.Ord.MAP.clm.cauchit.equidistant[[i]],test[[i]]$poverty),categories)
}

###

# MAE 

MAE.Ord.MAP.clm.logit.flexible<-vector()
for (i in 1:10)
{
  MAE.Ord.MAP.clm.logit.flexible[i]<-MAE(conf.mat.Ord.MAP.clm.logit.flexible[[i]])
}

MAE.Ord.MAP.clm.logit.symmetric<-vector()
for (i in 1:10)
{
  MAE.Ord.MAP.clm.logit.symmetric[i]<-MAE(conf.mat.Ord.MAP.clm.logit.symmetric[[i]])
}

MAE.Ord.MAP.clm.logit.symmetric2<-vector()
for (i in 1:10)
{
  MAE.Ord.MAP.clm.logit.symmetric2[i]<-MAE(conf.mat.Ord.MAP.clm.logit.symmetric2[[i]])
}

MAE.Ord.MAP.clm.logit.equidistant<-vector()
for (i in 1:10)
{
  MAE.Ord.MAP.clm.logit.equidistant[i]<-MAE(conf.mat.Ord.MAP.clm.logit.equidistant[[i]])
}

###

MAE.Ord.MAP.clm.probit.flexible<-vector()
for (i in 1:10)
{
  MAE.Ord.MAP.clm.probit.flexible[i]<-MAE(conf.mat.Ord.MAP.clm.probit.flexible[[i]])
}

MAE.Ord.MAP.clm.probit.symmetric<-vector()
for (i in 1:10)
{
  MAE.Ord.MAP.clm.probit.symmetric[i]<-MAE(conf.mat.Ord.MAP.clm.probit.symmetric[[i]])
}

MAE.Ord.MAP.clm.probit.symmetric2<-vector()
for (i in 1:10)
{
  MAE.Ord.MAP.clm.probit.symmetric2[i]<-MAE(conf.mat.Ord.MAP.clm.probit.symmetric2[[i]])
}

MAE.Ord.MAP.clm.probit.equidistant<-vector()
for (i in 1:10)
{
  MAE.Ord.MAP.clm.probit.equidistant[i]<-MAE(conf.mat.Ord.MAP.clm.probit.equidistant[[i]])
}

###


MAE.Ord.MAP.clm.loglog.flexible<-vector()
for (i in 1:10)
{
  MAE.Ord.MAP.clm.loglog.flexible[i]<-MAE(conf.mat.Ord.MAP.clm.loglog.flexible[[i]])
}

MAE.Ord.MAP.clm.loglog.symmetric<-vector()
for (i in 1:10)
{
  MAE.Ord.MAP.clm.loglog.symmetric[i]<-MAE(conf.mat.Ord.MAP.clm.loglog.symmetric[[i]])
}

MAE.Ord.MAP.clm.loglog.symmetric2<-vector()
for (i in 1:10)
{
  MAE.Ord.MAP.clm.loglog.symmetric2[i]<-MAE(conf.mat.Ord.MAP.clm.loglog.symmetric2[[i]])
}

MAE.Ord.MAP.clm.loglog.equidistant<-vector()
for (i in 1:10)
{
  MAE.Ord.MAP.clm.loglog.equidistant[i]<-MAE(conf.mat.Ord.MAP.clm.loglog.equidistant[[i]])
}

###


MAE.Ord.MAP.clm.cloglog.flexible<-vector()
for (i in 1:10)
{
  MAE.Ord.MAP.clm.cloglog.flexible[i]<-MAE(conf.mat.Ord.MAP.clm.cloglog.flexible[[i]])
}

MAE.Ord.MAP.clm.cloglog.symmetric<-vector()
for (i in 1:10)
{
  MAE.Ord.MAP.clm.cloglog.symmetric[i]<-MAE(conf.mat.Ord.MAP.clm.cloglog.symmetric[[i]])
}

MAE.Ord.MAP.clm.cloglog.symmetric2<-vector()
for (i in 1:10)
{
  MAE.Ord.MAP.clm.cloglog.symmetric2[i]<-MAE(conf.mat.Ord.MAP.clm.cloglog.symmetric2[[i]])
}

MAE.Ord.MAP.clm.cloglog.equidistant<-vector()
for (i in 1:10)
{
  MAE.Ord.MAP.clm.cloglog.equidistant[i]<-MAE(conf.mat.Ord.MAP.clm.cloglog.equidistant[[i]])
}

###


MAE.Ord.MAP.clm.cauchit.flexible<-vector()
for (i in 1:10)
{
  MAE.Ord.MAP.clm.cauchit.flexible[i]<-MAE(conf.mat.Ord.MAP.clm.cauchit.flexible[[i]])
}

MAE.Ord.MAP.clm.cauchit.symmetric<-vector()
for (i in 1:10)
{
  MAE.Ord.MAP.clm.cauchit.symmetric[i]<-MAE(conf.mat.Ord.MAP.clm.cauchit.symmetric[[i]])
}

MAE.Ord.MAP.clm.cauchit.symmetric2<-vector()
for (i in 1:10)
{
  MAE.Ord.MAP.clm.cauchit.symmetric2[i]<-MAE(conf.mat.Ord.MAP.clm.cauchit.symmetric2[[i]])
}

MAE.Ord.MAP.clm.cauchit.equidistant<-vector()
for (i in 1:10)
{
  MAE.Ord.MAP.clm.cauchit.equidistant[i]<-MAE(conf.mat.Ord.MAP.clm.cauchit.equidistant[[i]])
}


###
################################################################################
###. Comparison MAP vs Ord.MAP
###
###.  Tests MAE

p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.logit.flexible<-wilcox.test(MAE.MAP.clm.logit.flexible,MAE.Ord.MAP.clm.logit.flexible,paired = TRUE,alternative="greater")$p.value
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.logit.symmetric<-wilcox.test(MAE.MAP.clm.logit.symmetric,MAE.Ord.MAP.clm.logit.symmetric,paired = TRUE,alternative="greater")$p.value
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.logit.symmetric2<-wilcox.test(MAE.MAP.clm.logit.symmetric2,MAE.Ord.MAP.clm.logit.symmetric2,paired = TRUE,alternative="greater")$p.value
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.logit.equidistant<-wilcox.test(MAE.MAP.clm.logit.equidistant,MAE.Ord.MAP.clm.logit.equidistant,paired = TRUE,alternative="greater")$p.value

p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.logit.flexible.less<-wilcox.test(MAE.MAP.clm.logit.flexible,MAE.Ord.MAP.clm.logit.flexible,paired = TRUE,alternative="less")$p.value
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.logit.symmetric.less<-wilcox.test(MAE.MAP.clm.logit.symmetric,MAE.Ord.MAP.clm.logit.symmetric,paired = TRUE,alternative="less")$p.value
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.logit.symmetric2.less<-wilcox.test(MAE.MAP.clm.logit.symmetric2,MAE.Ord.MAP.clm.logit.symmetric2,paired = TRUE,alternative="less")$p.value
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.logit.equidistant.less<-wilcox.test(MAE.MAP.clm.logit.equidistant,MAE.Ord.MAP.clm.logit.equidistant,paired = TRUE,alternative="less")$p.value

###

p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.probit.flexible<-wilcox.test(MAE.MAP.clm.probit.flexible,MAE.Ord.MAP.clm.probit.flexible,paired = TRUE,alternative="greater")$p.value
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.probit.symmetric<-wilcox.test(MAE.MAP.clm.probit.symmetric,MAE.Ord.MAP.clm.probit.symmetric,paired = TRUE,alternative="greater")$p.value
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.probit.symmetric2<-wilcox.test(MAE.MAP.clm.probit.symmetric2,MAE.Ord.MAP.clm.probit.symmetric2,paired = TRUE,alternative="greater")$p.value
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.probit.equidistant<-wilcox.test(MAE.MAP.clm.probit.equidistant,MAE.Ord.MAP.clm.probit.equidistant,paired = TRUE,alternative="greater")$p.value

p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.probit.flexible.less<-wilcox.test(MAE.MAP.clm.probit.flexible,MAE.Ord.MAP.clm.probit.flexible,paired = TRUE,alternative="less")$p.value
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.probit.symmetric.less<-wilcox.test(MAE.MAP.clm.probit.symmetric,MAE.Ord.MAP.clm.probit.symmetric,paired = TRUE,alternative="less")$p.value
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.probit.symmetric2.less<-wilcox.test(MAE.MAP.clm.probit.symmetric2,MAE.Ord.MAP.clm.probit.symmetric2,paired = TRUE,alternative="less")$p.value
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.probit.equidistant.less<-wilcox.test(MAE.MAP.clm.probit.equidistant,MAE.Ord.MAP.clm.probit.equidistant,paired = TRUE,alternative="less")$p.value


###

p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.loglog.flexible<-wilcox.test(MAE.MAP.clm.loglog.flexible,MAE.Ord.MAP.clm.loglog.flexible,paired = TRUE,alternative="greater")$p.value
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.loglog.symmetric<-wilcox.test(MAE.MAP.clm.loglog.symmetric,MAE.Ord.MAP.clm.loglog.symmetric,paired = TRUE,alternative="greater")$p.value
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.loglog.symmetric2<-wilcox.test(MAE.MAP.clm.loglog.symmetric2,MAE.Ord.MAP.clm.loglog.symmetric2,paired = TRUE,alternative="greater")$p.value
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.loglog.equidistant<-wilcox.test(MAE.MAP.clm.loglog.equidistant,MAE.Ord.MAP.clm.loglog.equidistant,paired = TRUE,alternative="greater")$p.value

p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.loglog.flexible.less<-wilcox.test(MAE.MAP.clm.loglog.flexible,MAE.Ord.MAP.clm.loglog.flexible,paired = TRUE,alternative="less")$p.value
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.loglog.symmetric.less<-wilcox.test(MAE.MAP.clm.loglog.symmetric,MAE.Ord.MAP.clm.loglog.symmetric,paired = TRUE,alternative="less")$p.value
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.loglog.symmetric2.less<-wilcox.test(MAE.MAP.clm.loglog.symmetric2,MAE.Ord.MAP.clm.loglog.symmetric2,paired = TRUE,alternative="less")$p.value
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.loglog.equidistant.less<-wilcox.test(MAE.MAP.clm.loglog.equidistant,MAE.Ord.MAP.clm.loglog.equidistant,paired = TRUE,alternative="less")$p.value

###

p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.cloglog.flexible<-wilcox.test(MAE.MAP.clm.cloglog.flexible,MAE.Ord.MAP.clm.cloglog.flexible,paired = TRUE,alternative="greater")$p.value
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.cloglog.symmetric<-wilcox.test(MAE.MAP.clm.cloglog.symmetric,MAE.Ord.MAP.clm.cloglog.symmetric,paired = TRUE,alternative="greater")$p.value
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.cloglog.symmetric2<-wilcox.test(MAE.MAP.clm.cloglog.symmetric2,MAE.Ord.MAP.clm.cloglog.symmetric2,paired = TRUE,alternative="greater")$p.value
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.cloglog.equidistant<-wilcox.test(MAE.MAP.clm.cloglog.equidistant,MAE.Ord.MAP.clm.cloglog.equidistant,paired = TRUE,alternative="greater")$p.value

p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.cloglog.flexible.less<-wilcox.test(MAE.MAP.clm.cloglog.flexible,MAE.Ord.MAP.clm.cloglog.flexible,paired = TRUE,alternative="less")$p.value
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.cloglog.symmetric.less<-wilcox.test(MAE.MAP.clm.cloglog.symmetric,MAE.Ord.MAP.clm.cloglog.symmetric,paired = TRUE,alternative="less")$p.value
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.cloglog.symmetric2.less<-wilcox.test(MAE.MAP.clm.cloglog.symmetric2,MAE.Ord.MAP.clm.cloglog.symmetric2,paired = TRUE,alternative="less")$p.value
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.cloglog.equidistant.less<-wilcox.test(MAE.MAP.clm.cloglog.equidistant,MAE.Ord.MAP.clm.cloglog.equidistant,paired = TRUE,alternative="les")$p.value

###

p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.cauchit.flexible<-wilcox.test(MAE.MAP.clm.cauchit.flexible,MAE.Ord.MAP.clm.cauchit.flexible,paired = TRUE,alternative="greater")$p.value
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.cauchit.symmetric<-wilcox.test(MAE.MAP.clm.cauchit.symmetric,MAE.Ord.MAP.clm.cauchit.symmetric,paired = TRUE,alternative="greater")$p.value
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.cauchit.symmetric2<-wilcox.test(MAE.MAP.clm.cauchit.symmetric2,MAE.Ord.MAP.clm.cauchit.symmetric2,paired = TRUE,alternative="greater")$p.value
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.cauchit.equidistant<-wilcox.test(MAE.MAP.clm.cauchit.equidistant,MAE.Ord.MAP.clm.cauchit.equidistant,paired = TRUE,alternative="greater")$p.value

p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.cauchit.flexible.less<-wilcox.test(MAE.MAP.clm.cauchit.flexible,MAE.Ord.MAP.clm.cauchit.flexible,paired = TRUE,alternative="less")$p.value
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.cauchit.symmetric.less<-wilcox.test(MAE.MAP.clm.cauchit.symmetric,MAE.Ord.MAP.clm.cauchit.symmetric,paired = TRUE,alternative="less")$p.value
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.cauchit.symmetric2.less<-wilcox.test(MAE.MAP.clm.cauchit.symmetric2,MAE.Ord.MAP.clm.cauchit.symmetric2,paired = TRUE,alternative="less")$p.value
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.cauchit.equidistant.less<-wilcox.test(MAE.MAP.clm.cauchit.equidistant,MAE.Ord.MAP.clm.cauchit.equidistant,paired = TRUE,alternative="less")$p.value

###

p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.logit.flexible
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.logit.symmetric
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.logit.symmetric2
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.logit.equidistant

p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.logit.flexible.less
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.logit.symmetric.less
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.logit.symmetric2.less
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.logit.equidistant.less

# 

p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.probit.flexible
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.probit.symmetric
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.probit.symmetric2
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.probit.equidistant

p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.probit.flexible.less
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.probit.symmetric.less
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.probit.symmetric2.less
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.probit.equidistant.less

#

p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.loglog.flexible
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.loglog.symmetric
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.loglog.symmetric2
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.loglog.equidistant

p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.loglog.flexible.less
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.loglog.symmetric.less
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.loglog.symmetric2.less
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.loglog.equidistant.less

#

p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.cloglog.flexible
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.cloglog.symmetric
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.cloglog.symmetric2
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.cloglog.equidistant

p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.cloglog.flexible.less
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.cloglog.symmetric.less
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.cloglog.symmetric2.less
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.cloglog.equidistant.less

#

p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.cauchit.flexible
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.cauchit.symmetric
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.cauchit.symmetric2
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.cauchit.equidistant

p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.cauchit.flexible.less
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.cauchit.symmetric.less
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.cauchit.symmetric2.less
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.cauchit.equidistant.less


###############################
###############################
###   Boxplots   MAE

MAE.clm.flexible.results<-as.data.frame(cbind(MAE.MAP.clm.logit.flexible-MAE.Ord.MAP.clm.logit.flexible,
                           MAE.MAP.clm.probit.flexible-MAE.Ord.MAP.clm.probit.flexible,
                           MAE.MAP.clm.loglog.flexible-MAE.Ord.MAP.clm.loglog.flexible,
                           MAE.MAP.clm.cloglog.flexible-MAE.Ord.MAP.clm.cloglog.flexible,
                           MAE.MAP.clm.cauchit.flexible-MAE.Ord.MAP.clm.cauchit.flexible))

boxplot(MAE.clm.flexible.results,
        main=paste("(a) WVS dataset. Threshold flexible"), 
        xlab=" ", ylab="MAE increment (MAP - Ord.MAP)",names=c("logistic", "probit","loglog","cloglog","cauchit"),
        col="white")
abline(h = 0, col = 2, lwd = 2, lty=2)
# Dots
stripchart(MAE.clm.flexible.results,
           #method = "jitter",
           pch = 19,
           col = 1:5,
           vertical = TRUE,
           add = TRUE)

###

MAE.clm.symmetric.results<-as.data.frame(cbind(MAE.MAP.clm.logit.symmetric-MAE.Ord.MAP.clm.logit.symmetric,
                                              MAE.MAP.clm.probit.symmetric-MAE.Ord.MAP.clm.probit.symmetric,
                                              MAE.MAP.clm.loglog.symmetric-MAE.Ord.MAP.clm.loglog.symmetric,
                                              MAE.MAP.clm.cloglog.symmetric-MAE.Ord.MAP.clm.cloglog.symmetric,
                                              MAE.MAP.clm.cauchit.symmetric-MAE.Ord.MAP.clm.cauchit.symmetric))

boxplot(MAE.clm.symmetric.results,
        main=paste("(a) WVS dataset. Threshold symmetric"), 
        xlab=" ", ylab="MAE increment (MAP - Ord.MAP)",names=c("logistic", "probit","loglog","cloglog","cauchit"),
        col="white")
abline(h = 0, col = 2, lwd = 2, lty=2)
# Dots
stripchart(MAE.clm.symmetric.results,
           #method = "jitter",
           pch = 19,
           col = 1:5,
           vertical = TRUE,
           add = TRUE)

###

MAE.clm.symmetric2.results<-as.data.frame(cbind(MAE.MAP.clm.logit.symmetric2-MAE.Ord.MAP.clm.logit.symmetric2,
                                               MAE.MAP.clm.probit.symmetric2-MAE.Ord.MAP.clm.probit.symmetric2,
                                               MAE.MAP.clm.loglog.symmetric2-MAE.Ord.MAP.clm.loglog.symmetric2,
                                               MAE.MAP.clm.cloglog.symmetric2-MAE.Ord.MAP.clm.cloglog.symmetric2,
                                               MAE.MAP.clm.cauchit.symmetric2-MAE.Ord.MAP.clm.cauchit.symmetric2))

boxplot(MAE.clm.symmetric2.results,
        main=paste("(a) WVS dataset. Threshold symmetric2"), 
        xlab=" ", ylab="MAE increment (MAP - Ord.MAP)",names=c("logistic", "probit","loglog","cloglog","cauchit"),
        col="white")
abline(h = 0, col = 2, lwd = 2, lty=2)
# Dots
stripchart(MAE.clm.symmetric2.results,
           #method = "jitter",
           pch = 19,
           col = 1:5,
           vertical = TRUE,
           add = TRUE)

###

MAE.clm.equidistant.results<-as.data.frame(cbind(MAE.MAP.clm.logit.equidistant-MAE.Ord.MAP.clm.logit.equidistant,
                                                MAE.MAP.clm.probit.equidistant-MAE.Ord.MAP.clm.probit.equidistant,
                                                MAE.MAP.clm.loglog.equidistant-MAE.Ord.MAP.clm.loglog.equidistant,
                                                MAE.MAP.clm.cloglog.equidistant-MAE.Ord.MAP.clm.cloglog.equidistant,
                                                MAE.MAP.clm.cauchit.equidistant-MAE.Ord.MAP.clm.cauchit.equidistant))

boxplot(MAE.clm.equidistant.results,
        main=paste("(a) WVS dataset. Threshold equidistant"), 
        xlab=" ", ylab="MAE increment (MAP - Ord.MAP)",names=c("logistic", "probit","loglog","cloglog","cauchit"),
        col="white")
abline(h = 0, col = 2, lwd = 2, lty=2)
# Dots
stripchart(MAE.clm.equidistant.results,
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
# clm, from "ordinal" package, rfits cumulative link models (CLMs) such as the propotional odds model. 
# The model allows for various link functions and structured thresholds that restricts 
# the thresholds or cut-points to be e.g.,
# equidistant or symmetrically arranged around the central threshold(s). 
#
# clm(formula, scale, nominal, data, weights, start, subset, doFit = TRUE,
#     na.action, contrasts, model = TRUE, control=list(),
#     link = c("logit", "probit", "cloglog", "loglog", "cauchit",
#              "Aranda-Ordaz", "log-gamma"),
#     threshold = c("flexible", "symmetric", "symmetric2", "equidistant"), ...)

model.clm.logit.flexible<-list()
for (i in 1:10)
{model.clm.logit.flexible[[i]] <- clm(rating ~ temp + contact, data = training[[i]],
                                      link="logit", threshold="flexible")}


model.clm.logit.symmetric<-list()
for (i in 1:10)
{model.clm.logit.symmetric[[i]] <- clm(rating ~ temp + contact, data = training[[i]],
                                       link="logit", threshold="symmetric")}


model.clm.logit.symmetric2<-list()
for (i in 1:10)
{model.clm.logit.symmetric2[[i]] <- clm(rating ~ temp + contact, data = training[[i]],
                                        link="logit", threshold="symmetric2")}

model.clm.logit.equidistant<-list()
for (i in 1:10)
{model.clm.logit.equidistant[[i]] <- clm(rating ~ temp + contact, data = training[[i]],
                                         link="logit", threshold="equidistant")}


###

model.clm.probit.flexible<-list()
for (i in 1:10)
{model.clm.probit.flexible[[i]] <- clm(rating ~ temp + contact, data = training[[i]],
                                       link="probit", threshold="flexible")}


model.clm.probit.symmetric<-list()
for (i in 1:10)
{model.clm.probit.symmetric[[i]] <- clm(rating ~ temp + contact, data = training[[i]],
                                        link="probit", threshold="symmetric")}


model.clm.probit.symmetric2<-list()
for (i in 1:10)
{model.clm.probit.symmetric2[[i]] <- clm(rating ~ temp + contact, data = training[[i]],
                                         link="probit", threshold="symmetric2")}

model.clm.probit.equidistant<-list()
for (i in 1:10)
{model.clm.probit.equidistant[[i]] <- clm(rating ~ temp + contact, data = training[[i]],
                                          link="probit", threshold="equidistant")}


###

model.clm.loglog.flexible<-list()
for (i in 1:10)
{model.clm.loglog.flexible[[i]] <- clm(rating ~ temp + contact, data = training[[i]],
                                       link="loglog", threshold="flexible")}


model.clm.loglog.symmetric<-list()
for (i in 1:10)
{model.clm.loglog.symmetric[[i]] <- clm(rating ~ temp + contact, data = training[[i]],
                                        link="loglog", threshold="symmetric")}


model.clm.loglog.symmetric2<-list()
for (i in 1:10)
{model.clm.loglog.symmetric2[[i]] <- clm(rating ~ temp + contact, data = training[[i]],
                                         link="loglog", threshold="symmetric2")}

model.clm.loglog.equidistant<-list()
for (i in 1:10)
{model.clm.loglog.equidistant[[i]] <- clm(rating ~ temp + contact, data = training[[i]],
                                          link="loglog", threshold="equidistant")}


###

model.clm.cloglog.flexible<-list()
for (i in 1:10)
{model.clm.cloglog.flexible[[i]] <- clm(rating ~ temp + contact, data = training[[i]],
                                        link="cloglog", threshold="flexible")}


model.clm.cloglog.symmetric<-list()
for (i in 1:10)
{model.clm.cloglog.symmetric[[i]] <- clm(rating ~ temp + contact, data = training[[i]],
                                         link="cloglog", threshold="symmetric")}


model.clm.cloglog.symmetric2<-list()
for (i in 1:10)
{model.clm.cloglog.symmetric2[[i]] <- clm(rating ~ temp + contact, data = training[[i]],
                                          link="cloglog", threshold="symmetric2")}

model.clm.cloglog.equidistant<-list()
for (i in 1:10)
{model.clm.cloglog.equidistant[[i]] <- clm(rating ~ temp + contact, data = training[[i]],
                                           link="cloglog", threshold="equidistant")}


###

model.clm.cauchit.flexible<-list()
for (i in 1:10)
{model.clm.cauchit.flexible[[i]] <- clm(rating ~ temp + contact, data = training[[i]],
                                        link="cauchit", threshold="flexible")}


model.clm.cauchit.symmetric<-list()
for (i in 1:10)
{model.clm.cauchit.symmetric[[i]] <- clm(rating ~ temp + contact, data = training[[i]],
                                         link="cauchit", threshold="symmetric")}


model.clm.cauchit.symmetric2<-list()
for (i in 1:10)
{model.clm.cauchit.symmetric2[[i]] <- clm(rating ~ temp + contact, data = training[[i]],
                                          link="cauchit", threshold="symmetric2")}

model.clm.cauchit.equidistant<-list()
for (i in 1:10)
{model.clm.cauchit.equidistant[[i]] <- clm(rating ~ temp + contact, data = training[[i]],
                                           link="cauchit", threshold="equidistant")}



################################################################################
#
# Predictions of any of the test set
# 
#

pred.test.clm.logit.flexible<-list()
for (i in 1:10)
{pred.test.clm.logit.flexible[[i]] <- predict(model.clm.logit.flexible[[i]],test[[i]][,-2],type = "p")$fit
}

pred.test.clm.logit.symmetric<-list()
for (i in 1:10)
{pred.test.clm.logit.symmetric[[i]] <- predict(model.clm.logit.symmetric[[i]],test[[i]][,-2],type = "p")$fit
}

pred.test.clm.logit.symmetric2<-list()
for (i in 1:10)
{pred.test.clm.logit.symmetric2[[i]] <- predict(model.clm.logit.symmetric2[[i]],test[[i]][,-2],type = "p")$fit
}

pred.test.clm.logit.equidistant<-list()
for (i in 1:10)
{pred.test.clm.logit.equidistant[[i]] <- predict(model.clm.logit.equidistant[[i]],test[[i]][,-2],type = "p")$fit
}

###

pred.test.clm.probit.flexible<-list()
for (i in 1:10)
{pred.test.clm.probit.flexible[[i]] <- predict(model.clm.probit.flexible[[i]],test[[i]][,-2],type = "p")$fit
}

pred.test.clm.probit.symmetric<-list()
for (i in 1:10)
{pred.test.clm.probit.symmetric[[i]] <- predict(model.clm.probit.symmetric[[i]],test[[i]][,-2],type = "p")$fit
}

pred.test.clm.probit.symmetric2<-list()
for (i in 1:10)
{pred.test.clm.probit.symmetric2[[i]] <- predict(model.clm.probit.symmetric2[[i]],test[[i]][,-2],type = "p")$fit
}

pred.test.clm.probit.equidistant<-list()
for (i in 1:10)
{pred.test.clm.probit.equidistant[[i]] <- predict(model.clm.probit.equidistant[[i]],test[[i]][,-2],type = "p")$fit
}

###

pred.test.clm.loglog.flexible<-list()
for (i in 1:10)
{pred.test.clm.loglog.flexible[[i]] <- predict(model.clm.loglog.flexible[[i]],test[[i]][,-2],type = "p")$fit
}

pred.test.clm.loglog.symmetric<-list()
for (i in 1:10)
{pred.test.clm.loglog.symmetric[[i]] <- predict(model.clm.loglog.symmetric[[i]],test[[i]][,-2],type = "p")$fit
}

pred.test.clm.loglog.symmetric2<-list()
for (i in 1:10)
{pred.test.clm.loglog.symmetric2[[i]] <- predict(model.clm.loglog.symmetric2[[i]],test[[i]][,-2],type = "p")$fit
}

pred.test.clm.loglog.equidistant<-list()
for (i in 1:10)
{pred.test.clm.loglog.equidistant[[i]] <- predict(model.clm.loglog.equidistant[[i]],test[[i]][,-2],type = "p")$fit
}

###


pred.test.clm.cloglog.flexible<-list()
for (i in 1:10)
{pred.test.clm.cloglog.flexible[[i]] <- predict(model.clm.cloglog.flexible[[i]],test[[i]][,-2],type = "p")$fit
}

pred.test.clm.cloglog.symmetric<-list()
for (i in 1:10)
{pred.test.clm.cloglog.symmetric[[i]] <- predict(model.clm.cloglog.symmetric[[i]],test[[i]][,-2],type = "p")$fit
}

pred.test.clm.cloglog.symmetric2<-list()
for (i in 1:10)
{pred.test.clm.cloglog.symmetric2[[i]] <- predict(model.clm.cloglog.symmetric2[[i]],test[[i]][,-2],type = "p")$fit
}

pred.test.clm.cloglog.equidistant<-list()
for (i in 1:10)
{pred.test.clm.cloglog.equidistant[[i]] <- predict(model.clm.cloglog.equidistant[[i]],test[[i]][,-2],type = "p")$fit
}

###


pred.test.clm.cauchit.flexible<-list()
for (i in 1:10)
{pred.test.clm.cauchit.flexible[[i]] <- predict(model.clm.cauchit.flexible[[i]],test[[i]][,-2],type = "p")$fit
}

pred.test.clm.cauchit.symmetric<-list()
for (i in 1:10)
{pred.test.clm.cauchit.symmetric[[i]] <- predict(model.clm.cauchit.symmetric[[i]],test[[i]][,-2],type = "p")$fit
}

pred.test.clm.cauchit.symmetric2<-list()
for (i in 1:10)
{pred.test.clm.cauchit.symmetric2[[i]] <- predict(model.clm.cauchit.symmetric2[[i]],test[[i]][,-2],type = "p")$fit
}

pred.test.clm.cauchit.equidistant<-list()
for (i in 1:10)
{pred.test.clm.cauchit.equidistant[[i]] <- predict(model.clm.cauchit.equidistant[[i]],test[[i]][,-2],type = "p")$fit
}

###

################################################################################
### predictions with standard MAP criterion:
###

pred.MAP.clm.logit.flexible<-list()
for (i in 1:10)
{pred.MAP.clm.logit.flexible[[i]]<-vector()
for (j in 1:dim(pred.test.clm.logit.flexible[[i]])[1])
{pred.MAP.clm.logit.flexible[[i]][j]<-categories[which.max(pred.test.clm.logit.flexible[[i]][j,])]
}
}

pred.MAP.clm.logit.symmetric<-list()
for (i in 1:10)
{pred.MAP.clm.logit.symmetric[[i]]<-vector()
for (j in 1:dim(pred.test.clm.logit.symmetric[[i]])[1])
{pred.MAP.clm.logit.symmetric[[i]][j]<-categories[which.max(pred.test.clm.logit.symmetric[[i]][j,])]
}
}

pred.MAP.clm.logit.symmetric2<-list()
for (i in 1:10)
{pred.MAP.clm.logit.symmetric2[[i]]<-vector()
for (j in 1:dim(pred.test.clm.logit.symmetric2[[i]])[1])
{pred.MAP.clm.logit.symmetric2[[i]][j]<-categories[which.max(pred.test.clm.logit.symmetric2[[i]][j,])]
}
}

pred.MAP.clm.logit.equidistant<-list()
for (i in 1:10)
{pred.MAP.clm.logit.equidistant[[i]]<-vector()
for (j in 1:dim(pred.test.clm.logit.equidistant[[i]])[1])
{pred.MAP.clm.logit.equidistant[[i]][j]<-categories[which.max(pred.test.clm.logit.equidistant[[i]][j,])]
}
}

###


pred.MAP.clm.probit.flexible<-list()
for (i in 1:10)
{pred.MAP.clm.probit.flexible[[i]]<-vector()
for (j in 1:dim(pred.test.clm.probit.flexible[[i]])[1])
{pred.MAP.clm.probit.flexible[[i]][j]<-categories[which.max(pred.test.clm.probit.flexible[[i]][j,])]
}
}

pred.MAP.clm.probit.symmetric<-list()
for (i in 1:10)
{pred.MAP.clm.probit.symmetric[[i]]<-vector()
for (j in 1:dim(pred.test.clm.probit.symmetric[[i]])[1])
{pred.MAP.clm.probit.symmetric[[i]][j]<-categories[which.max(pred.test.clm.probit.symmetric[[i]][j,])]
}
}

pred.MAP.clm.probit.symmetric2<-list()
for (i in 1:10)
{pred.MAP.clm.probit.symmetric2[[i]]<-vector()
for (j in 1:dim(pred.test.clm.probit.symmetric2[[i]])[1])
{pred.MAP.clm.probit.symmetric2[[i]][j]<-categories[which.max(pred.test.clm.probit.symmetric2[[i]][j,])]
}
}

pred.MAP.clm.probit.equidistant<-list()
for (i in 1:10)
{pred.MAP.clm.probit.equidistant[[i]]<-vector()
for (j in 1:dim(pred.test.clm.probit.equidistant[[i]])[1])
{pred.MAP.clm.probit.equidistant[[i]][j]<-categories[which.max(pred.test.clm.probit.equidistant[[i]][j,])]
}
}

###


pred.MAP.clm.loglog.flexible<-list()
for (i in 1:10)
{pred.MAP.clm.loglog.flexible[[i]]<-vector()
for (j in 1:dim(pred.test.clm.loglog.flexible[[i]])[1])
{pred.MAP.clm.loglog.flexible[[i]][j]<-categories[which.max(pred.test.clm.loglog.flexible[[i]][j,])]
}
}

pred.MAP.clm.loglog.symmetric<-list()
for (i in 1:10)
{pred.MAP.clm.loglog.symmetric[[i]]<-vector()
for (j in 1:dim(pred.test.clm.loglog.symmetric[[i]])[1])
{pred.MAP.clm.loglog.symmetric[[i]][j]<-categories[which.max(pred.test.clm.loglog.symmetric[[i]][j,])]
}
}

pred.MAP.clm.loglog.symmetric2<-list()
for (i in 1:10)
{pred.MAP.clm.loglog.symmetric2[[i]]<-vector()
for (j in 1:dim(pred.test.clm.loglog.symmetric2[[i]])[1])
{pred.MAP.clm.loglog.symmetric2[[i]][j]<-categories[which.max(pred.test.clm.loglog.symmetric2[[i]][j,])]
}
}

pred.MAP.clm.loglog.equidistant<-list()
for (i in 1:10)
{pred.MAP.clm.loglog.equidistant[[i]]<-vector()
for (j in 1:dim(pred.test.clm.loglog.equidistant[[i]])[1])
{pred.MAP.clm.loglog.equidistant[[i]][j]<-categories[which.max(pred.test.clm.loglog.equidistant[[i]][j,])]
}
}

###


pred.MAP.clm.cloglog.flexible<-list()
for (i in 1:10)
{pred.MAP.clm.cloglog.flexible[[i]]<-vector()
for (j in 1:dim(pred.test.clm.cloglog.flexible[[i]])[1])
{pred.MAP.clm.cloglog.flexible[[i]][j]<-categories[which.max(pred.test.clm.cloglog.flexible[[i]][j,])]
}
}

pred.MAP.clm.cloglog.symmetric<-list()
for (i in 1:10)
{pred.MAP.clm.cloglog.symmetric[[i]]<-vector()
for (j in 1:dim(pred.test.clm.cloglog.symmetric[[i]])[1])
{pred.MAP.clm.cloglog.symmetric[[i]][j]<-categories[which.max(pred.test.clm.cloglog.symmetric[[i]][j,])]
}
}

pred.MAP.clm.cloglog.symmetric2<-list()
for (i in 1:10)
{pred.MAP.clm.cloglog.symmetric2[[i]]<-vector()
for (j in 1:dim(pred.test.clm.cloglog.symmetric2[[i]])[1])
{pred.MAP.clm.cloglog.symmetric2[[i]][j]<-categories[which.max(pred.test.clm.cloglog.symmetric2[[i]][j,])]
}
}

pred.MAP.clm.cloglog.equidistant<-list()
for (i in 1:10)
{pred.MAP.clm.cloglog.equidistant[[i]]<-vector()
for (j in 1:dim(pred.test.clm.cloglog.equidistant[[i]])[1])
{pred.MAP.clm.cloglog.equidistant[[i]][j]<-categories[which.max(pred.test.clm.cloglog.equidistant[[i]][j,])]
}
}

###


pred.MAP.clm.cauchit.flexible<-list()
for (i in 1:10)
{pred.MAP.clm.cauchit.flexible[[i]]<-vector()
for (j in 1:dim(pred.test.clm.cauchit.flexible[[i]])[1])
{pred.MAP.clm.cauchit.flexible[[i]][j]<-categories[which.max(pred.test.clm.cauchit.flexible[[i]][j,])]
}
}

pred.MAP.clm.cauchit.symmetric<-list()
for (i in 1:10)
{pred.MAP.clm.cauchit.symmetric[[i]]<-vector()
for (j in 1:dim(pred.test.clm.cauchit.symmetric[[i]])[1])
{pred.MAP.clm.cauchit.symmetric[[i]][j]<-categories[which.max(pred.test.clm.cauchit.symmetric[[i]][j,])]
}
}

pred.MAP.clm.cauchit.symmetric2<-list()
for (i in 1:10)
{pred.MAP.clm.cauchit.symmetric2[[i]]<-vector()
for (j in 1:dim(pred.test.clm.cauchit.symmetric2[[i]])[1])
{pred.MAP.clm.cauchit.symmetric2[[i]][j]<-categories[which.max(pred.test.clm.cauchit.symmetric2[[i]][j,])]
}
}

pred.MAP.clm.cauchit.equidistant<-list()
for (i in 1:10)
{pred.MAP.clm.cauchit.equidistant[[i]]<-vector()
for (j in 1:dim(pred.test.clm.cauchit.equidistant[[i]])[1])
{pred.MAP.clm.cauchit.equidistant[[i]][j]<-categories[which.max(pred.test.clm.cauchit.equidistant[[i]][j,])]
}
}

# confusion matrices


conf.mat.MAP.clm.logit.flexible<-list()
for(i in 1:10)
{
  conf.mat.MAP.clm.logit.flexible[[i]]<-mat.square(table(pred.MAP.clm.logit.flexible[[i]],test[[i]]$rating),categories)
}


conf.mat.MAP.clm.logit.symmetric<-list()
for(i in 1:10)
{
  conf.mat.MAP.clm.logit.symmetric[[i]]<-mat.square(table(pred.MAP.clm.logit.symmetric[[i]],test[[i]]$rating),categories)
}


conf.mat.MAP.clm.logit.symmetric2<-list()
for(i in 1:10)
{
  conf.mat.MAP.clm.logit.symmetric2[[i]]<-mat.square(table(pred.MAP.clm.logit.symmetric2[[i]],test[[i]]$rating),categories)
}


conf.mat.MAP.clm.logit.equidistant<-list()
for(i in 1:10)
{
  conf.mat.MAP.clm.logit.equidistant[[i]]<-mat.square(table(pred.MAP.clm.logit.equidistant[[i]],test[[i]]$rating),categories)
}

###

conf.mat.MAP.clm.probit.flexible<-list()
for(i in 1:10)
{
  conf.mat.MAP.clm.probit.flexible[[i]]<-mat.square(table(pred.MAP.clm.probit.flexible[[i]],test[[i]]$rating),categories)
}


conf.mat.MAP.clm.probit.symmetric<-list()
for(i in 1:10)
{
  conf.mat.MAP.clm.probit.symmetric[[i]]<-mat.square(table(pred.MAP.clm.probit.symmetric[[i]],test[[i]]$rating),categories)
}


conf.mat.MAP.clm.probit.symmetric2<-list()
for(i in 1:10)
{
  conf.mat.MAP.clm.probit.symmetric2[[i]]<-mat.square(table(pred.MAP.clm.probit.symmetric2[[i]],test[[i]]$rating),categories)
}


conf.mat.MAP.clm.probit.equidistant<-list()
for(i in 1:10)
{
  conf.mat.MAP.clm.probit.equidistant[[i]]<-mat.square(table(pred.MAP.clm.probit.equidistant[[i]],test[[i]]$rating),categories)
}

###

conf.mat.MAP.clm.loglog.flexible<-list()
for(i in 1:10)
{
  conf.mat.MAP.clm.loglog.flexible[[i]]<-mat.square(table(pred.MAP.clm.loglog.flexible[[i]],test[[i]]$rating),categories)
}


conf.mat.MAP.clm.loglog.symmetric<-list()
for(i in 1:10)
{
  conf.mat.MAP.clm.loglog.symmetric[[i]]<-mat.square(table(pred.MAP.clm.loglog.symmetric[[i]],test[[i]]$rating),categories)
}


conf.mat.MAP.clm.loglog.symmetric2<-list()
for(i in 1:10)
{
  conf.mat.MAP.clm.loglog.symmetric2[[i]]<-mat.square(table(pred.MAP.clm.loglog.symmetric2[[i]],test[[i]]$rating),categories)
}


conf.mat.MAP.clm.loglog.equidistant<-list()
for(i in 1:10)
{
  conf.mat.MAP.clm.loglog.equidistant[[i]]<-mat.square(table(pred.MAP.clm.loglog.equidistant[[i]],test[[i]]$rating),categories)
}

###

conf.mat.MAP.clm.cloglog.flexible<-list()
for(i in 1:10)
{
  conf.mat.MAP.clm.cloglog.flexible[[i]]<-mat.square(table(pred.MAP.clm.cloglog.flexible[[i]],test[[i]]$rating),categories)
}


conf.mat.MAP.clm.cloglog.symmetric<-list()
for(i in 1:10)
{
  conf.mat.MAP.clm.cloglog.symmetric[[i]]<-mat.square(table(pred.MAP.clm.cloglog.symmetric[[i]],test[[i]]$rating),categories)
}


conf.mat.MAP.clm.cloglog.symmetric2<-list()
for(i in 1:10)
{
  conf.mat.MAP.clm.cloglog.symmetric2[[i]]<-mat.square(table(pred.MAP.clm.cloglog.symmetric2[[i]],test[[i]]$rating),categories)
}


conf.mat.MAP.clm.cloglog.equidistant<-list()
for(i in 1:10)
{
  conf.mat.MAP.clm.cloglog.equidistant[[i]]<-mat.square(table(pred.MAP.clm.cloglog.equidistant[[i]],test[[i]]$rating),categories)
}

###

conf.mat.MAP.clm.cauchit.flexible<-list()
for(i in 1:10)
{
  conf.mat.MAP.clm.cauchit.flexible[[i]]<-mat.square(table(pred.MAP.clm.cauchit.flexible[[i]],test[[i]]$rating),categories)
}


conf.mat.MAP.clm.cauchit.symmetric<-list()
for(i in 1:10)
{
  conf.mat.MAP.clm.cauchit.symmetric[[i]]<-mat.square(table(pred.MAP.clm.cauchit.symmetric[[i]],test[[i]]$rating),categories)
}


conf.mat.MAP.clm.cauchit.symmetric2<-list()
for(i in 1:10)
{
  conf.mat.MAP.clm.cauchit.symmetric2[[i]]<-mat.square(table(pred.MAP.clm.cauchit.symmetric2[[i]],test[[i]]$rating),categories)
}


conf.mat.MAP.clm.cauchit.equidistant<-list()
for(i in 1:10)
{
  conf.mat.MAP.clm.cauchit.equidistant[[i]]<-mat.square(table(pred.MAP.clm.cauchit.equidistant[[i]],test[[i]]$rating),categories)
}


# MAE 

MAE.MAP.clm.logit.flexible<-vector()
for (i in 1:10)
{
  MAE.MAP.clm.logit.flexible[i]<-MAE(conf.mat.MAP.clm.logit.flexible[[i]])
}

MAE.MAP.clm.probit.flexible<-vector()
for (i in 1:10)
{
  MAE.MAP.clm.probit.flexible[i]<-MAE(conf.mat.MAP.clm.probit.flexible[[i]])
}

MAE.MAP.clm.loglog.flexible<-vector()
for (i in 1:10)
{
  MAE.MAP.clm.loglog.flexible[i]<-MAE(conf.mat.MAP.clm.loglog.flexible[[i]])
}

MAE.MAP.clm.cloglog.flexible<-vector()
for (i in 1:10)
{
  MAE.MAP.clm.cloglog.flexible[i]<-MAE(conf.mat.MAP.clm.cloglog.flexible[[i]])
}

MAE.MAP.clm.cauchit.flexible<-vector()
for (i in 1:10)
{
  MAE.MAP.clm.cauchit.flexible[i]<-MAE(conf.mat.MAP.clm.cauchit.flexible[[i]])
}

###

MAE.MAP.clm.logit.symmetric<-vector()
for (i in 1:10)
{
  MAE.MAP.clm.logit.symmetric[i]<-MAE(conf.mat.MAP.clm.logit.symmetric[[i]])
}

MAE.MAP.clm.probit.symmetric<-vector()
for (i in 1:10)
{
  MAE.MAP.clm.probit.symmetric[i]<-MAE(conf.mat.MAP.clm.probit.symmetric[[i]])
}

MAE.MAP.clm.loglog.symmetric<-vector()
for (i in 1:10)
{
  MAE.MAP.clm.loglog.symmetric[i]<-MAE(conf.mat.MAP.clm.loglog.symmetric[[i]])
}

MAE.MAP.clm.cloglog.symmetric<-vector()
for (i in 1:10)
{
  MAE.MAP.clm.cloglog.symmetric[i]<-MAE(conf.mat.MAP.clm.cloglog.symmetric[[i]])
}

MAE.MAP.clm.cauchit.symmetric<-vector()
for (i in 1:10)
{
  MAE.MAP.clm.cauchit.symmetric[i]<-MAE(conf.mat.MAP.clm.cauchit.symmetric[[i]])
}

###

MAE.MAP.clm.logit.symmetric2<-vector()
for (i in 1:10)
{
  MAE.MAP.clm.logit.symmetric2[i]<-MAE(conf.mat.MAP.clm.logit.symmetric2[[i]])
}

MAE.MAP.clm.probit.symmetric2<-vector()
for (i in 1:10)
{
  MAE.MAP.clm.probit.symmetric2[i]<-MAE(conf.mat.MAP.clm.probit.symmetric2[[i]])
}

MAE.MAP.clm.loglog.symmetric2<-vector()
for (i in 1:10)
{
  MAE.MAP.clm.loglog.symmetric2[i]<-MAE(conf.mat.MAP.clm.loglog.symmetric2[[i]])
}

MAE.MAP.clm.cloglog.symmetric2<-vector()
for (i in 1:10)
{
  MAE.MAP.clm.cloglog.symmetric2[i]<-MAE(conf.mat.MAP.clm.cloglog.symmetric2[[i]])
}

MAE.MAP.clm.cauchit.symmetric2<-vector()
for (i in 1:10)
{
  MAE.MAP.clm.cauchit.symmetric2[i]<-MAE(conf.mat.MAP.clm.cauchit.symmetric2[[i]])
}

###

MAE.MAP.clm.logit.equidistant<-vector()
for (i in 1:10)
{
  MAE.MAP.clm.logit.equidistant[i]<-MAE(conf.mat.MAP.clm.logit.equidistant[[i]])
}

MAE.MAP.clm.probit.equidistant<-vector()
for (i in 1:10)
{
  MAE.MAP.clm.probit.equidistant[i]<-MAE(conf.mat.MAP.clm.probit.equidistant[[i]])
}

MAE.MAP.clm.loglog.equidistant<-vector()
for (i in 1:10)
{
  MAE.MAP.clm.loglog.equidistant[i]<-MAE(conf.mat.MAP.clm.loglog.equidistant[[i]])
}

MAE.MAP.clm.cloglog.equidistant<-vector()
for (i in 1:10)
{
  MAE.MAP.clm.cloglog.equidistant[i]<-MAE(conf.mat.MAP.clm.cloglog.equidistant[[i]])
}

MAE.MAP.clm.cauchit.equidistant<-vector()
for (i in 1:10)
{
  MAE.MAP.clm.cauchit.equidistant[i]<-MAE(conf.mat.MAP.clm.cauchit.equidistant[[i]])
}



################################################################################
### predictions with Ord.MAP criterion:
###

pred.Ord.MAP.clm.logit.flexible<-list()
for (i in 1:10)
{ pred.Ord.MAP.clm.logit.flexible[[i]]<-vector()
for (j in 1:dim(pred.test.clm.logit.flexible[[i]])[1])
{h<-min(which(cumsum(pred.test.clm.logit.flexible[[i]][j,])>=0.5))
pred.Ord.MAP.clm.logit.flexible[[i]][j]<-categories[h]
}
}

pred.Ord.MAP.clm.logit.symmetric<-list()
for (i in 1:10)
{ pred.Ord.MAP.clm.logit.symmetric[[i]]<-vector()
for (j in 1:dim(pred.test.clm.logit.symmetric[[i]])[1])
{h<-min(which(cumsum(pred.test.clm.logit.symmetric[[i]][j,])>=0.5))
pred.Ord.MAP.clm.logit.symmetric[[i]][j]<-categories[h]
}
}

pred.Ord.MAP.clm.logit.symmetric2<-list()
for (i in 1:10)
{ pred.Ord.MAP.clm.logit.symmetric2[[i]]<-vector()
for (j in 1:dim(pred.test.clm.logit.symmetric2[[i]])[1])
{h<-min(which(cumsum(pred.test.clm.logit.symmetric2[[i]][j,])>=0.5))
pred.Ord.MAP.clm.logit.symmetric2[[i]][j]<-categories[h]
}
}

pred.Ord.MAP.clm.logit.equidistant<-list()
for (i in 1:10)
{ pred.Ord.MAP.clm.logit.equidistant[[i]]<-vector()
for (j in 1:dim(pred.test.clm.logit.equidistant[[i]])[1])
{h<-min(which(cumsum(pred.test.clm.logit.equidistant[[i]][j,])>=0.5))
pred.Ord.MAP.clm.logit.equidistant[[i]][j]<-categories[h]
}
}

###

pred.Ord.MAP.clm.probit.flexible<-list()
for (i in 1:10)
{ pred.Ord.MAP.clm.probit.flexible[[i]]<-vector()
for (j in 1:dim(pred.test.clm.probit.flexible[[i]])[1])
{h<-min(which(cumsum(pred.test.clm.probit.flexible[[i]][j,])>=0.5))
pred.Ord.MAP.clm.probit.flexible[[i]][j]<-categories[h]
}
}

pred.Ord.MAP.clm.probit.symmetric<-list()
for (i in 1:10)
{ pred.Ord.MAP.clm.probit.symmetric[[i]]<-vector()
for (j in 1:dim(pred.test.clm.probit.symmetric[[i]])[1])
{h<-min(which(cumsum(pred.test.clm.probit.symmetric[[i]][j,])>=0.5))
pred.Ord.MAP.clm.probit.symmetric[[i]][j]<-categories[h]
}
}

pred.Ord.MAP.clm.probit.symmetric2<-list()
for (i in 1:10)
{ pred.Ord.MAP.clm.probit.symmetric2[[i]]<-vector()
for (j in 1:dim(pred.test.clm.probit.symmetric2[[i]])[1])
{h<-min(which(cumsum(pred.test.clm.probit.symmetric2[[i]][j,])>=0.5))
pred.Ord.MAP.clm.probit.symmetric2[[i]][j]<-categories[h]
}
}

pred.Ord.MAP.clm.probit.equidistant<-list()
for (i in 1:10)
{ pred.Ord.MAP.clm.probit.equidistant[[i]]<-vector()
for (j in 1:dim(pred.test.clm.probit.equidistant[[i]])[1])
{h<-min(which(cumsum(pred.test.clm.probit.equidistant[[i]][j,])>=0.5))
pred.Ord.MAP.clm.probit.equidistant[[i]][j]<-categories[h]
}
}

###

pred.Ord.MAP.clm.loglog.flexible<-list()
for (i in 1:10)
{ pred.Ord.MAP.clm.loglog.flexible[[i]]<-vector()
for (j in 1:dim(pred.test.clm.loglog.flexible[[i]])[1])
{h<-min(which(cumsum(pred.test.clm.loglog.flexible[[i]][j,])>=0.5))
pred.Ord.MAP.clm.loglog.flexible[[i]][j]<-categories[h]
}
}

pred.Ord.MAP.clm.loglog.symmetric<-list()
for (i in 1:10)
{ pred.Ord.MAP.clm.loglog.symmetric[[i]]<-vector()
for (j in 1:dim(pred.test.clm.loglog.symmetric[[i]])[1])
{h<-min(which(cumsum(pred.test.clm.loglog.symmetric[[i]][j,])>=0.5))
pred.Ord.MAP.clm.loglog.symmetric[[i]][j]<-categories[h]
}
}

pred.Ord.MAP.clm.loglog.symmetric2<-list()
for (i in 1:10)
{ pred.Ord.MAP.clm.loglog.symmetric2[[i]]<-vector()
for (j in 1:dim(pred.test.clm.loglog.symmetric2[[i]])[1])
{h<-min(which(cumsum(pred.test.clm.loglog.symmetric2[[i]][j,])>=0.5))
pred.Ord.MAP.clm.loglog.symmetric2[[i]][j]<-categories[h]
}
}

pred.Ord.MAP.clm.loglog.equidistant<-list()
for (i in 1:10)
{ pred.Ord.MAP.clm.loglog.equidistant[[i]]<-vector()
for (j in 1:dim(pred.test.clm.loglog.equidistant[[i]])[1])
{h<-min(which(cumsum(pred.test.clm.loglog.equidistant[[i]][j,])>=0.5))
pred.Ord.MAP.clm.loglog.equidistant[[i]][j]<-categories[h]
}
}

###

pred.Ord.MAP.clm.cloglog.flexible<-list()
for (i in 1:10)
{ pred.Ord.MAP.clm.cloglog.flexible[[i]]<-vector()
for (j in 1:dim(pred.test.clm.cloglog.flexible[[i]])[1])
{h<-min(which(cumsum(pred.test.clm.cloglog.flexible[[i]][j,])>=0.5))
pred.Ord.MAP.clm.cloglog.flexible[[i]][j]<-categories[h]
}
}

pred.Ord.MAP.clm.cloglog.symmetric<-list()
for (i in 1:10)
{ pred.Ord.MAP.clm.cloglog.symmetric[[i]]<-vector()
for (j in 1:dim(pred.test.clm.cloglog.symmetric[[i]])[1])
{h<-min(which(cumsum(pred.test.clm.cloglog.symmetric[[i]][j,])>=0.5))
pred.Ord.MAP.clm.cloglog.symmetric[[i]][j]<-categories[h]
}
}

pred.Ord.MAP.clm.cloglog.symmetric2<-list()
for (i in 1:10)
{ pred.Ord.MAP.clm.cloglog.symmetric2[[i]]<-vector()
for (j in 1:dim(pred.test.clm.cloglog.symmetric2[[i]])[1])
{h<-min(which(cumsum(pred.test.clm.cloglog.symmetric2[[i]][j,])>=0.5))
pred.Ord.MAP.clm.cloglog.symmetric2[[i]][j]<-categories[h]
}
}

pred.Ord.MAP.clm.cloglog.equidistant<-list()
for (i in 1:10)
{ pred.Ord.MAP.clm.cloglog.equidistant[[i]]<-vector()
for (j in 1:dim(pred.test.clm.cloglog.equidistant[[i]])[1])
{h<-min(which(cumsum(pred.test.clm.cloglog.equidistant[[i]][j,])>=0.5))
pred.Ord.MAP.clm.cloglog.equidistant[[i]][j]<-categories[h]
}
}

###

pred.Ord.MAP.clm.cauchit.flexible<-list()
for (i in 1:10)
{ pred.Ord.MAP.clm.cauchit.flexible[[i]]<-vector()
for (j in 1:dim(pred.test.clm.cauchit.flexible[[i]])[1])
{h<-min(which(cumsum(pred.test.clm.cauchit.flexible[[i]][j,])>=0.5))
pred.Ord.MAP.clm.cauchit.flexible[[i]][j]<-categories[h]
}
}

pred.Ord.MAP.clm.cauchit.symmetric<-list()
for (i in 1:10)
{ pred.Ord.MAP.clm.cauchit.symmetric[[i]]<-vector()
for (j in 1:dim(pred.test.clm.cauchit.symmetric[[i]])[1])
{h<-min(which(cumsum(pred.test.clm.cauchit.symmetric[[i]][j,])>=0.5))
pred.Ord.MAP.clm.cauchit.symmetric[[i]][j]<-categories[h]
}
}

pred.Ord.MAP.clm.cauchit.symmetric2<-list()
for (i in 1:10)
{ pred.Ord.MAP.clm.cauchit.symmetric2[[i]]<-vector()
for (j in 1:dim(pred.test.clm.cauchit.symmetric2[[i]])[1])
{h<-min(which(cumsum(pred.test.clm.cauchit.symmetric2[[i]][j,])>=0.5))
pred.Ord.MAP.clm.cauchit.symmetric2[[i]][j]<-categories[h]
}
}

pred.Ord.MAP.clm.cauchit.equidistant<-list()
for (i in 1:10)
{ pred.Ord.MAP.clm.cauchit.equidistant[[i]]<-vector()
for (j in 1:dim(pred.test.clm.cauchit.equidistant[[i]])[1])
{h<-min(which(cumsum(pred.test.clm.cauchit.equidistant[[i]][j,])>=0.5))
pred.Ord.MAP.clm.cauchit.equidistant[[i]][j]<-categories[h]
}
}

###
# confusion matrices

conf.mat.Ord.MAP.clm.logit.flexible<-list()
for(i in 1:10)
{
  conf.mat.Ord.MAP.clm.logit.flexible[[i]]<-mat.square(table(pred.Ord.MAP.clm.logit.flexible[[i]],test[[i]]$rating),categories)
}

conf.mat.Ord.MAP.clm.logit.symmetric<-list()
for(i in 1:10)
{
  conf.mat.Ord.MAP.clm.logit.symmetric[[i]]<-mat.square(table(pred.Ord.MAP.clm.logit.symmetric[[i]],test[[i]]$rating),categories)
}

conf.mat.Ord.MAP.clm.logit.symmetric2<-list()
for(i in 1:10)
{
  conf.mat.Ord.MAP.clm.logit.symmetric2[[i]]<-mat.square(table(pred.Ord.MAP.clm.logit.symmetric2[[i]],test[[i]]$rating),categories)
}

conf.mat.Ord.MAP.clm.logit.equidistant<-list()
for(i in 1:10)
{
  conf.mat.Ord.MAP.clm.logit.equidistant[[i]]<-mat.square(table(pred.Ord.MAP.clm.logit.equidistant[[i]],test[[i]]$rating),categories)
}

###

conf.mat.Ord.MAP.clm.probit.flexible<-list()
for(i in 1:10)
{
  conf.mat.Ord.MAP.clm.probit.flexible[[i]]<-mat.square(table(pred.Ord.MAP.clm.probit.flexible[[i]],test[[i]]$rating),categories)
}

conf.mat.Ord.MAP.clm.probit.symmetric<-list()
for(i in 1:10)
{
  conf.mat.Ord.MAP.clm.probit.symmetric[[i]]<-mat.square(table(pred.Ord.MAP.clm.probit.symmetric[[i]],test[[i]]$rating),categories)
}

conf.mat.Ord.MAP.clm.probit.symmetric2<-list()
for(i in 1:10)
{
  conf.mat.Ord.MAP.clm.probit.symmetric2[[i]]<-mat.square(table(pred.Ord.MAP.clm.probit.symmetric2[[i]],test[[i]]$rating),categories)
}

conf.mat.Ord.MAP.clm.probit.equidistant<-list()
for(i in 1:10)
{
  conf.mat.Ord.MAP.clm.probit.equidistant[[i]]<-mat.square(table(pred.Ord.MAP.clm.probit.equidistant[[i]],test[[i]]$rating),categories)
}

###

conf.mat.Ord.MAP.clm.loglog.flexible<-list()
for(i in 1:10)
{
  conf.mat.Ord.MAP.clm.loglog.flexible[[i]]<-mat.square(table(pred.Ord.MAP.clm.loglog.flexible[[i]],test[[i]]$rating),categories)
}

conf.mat.Ord.MAP.clm.loglog.symmetric<-list()
for(i in 1:10)
{
  conf.mat.Ord.MAP.clm.loglog.symmetric[[i]]<-mat.square(table(pred.Ord.MAP.clm.loglog.symmetric[[i]],test[[i]]$rating),categories)
}

conf.mat.Ord.MAP.clm.loglog.symmetric2<-list()
for(i in 1:10)
{
  conf.mat.Ord.MAP.clm.loglog.symmetric2[[i]]<-mat.square(table(pred.Ord.MAP.clm.loglog.symmetric2[[i]],test[[i]]$rating),categories)
}

conf.mat.Ord.MAP.clm.loglog.equidistant<-list()
for(i in 1:10)
{
  conf.mat.Ord.MAP.clm.loglog.equidistant[[i]]<-mat.square(table(pred.Ord.MAP.clm.loglog.equidistant[[i]],test[[i]]$rating),categories)
}

###


conf.mat.Ord.MAP.clm.cloglog.flexible<-list()
for(i in 1:10)
{
  conf.mat.Ord.MAP.clm.cloglog.flexible[[i]]<-mat.square(table(pred.Ord.MAP.clm.cloglog.flexible[[i]],test[[i]]$rating),categories)
}

conf.mat.Ord.MAP.clm.cloglog.symmetric<-list()
for(i in 1:10)
{
  conf.mat.Ord.MAP.clm.cloglog.symmetric[[i]]<-mat.square(table(pred.Ord.MAP.clm.cloglog.symmetric[[i]],test[[i]]$rating),categories)
}

conf.mat.Ord.MAP.clm.cloglog.symmetric2<-list()
for(i in 1:10)
{
  conf.mat.Ord.MAP.clm.cloglog.symmetric2[[i]]<-mat.square(table(pred.Ord.MAP.clm.cloglog.symmetric2[[i]],test[[i]]$rating),categories)
}

conf.mat.Ord.MAP.clm.cloglog.equidistant<-list()
for(i in 1:10)
{
  conf.mat.Ord.MAP.clm.cloglog.equidistant[[i]]<-mat.square(table(pred.Ord.MAP.clm.cloglog.equidistant[[i]],test[[i]]$rating),categories)
}

###

conf.mat.Ord.MAP.clm.cauchit.flexible<-list()
for(i in 1:10)
{
  conf.mat.Ord.MAP.clm.cauchit.flexible[[i]]<-mat.square(table(pred.Ord.MAP.clm.cauchit.flexible[[i]],test[[i]]$rating),categories)
}

conf.mat.Ord.MAP.clm.cauchit.symmetric<-list()
for(i in 1:10)
{
  conf.mat.Ord.MAP.clm.cauchit.symmetric[[i]]<-mat.square(table(pred.Ord.MAP.clm.cauchit.symmetric[[i]],test[[i]]$rating),categories)
}

conf.mat.Ord.MAP.clm.cauchit.symmetric2<-list()
for(i in 1:10)
{
  conf.mat.Ord.MAP.clm.cauchit.symmetric2[[i]]<-mat.square(table(pred.Ord.MAP.clm.cauchit.symmetric2[[i]],test[[i]]$rating),categories)
}

conf.mat.Ord.MAP.clm.cauchit.equidistant<-list()
for(i in 1:10)
{
  conf.mat.Ord.MAP.clm.cauchit.equidistant[[i]]<-mat.square(table(pred.Ord.MAP.clm.cauchit.equidistant[[i]],test[[i]]$rating),categories)
}

###

# MAE 

MAE.Ord.MAP.clm.logit.flexible<-vector()
for (i in 1:10)
{
  MAE.Ord.MAP.clm.logit.flexible[i]<-MAE(conf.mat.Ord.MAP.clm.logit.flexible[[i]])
}

MAE.Ord.MAP.clm.logit.symmetric<-vector()
for (i in 1:10)
{
  MAE.Ord.MAP.clm.logit.symmetric[i]<-MAE(conf.mat.Ord.MAP.clm.logit.symmetric[[i]])
}

MAE.Ord.MAP.clm.logit.symmetric2<-vector()
for (i in 1:10)
{
  MAE.Ord.MAP.clm.logit.symmetric2[i]<-MAE(conf.mat.Ord.MAP.clm.logit.symmetric2[[i]])
}

MAE.Ord.MAP.clm.logit.equidistant<-vector()
for (i in 1:10)
{
  MAE.Ord.MAP.clm.logit.equidistant[i]<-MAE(conf.mat.Ord.MAP.clm.logit.equidistant[[i]])
}

###

MAE.Ord.MAP.clm.probit.flexible<-vector()
for (i in 1:10)
{
  MAE.Ord.MAP.clm.probit.flexible[i]<-MAE(conf.mat.Ord.MAP.clm.probit.flexible[[i]])
}

MAE.Ord.MAP.clm.probit.symmetric<-vector()
for (i in 1:10)
{
  MAE.Ord.MAP.clm.probit.symmetric[i]<-MAE(conf.mat.Ord.MAP.clm.probit.symmetric[[i]])
}

MAE.Ord.MAP.clm.probit.symmetric2<-vector()
for (i in 1:10)
{
  MAE.Ord.MAP.clm.probit.symmetric2[i]<-MAE(conf.mat.Ord.MAP.clm.probit.symmetric2[[i]])
}

MAE.Ord.MAP.clm.probit.equidistant<-vector()
for (i in 1:10)
{
  MAE.Ord.MAP.clm.probit.equidistant[i]<-MAE(conf.mat.Ord.MAP.clm.probit.equidistant[[i]])
}

###


MAE.Ord.MAP.clm.loglog.flexible<-vector()
for (i in 1:10)
{
  MAE.Ord.MAP.clm.loglog.flexible[i]<-MAE(conf.mat.Ord.MAP.clm.loglog.flexible[[i]])
}

MAE.Ord.MAP.clm.loglog.symmetric<-vector()
for (i in 1:10)
{
  MAE.Ord.MAP.clm.loglog.symmetric[i]<-MAE(conf.mat.Ord.MAP.clm.loglog.symmetric[[i]])
}

MAE.Ord.MAP.clm.loglog.symmetric2<-vector()
for (i in 1:10)
{
  MAE.Ord.MAP.clm.loglog.symmetric2[i]<-MAE(conf.mat.Ord.MAP.clm.loglog.symmetric2[[i]])
}

MAE.Ord.MAP.clm.loglog.equidistant<-vector()
for (i in 1:10)
{
  MAE.Ord.MAP.clm.loglog.equidistant[i]<-MAE(conf.mat.Ord.MAP.clm.loglog.equidistant[[i]])
}

###


MAE.Ord.MAP.clm.cloglog.flexible<-vector()
for (i in 1:10)
{
  MAE.Ord.MAP.clm.cloglog.flexible[i]<-MAE(conf.mat.Ord.MAP.clm.cloglog.flexible[[i]])
}

MAE.Ord.MAP.clm.cloglog.symmetric<-vector()
for (i in 1:10)
{
  MAE.Ord.MAP.clm.cloglog.symmetric[i]<-MAE(conf.mat.Ord.MAP.clm.cloglog.symmetric[[i]])
}

MAE.Ord.MAP.clm.cloglog.symmetric2<-vector()
for (i in 1:10)
{
  MAE.Ord.MAP.clm.cloglog.symmetric2[i]<-MAE(conf.mat.Ord.MAP.clm.cloglog.symmetric2[[i]])
}

MAE.Ord.MAP.clm.cloglog.equidistant<-vector()
for (i in 1:10)
{
  MAE.Ord.MAP.clm.cloglog.equidistant[i]<-MAE(conf.mat.Ord.MAP.clm.cloglog.equidistant[[i]])
}

###


MAE.Ord.MAP.clm.cauchit.flexible<-vector()
for (i in 1:10)
{
  MAE.Ord.MAP.clm.cauchit.flexible[i]<-MAE(conf.mat.Ord.MAP.clm.cauchit.flexible[[i]])
}

MAE.Ord.MAP.clm.cauchit.symmetric<-vector()
for (i in 1:10)
{
  MAE.Ord.MAP.clm.cauchit.symmetric[i]<-MAE(conf.mat.Ord.MAP.clm.cauchit.symmetric[[i]])
}

MAE.Ord.MAP.clm.cauchit.symmetric2<-vector()
for (i in 1:10)
{
  MAE.Ord.MAP.clm.cauchit.symmetric2[i]<-MAE(conf.mat.Ord.MAP.clm.cauchit.symmetric2[[i]])
}

MAE.Ord.MAP.clm.cauchit.equidistant<-vector()
for (i in 1:10)
{
  MAE.Ord.MAP.clm.cauchit.equidistant[i]<-MAE(conf.mat.Ord.MAP.clm.cauchit.equidistant[[i]])
}


###
################################################################################
###. Comparison MAP vs Ord.MAP
###
###.  Tests MAE

p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.logit.flexible<-wilcox.test(MAE.MAP.clm.logit.flexible,MAE.Ord.MAP.clm.logit.flexible,paired = TRUE,alternative="greater")$p.value
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.logit.symmetric<-wilcox.test(MAE.MAP.clm.logit.symmetric,MAE.Ord.MAP.clm.logit.symmetric,paired = TRUE,alternative="greater")$p.value
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.logit.symmetric2<-wilcox.test(MAE.MAP.clm.logit.symmetric2,MAE.Ord.MAP.clm.logit.symmetric2,paired = TRUE,alternative="greater")$p.value
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.logit.equidistant<-wilcox.test(MAE.MAP.clm.logit.equidistant,MAE.Ord.MAP.clm.logit.equidistant,paired = TRUE,alternative="greater")$p.value

p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.logit.flexible.less<-wilcox.test(MAE.MAP.clm.logit.flexible,MAE.Ord.MAP.clm.logit.flexible,paired = TRUE,alternative="less")$p.value
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.logit.symmetric.less<-wilcox.test(MAE.MAP.clm.logit.symmetric,MAE.Ord.MAP.clm.logit.symmetric,paired = TRUE,alternative="less")$p.value
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.logit.symmetric2.less<-wilcox.test(MAE.MAP.clm.logit.symmetric2,MAE.Ord.MAP.clm.logit.symmetric2,paired = TRUE,alternative="less")$p.value
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.logit.equidistant.less<-wilcox.test(MAE.MAP.clm.logit.equidistant,MAE.Ord.MAP.clm.logit.equidistant,paired = TRUE,alternative="less")$p.value

###

p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.probit.flexible<-wilcox.test(MAE.MAP.clm.probit.flexible,MAE.Ord.MAP.clm.probit.flexible,paired = TRUE,alternative="greater")$p.value
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.probit.symmetric<-wilcox.test(MAE.MAP.clm.probit.symmetric,MAE.Ord.MAP.clm.probit.symmetric,paired = TRUE,alternative="greater")$p.value
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.probit.symmetric2<-wilcox.test(MAE.MAP.clm.probit.symmetric2,MAE.Ord.MAP.clm.probit.symmetric2,paired = TRUE,alternative="greater")$p.value
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.probit.equidistant<-wilcox.test(MAE.MAP.clm.probit.equidistant,MAE.Ord.MAP.clm.probit.equidistant,paired = TRUE,alternative="greater")$p.value

p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.probit.flexible.less<-wilcox.test(MAE.MAP.clm.probit.flexible,MAE.Ord.MAP.clm.probit.flexible,paired = TRUE,alternative="less")$p.value
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.probit.symmetric.less<-wilcox.test(MAE.MAP.clm.probit.symmetric,MAE.Ord.MAP.clm.probit.symmetric,paired = TRUE,alternative="less")$p.value
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.probit.symmetric2.less<-wilcox.test(MAE.MAP.clm.probit.symmetric2,MAE.Ord.MAP.clm.probit.symmetric2,paired = TRUE,alternative="less")$p.value
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.probit.equidistant.less<-wilcox.test(MAE.MAP.clm.probit.equidistant,MAE.Ord.MAP.clm.probit.equidistant,paired = TRUE,alternative="less")$p.value


###

p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.loglog.flexible<-wilcox.test(MAE.MAP.clm.loglog.flexible,MAE.Ord.MAP.clm.loglog.flexible,paired = TRUE,alternative="greater")$p.value
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.loglog.symmetric<-wilcox.test(MAE.MAP.clm.loglog.symmetric,MAE.Ord.MAP.clm.loglog.symmetric,paired = TRUE,alternative="greater")$p.value
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.loglog.symmetric2<-wilcox.test(MAE.MAP.clm.loglog.symmetric2,MAE.Ord.MAP.clm.loglog.symmetric2,paired = TRUE,alternative="greater")$p.value
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.loglog.equidistant<-wilcox.test(MAE.MAP.clm.loglog.equidistant,MAE.Ord.MAP.clm.loglog.equidistant,paired = TRUE,alternative="greater")$p.value

p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.loglog.flexible.less<-wilcox.test(MAE.MAP.clm.loglog.flexible,MAE.Ord.MAP.clm.loglog.flexible,paired = TRUE,alternative="less")$p.value
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.loglog.symmetric.less<-wilcox.test(MAE.MAP.clm.loglog.symmetric,MAE.Ord.MAP.clm.loglog.symmetric,paired = TRUE,alternative="less")$p.value
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.loglog.symmetric2.less<-wilcox.test(MAE.MAP.clm.loglog.symmetric2,MAE.Ord.MAP.clm.loglog.symmetric2,paired = TRUE,alternative="less")$p.value
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.loglog.equidistant.less<-wilcox.test(MAE.MAP.clm.loglog.equidistant,MAE.Ord.MAP.clm.loglog.equidistant,paired = TRUE,alternative="less")$p.value

###

p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.cloglog.flexible<-wilcox.test(MAE.MAP.clm.cloglog.flexible,MAE.Ord.MAP.clm.cloglog.flexible,paired = TRUE,alternative="greater")$p.value
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.cloglog.symmetric<-wilcox.test(MAE.MAP.clm.cloglog.symmetric,MAE.Ord.MAP.clm.cloglog.symmetric,paired = TRUE,alternative="greater")$p.value
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.cloglog.symmetric2<-wilcox.test(MAE.MAP.clm.cloglog.symmetric2,MAE.Ord.MAP.clm.cloglog.symmetric2,paired = TRUE,alternative="greater")$p.value
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.cloglog.equidistant<-wilcox.test(MAE.MAP.clm.cloglog.equidistant,MAE.Ord.MAP.clm.cloglog.equidistant,paired = TRUE,alternative="greater")$p.value

p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.cloglog.flexible.less<-wilcox.test(MAE.MAP.clm.cloglog.flexible,MAE.Ord.MAP.clm.cloglog.flexible,paired = TRUE,alternative="less")$p.value
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.cloglog.symmetric.less<-wilcox.test(MAE.MAP.clm.cloglog.symmetric,MAE.Ord.MAP.clm.cloglog.symmetric,paired = TRUE,alternative="less")$p.value
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.cloglog.symmetric2.less<-wilcox.test(MAE.MAP.clm.cloglog.symmetric2,MAE.Ord.MAP.clm.cloglog.symmetric2,paired = TRUE,alternative="less")$p.value
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.cloglog.equidistant.less<-wilcox.test(MAE.MAP.clm.cloglog.equidistant,MAE.Ord.MAP.clm.cloglog.equidistant,paired = TRUE,alternative="les")$p.value

###

p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.cauchit.flexible<-wilcox.test(MAE.MAP.clm.cauchit.flexible,MAE.Ord.MAP.clm.cauchit.flexible,paired = TRUE,alternative="greater")$p.value
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.cauchit.symmetric<-wilcox.test(MAE.MAP.clm.cauchit.symmetric,MAE.Ord.MAP.clm.cauchit.symmetric,paired = TRUE,alternative="greater")$p.value
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.cauchit.symmetric2<-wilcox.test(MAE.MAP.clm.cauchit.symmetric2,MAE.Ord.MAP.clm.cauchit.symmetric2,paired = TRUE,alternative="greater")$p.value
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.cauchit.equidistant<-wilcox.test(MAE.MAP.clm.cauchit.equidistant,MAE.Ord.MAP.clm.cauchit.equidistant,paired = TRUE,alternative="greater")$p.value

p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.cauchit.flexible.less<-wilcox.test(MAE.MAP.clm.cauchit.flexible,MAE.Ord.MAP.clm.cauchit.flexible,paired = TRUE,alternative="less")$p.value
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.cauchit.symmetric.less<-wilcox.test(MAE.MAP.clm.cauchit.symmetric,MAE.Ord.MAP.clm.cauchit.symmetric,paired = TRUE,alternative="less")$p.value
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.cauchit.symmetric2.less<-wilcox.test(MAE.MAP.clm.cauchit.symmetric2,MAE.Ord.MAP.clm.cauchit.symmetric2,paired = TRUE,alternative="less")$p.value
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.cauchit.equidistant.less<-wilcox.test(MAE.MAP.clm.cauchit.equidistant,MAE.Ord.MAP.clm.cauchit.equidistant,paired = TRUE,alternative="less")$p.value

###

p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.logit.flexible
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.logit.symmetric
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.logit.symmetric2
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.logit.equidistant

p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.logit.flexible.less
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.logit.symmetric.less
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.logit.symmetric2.less
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.logit.equidistant.less

# 

p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.probit.flexible
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.probit.symmetric
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.probit.symmetric2
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.probit.equidistant

p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.probit.flexible.less
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.probit.symmetric.less
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.probit.symmetric2.less
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.probit.equidistant.less

#

p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.loglog.flexible
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.loglog.symmetric
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.loglog.symmetric2
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.loglog.equidistant

p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.loglog.flexible.less
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.loglog.symmetric.less
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.loglog.symmetric2.less
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.loglog.equidistant.less

#

p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.cloglog.flexible
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.cloglog.symmetric
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.cloglog.symmetric2
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.cloglog.equidistant

p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.cloglog.flexible.less
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.cloglog.symmetric.less
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.cloglog.symmetric2.less
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.cloglog.equidistant.less

#

p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.cauchit.flexible
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.cauchit.symmetric
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.cauchit.symmetric2
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.cauchit.equidistant

p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.cauchit.flexible.less
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.cauchit.symmetric.less
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.cauchit.symmetric2.less
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.cauchit.equidistant.less


###############################
###############################
###   Boxplots   MAE

MAE.clm.flexible.results<-as.data.frame(cbind(MAE.MAP.clm.logit.flexible-MAE.Ord.MAP.clm.logit.flexible,
                                              MAE.MAP.clm.probit.flexible-MAE.Ord.MAP.clm.probit.flexible,
                                              MAE.MAP.clm.loglog.flexible-MAE.Ord.MAP.clm.loglog.flexible,
                                              MAE.MAP.clm.cloglog.flexible-MAE.Ord.MAP.clm.cloglog.flexible,
                                              MAE.MAP.clm.cauchit.flexible-MAE.Ord.MAP.clm.cauchit.flexible))

boxplot(MAE.clm.flexible.results,
        main=paste("(b) Wine dataset. Threshold flexible"), 
        xlab=" ", ylab="MAE increment (MAP - Ord.MAP)",names=c("logistic", "probit","loglog","cloglog","cauchit"),
        col="white")
abline(h = 0, col = 2, lwd = 2, lty=2)
# Dots
stripchart(MAE.clm.flexible.results,
           #method = "jitter",
           pch = 19,
           col = 1:5,
           vertical = TRUE,
           add = TRUE)

###

MAE.clm.symmetric.results<-as.data.frame(cbind(MAE.MAP.clm.logit.symmetric-MAE.Ord.MAP.clm.logit.symmetric,
                                               MAE.MAP.clm.probit.symmetric-MAE.Ord.MAP.clm.probit.symmetric,
                                               MAE.MAP.clm.loglog.symmetric-MAE.Ord.MAP.clm.loglog.symmetric,
                                               MAE.MAP.clm.cloglog.symmetric-MAE.Ord.MAP.clm.cloglog.symmetric,
                                               MAE.MAP.clm.cauchit.symmetric-MAE.Ord.MAP.clm.cauchit.symmetric))

boxplot(MAE.clm.symmetric.results,
        main=paste("(b) Wine dataset. Threshold symmetric"), 
        xlab=" ", ylab="MAE increment (MAP - Ord.MAP)",names=c("logistic", "probit","loglog","cloglog","cauchit"),
        col="white")
abline(h = 0, col = 2, lwd = 2, lty=2)
# Dots
stripchart(MAE.clm.symmetric.results,
           #method = "jitter",
           pch = 19,
           col = 1:5,
           vertical = TRUE,
           add = TRUE)

###

MAE.clm.symmetric2.results<-as.data.frame(cbind(MAE.MAP.clm.logit.symmetric2-MAE.Ord.MAP.clm.logit.symmetric2,
                                                MAE.MAP.clm.probit.symmetric2-MAE.Ord.MAP.clm.probit.symmetric2,
                                                MAE.MAP.clm.loglog.symmetric2-MAE.Ord.MAP.clm.loglog.symmetric2,
                                                MAE.MAP.clm.cloglog.symmetric2-MAE.Ord.MAP.clm.cloglog.symmetric2,
                                                MAE.MAP.clm.cauchit.symmetric2-MAE.Ord.MAP.clm.cauchit.symmetric2))

boxplot(MAE.clm.symmetric2.results,
        main=paste("(b) Wine dataset. Threshold symmetric2"), 
        xlab=" ", ylab="MAE increment (MAP - Ord.MAP)",names=c("logistic", "probit","loglog","cloglog","cauchit"),
        col="white")
abline(h = 0, col = 2, lwd = 2, lty=2)
# Dots
stripchart(MAE.clm.symmetric2.results,
           #method = "jitter",
           pch = 19,
           col = 1:5,
           vertical = TRUE,
           add = TRUE)

###

MAE.clm.equidistant.results<-as.data.frame(cbind(MAE.MAP.clm.logit.equidistant-MAE.Ord.MAP.clm.logit.equidistant,
                                                 MAE.MAP.clm.probit.equidistant-MAE.Ord.MAP.clm.probit.equidistant,
                                                 MAE.MAP.clm.loglog.equidistant-MAE.Ord.MAP.clm.loglog.equidistant,
                                                 MAE.MAP.clm.cloglog.equidistant-MAE.Ord.MAP.clm.cloglog.equidistant,
                                                 MAE.MAP.clm.cauchit.equidistant-MAE.Ord.MAP.clm.cauchit.equidistant))

boxplot(MAE.clm.equidistant.results,
        main=paste("(b) Wine dataset. Threshold equidistant"), 
        xlab=" ", ylab="MAE increment (MAP - Ord.MAP)",names=c("logistic", "probit","loglog","cloglog","cauchit"),
        col="white")
abline(h = 0, col = 2, lwd = 2, lty=2)
# Dots
stripchart(MAE.clm.equidistant.results,
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
# clm, from "ordinal" package, rfits cumulative link models (CLMs) such as the propotional odds model. 
# The model allows for various link functions and structured thresholds that restricts 
# the thresholds or cut-points to be e.g.,
# equidistant or symmetrically arranged around the central threshold(s). 
#
# clm(formula, scale, nominal, data, weights, start, subset, doFit = TRUE,
#     na.action, contrasts, model = TRUE, control=list(),
#     link = c("logit", "probit", "cloglog", "loglog", "cauchit",
#              "Aranda-Ordaz", "log-gamma"),
#     threshold = c("flexible", "symmetric", "symmetric2", "equidistant"), ...)

model.clm.logit.flexible<-list()
for (i in 1:10)
{model.clm.logit.flexible[[i]] <- clm(Class ~ age + sex + chest_pain + trestbps + chol +
                                        fbs + restecg + thalach + exang + oldpeak, data = training[[i]],
                                      link="logit", threshold="flexible")}


model.clm.logit.symmetric<-list()
for (i in 1:10)
{model.clm.logit.symmetric[[i]] <- clm(Class ~ age + sex + chest_pain + trestbps + chol +
                                         fbs + restecg + thalach + exang + oldpeak, data = training[[i]],
                                       link="logit", threshold="symmetric")}


model.clm.logit.symmetric2<-list()
for (i in 1:10)
{model.clm.logit.symmetric2[[i]] <- clm(Class ~ age + sex + chest_pain + trestbps + chol +
                                          fbs + restecg + thalach + exang + oldpeak, data = training[[i]],
                                        link="logit", threshold="symmetric2")}

model.clm.logit.equidistant<-list()
for (i in 1:10)
{model.clm.logit.equidistant[[i]] <- clm(Class ~ age + sex + chest_pain + trestbps + chol +
                                           fbs + restecg + thalach + exang + oldpeak, data = training[[i]],
                                         link="logit", threshold="equidistant")}


###

model.clm.probit.flexible<-list()
for (i in 1:10)
{model.clm.probit.flexible[[i]] <- clm(Class ~ age + sex + chest_pain + trestbps + chol +
                                         fbs + restecg + thalach + exang + oldpeak, data = training[[i]],
                                       link="probit", threshold="flexible")}


model.clm.probit.symmetric<-list()
for (i in 1:10)
{model.clm.probit.symmetric[[i]] <- clm(Class ~ age + sex + chest_pain + trestbps + chol +
                                          fbs + restecg + thalach + exang + oldpeak, data = training[[i]],
                                        link="probit", threshold="symmetric")}


model.clm.probit.symmetric2<-list()
for (i in 1:10)
{model.clm.probit.symmetric2[[i]] <- clm(Class ~ age + sex + chest_pain + trestbps + chol +
                                           fbs + restecg + thalach + exang + oldpeak, data = training[[i]],
                                         link="probit", threshold="symmetric2")}

model.clm.probit.equidistant<-list()
for (i in 1:10)
{model.clm.probit.equidistant[[i]] <- clm(Class ~ age + sex + chest_pain + trestbps + chol +
                                            fbs + restecg + thalach + exang + oldpeak, data = training[[i]],
                                          link="probit", threshold="equidistant")}


###

model.clm.loglog.flexible<-list()
for (i in 1:10)
{model.clm.loglog.flexible[[i]] <- clm(Class ~ age + sex + chest_pain + trestbps + chol +
                                         fbs + restecg + thalach + exang + oldpeak, data = training[[i]],
                                       link="loglog", threshold="flexible")}


model.clm.loglog.symmetric<-list()
for (i in 1:10)
{model.clm.loglog.symmetric[[i]] <- clm(Class ~ age + sex + chest_pain + trestbps + chol +
                                          fbs + restecg + thalach + exang + oldpeak, data = training[[i]],
                                        link="loglog", threshold="symmetric")}


model.clm.loglog.symmetric2<-list()
for (i in 1:10)
{model.clm.loglog.symmetric2[[i]] <- clm(Class ~ age + sex + chest_pain + trestbps + chol +
                                           fbs + restecg + thalach + exang + oldpeak, data = training[[i]],
                                         link="loglog", threshold="symmetric2")}

model.clm.loglog.equidistant<-list()
for (i in 1:10)
{model.clm.loglog.equidistant[[i]] <- clm(Class ~ age + sex + chest_pain + trestbps + chol +
                                            fbs + restecg + thalach + exang + oldpeak, data = training[[i]],
                                          link="loglog", threshold="equidistant")}


###

model.clm.cloglog.flexible<-list()
for (i in 1:10)
{model.clm.cloglog.flexible[[i]] <- clm(Class ~ age + sex + chest_pain + trestbps + chol +
                                          fbs + restecg + thalach + exang + oldpeak, data = training[[i]],
                                        link="cloglog", threshold="flexible")}


model.clm.cloglog.symmetric<-list()
for (i in 1:10)
{model.clm.cloglog.symmetric[[i]] <- clm(Class ~ age + sex + chest_pain + trestbps + chol +
                                           fbs + restecg + thalach + exang + oldpeak, data = training[[i]],
                                         link="cloglog", threshold="symmetric")}


model.clm.cloglog.symmetric2<-list()
for (i in 1:10)
{model.clm.cloglog.symmetric2[[i]] <- clm(Class ~ age + sex + chest_pain + trestbps + chol +
                                            fbs + restecg + thalach + exang + oldpeak, data = training[[i]],
                                          link="cloglog", threshold="symmetric2")}

model.clm.cloglog.equidistant<-list()
for (i in 1:10)
{model.clm.cloglog.equidistant[[i]] <- clm(Class ~ age + sex + chest_pain + trestbps + chol +
                                             fbs + restecg + thalach + exang + oldpeak, data = training[[i]],
                                           link="cloglog", threshold="equidistant")}


###

model.clm.cauchit.flexible<-list()
for (i in 1:10)
{model.clm.cauchit.flexible[[i]] <- clm(Class ~ age + sex + chest_pain + trestbps + chol +
                                          fbs + restecg + thalach + exang + oldpeak, data = training[[i]],
                                        link="cauchit", threshold="flexible")}


model.clm.cauchit.symmetric<-list()
for (i in 1:10)
{model.clm.cauchit.symmetric[[i]] <- clm(Class ~ age + sex + chest_pain + trestbps + chol +
                                           fbs + restecg + thalach + exang + oldpeak, data = training[[i]],
                                         link="cauchit", threshold="symmetric")}


model.clm.cauchit.symmetric2<-list()
for (i in 1:10)
{model.clm.cauchit.symmetric2[[i]] <- clm(Class ~ age + sex + chest_pain + trestbps + chol +
                                            fbs + restecg + thalach + exang + oldpeak, data = training[[i]],
                                          link="cauchit", threshold="symmetric2")}

model.clm.cauchit.equidistant<-list()
for (i in 1:10)
{model.clm.cauchit.equidistant[[i]] <- clm(Class ~ age + sex + chest_pain + trestbps + chol +
                                             fbs + restecg + thalach + exang + oldpeak, data = training[[i]],
                                           link="cauchit", threshold="equidistant")}



################################################################################
#
# Predictions of any of the test set
# 
#

pred.test.clm.logit.flexible<-list()
for (i in 1:10)
{pred.test.clm.logit.flexible[[i]] <- predict(model.clm.logit.flexible[[i]],test[[i]][,-11],type = "p")$fit
}

pred.test.clm.logit.symmetric<-list()
for (i in 1:10)
{pred.test.clm.logit.symmetric[[i]] <- predict(model.clm.logit.symmetric[[i]],test[[i]][,-11],type = "p")$fit
}

pred.test.clm.logit.symmetric2<-list()
for (i in 1:10)
{pred.test.clm.logit.symmetric2[[i]] <- predict(model.clm.logit.symmetric2[[i]],test[[i]][,-11],type = "p")$fit
}

pred.test.clm.logit.equidistant<-list()
for (i in 1:10)
{pred.test.clm.logit.equidistant[[i]] <- predict(model.clm.logit.equidistant[[i]],test[[i]][,-11],type = "p")$fit
}

###

pred.test.clm.probit.flexible<-list()
for (i in 1:10)
{pred.test.clm.probit.flexible[[i]] <- predict(model.clm.probit.flexible[[i]],test[[i]][,-11],type = "p")$fit
}

pred.test.clm.probit.symmetric<-list()
for (i in 1:10)
{pred.test.clm.probit.symmetric[[i]] <- predict(model.clm.probit.symmetric[[i]],test[[i]][,-11],type = "p")$fit
}

pred.test.clm.probit.symmetric2<-list()
for (i in 1:10)
{pred.test.clm.probit.symmetric2[[i]] <- predict(model.clm.probit.symmetric2[[i]],test[[i]][,-11],type = "p")$fit
}

pred.test.clm.probit.equidistant<-list()
for (i in 1:10)
{pred.test.clm.probit.equidistant[[i]] <- predict(model.clm.probit.equidistant[[i]],test[[i]][,-11],type = "p")$fit
}

###

pred.test.clm.loglog.flexible<-list()
for (i in 1:10)
{pred.test.clm.loglog.flexible[[i]] <- predict(model.clm.loglog.flexible[[i]],test[[i]][,-11],type = "p")$fit
}

pred.test.clm.loglog.symmetric<-list()
for (i in 1:10)
{pred.test.clm.loglog.symmetric[[i]] <- predict(model.clm.loglog.symmetric[[i]],test[[i]][,-11],type = "p")$fit
}

pred.test.clm.loglog.symmetric2<-list()
for (i in 1:10)
{pred.test.clm.loglog.symmetric2[[i]] <- predict(model.clm.loglog.symmetric2[[i]],test[[i]][,-11],type = "p")$fit
}

pred.test.clm.loglog.equidistant<-list()
for (i in 1:10)
{pred.test.clm.loglog.equidistant[[i]] <- predict(model.clm.loglog.equidistant[[i]],test[[i]][,-11],type = "p")$fit
}

###


pred.test.clm.cloglog.flexible<-list()
for (i in 1:10)
{pred.test.clm.cloglog.flexible[[i]] <- predict(model.clm.cloglog.flexible[[i]],test[[i]][,-11],type = "p")$fit
}

pred.test.clm.cloglog.symmetric<-list()
for (i in 1:10)
{pred.test.clm.cloglog.symmetric[[i]] <- predict(model.clm.cloglog.symmetric[[i]],test[[i]][,-11],type = "p")$fit
}

pred.test.clm.cloglog.symmetric2<-list()
for (i in 1:10)
{pred.test.clm.cloglog.symmetric2[[i]] <- predict(model.clm.cloglog.symmetric2[[i]],test[[i]][,-11],type = "p")$fit
}

pred.test.clm.cloglog.equidistant<-list()
for (i in 1:10)
{pred.test.clm.cloglog.equidistant[[i]] <- predict(model.clm.cloglog.equidistant[[i]],test[[i]][,-11],type = "p")$fit
}

###


pred.test.clm.cauchit.flexible<-list()
for (i in 1:10)
{pred.test.clm.cauchit.flexible[[i]] <- predict(model.clm.cauchit.flexible[[i]],test[[i]][,-11],type = "p")$fit
}

pred.test.clm.cauchit.symmetric<-list()
for (i in 1:10)
{pred.test.clm.cauchit.symmetric[[i]] <- predict(model.clm.cauchit.symmetric[[i]],test[[i]][,-11],type = "p")$fit
}

pred.test.clm.cauchit.symmetric2<-list()
for (i in 1:10)
{pred.test.clm.cauchit.symmetric2[[i]] <- predict(model.clm.cauchit.symmetric2[[i]],test[[i]][,-11],type = "p")$fit
}

pred.test.clm.cauchit.equidistant<-list()
for (i in 1:10)
{pred.test.clm.cauchit.equidistant[[i]] <- predict(model.clm.cauchit.equidistant[[i]],test[[i]][,-11],type = "p")$fit
}

###

################################################################################
### predictions with standard MAP criterion:
###

pred.MAP.clm.logit.flexible<-list()
for (i in 1:10)
{pred.MAP.clm.logit.flexible[[i]]<-vector()
for (j in 1:dim(pred.test.clm.logit.flexible[[i]])[1])
{pred.MAP.clm.logit.flexible[[i]][j]<-categories[which.max(pred.test.clm.logit.flexible[[i]][j,])]
}
}

pred.MAP.clm.logit.symmetric<-list()
for (i in 1:10)
{pred.MAP.clm.logit.symmetric[[i]]<-vector()
for (j in 1:dim(pred.test.clm.logit.symmetric[[i]])[1])
{pred.MAP.clm.logit.symmetric[[i]][j]<-categories[which.max(pred.test.clm.logit.symmetric[[i]][j,])]
}
}

pred.MAP.clm.logit.symmetric2<-list()
for (i in 1:10)
{pred.MAP.clm.logit.symmetric2[[i]]<-vector()
for (j in 1:dim(pred.test.clm.logit.symmetric2[[i]])[1])
{pred.MAP.clm.logit.symmetric2[[i]][j]<-categories[which.max(pred.test.clm.logit.symmetric2[[i]][j,])]
}
}

pred.MAP.clm.logit.equidistant<-list()
for (i in 1:10)
{pred.MAP.clm.logit.equidistant[[i]]<-vector()
for (j in 1:dim(pred.test.clm.logit.equidistant[[i]])[1])
{pred.MAP.clm.logit.equidistant[[i]][j]<-categories[which.max(pred.test.clm.logit.equidistant[[i]][j,])]
}
}

###


pred.MAP.clm.probit.flexible<-list()
for (i in 1:10)
{pred.MAP.clm.probit.flexible[[i]]<-vector()
for (j in 1:dim(pred.test.clm.probit.flexible[[i]])[1])
{pred.MAP.clm.probit.flexible[[i]][j]<-categories[which.max(pred.test.clm.probit.flexible[[i]][j,])]
}
}

pred.MAP.clm.probit.symmetric<-list()
for (i in 1:10)
{pred.MAP.clm.probit.symmetric[[i]]<-vector()
for (j in 1:dim(pred.test.clm.probit.symmetric[[i]])[1])
{pred.MAP.clm.probit.symmetric[[i]][j]<-categories[which.max(pred.test.clm.probit.symmetric[[i]][j,])]
}
}

pred.MAP.clm.probit.symmetric2<-list()
for (i in 1:10)
{pred.MAP.clm.probit.symmetric2[[i]]<-vector()
for (j in 1:dim(pred.test.clm.probit.symmetric2[[i]])[1])
{pred.MAP.clm.probit.symmetric2[[i]][j]<-categories[which.max(pred.test.clm.probit.symmetric2[[i]][j,])]
}
}

pred.MAP.clm.probit.equidistant<-list()
for (i in 1:10)
{pred.MAP.clm.probit.equidistant[[i]]<-vector()
for (j in 1:dim(pred.test.clm.probit.equidistant[[i]])[1])
{pred.MAP.clm.probit.equidistant[[i]][j]<-categories[which.max(pred.test.clm.probit.equidistant[[i]][j,])]
}
}

###


pred.MAP.clm.loglog.flexible<-list()
for (i in 1:10)
{pred.MAP.clm.loglog.flexible[[i]]<-vector()
for (j in 1:dim(pred.test.clm.loglog.flexible[[i]])[1])
{pred.MAP.clm.loglog.flexible[[i]][j]<-categories[which.max(pred.test.clm.loglog.flexible[[i]][j,])]
}
}

pred.MAP.clm.loglog.symmetric<-list()
for (i in 1:10)
{pred.MAP.clm.loglog.symmetric[[i]]<-vector()
for (j in 1:dim(pred.test.clm.loglog.symmetric[[i]])[1])
{pred.MAP.clm.loglog.symmetric[[i]][j]<-categories[which.max(pred.test.clm.loglog.symmetric[[i]][j,])]
}
}

pred.MAP.clm.loglog.symmetric2<-list()
for (i in 1:10)
{pred.MAP.clm.loglog.symmetric2[[i]]<-vector()
for (j in 1:dim(pred.test.clm.loglog.symmetric2[[i]])[1])
{pred.MAP.clm.loglog.symmetric2[[i]][j]<-categories[which.max(pred.test.clm.loglog.symmetric2[[i]][j,])]
}
}

pred.MAP.clm.loglog.equidistant<-list()
for (i in 1:10)
{pred.MAP.clm.loglog.equidistant[[i]]<-vector()
for (j in 1:dim(pred.test.clm.loglog.equidistant[[i]])[1])
{pred.MAP.clm.loglog.equidistant[[i]][j]<-categories[which.max(pred.test.clm.loglog.equidistant[[i]][j,])]
}
}

###


pred.MAP.clm.cloglog.flexible<-list()
for (i in 1:10)
{pred.MAP.clm.cloglog.flexible[[i]]<-vector()
for (j in 1:dim(pred.test.clm.cloglog.flexible[[i]])[1])
{pred.MAP.clm.cloglog.flexible[[i]][j]<-categories[which.max(pred.test.clm.cloglog.flexible[[i]][j,])]
}
}

pred.MAP.clm.cloglog.symmetric<-list()
for (i in 1:10)
{pred.MAP.clm.cloglog.symmetric[[i]]<-vector()
for (j in 1:dim(pred.test.clm.cloglog.symmetric[[i]])[1])
{pred.MAP.clm.cloglog.symmetric[[i]][j]<-categories[which.max(pred.test.clm.cloglog.symmetric[[i]][j,])]
}
}

pred.MAP.clm.cloglog.symmetric2<-list()
for (i in 1:10)
{pred.MAP.clm.cloglog.symmetric2[[i]]<-vector()
for (j in 1:dim(pred.test.clm.cloglog.symmetric2[[i]])[1])
{pred.MAP.clm.cloglog.symmetric2[[i]][j]<-categories[which.max(pred.test.clm.cloglog.symmetric2[[i]][j,])]
}
}

pred.MAP.clm.cloglog.equidistant<-list()
for (i in 1:10)
{pred.MAP.clm.cloglog.equidistant[[i]]<-vector()
for (j in 1:dim(pred.test.clm.cloglog.equidistant[[i]])[1])
{pred.MAP.clm.cloglog.equidistant[[i]][j]<-categories[which.max(pred.test.clm.cloglog.equidistant[[i]][j,])]
}
}

###


pred.MAP.clm.cauchit.flexible<-list()
for (i in 1:10)
{pred.MAP.clm.cauchit.flexible[[i]]<-vector()
for (j in 1:dim(pred.test.clm.cauchit.flexible[[i]])[1])
{pred.MAP.clm.cauchit.flexible[[i]][j]<-categories[which.max(pred.test.clm.cauchit.flexible[[i]][j,])]
}
}

pred.MAP.clm.cauchit.symmetric<-list()
for (i in 1:10)
{pred.MAP.clm.cauchit.symmetric[[i]]<-vector()
for (j in 1:dim(pred.test.clm.cauchit.symmetric[[i]])[1])
{pred.MAP.clm.cauchit.symmetric[[i]][j]<-categories[which.max(pred.test.clm.cauchit.symmetric[[i]][j,])]
}
}

pred.MAP.clm.cauchit.symmetric2<-list()
for (i in 1:10)
{pred.MAP.clm.cauchit.symmetric2[[i]]<-vector()
for (j in 1:dim(pred.test.clm.cauchit.symmetric2[[i]])[1])
{pred.MAP.clm.cauchit.symmetric2[[i]][j]<-categories[which.max(pred.test.clm.cauchit.symmetric2[[i]][j,])]
}
}

pred.MAP.clm.cauchit.equidistant<-list()
for (i in 1:10)
{pred.MAP.clm.cauchit.equidistant[[i]]<-vector()
for (j in 1:dim(pred.test.clm.cauchit.equidistant[[i]])[1])
{pred.MAP.clm.cauchit.equidistant[[i]][j]<-categories[which.max(pred.test.clm.cauchit.equidistant[[i]][j,])]
}
}

# confusion matrices


conf.mat.MAP.clm.logit.flexible<-list()
for(i in 1:10)
{
  conf.mat.MAP.clm.logit.flexible[[i]]<-mat.square(table(pred.MAP.clm.logit.flexible[[i]],test[[i]]$Class),categories)
}


conf.mat.MAP.clm.logit.symmetric<-list()
for(i in 1:10)
{
  conf.mat.MAP.clm.logit.symmetric[[i]]<-mat.square(table(pred.MAP.clm.logit.symmetric[[i]],test[[i]]$Class),categories)
}


conf.mat.MAP.clm.logit.symmetric2<-list()
for(i in 1:10)
{
  conf.mat.MAP.clm.logit.symmetric2[[i]]<-mat.square(table(pred.MAP.clm.logit.symmetric2[[i]],test[[i]]$Class),categories)
}


conf.mat.MAP.clm.logit.equidistant<-list()
for(i in 1:10)
{
  conf.mat.MAP.clm.logit.equidistant[[i]]<-mat.square(table(pred.MAP.clm.logit.equidistant[[i]],test[[i]]$Class),categories)
}

###

conf.mat.MAP.clm.probit.flexible<-list()
for(i in 1:10)
{
  conf.mat.MAP.clm.probit.flexible[[i]]<-mat.square(table(pred.MAP.clm.probit.flexible[[i]],test[[i]]$Class),categories)
}


conf.mat.MAP.clm.probit.symmetric<-list()
for(i in 1:10)
{
  conf.mat.MAP.clm.probit.symmetric[[i]]<-mat.square(table(pred.MAP.clm.probit.symmetric[[i]],test[[i]]$Class),categories)
}


conf.mat.MAP.clm.probit.symmetric2<-list()
for(i in 1:10)
{
  conf.mat.MAP.clm.probit.symmetric2[[i]]<-mat.square(table(pred.MAP.clm.probit.symmetric2[[i]],test[[i]]$Class),categories)
}


conf.mat.MAP.clm.probit.equidistant<-list()
for(i in 1:10)
{
  conf.mat.MAP.clm.probit.equidistant[[i]]<-mat.square(table(pred.MAP.clm.probit.equidistant[[i]],test[[i]]$Class),categories)
}

###

conf.mat.MAP.clm.loglog.flexible<-list()
for(i in 1:10)
{
  conf.mat.MAP.clm.loglog.flexible[[i]]<-mat.square(table(pred.MAP.clm.loglog.flexible[[i]],test[[i]]$Class),categories)
}


conf.mat.MAP.clm.loglog.symmetric<-list()
for(i in 1:10)
{
  conf.mat.MAP.clm.loglog.symmetric[[i]]<-mat.square(table(pred.MAP.clm.loglog.symmetric[[i]],test[[i]]$Class),categories)
}


conf.mat.MAP.clm.loglog.symmetric2<-list()
for(i in 1:10)
{
  conf.mat.MAP.clm.loglog.symmetric2[[i]]<-mat.square(table(pred.MAP.clm.loglog.symmetric2[[i]],test[[i]]$Class),categories)
}


conf.mat.MAP.clm.loglog.equidistant<-list()
for(i in 1:10)
{
  conf.mat.MAP.clm.loglog.equidistant[[i]]<-mat.square(table(pred.MAP.clm.loglog.equidistant[[i]],test[[i]]$Class),categories)
}

###

conf.mat.MAP.clm.cloglog.flexible<-list()
for(i in 1:10)
{
  conf.mat.MAP.clm.cloglog.flexible[[i]]<-mat.square(table(pred.MAP.clm.cloglog.flexible[[i]],test[[i]]$Class),categories)
}


conf.mat.MAP.clm.cloglog.symmetric<-list()
for(i in 1:10)
{
  conf.mat.MAP.clm.cloglog.symmetric[[i]]<-mat.square(table(pred.MAP.clm.cloglog.symmetric[[i]],test[[i]]$Class),categories)
}


conf.mat.MAP.clm.cloglog.symmetric2<-list()
for(i in 1:10)
{
  conf.mat.MAP.clm.cloglog.symmetric2[[i]]<-mat.square(table(pred.MAP.clm.cloglog.symmetric2[[i]],test[[i]]$Class),categories)
}


conf.mat.MAP.clm.cloglog.equidistant<-list()
for(i in 1:10)
{
  conf.mat.MAP.clm.cloglog.equidistant[[i]]<-mat.square(table(pred.MAP.clm.cloglog.equidistant[[i]],test[[i]]$Class),categories)
}

###

conf.mat.MAP.clm.cauchit.flexible<-list()
for(i in 1:10)
{
  conf.mat.MAP.clm.cauchit.flexible[[i]]<-mat.square(table(pred.MAP.clm.cauchit.flexible[[i]],test[[i]]$Class),categories)
}


conf.mat.MAP.clm.cauchit.symmetric<-list()
for(i in 1:10)
{
  conf.mat.MAP.clm.cauchit.symmetric[[i]]<-mat.square(table(pred.MAP.clm.cauchit.symmetric[[i]],test[[i]]$Class),categories)
}


conf.mat.MAP.clm.cauchit.symmetric2<-list()
for(i in 1:10)
{
  conf.mat.MAP.clm.cauchit.symmetric2[[i]]<-mat.square(table(pred.MAP.clm.cauchit.symmetric2[[i]],test[[i]]$Class),categories)
}


conf.mat.MAP.clm.cauchit.equidistant<-list()
for(i in 1:10)
{
  conf.mat.MAP.clm.cauchit.equidistant[[i]]<-mat.square(table(pred.MAP.clm.cauchit.equidistant[[i]],test[[i]]$Class),categories)
}


# MAE 

MAE.MAP.clm.logit.flexible<-vector()
for (i in 1:10)
{
  MAE.MAP.clm.logit.flexible[i]<-MAE(conf.mat.MAP.clm.logit.flexible[[i]])
}

MAE.MAP.clm.probit.flexible<-vector()
for (i in 1:10)
{
  MAE.MAP.clm.probit.flexible[i]<-MAE(conf.mat.MAP.clm.probit.flexible[[i]])
}

MAE.MAP.clm.loglog.flexible<-vector()
for (i in 1:10)
{
  MAE.MAP.clm.loglog.flexible[i]<-MAE(conf.mat.MAP.clm.loglog.flexible[[i]])
}

MAE.MAP.clm.cloglog.flexible<-vector()
for (i in 1:10)
{
  MAE.MAP.clm.cloglog.flexible[i]<-MAE(conf.mat.MAP.clm.cloglog.flexible[[i]])
}

MAE.MAP.clm.cauchit.flexible<-vector()
for (i in 1:10)
{
  MAE.MAP.clm.cauchit.flexible[i]<-MAE(conf.mat.MAP.clm.cauchit.flexible[[i]])
}

###

MAE.MAP.clm.logit.symmetric<-vector()
for (i in 1:10)
{
  MAE.MAP.clm.logit.symmetric[i]<-MAE(conf.mat.MAP.clm.logit.symmetric[[i]])
}

MAE.MAP.clm.probit.symmetric<-vector()
for (i in 1:10)
{
  MAE.MAP.clm.probit.symmetric[i]<-MAE(conf.mat.MAP.clm.probit.symmetric[[i]])
}

MAE.MAP.clm.loglog.symmetric<-vector()
for (i in 1:10)
{
  MAE.MAP.clm.loglog.symmetric[i]<-MAE(conf.mat.MAP.clm.loglog.symmetric[[i]])
}

MAE.MAP.clm.cloglog.symmetric<-vector()
for (i in 1:10)
{
  MAE.MAP.clm.cloglog.symmetric[i]<-MAE(conf.mat.MAP.clm.cloglog.symmetric[[i]])
}

MAE.MAP.clm.cauchit.symmetric<-vector()
for (i in 1:10)
{
  MAE.MAP.clm.cauchit.symmetric[i]<-MAE(conf.mat.MAP.clm.cauchit.symmetric[[i]])
}

###

MAE.MAP.clm.logit.symmetric2<-vector()
for (i in 1:10)
{
  MAE.MAP.clm.logit.symmetric2[i]<-MAE(conf.mat.MAP.clm.logit.symmetric2[[i]])
}

MAE.MAP.clm.probit.symmetric2<-vector()
for (i in 1:10)
{
  MAE.MAP.clm.probit.symmetric2[i]<-MAE(conf.mat.MAP.clm.probit.symmetric2[[i]])
}

MAE.MAP.clm.loglog.symmetric2<-vector()
for (i in 1:10)
{
  MAE.MAP.clm.loglog.symmetric2[i]<-MAE(conf.mat.MAP.clm.loglog.symmetric2[[i]])
}

MAE.MAP.clm.cloglog.symmetric2<-vector()
for (i in 1:10)
{
  MAE.MAP.clm.cloglog.symmetric2[i]<-MAE(conf.mat.MAP.clm.cloglog.symmetric2[[i]])
}

MAE.MAP.clm.cauchit.symmetric2<-vector()
for (i in 1:10)
{
  MAE.MAP.clm.cauchit.symmetric2[i]<-MAE(conf.mat.MAP.clm.cauchit.symmetric2[[i]])
}

###

MAE.MAP.clm.logit.equidistant<-vector()
for (i in 1:10)
{
  MAE.MAP.clm.logit.equidistant[i]<-MAE(conf.mat.MAP.clm.logit.equidistant[[i]])
}

MAE.MAP.clm.probit.equidistant<-vector()
for (i in 1:10)
{
  MAE.MAP.clm.probit.equidistant[i]<-MAE(conf.mat.MAP.clm.probit.equidistant[[i]])
}

MAE.MAP.clm.loglog.equidistant<-vector()
for (i in 1:10)
{
  MAE.MAP.clm.loglog.equidistant[i]<-MAE(conf.mat.MAP.clm.loglog.equidistant[[i]])
}

MAE.MAP.clm.cloglog.equidistant<-vector()
for (i in 1:10)
{
  MAE.MAP.clm.cloglog.equidistant[i]<-MAE(conf.mat.MAP.clm.cloglog.equidistant[[i]])
}

MAE.MAP.clm.cauchit.equidistant<-vector()
for (i in 1:10)
{
  MAE.MAP.clm.cauchit.equidistant[i]<-MAE(conf.mat.MAP.clm.cauchit.equidistant[[i]])
}



################################################################################
### predictions with Ord.MAP criterion:
###

pred.Ord.MAP.clm.logit.flexible<-list()
for (i in 1:10)
{ pred.Ord.MAP.clm.logit.flexible[[i]]<-vector()
for (j in 1:dim(pred.test.clm.logit.flexible[[i]])[1])
{h<-min(which(cumsum(pred.test.clm.logit.flexible[[i]][j,])>=0.5))
pred.Ord.MAP.clm.logit.flexible[[i]][j]<-categories[h]
}
}

pred.Ord.MAP.clm.logit.symmetric<-list()
for (i in 1:10)
{ pred.Ord.MAP.clm.logit.symmetric[[i]]<-vector()
for (j in 1:dim(pred.test.clm.logit.symmetric[[i]])[1])
{h<-min(which(cumsum(pred.test.clm.logit.symmetric[[i]][j,])>=0.5))
pred.Ord.MAP.clm.logit.symmetric[[i]][j]<-categories[h]
}
}

pred.Ord.MAP.clm.logit.symmetric2<-list()
for (i in 1:10)
{ pred.Ord.MAP.clm.logit.symmetric2[[i]]<-vector()
for (j in 1:dim(pred.test.clm.logit.symmetric2[[i]])[1])
{h<-min(which(cumsum(pred.test.clm.logit.symmetric2[[i]][j,])>=0.5))
pred.Ord.MAP.clm.logit.symmetric2[[i]][j]<-categories[h]
}
}

pred.Ord.MAP.clm.logit.equidistant<-list()
for (i in 1:10)
{ pred.Ord.MAP.clm.logit.equidistant[[i]]<-vector()
for (j in 1:dim(pred.test.clm.logit.equidistant[[i]])[1])
{h<-min(which(cumsum(pred.test.clm.logit.equidistant[[i]][j,])>=0.5))
pred.Ord.MAP.clm.logit.equidistant[[i]][j]<-categories[h]
}
}

###

pred.Ord.MAP.clm.probit.flexible<-list()
for (i in 1:10)
{ pred.Ord.MAP.clm.probit.flexible[[i]]<-vector()
for (j in 1:dim(pred.test.clm.probit.flexible[[i]])[1])
{h<-min(which(cumsum(pred.test.clm.probit.flexible[[i]][j,])>=0.5))
pred.Ord.MAP.clm.probit.flexible[[i]][j]<-categories[h]
}
}

pred.Ord.MAP.clm.probit.symmetric<-list()
for (i in 1:10)
{ pred.Ord.MAP.clm.probit.symmetric[[i]]<-vector()
for (j in 1:dim(pred.test.clm.probit.symmetric[[i]])[1])
{h<-min(which(cumsum(pred.test.clm.probit.symmetric[[i]][j,])>=0.5))
pred.Ord.MAP.clm.probit.symmetric[[i]][j]<-categories[h]
}
}

pred.Ord.MAP.clm.probit.symmetric2<-list()
for (i in 1:10)
{ pred.Ord.MAP.clm.probit.symmetric2[[i]]<-vector()
for (j in 1:dim(pred.test.clm.probit.symmetric2[[i]])[1])
{h<-min(which(cumsum(pred.test.clm.probit.symmetric2[[i]][j,])>=0.5))
pred.Ord.MAP.clm.probit.symmetric2[[i]][j]<-categories[h]
}
}

pred.Ord.MAP.clm.probit.equidistant<-list()
for (i in 1:10)
{ pred.Ord.MAP.clm.probit.equidistant[[i]]<-vector()
for (j in 1:dim(pred.test.clm.probit.equidistant[[i]])[1])
{h<-min(which(cumsum(pred.test.clm.probit.equidistant[[i]][j,])>=0.5))
pred.Ord.MAP.clm.probit.equidistant[[i]][j]<-categories[h]
}
}

###

pred.Ord.MAP.clm.loglog.flexible<-list()
for (i in 1:10)
{ pred.Ord.MAP.clm.loglog.flexible[[i]]<-vector()
for (j in 1:dim(pred.test.clm.loglog.flexible[[i]])[1])
{h<-min(which(cumsum(pred.test.clm.loglog.flexible[[i]][j,])>=0.5))
pred.Ord.MAP.clm.loglog.flexible[[i]][j]<-categories[h]
}
}

pred.Ord.MAP.clm.loglog.symmetric<-list()
for (i in 1:10)
{ pred.Ord.MAP.clm.loglog.symmetric[[i]]<-vector()
for (j in 1:dim(pred.test.clm.loglog.symmetric[[i]])[1])
{h<-min(which(cumsum(pred.test.clm.loglog.symmetric[[i]][j,])>=0.5))
pred.Ord.MAP.clm.loglog.symmetric[[i]][j]<-categories[h]
}
}

pred.Ord.MAP.clm.loglog.symmetric2<-list()
for (i in 1:10)
{ pred.Ord.MAP.clm.loglog.symmetric2[[i]]<-vector()
for (j in 1:dim(pred.test.clm.loglog.symmetric2[[i]])[1])
{h<-min(which(cumsum(pred.test.clm.loglog.symmetric2[[i]][j,])>=0.5))
pred.Ord.MAP.clm.loglog.symmetric2[[i]][j]<-categories[h]
}
}

pred.Ord.MAP.clm.loglog.equidistant<-list()
for (i in 1:10)
{ pred.Ord.MAP.clm.loglog.equidistant[[i]]<-vector()
for (j in 1:dim(pred.test.clm.loglog.equidistant[[i]])[1])
{h<-min(which(cumsum(pred.test.clm.loglog.equidistant[[i]][j,])>=0.5))
pred.Ord.MAP.clm.loglog.equidistant[[i]][j]<-categories[h]
}
}

###

pred.Ord.MAP.clm.cloglog.flexible<-list()
for (i in 1:10)
{ pred.Ord.MAP.clm.cloglog.flexible[[i]]<-vector()
for (j in 1:dim(pred.test.clm.cloglog.flexible[[i]])[1])
{h<-min(which(cumsum(pred.test.clm.cloglog.flexible[[i]][j,])>=0.5))
pred.Ord.MAP.clm.cloglog.flexible[[i]][j]<-categories[h]
}
}

pred.Ord.MAP.clm.cloglog.symmetric<-list()
for (i in 1:10)
{ pred.Ord.MAP.clm.cloglog.symmetric[[i]]<-vector()
for (j in 1:dim(pred.test.clm.cloglog.symmetric[[i]])[1])
{h<-min(which(cumsum(pred.test.clm.cloglog.symmetric[[i]][j,])>=0.5))
pred.Ord.MAP.clm.cloglog.symmetric[[i]][j]<-categories[h]
}
}

pred.Ord.MAP.clm.cloglog.symmetric2<-list()
for (i in 1:10)
{ pred.Ord.MAP.clm.cloglog.symmetric2[[i]]<-vector()
for (j in 1:dim(pred.test.clm.cloglog.symmetric2[[i]])[1])
{h<-min(which(cumsum(pred.test.clm.cloglog.symmetric2[[i]][j,])>=0.5))
pred.Ord.MAP.clm.cloglog.symmetric2[[i]][j]<-categories[h]
}
}

pred.Ord.MAP.clm.cloglog.equidistant<-list()
for (i in 1:10)
{ pred.Ord.MAP.clm.cloglog.equidistant[[i]]<-vector()
for (j in 1:dim(pred.test.clm.cloglog.equidistant[[i]])[1])
{h<-min(which(cumsum(pred.test.clm.cloglog.equidistant[[i]][j,])>=0.5))
pred.Ord.MAP.clm.cloglog.equidistant[[i]][j]<-categories[h]
}
}

###

pred.Ord.MAP.clm.cauchit.flexible<-list()
for (i in 1:10)
{ pred.Ord.MAP.clm.cauchit.flexible[[i]]<-vector()
for (j in 1:dim(pred.test.clm.cauchit.flexible[[i]])[1])
{h<-min(which(cumsum(pred.test.clm.cauchit.flexible[[i]][j,])>=0.5))
pred.Ord.MAP.clm.cauchit.flexible[[i]][j]<-categories[h]
}
}

pred.Ord.MAP.clm.cauchit.symmetric<-list()
for (i in 1:10)
{ pred.Ord.MAP.clm.cauchit.symmetric[[i]]<-vector()
for (j in 1:dim(pred.test.clm.cauchit.symmetric[[i]])[1])
{h<-min(which(cumsum(pred.test.clm.cauchit.symmetric[[i]][j,])>=0.5))
pred.Ord.MAP.clm.cauchit.symmetric[[i]][j]<-categories[h]
}
}

pred.Ord.MAP.clm.cauchit.symmetric2<-list()
for (i in 1:10)
{ pred.Ord.MAP.clm.cauchit.symmetric2[[i]]<-vector()
for (j in 1:dim(pred.test.clm.cauchit.symmetric2[[i]])[1])
{h<-min(which(cumsum(pred.test.clm.cauchit.symmetric2[[i]][j,])>=0.5))
pred.Ord.MAP.clm.cauchit.symmetric2[[i]][j]<-categories[h]
}
}

pred.Ord.MAP.clm.cauchit.equidistant<-list()
for (i in 1:10)
{ pred.Ord.MAP.clm.cauchit.equidistant[[i]]<-vector()
for (j in 1:dim(pred.test.clm.cauchit.equidistant[[i]])[1])
{h<-min(which(cumsum(pred.test.clm.cauchit.equidistant[[i]][j,])>=0.5))
pred.Ord.MAP.clm.cauchit.equidistant[[i]][j]<-categories[h]
}
}

###
# confusion matrices

conf.mat.Ord.MAP.clm.logit.flexible<-list()
for(i in 1:10)
{
  conf.mat.Ord.MAP.clm.logit.flexible[[i]]<-mat.square(table(pred.Ord.MAP.clm.logit.flexible[[i]],test[[i]]$Class),categories)
}

conf.mat.Ord.MAP.clm.logit.symmetric<-list()
for(i in 1:10)
{
  conf.mat.Ord.MAP.clm.logit.symmetric[[i]]<-mat.square(table(pred.Ord.MAP.clm.logit.symmetric[[i]],test[[i]]$Class),categories)
}

conf.mat.Ord.MAP.clm.logit.symmetric2<-list()
for(i in 1:10)
{
  conf.mat.Ord.MAP.clm.logit.symmetric2[[i]]<-mat.square(table(pred.Ord.MAP.clm.logit.symmetric2[[i]],test[[i]]$Class),categories)
}

conf.mat.Ord.MAP.clm.logit.equidistant<-list()
for(i in 1:10)
{
  conf.mat.Ord.MAP.clm.logit.equidistant[[i]]<-mat.square(table(pred.Ord.MAP.clm.logit.equidistant[[i]],test[[i]]$Class),categories)
}

###

conf.mat.Ord.MAP.clm.probit.flexible<-list()
for(i in 1:10)
{
  conf.mat.Ord.MAP.clm.probit.flexible[[i]]<-mat.square(table(pred.Ord.MAP.clm.probit.flexible[[i]],test[[i]]$Class),categories)
}

conf.mat.Ord.MAP.clm.probit.symmetric<-list()
for(i in 1:10)
{
  conf.mat.Ord.MAP.clm.probit.symmetric[[i]]<-mat.square(table(pred.Ord.MAP.clm.probit.symmetric[[i]],test[[i]]$Class),categories)
}

conf.mat.Ord.MAP.clm.probit.symmetric2<-list()
for(i in 1:10)
{
  conf.mat.Ord.MAP.clm.probit.symmetric2[[i]]<-mat.square(table(pred.Ord.MAP.clm.probit.symmetric2[[i]],test[[i]]$Class),categories)
}

conf.mat.Ord.MAP.clm.probit.equidistant<-list()
for(i in 1:10)
{
  conf.mat.Ord.MAP.clm.probit.equidistant[[i]]<-mat.square(table(pred.Ord.MAP.clm.probit.equidistant[[i]],test[[i]]$Class),categories)
}

###

conf.mat.Ord.MAP.clm.loglog.flexible<-list()
for(i in 1:10)
{
  conf.mat.Ord.MAP.clm.loglog.flexible[[i]]<-mat.square(table(pred.Ord.MAP.clm.loglog.flexible[[i]],test[[i]]$Class),categories)
}

conf.mat.Ord.MAP.clm.loglog.symmetric<-list()
for(i in 1:10)
{
  conf.mat.Ord.MAP.clm.loglog.symmetric[[i]]<-mat.square(table(pred.Ord.MAP.clm.loglog.symmetric[[i]],test[[i]]$Class),categories)
}

conf.mat.Ord.MAP.clm.loglog.symmetric2<-list()
for(i in 1:10)
{
  conf.mat.Ord.MAP.clm.loglog.symmetric2[[i]]<-mat.square(table(pred.Ord.MAP.clm.loglog.symmetric2[[i]],test[[i]]$Class),categories)
}

conf.mat.Ord.MAP.clm.loglog.equidistant<-list()
for(i in 1:10)
{
  conf.mat.Ord.MAP.clm.loglog.equidistant[[i]]<-mat.square(table(pred.Ord.MAP.clm.loglog.equidistant[[i]],test[[i]]$Class),categories)
}

###


conf.mat.Ord.MAP.clm.cloglog.flexible<-list()
for(i in 1:10)
{
  conf.mat.Ord.MAP.clm.cloglog.flexible[[i]]<-mat.square(table(pred.Ord.MAP.clm.cloglog.flexible[[i]],test[[i]]$Class),categories)
}

conf.mat.Ord.MAP.clm.cloglog.symmetric<-list()
for(i in 1:10)
{
  conf.mat.Ord.MAP.clm.cloglog.symmetric[[i]]<-mat.square(table(pred.Ord.MAP.clm.cloglog.symmetric[[i]],test[[i]]$Class),categories)
}

conf.mat.Ord.MAP.clm.cloglog.symmetric2<-list()
for(i in 1:10)
{
  conf.mat.Ord.MAP.clm.cloglog.symmetric2[[i]]<-mat.square(table(pred.Ord.MAP.clm.cloglog.symmetric2[[i]],test[[i]]$Class),categories)
}

conf.mat.Ord.MAP.clm.cloglog.equidistant<-list()
for(i in 1:10)
{
  conf.mat.Ord.MAP.clm.cloglog.equidistant[[i]]<-mat.square(table(pred.Ord.MAP.clm.cloglog.equidistant[[i]],test[[i]]$Class),categories)
}

###

conf.mat.Ord.MAP.clm.cauchit.flexible<-list()
for(i in 1:10)
{
  conf.mat.Ord.MAP.clm.cauchit.flexible[[i]]<-mat.square(table(pred.Ord.MAP.clm.cauchit.flexible[[i]],test[[i]]$Class),categories)
}

conf.mat.Ord.MAP.clm.cauchit.symmetric<-list()
for(i in 1:10)
{
  conf.mat.Ord.MAP.clm.cauchit.symmetric[[i]]<-mat.square(table(pred.Ord.MAP.clm.cauchit.symmetric[[i]],test[[i]]$Class),categories)
}

conf.mat.Ord.MAP.clm.cauchit.symmetric2<-list()
for(i in 1:10)
{
  conf.mat.Ord.MAP.clm.cauchit.symmetric2[[i]]<-mat.square(table(pred.Ord.MAP.clm.cauchit.symmetric2[[i]],test[[i]]$Class),categories)
}

conf.mat.Ord.MAP.clm.cauchit.equidistant<-list()
for(i in 1:10)
{
  conf.mat.Ord.MAP.clm.cauchit.equidistant[[i]]<-mat.square(table(pred.Ord.MAP.clm.cauchit.equidistant[[i]],test[[i]]$Class),categories)
}

###

# MAE 

MAE.Ord.MAP.clm.logit.flexible<-vector()
for (i in 1:10)
{
  MAE.Ord.MAP.clm.logit.flexible[i]<-MAE(conf.mat.Ord.MAP.clm.logit.flexible[[i]])
}

MAE.Ord.MAP.clm.logit.symmetric<-vector()
for (i in 1:10)
{
  MAE.Ord.MAP.clm.logit.symmetric[i]<-MAE(conf.mat.Ord.MAP.clm.logit.symmetric[[i]])
}

MAE.Ord.MAP.clm.logit.symmetric2<-vector()
for (i in 1:10)
{
  MAE.Ord.MAP.clm.logit.symmetric2[i]<-MAE(conf.mat.Ord.MAP.clm.logit.symmetric2[[i]])
}

MAE.Ord.MAP.clm.logit.equidistant<-vector()
for (i in 1:10)
{
  MAE.Ord.MAP.clm.logit.equidistant[i]<-MAE(conf.mat.Ord.MAP.clm.logit.equidistant[[i]])
}

###

MAE.Ord.MAP.clm.probit.flexible<-vector()
for (i in 1:10)
{
  MAE.Ord.MAP.clm.probit.flexible[i]<-MAE(conf.mat.Ord.MAP.clm.probit.flexible[[i]])
}

MAE.Ord.MAP.clm.probit.symmetric<-vector()
for (i in 1:10)
{
  MAE.Ord.MAP.clm.probit.symmetric[i]<-MAE(conf.mat.Ord.MAP.clm.probit.symmetric[[i]])
}

MAE.Ord.MAP.clm.probit.symmetric2<-vector()
for (i in 1:10)
{
  MAE.Ord.MAP.clm.probit.symmetric2[i]<-MAE(conf.mat.Ord.MAP.clm.probit.symmetric2[[i]])
}

MAE.Ord.MAP.clm.probit.equidistant<-vector()
for (i in 1:10)
{
  MAE.Ord.MAP.clm.probit.equidistant[i]<-MAE(conf.mat.Ord.MAP.clm.probit.equidistant[[i]])
}

###


MAE.Ord.MAP.clm.loglog.flexible<-vector()
for (i in 1:10)
{
  MAE.Ord.MAP.clm.loglog.flexible[i]<-MAE(conf.mat.Ord.MAP.clm.loglog.flexible[[i]])
}

MAE.Ord.MAP.clm.loglog.symmetric<-vector()
for (i in 1:10)
{
  MAE.Ord.MAP.clm.loglog.symmetric[i]<-MAE(conf.mat.Ord.MAP.clm.loglog.symmetric[[i]])
}

MAE.Ord.MAP.clm.loglog.symmetric2<-vector()
for (i in 1:10)
{
  MAE.Ord.MAP.clm.loglog.symmetric2[i]<-MAE(conf.mat.Ord.MAP.clm.loglog.symmetric2[[i]])
}

MAE.Ord.MAP.clm.loglog.equidistant<-vector()
for (i in 1:10)
{
  MAE.Ord.MAP.clm.loglog.equidistant[i]<-MAE(conf.mat.Ord.MAP.clm.loglog.equidistant[[i]])
}

###


MAE.Ord.MAP.clm.cloglog.flexible<-vector()
for (i in 1:10)
{
  MAE.Ord.MAP.clm.cloglog.flexible[i]<-MAE(conf.mat.Ord.MAP.clm.cloglog.flexible[[i]])
}

MAE.Ord.MAP.clm.cloglog.symmetric<-vector()
for (i in 1:10)
{
  MAE.Ord.MAP.clm.cloglog.symmetric[i]<-MAE(conf.mat.Ord.MAP.clm.cloglog.symmetric[[i]])
}

MAE.Ord.MAP.clm.cloglog.symmetric2<-vector()
for (i in 1:10)
{
  MAE.Ord.MAP.clm.cloglog.symmetric2[i]<-MAE(conf.mat.Ord.MAP.clm.cloglog.symmetric2[[i]])
}

MAE.Ord.MAP.clm.cloglog.equidistant<-vector()
for (i in 1:10)
{
  MAE.Ord.MAP.clm.cloglog.equidistant[i]<-MAE(conf.mat.Ord.MAP.clm.cloglog.equidistant[[i]])
}

###


MAE.Ord.MAP.clm.cauchit.flexible<-vector()
for (i in 1:10)
{
  MAE.Ord.MAP.clm.cauchit.flexible[i]<-MAE(conf.mat.Ord.MAP.clm.cauchit.flexible[[i]])
}

MAE.Ord.MAP.clm.cauchit.symmetric<-vector()
for (i in 1:10)
{
  MAE.Ord.MAP.clm.cauchit.symmetric[i]<-MAE(conf.mat.Ord.MAP.clm.cauchit.symmetric[[i]])
}

MAE.Ord.MAP.clm.cauchit.symmetric2<-vector()
for (i in 1:10)
{
  MAE.Ord.MAP.clm.cauchit.symmetric2[i]<-MAE(conf.mat.Ord.MAP.clm.cauchit.symmetric2[[i]])
}

MAE.Ord.MAP.clm.cauchit.equidistant<-vector()
for (i in 1:10)
{
  MAE.Ord.MAP.clm.cauchit.equidistant[i]<-MAE(conf.mat.Ord.MAP.clm.cauchit.equidistant[[i]])
}


###
################################################################################
###. Comparison MAP vs Ord.MAP
###
###.  Tests MAE

p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.logit.flexible<-wilcox.test(MAE.MAP.clm.logit.flexible,MAE.Ord.MAP.clm.logit.flexible,paired = TRUE,alternative="greater")$p.value
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.logit.symmetric<-wilcox.test(MAE.MAP.clm.logit.symmetric,MAE.Ord.MAP.clm.logit.symmetric,paired = TRUE,alternative="greater")$p.value
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.logit.symmetric2<-wilcox.test(MAE.MAP.clm.logit.symmetric2,MAE.Ord.MAP.clm.logit.symmetric2,paired = TRUE,alternative="greater")$p.value
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.logit.equidistant<-wilcox.test(MAE.MAP.clm.logit.equidistant,MAE.Ord.MAP.clm.logit.equidistant,paired = TRUE,alternative="greater")$p.value

p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.logit.flexible.less<-wilcox.test(MAE.MAP.clm.logit.flexible,MAE.Ord.MAP.clm.logit.flexible,paired = TRUE,alternative="less")$p.value
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.logit.symmetric.less<-wilcox.test(MAE.MAP.clm.logit.symmetric,MAE.Ord.MAP.clm.logit.symmetric,paired = TRUE,alternative="less")$p.value
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.logit.symmetric2.less<-wilcox.test(MAE.MAP.clm.logit.symmetric2,MAE.Ord.MAP.clm.logit.symmetric2,paired = TRUE,alternative="less")$p.value
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.logit.equidistant.less<-wilcox.test(MAE.MAP.clm.logit.equidistant,MAE.Ord.MAP.clm.logit.equidistant,paired = TRUE,alternative="less")$p.value

###

p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.probit.flexible<-wilcox.test(MAE.MAP.clm.probit.flexible,MAE.Ord.MAP.clm.probit.flexible,paired = TRUE,alternative="greater")$p.value
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.probit.symmetric<-wilcox.test(MAE.MAP.clm.probit.symmetric,MAE.Ord.MAP.clm.probit.symmetric,paired = TRUE,alternative="greater")$p.value
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.probit.symmetric2<-wilcox.test(MAE.MAP.clm.probit.symmetric2,MAE.Ord.MAP.clm.probit.symmetric2,paired = TRUE,alternative="greater")$p.value
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.probit.equidistant<-wilcox.test(MAE.MAP.clm.probit.equidistant,MAE.Ord.MAP.clm.probit.equidistant,paired = TRUE,alternative="greater")$p.value

p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.probit.flexible.less<-wilcox.test(MAE.MAP.clm.probit.flexible,MAE.Ord.MAP.clm.probit.flexible,paired = TRUE,alternative="less")$p.value
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.probit.symmetric.less<-wilcox.test(MAE.MAP.clm.probit.symmetric,MAE.Ord.MAP.clm.probit.symmetric,paired = TRUE,alternative="less")$p.value
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.probit.symmetric2.less<-wilcox.test(MAE.MAP.clm.probit.symmetric2,MAE.Ord.MAP.clm.probit.symmetric2,paired = TRUE,alternative="less")$p.value
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.probit.equidistant.less<-wilcox.test(MAE.MAP.clm.probit.equidistant,MAE.Ord.MAP.clm.probit.equidistant,paired = TRUE,alternative="less")$p.value


###

p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.loglog.flexible<-wilcox.test(MAE.MAP.clm.loglog.flexible,MAE.Ord.MAP.clm.loglog.flexible,paired = TRUE,alternative="greater")$p.value
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.loglog.symmetric<-wilcox.test(MAE.MAP.clm.loglog.symmetric,MAE.Ord.MAP.clm.loglog.symmetric,paired = TRUE,alternative="greater")$p.value
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.loglog.symmetric2<-wilcox.test(MAE.MAP.clm.loglog.symmetric2,MAE.Ord.MAP.clm.loglog.symmetric2,paired = TRUE,alternative="greater")$p.value
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.loglog.equidistant<-wilcox.test(MAE.MAP.clm.loglog.equidistant,MAE.Ord.MAP.clm.loglog.equidistant,paired = TRUE,alternative="greater")$p.value

p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.loglog.flexible.less<-wilcox.test(MAE.MAP.clm.loglog.flexible,MAE.Ord.MAP.clm.loglog.flexible,paired = TRUE,alternative="less")$p.value
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.loglog.symmetric.less<-wilcox.test(MAE.MAP.clm.loglog.symmetric,MAE.Ord.MAP.clm.loglog.symmetric,paired = TRUE,alternative="less")$p.value
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.loglog.symmetric2.less<-wilcox.test(MAE.MAP.clm.loglog.symmetric2,MAE.Ord.MAP.clm.loglog.symmetric2,paired = TRUE,alternative="less")$p.value
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.loglog.equidistant.less<-wilcox.test(MAE.MAP.clm.loglog.equidistant,MAE.Ord.MAP.clm.loglog.equidistant,paired = TRUE,alternative="less")$p.value

###

p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.cloglog.flexible<-wilcox.test(MAE.MAP.clm.cloglog.flexible,MAE.Ord.MAP.clm.cloglog.flexible,paired = TRUE,alternative="greater")$p.value
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.cloglog.symmetric<-wilcox.test(MAE.MAP.clm.cloglog.symmetric,MAE.Ord.MAP.clm.cloglog.symmetric,paired = TRUE,alternative="greater")$p.value
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.cloglog.symmetric2<-wilcox.test(MAE.MAP.clm.cloglog.symmetric2,MAE.Ord.MAP.clm.cloglog.symmetric2,paired = TRUE,alternative="greater")$p.value
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.cloglog.equidistant<-wilcox.test(MAE.MAP.clm.cloglog.equidistant,MAE.Ord.MAP.clm.cloglog.equidistant,paired = TRUE,alternative="greater")$p.value

p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.cloglog.flexible.less<-wilcox.test(MAE.MAP.clm.cloglog.flexible,MAE.Ord.MAP.clm.cloglog.flexible,paired = TRUE,alternative="less")$p.value
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.cloglog.symmetric.less<-wilcox.test(MAE.MAP.clm.cloglog.symmetric,MAE.Ord.MAP.clm.cloglog.symmetric,paired = TRUE,alternative="less")$p.value
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.cloglog.symmetric2.less<-wilcox.test(MAE.MAP.clm.cloglog.symmetric2,MAE.Ord.MAP.clm.cloglog.symmetric2,paired = TRUE,alternative="less")$p.value
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.cloglog.equidistant.less<-wilcox.test(MAE.MAP.clm.cloglog.equidistant,MAE.Ord.MAP.clm.cloglog.equidistant,paired = TRUE,alternative="les")$p.value

###

p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.cauchit.flexible<-wilcox.test(MAE.MAP.clm.cauchit.flexible,MAE.Ord.MAP.clm.cauchit.flexible,paired = TRUE,alternative="greater")$p.value
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.cauchit.symmetric<-wilcox.test(MAE.MAP.clm.cauchit.symmetric,MAE.Ord.MAP.clm.cauchit.symmetric,paired = TRUE,alternative="greater")$p.value
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.cauchit.symmetric2<-wilcox.test(MAE.MAP.clm.cauchit.symmetric2,MAE.Ord.MAP.clm.cauchit.symmetric2,paired = TRUE,alternative="greater")$p.value
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.cauchit.equidistant<-wilcox.test(MAE.MAP.clm.cauchit.equidistant,MAE.Ord.MAP.clm.cauchit.equidistant,paired = TRUE,alternative="greater")$p.value

p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.cauchit.flexible.less<-wilcox.test(MAE.MAP.clm.cauchit.flexible,MAE.Ord.MAP.clm.cauchit.flexible,paired = TRUE,alternative="less")$p.value
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.cauchit.symmetric.less<-wilcox.test(MAE.MAP.clm.cauchit.symmetric,MAE.Ord.MAP.clm.cauchit.symmetric,paired = TRUE,alternative="less")$p.value
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.cauchit.symmetric2.less<-wilcox.test(MAE.MAP.clm.cauchit.symmetric2,MAE.Ord.MAP.clm.cauchit.symmetric2,paired = TRUE,alternative="less")$p.value
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.cauchit.equidistant.less<-wilcox.test(MAE.MAP.clm.cauchit.equidistant,MAE.Ord.MAP.clm.cauchit.equidistant,paired = TRUE,alternative="less")$p.value

###

p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.logit.flexible
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.logit.symmetric
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.logit.symmetric2
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.logit.equidistant

p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.logit.flexible.less
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.logit.symmetric.less
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.logit.symmetric2.less
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.logit.equidistant.less

# 

p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.probit.flexible
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.probit.symmetric
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.probit.symmetric2
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.probit.equidistant

p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.probit.flexible.less
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.probit.symmetric.less
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.probit.symmetric2.less
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.probit.equidistant.less

#

p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.loglog.flexible
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.loglog.symmetric
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.loglog.symmetric2
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.loglog.equidistant

p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.loglog.flexible.less
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.loglog.symmetric.less
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.loglog.symmetric2.less
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.loglog.equidistant.less

#

p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.cloglog.flexible
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.cloglog.symmetric
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.cloglog.symmetric2
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.cloglog.equidistant

p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.cloglog.flexible.less
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.cloglog.symmetric.less
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.cloglog.symmetric2.less
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.cloglog.equidistant.less

#

p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.cauchit.flexible
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.cauchit.symmetric
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.cauchit.symmetric2
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.cauchit.equidistant

p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.cauchit.flexible.less
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.cauchit.symmetric.less
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.cauchit.symmetric2.less
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.cauchit.equidistant.less


###############################
###############################
###   Boxplots   MAE

MAE.clm.flexible.results<-as.data.frame(cbind(MAE.MAP.clm.logit.flexible-MAE.Ord.MAP.clm.logit.flexible,
                                              MAE.MAP.clm.probit.flexible-MAE.Ord.MAP.clm.probit.flexible,
                                              MAE.MAP.clm.loglog.flexible-MAE.Ord.MAP.clm.loglog.flexible,
                                              MAE.MAP.clm.cloglog.flexible-MAE.Ord.MAP.clm.cloglog.flexible,
                                              MAE.MAP.clm.cauchit.flexible-MAE.Ord.MAP.clm.cauchit.flexible))

boxplot(MAE.clm.flexible.results,
        main=paste("(c) Hearth dataset. Threshold flexible"), 
        xlab=" ", ylab="MAE increment (MAP - Ord.MAP)",names=c("logistic", "probit","loglog","cloglog","cauchit"),
        col="white")
abline(h = 0, col = 2, lwd = 2, lty=2)
# Dots
stripchart(MAE.clm.flexible.results,
           #method = "jitter",
           pch = 19,
           col = 1:5,
           vertical = TRUE,
           add = TRUE)

###

MAE.clm.symmetric.results<-as.data.frame(cbind(MAE.MAP.clm.logit.symmetric-MAE.Ord.MAP.clm.logit.symmetric,
                                               MAE.MAP.clm.probit.symmetric-MAE.Ord.MAP.clm.probit.symmetric,
                                               MAE.MAP.clm.loglog.symmetric-MAE.Ord.MAP.clm.loglog.symmetric,
                                               MAE.MAP.clm.cloglog.symmetric-MAE.Ord.MAP.clm.cloglog.symmetric,
                                               MAE.MAP.clm.cauchit.symmetric-MAE.Ord.MAP.clm.cauchit.symmetric))

boxplot(MAE.clm.symmetric.results,
        main=paste("(c) Hearth dataset. Threshold symmetric"), 
        xlab=" ", ylab="MAE increment (MAP - Ord.MAP)",names=c("logistic", "probit","loglog","cloglog","cauchit"),
        col="white")
abline(h = 0, col = 2, lwd = 2, lty=2)
# Dots
stripchart(MAE.clm.symmetric.results,
           #method = "jitter",
           pch = 19,
           col = 1:5,
           vertical = TRUE,
           add = TRUE)

###

MAE.clm.symmetric2.results<-as.data.frame(cbind(MAE.MAP.clm.logit.symmetric2-MAE.Ord.MAP.clm.logit.symmetric2,
                                                MAE.MAP.clm.probit.symmetric2-MAE.Ord.MAP.clm.probit.symmetric2,
                                                MAE.MAP.clm.loglog.symmetric2-MAE.Ord.MAP.clm.loglog.symmetric2,
                                                MAE.MAP.clm.cloglog.symmetric2-MAE.Ord.MAP.clm.cloglog.symmetric2,
                                                MAE.MAP.clm.cauchit.symmetric2-MAE.Ord.MAP.clm.cauchit.symmetric2))

boxplot(MAE.clm.symmetric2.results,
        main=paste("(c) Hearth dataset. Threshold symmetric2"), 
        xlab=" ", ylab="MAE increment (MAP - Ord.MAP)",names=c("logistic", "probit","loglog","cloglog","cauchit"),
        col="white")
abline(h = 0, col = 2, lwd = 2, lty=2)
# Dots
stripchart(MAE.clm.symmetric2.results,
           #method = "jitter",
           pch = 19,
           col = 1:5,
           vertical = TRUE,
           add = TRUE)

###

MAE.clm.equidistant.results<-as.data.frame(cbind(MAE.MAP.clm.logit.equidistant-MAE.Ord.MAP.clm.logit.equidistant,
                                                 MAE.MAP.clm.probit.equidistant-MAE.Ord.MAP.clm.probit.equidistant,
                                                 MAE.MAP.clm.loglog.equidistant-MAE.Ord.MAP.clm.loglog.equidistant,
                                                 MAE.MAP.clm.cloglog.equidistant-MAE.Ord.MAP.clm.cloglog.equidistant,
                                                 MAE.MAP.clm.cauchit.equidistant-MAE.Ord.MAP.clm.cauchit.equidistant))

boxplot(MAE.clm.equidistant.results,
        main=paste("(c) Hearth dataset. Threshold equidistant"), 
        xlab=" ", ylab="MAE increment (MAP - Ord.MAP)",names=c("logistic", "probit","loglog","cloglog","cauchit"),
        col="white")
abline(h = 0, col = 2, lwd = 2, lty=2)
# Dots
stripchart(MAE.clm.equidistant.results,
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
# clm, from "ordinal" package, rfits cumulative link models (CLMs) such as the propotional odds model. 
# The model allows for various link functions and structured thresholds that restricts 
# the thresholds or cut-points to be e.g.,
# equidistant or symmetrically arranged around the central threshold(s). 
#
# clm(formula, scale, nominal, data, weights, start, subset, doFit = TRUE,
#     na.action, contrasts, model = TRUE, control=list(),
#     link = c("logit", "probit", "cloglog", "loglog", "cauchit",
#              "Aranda-Ordaz", "log-gamma"),
#     threshold = c("flexible", "symmetric", "symmetric2", "equidistant"), ...)

model.clm.logit.flexible<-list()
for (i in 1:10)
{model.clm.logit.flexible[[i]] <- clm(V5.bin.num.factor ~ V2 + V3 + V4 + V7 + V8 + V9 + V10 +
                                        V11 + V12 + V13 + V14 + V15 + V16 + V17 + V18 + V19 + 
                                        V20 + V21 +V22, data = training[[i]],
                                      link="logit", threshold="flexible")}


model.clm.logit.symmetric<-list()
for (i in 1:10)
{model.clm.logit.symmetric[[i]] <- clm(V5.bin.num.factor ~ V2 + V3 + V4 + V7 + V8 + V9 + V10 +
                                         V11 + V12 + V13 + V14 + V15 + V16 + V17 + V18 + V19 + 
                                         V20 + V21 +V22, data = training[[i]],
                                       link="logit", threshold="symmetric")}


model.clm.logit.symmetric2<-list()
for (i in 1:10)
{model.clm.logit.symmetric2[[i]] <- clm(V5.bin.num.factor ~ V2 + V3 + V4 + V7 + V8 + V9 + V10 +
                                          V11 + V12 + V13 + V14 + V15 + V16 + V17 + V18 + V19 + 
                                          V20 + V21 +V22, data = training[[i]],
                                        link="logit", threshold="symmetric2")}

model.clm.logit.equidistant<-list()
for (i in 1:10)
{model.clm.logit.equidistant[[i]] <- clm(V5.bin.num.factor ~ V2 + V3 + V4 + V7 + V8 + V9 + V10 +
                                           V11 + V12 + V13 + V14 + V15 + V16 + V17 + V18 + V19 + 
                                           V20 + V21 +V22, data = training[[i]],
                                         link="logit", threshold="equidistant")}


###

model.clm.probit.flexible<-list()
for (i in 1:10)
{model.clm.probit.flexible[[i]] <- clm(V5.bin.num.factor ~ V2 + V3 + V4 + V7 + V8 + V9 + V10 +
                                         V11 + V12 + V13 + V14 + V15 + V16 + V17 + V18 + V19 + 
                                         V20 + V21 +V22, data = training[[i]],
                                       link="probit", threshold="flexible")}


model.clm.probit.symmetric<-list()
for (i in 1:10)
{model.clm.probit.symmetric[[i]] <- clm(V5.bin.num.factor ~ V2 + V3 + V4 + V7 + V8 + V9 + V10 +
                                          V11 + V12 + V13 + V14 + V15 + V16 + V17 + V18 + V19 + 
                                          V20 + V21 +V22, data = training[[i]],
                                        link="probit", threshold="symmetric")}


model.clm.probit.symmetric2<-list()
for (i in 1:10)
{model.clm.probit.symmetric2[[i]] <- clm(V5.bin.num.factor ~ V2 + V3 + V4 + V7 + V8 + V9 + V10 +
                                           V11 + V12 + V13 + V14 + V15 + V16 + V17 + V18 + V19 + 
                                           V20 + V21 +V22, data = training[[i]],
                                         link="probit", threshold="symmetric2")}

model.clm.probit.equidistant<-list()
for (i in 1:10)
{model.clm.probit.equidistant[[i]] <- clm(V5.bin.num.factor ~ V2 + V3 + V4 + V7 + V8 + V9 + V10 +
                                            V11 + V12 + V13 + V14 + V15 + V16 + V17 + V18 + V19 + 
                                            V20 + V21 +V22, data = training[[i]],
                                          link="probit", threshold="equidistant")}


###

model.clm.loglog.flexible<-list()
for (i in 1:10)
{model.clm.loglog.flexible[[i]] <- clm(V5.bin.num.factor ~ V2 + V3 + V4 + V7 + V8 + V9 + V10 +
                                         V11 + V12 + V13 + V14 + V15 + V16 + V17 + V18 + V19 + 
                                         V20 + V21 +V22, data = training[[i]],
                                       link="loglog", threshold="flexible")}


model.clm.loglog.symmetric<-list()
for (i in 1:10)
{model.clm.loglog.symmetric[[i]] <- clm(V5.bin.num.factor ~ V2 + V3 + V4 + V7 + V8 + V9 + V10 +
                                          V11 + V12 + V13 + V14 + V15 + V16 + V17 + V18 + V19 + 
                                          V20 + V21 +V22, data = training[[i]],
                                        link="loglog", threshold="symmetric")}


model.clm.loglog.symmetric2<-list()
for (i in 1:10)
{model.clm.loglog.symmetric2[[i]] <- clm(V5.bin.num.factor ~ V2 + V3 + V4 + V7 + V8 + V9 + V10 +
                                           V11 + V12 + V13 + V14 + V15 + V16 + V17 + V18 + V19 + 
                                           V20 + V21 +V22, data = training[[i]],
                                         link="loglog", threshold="symmetric2")}

model.clm.loglog.equidistant<-list()
for (i in 1:10)
{model.clm.loglog.equidistant[[i]] <- clm(V5.bin.num.factor ~ V2 + V3 + V4 + V7 + V8 + V9 + V10 +
                                            V11 + V12 + V13 + V14 + V15 + V16 + V17 + V18 + V19 + 
                                            V20 + V21 +V22, data = training[[i]],
                                          link="loglog", threshold="equidistant")}


###

model.clm.cloglog.flexible<-list()
for (i in 1:10)
{model.clm.cloglog.flexible[[i]] <- clm(V5.bin.num.factor ~ V2 + V3 + V4 + V7 + V8 + V9 + V10 +
                                          V11 + V12 + V13 + V14 + V15 + V16 + V17 + V18 + V19 + 
                                          V20 + V21 +V22, data = training[[i]],
                                        link="cloglog", threshold="flexible")}


model.clm.cloglog.symmetric<-list()
for (i in 1:10)
{model.clm.cloglog.symmetric[[i]] <- clm(V5.bin.num.factor ~ V2 + V3 + V4 + V7 + V8 + V9 + V10 +
                                           V11 + V12 + V13 + V14 + V15 + V16 + V17 + V18 + V19 + 
                                           V20 + V21 +V22, data = training[[i]],
                                         link="cloglog", threshold="symmetric")}


model.clm.cloglog.symmetric2<-list()
for (i in 1:10)
{model.clm.cloglog.symmetric2[[i]] <- clm(V5.bin.num.factor ~ V2 + V3 + V4 + V7 + V8 + V9 + V10 +
                                            V11 + V12 + V13 + V14 + V15 + V16 + V17 + V18 + V19 + 
                                            V20 + V21 +V22, data = training[[i]],
                                          link="cloglog", threshold="symmetric2")}

model.clm.cloglog.equidistant<-list()
for (i in 1:10)
{model.clm.cloglog.equidistant[[i]] <- clm(V5.bin.num.factor ~ V2 + V3 + V4 + V7 + V8 + V9 + V10 +
                                             V11 + V12 + V13 + V14 + V15 + V16 + V17 + V18 + V19 + 
                                             V20 + V21 +V22, data = training[[i]],
                                           link="cloglog", threshold="equidistant")}


###

model.clm.cauchit.flexible<-list()
for (i in 1:10)
{model.clm.cauchit.flexible[[i]] <- clm(V5.bin.num.factor ~ V2 + V3 + V4 + V7 + V8 + V9 + V10 +
                                          V11 + V12 + V13 + V14 + V15 + V16 + V17 + V18 + V19 + 
                                          V20 + V21 +V22, data = training[[i]],
                                        link="cauchit", threshold="flexible")}


model.clm.cauchit.symmetric<-list()
for (i in 1:10)
{model.clm.cauchit.symmetric[[i]] <- clm(V5.bin.num.factor ~ V2 + V3 + V4 + V7 + V8 + V9 + V10 +
                                           V11 + V12 + V13 + V14 + V15 + V16 + V17 + V18 + V19 + 
                                           V20 + V21 +V22, data = training[[i]],
                                         link="cauchit", threshold="symmetric")}


model.clm.cauchit.symmetric2<-list()
for (i in 1:10)
{model.clm.cauchit.symmetric2[[i]] <- clm(V5.bin.num.factor ~ V2 + V3 + V4 + V7 + V8 + V9 + V10 +
                                            V11 + V12 + V13 + V14 + V15 + V16 + V17 + V18 + V19 + 
                                            V20 + V21 +V22, data = training[[i]],
                                          link="cauchit", threshold="symmetric2")}

model.clm.cauchit.equidistant<-list()
for (i in 1:10)
{model.clm.cauchit.equidistant[[i]] <- clm(V5.bin.num.factor ~ V2 + V3 + V4 + V7 + V8 + V9 + V10 +
                                             V11 + V12 + V13 + V14 + V15 + V16 + V17 + V18 + V19 + 
                                             V20 + V21 +V22, data = training[[i]],
                                           link="cauchit", threshold="equidistant")}



################################################################################
#
# Predictions of any of the test set
# 
#

vector.no<-c(1,5,6,23:28)

pred.test.clm.logit.flexible<-list()
for (i in 1:10)
{pred.test.clm.logit.flexible[[i]] <- predict(model.clm.logit.flexible[[i]],test[[i]][,-vector.no],type = "p")$fit
}

pred.test.clm.logit.symmetric<-list()
for (i in 1:10)
{pred.test.clm.logit.symmetric[[i]] <- predict(model.clm.logit.symmetric[[i]],test[[i]][,-vector.no],type = "p")$fit
}

pred.test.clm.logit.symmetric2<-list()
for (i in 1:10)
{pred.test.clm.logit.symmetric2[[i]] <- predict(model.clm.logit.symmetric2[[i]],test[[i]][,-vector.no],type = "p")$fit
}

pred.test.clm.logit.equidistant<-list()
for (i in 1:10)
{pred.test.clm.logit.equidistant[[i]] <- predict(model.clm.logit.equidistant[[i]],test[[i]][,-vector.no],type = "p")$fit
}

###

pred.test.clm.probit.flexible<-list()
for (i in 1:10)
{pred.test.clm.probit.flexible[[i]] <- predict(model.clm.probit.flexible[[i]],test[[i]][,-vector.no],type = "p")$fit
}

pred.test.clm.probit.symmetric<-list()
for (i in 1:10)
{pred.test.clm.probit.symmetric[[i]] <- predict(model.clm.probit.symmetric[[i]],test[[i]][,-vector.no],type = "p")$fit
}

pred.test.clm.probit.symmetric2<-list()
for (i in 1:10)
{pred.test.clm.probit.symmetric2[[i]] <- predict(model.clm.probit.symmetric2[[i]],test[[i]][,-vector.no],type = "p")$fit
}

pred.test.clm.probit.equidistant<-list()
for (i in 1:10)
{pred.test.clm.probit.equidistant[[i]] <- predict(model.clm.probit.equidistant[[i]],test[[i]][,-vector.no],type = "p")$fit
}

###

pred.test.clm.loglog.flexible<-list()
for (i in 1:10)
{pred.test.clm.loglog.flexible[[i]] <- predict(model.clm.loglog.flexible[[i]],test[[i]][,-vector.no],type = "p")$fit
}

pred.test.clm.loglog.symmetric<-list()
for (i in 1:10)
{pred.test.clm.loglog.symmetric[[i]] <- predict(model.clm.loglog.symmetric[[i]],test[[i]][,-vector.no],type = "p")$fit
}

pred.test.clm.loglog.symmetric2<-list()
for (i in 1:10)
{pred.test.clm.loglog.symmetric2[[i]] <- predict(model.clm.loglog.symmetric2[[i]],test[[i]][,-vector.no],type = "p")$fit
}

pred.test.clm.loglog.equidistant<-list()
for (i in 1:10)
{pred.test.clm.loglog.equidistant[[i]] <- predict(model.clm.loglog.equidistant[[i]],test[[i]][,-vector.no],type = "p")$fit
}

###


pred.test.clm.cloglog.flexible<-list()
for (i in 1:10)
{pred.test.clm.cloglog.flexible[[i]] <- predict(model.clm.cloglog.flexible[[i]],test[[i]][,-vector.no],type = "p")$fit
}

pred.test.clm.cloglog.symmetric<-list()
for (i in 1:10)
{pred.test.clm.cloglog.symmetric[[i]] <- predict(model.clm.cloglog.symmetric[[i]],test[[i]][,-vector.no],type = "p")$fit
}

pred.test.clm.cloglog.symmetric2<-list()
for (i in 1:10)
{pred.test.clm.cloglog.symmetric2[[i]] <- predict(model.clm.cloglog.symmetric2[[i]],test[[i]][,-vector.no],type = "p")$fit
}

pred.test.clm.cloglog.equidistant<-list()
for (i in 1:10)
{pred.test.clm.cloglog.equidistant[[i]] <- predict(model.clm.cloglog.equidistant[[i]],test[[i]][,-vector.no],type = "p")$fit
}

###


pred.test.clm.cauchit.flexible<-list()
for (i in 1:10)
{pred.test.clm.cauchit.flexible[[i]] <- predict(model.clm.cauchit.flexible[[i]],test[[i]][,-vector.no],type = "p")$fit
}

pred.test.clm.cauchit.symmetric<-list()
for (i in 1:10)
{pred.test.clm.cauchit.symmetric[[i]] <- predict(model.clm.cauchit.symmetric[[i]],test[[i]][,-vector.no],type = "p")$fit
}

pred.test.clm.cauchit.symmetric2<-list()
for (i in 1:10)
{pred.test.clm.cauchit.symmetric2[[i]] <- predict(model.clm.cauchit.symmetric2[[i]],test[[i]][,-vector.no],type = "p")$fit
}

pred.test.clm.cauchit.equidistant<-list()
for (i in 1:10)
{pred.test.clm.cauchit.equidistant[[i]] <- predict(model.clm.cauchit.equidistant[[i]],test[[i]][,-vector.no],type = "p")$fit
}

###

################################################################################
### predictions with standard MAP criterion:
###

pred.MAP.clm.logit.flexible<-list()
for (i in 1:10)
{pred.MAP.clm.logit.flexible[[i]]<-vector()
for (j in 1:dim(pred.test.clm.logit.flexible[[i]])[1])
{pred.MAP.clm.logit.flexible[[i]][j]<-categories[which.max(pred.test.clm.logit.flexible[[i]][j,])]
}
}

pred.MAP.clm.logit.symmetric<-list()
for (i in 1:10)
{pred.MAP.clm.logit.symmetric[[i]]<-vector()
for (j in 1:dim(pred.test.clm.logit.symmetric[[i]])[1])
{pred.MAP.clm.logit.symmetric[[i]][j]<-categories[which.max(pred.test.clm.logit.symmetric[[i]][j,])]
}
}

pred.MAP.clm.logit.symmetric2<-list()
for (i in 1:10)
{pred.MAP.clm.logit.symmetric2[[i]]<-vector()
for (j in 1:dim(pred.test.clm.logit.symmetric2[[i]])[1])
{pred.MAP.clm.logit.symmetric2[[i]][j]<-categories[which.max(pred.test.clm.logit.symmetric2[[i]][j,])]
}
}

pred.MAP.clm.logit.equidistant<-list()
for (i in 1:10)
{pred.MAP.clm.logit.equidistant[[i]]<-vector()
for (j in 1:dim(pred.test.clm.logit.equidistant[[i]])[1])
{pred.MAP.clm.logit.equidistant[[i]][j]<-categories[which.max(pred.test.clm.logit.equidistant[[i]][j,])]
}
}

###


pred.MAP.clm.probit.flexible<-list()
for (i in 1:10)
{pred.MAP.clm.probit.flexible[[i]]<-vector()
for (j in 1:dim(pred.test.clm.probit.flexible[[i]])[1])
{pred.MAP.clm.probit.flexible[[i]][j]<-categories[which.max(pred.test.clm.probit.flexible[[i]][j,])]
}
}

pred.MAP.clm.probit.symmetric<-list()
for (i in 1:10)
{pred.MAP.clm.probit.symmetric[[i]]<-vector()
for (j in 1:dim(pred.test.clm.probit.symmetric[[i]])[1])
{pred.MAP.clm.probit.symmetric[[i]][j]<-categories[which.max(pred.test.clm.probit.symmetric[[i]][j,])]
}
}

pred.MAP.clm.probit.symmetric2<-list()
for (i in 1:10)
{pred.MAP.clm.probit.symmetric2[[i]]<-vector()
for (j in 1:dim(pred.test.clm.probit.symmetric2[[i]])[1])
{pred.MAP.clm.probit.symmetric2[[i]][j]<-categories[which.max(pred.test.clm.probit.symmetric2[[i]][j,])]
}
}

pred.MAP.clm.probit.equidistant<-list()
for (i in 1:10)
{pred.MAP.clm.probit.equidistant[[i]]<-vector()
for (j in 1:dim(pred.test.clm.probit.equidistant[[i]])[1])
{pred.MAP.clm.probit.equidistant[[i]][j]<-categories[which.max(pred.test.clm.probit.equidistant[[i]][j,])]
}
}

###


pred.MAP.clm.loglog.flexible<-list()
for (i in 1:10)
{pred.MAP.clm.loglog.flexible[[i]]<-vector()
for (j in 1:dim(pred.test.clm.loglog.flexible[[i]])[1])
{pred.MAP.clm.loglog.flexible[[i]][j]<-categories[which.max(pred.test.clm.loglog.flexible[[i]][j,])]
}
}

pred.MAP.clm.loglog.symmetric<-list()
for (i in 1:10)
{pred.MAP.clm.loglog.symmetric[[i]]<-vector()
for (j in 1:dim(pred.test.clm.loglog.symmetric[[i]])[1])
{pred.MAP.clm.loglog.symmetric[[i]][j]<-categories[which.max(pred.test.clm.loglog.symmetric[[i]][j,])]
}
}

pred.MAP.clm.loglog.symmetric2<-list()
for (i in 1:10)
{pred.MAP.clm.loglog.symmetric2[[i]]<-vector()
for (j in 1:dim(pred.test.clm.loglog.symmetric2[[i]])[1])
{pred.MAP.clm.loglog.symmetric2[[i]][j]<-categories[which.max(pred.test.clm.loglog.symmetric2[[i]][j,])]
}
}

pred.MAP.clm.loglog.equidistant<-list()
for (i in 1:10)
{pred.MAP.clm.loglog.equidistant[[i]]<-vector()
for (j in 1:dim(pred.test.clm.loglog.equidistant[[i]])[1])
{pred.MAP.clm.loglog.equidistant[[i]][j]<-categories[which.max(pred.test.clm.loglog.equidistant[[i]][j,])]
}
}

###


pred.MAP.clm.cloglog.flexible<-list()
for (i in 1:10)
{pred.MAP.clm.cloglog.flexible[[i]]<-vector()
for (j in 1:dim(pred.test.clm.cloglog.flexible[[i]])[1])
{pred.MAP.clm.cloglog.flexible[[i]][j]<-categories[which.max(pred.test.clm.cloglog.flexible[[i]][j,])]
}
}

pred.MAP.clm.cloglog.symmetric<-list()
for (i in 1:10)
{pred.MAP.clm.cloglog.symmetric[[i]]<-vector()
for (j in 1:dim(pred.test.clm.cloglog.symmetric[[i]])[1])
{pred.MAP.clm.cloglog.symmetric[[i]][j]<-categories[which.max(pred.test.clm.cloglog.symmetric[[i]][j,])]
}
}

pred.MAP.clm.cloglog.symmetric2<-list()
for (i in 1:10)
{pred.MAP.clm.cloglog.symmetric2[[i]]<-vector()
for (j in 1:dim(pred.test.clm.cloglog.symmetric2[[i]])[1])
{pred.MAP.clm.cloglog.symmetric2[[i]][j]<-categories[which.max(pred.test.clm.cloglog.symmetric2[[i]][j,])]
}
}

pred.MAP.clm.cloglog.equidistant<-list()
for (i in 1:10)
{pred.MAP.clm.cloglog.equidistant[[i]]<-vector()
for (j in 1:dim(pred.test.clm.cloglog.equidistant[[i]])[1])
{pred.MAP.clm.cloglog.equidistant[[i]][j]<-categories[which.max(pred.test.clm.cloglog.equidistant[[i]][j,])]
}
}

###


pred.MAP.clm.cauchit.flexible<-list()
for (i in 1:10)
{pred.MAP.clm.cauchit.flexible[[i]]<-vector()
for (j in 1:dim(pred.test.clm.cauchit.flexible[[i]])[1])
{pred.MAP.clm.cauchit.flexible[[i]][j]<-categories[which.max(pred.test.clm.cauchit.flexible[[i]][j,])]
}
}

pred.MAP.clm.cauchit.symmetric<-list()
for (i in 1:10)
{pred.MAP.clm.cauchit.symmetric[[i]]<-vector()
for (j in 1:dim(pred.test.clm.cauchit.symmetric[[i]])[1])
{pred.MAP.clm.cauchit.symmetric[[i]][j]<-categories[which.max(pred.test.clm.cauchit.symmetric[[i]][j,])]
}
}

pred.MAP.clm.cauchit.symmetric2<-list()
for (i in 1:10)
{pred.MAP.clm.cauchit.symmetric2[[i]]<-vector()
for (j in 1:dim(pred.test.clm.cauchit.symmetric2[[i]])[1])
{pred.MAP.clm.cauchit.symmetric2[[i]][j]<-categories[which.max(pred.test.clm.cauchit.symmetric2[[i]][j,])]
}
}

pred.MAP.clm.cauchit.equidistant<-list()
for (i in 1:10)
{pred.MAP.clm.cauchit.equidistant[[i]]<-vector()
for (j in 1:dim(pred.test.clm.cauchit.equidistant[[i]])[1])
{pred.MAP.clm.cauchit.equidistant[[i]][j]<-categories[which.max(pred.test.clm.cauchit.equidistant[[i]][j,])]
}
}

# confusion matrices


conf.mat.MAP.clm.logit.flexible<-list()
for(i in 1:10)
{
  conf.mat.MAP.clm.logit.flexible[[i]]<-mat.square(table(pred.MAP.clm.logit.flexible[[i]],test[[i]]$V5.bin.num.factor),categories)
}


conf.mat.MAP.clm.logit.symmetric<-list()
for(i in 1:10)
{
  conf.mat.MAP.clm.logit.symmetric[[i]]<-mat.square(table(pred.MAP.clm.logit.symmetric[[i]],test[[i]]$V5.bin.num.factor),categories)
}


conf.mat.MAP.clm.logit.symmetric2<-list()
for(i in 1:10)
{
  conf.mat.MAP.clm.logit.symmetric2[[i]]<-mat.square(table(pred.MAP.clm.logit.symmetric2[[i]],test[[i]]$V5.bin.num.factor),categories)
}


conf.mat.MAP.clm.logit.equidistant<-list()
for(i in 1:10)
{
  conf.mat.MAP.clm.logit.equidistant[[i]]<-mat.square(table(pred.MAP.clm.logit.equidistant[[i]],test[[i]]$V5.bin.num.factor),categories)
}

###

conf.mat.MAP.clm.probit.flexible<-list()
for(i in 1:10)
{
  conf.mat.MAP.clm.probit.flexible[[i]]<-mat.square(table(pred.MAP.clm.probit.flexible[[i]],test[[i]]$V5.bin.num.factor),categories)
}


conf.mat.MAP.clm.probit.symmetric<-list()
for(i in 1:10)
{
  conf.mat.MAP.clm.probit.symmetric[[i]]<-mat.square(table(pred.MAP.clm.probit.symmetric[[i]],test[[i]]$V5.bin.num.factor),categories)
}


conf.mat.MAP.clm.probit.symmetric2<-list()
for(i in 1:10)
{
  conf.mat.MAP.clm.probit.symmetric2[[i]]<-mat.square(table(pred.MAP.clm.probit.symmetric2[[i]],test[[i]]$V5.bin.num.factor),categories)
}


conf.mat.MAP.clm.probit.equidistant<-list()
for(i in 1:10)
{
  conf.mat.MAP.clm.probit.equidistant[[i]]<-mat.square(table(pred.MAP.clm.probit.equidistant[[i]],test[[i]]$V5.bin.num.factor),categories)
}

###

conf.mat.MAP.clm.loglog.flexible<-list()
for(i in 1:10)
{
  conf.mat.MAP.clm.loglog.flexible[[i]]<-mat.square(table(pred.MAP.clm.loglog.flexible[[i]],test[[i]]$V5.bin.num.factor),categories)
}


conf.mat.MAP.clm.loglog.symmetric<-list()
for(i in 1:10)
{
  conf.mat.MAP.clm.loglog.symmetric[[i]]<-mat.square(table(pred.MAP.clm.loglog.symmetric[[i]],test[[i]]$V5.bin.num.factor),categories)
}


conf.mat.MAP.clm.loglog.symmetric2<-list()
for(i in 1:10)
{
  conf.mat.MAP.clm.loglog.symmetric2[[i]]<-mat.square(table(pred.MAP.clm.loglog.symmetric2[[i]],test[[i]]$V5.bin.num.factor),categories)
}


conf.mat.MAP.clm.loglog.equidistant<-list()
for(i in 1:10)
{
  conf.mat.MAP.clm.loglog.equidistant[[i]]<-mat.square(table(pred.MAP.clm.loglog.equidistant[[i]],test[[i]]$V5.bin.num.factor),categories)
}

###

conf.mat.MAP.clm.cloglog.flexible<-list()
for(i in 1:10)
{
  conf.mat.MAP.clm.cloglog.flexible[[i]]<-mat.square(table(pred.MAP.clm.cloglog.flexible[[i]],test[[i]]$V5.bin.num.factor),categories)
}


conf.mat.MAP.clm.cloglog.symmetric<-list()
for(i in 1:10)
{
  conf.mat.MAP.clm.cloglog.symmetric[[i]]<-mat.square(table(pred.MAP.clm.cloglog.symmetric[[i]],test[[i]]$V5.bin.num.factor),categories)
}


conf.mat.MAP.clm.cloglog.symmetric2<-list()
for(i in 1:10)
{
  conf.mat.MAP.clm.cloglog.symmetric2[[i]]<-mat.square(table(pred.MAP.clm.cloglog.symmetric2[[i]],test[[i]]$V5.bin.num.factor),categories)
}


conf.mat.MAP.clm.cloglog.equidistant<-list()
for(i in 1:10)
{
  conf.mat.MAP.clm.cloglog.equidistant[[i]]<-mat.square(table(pred.MAP.clm.cloglog.equidistant[[i]],test[[i]]$V5.bin.num.factor),categories)
}

###

conf.mat.MAP.clm.cauchit.flexible<-list()
for(i in 1:10)
{
  conf.mat.MAP.clm.cauchit.flexible[[i]]<-mat.square(table(pred.MAP.clm.cauchit.flexible[[i]],test[[i]]$V5.bin.num.factor),categories)
}


conf.mat.MAP.clm.cauchit.symmetric<-list()
for(i in 1:10)
{
  conf.mat.MAP.clm.cauchit.symmetric[[i]]<-mat.square(table(pred.MAP.clm.cauchit.symmetric[[i]],test[[i]]$V5.bin.num.factor),categories)
}


conf.mat.MAP.clm.cauchit.symmetric2<-list()
for(i in 1:10)
{
  conf.mat.MAP.clm.cauchit.symmetric2[[i]]<-mat.square(table(pred.MAP.clm.cauchit.symmetric2[[i]],test[[i]]$V5.bin.num.factor),categories)
}


conf.mat.MAP.clm.cauchit.equidistant<-list()
for(i in 1:10)
{
  conf.mat.MAP.clm.cauchit.equidistant[[i]]<-mat.square(table(pred.MAP.clm.cauchit.equidistant[[i]],test[[i]]$V5.bin.num.factor),categories)
}


# MAE 

MAE.MAP.clm.logit.flexible<-vector()
for (i in 1:10)
{
  MAE.MAP.clm.logit.flexible[i]<-MAE(conf.mat.MAP.clm.logit.flexible[[i]])
}

MAE.MAP.clm.probit.flexible<-vector()
for (i in 1:10)
{
  MAE.MAP.clm.probit.flexible[i]<-MAE(conf.mat.MAP.clm.probit.flexible[[i]])
}

MAE.MAP.clm.loglog.flexible<-vector()
for (i in 1:10)
{
  MAE.MAP.clm.loglog.flexible[i]<-MAE(conf.mat.MAP.clm.loglog.flexible[[i]])
}

MAE.MAP.clm.cloglog.flexible<-vector()
for (i in 1:10)
{
  MAE.MAP.clm.cloglog.flexible[i]<-MAE(conf.mat.MAP.clm.cloglog.flexible[[i]])
}

MAE.MAP.clm.cauchit.flexible<-vector()
for (i in 1:10)
{
  MAE.MAP.clm.cauchit.flexible[i]<-MAE(conf.mat.MAP.clm.cauchit.flexible[[i]])
}

###

MAE.MAP.clm.logit.symmetric<-vector()
for (i in 1:10)
{
  MAE.MAP.clm.logit.symmetric[i]<-MAE(conf.mat.MAP.clm.logit.symmetric[[i]])
}

MAE.MAP.clm.probit.symmetric<-vector()
for (i in 1:10)
{
  MAE.MAP.clm.probit.symmetric[i]<-MAE(conf.mat.MAP.clm.probit.symmetric[[i]])
}

MAE.MAP.clm.loglog.symmetric<-vector()
for (i in 1:10)
{
  MAE.MAP.clm.loglog.symmetric[i]<-MAE(conf.mat.MAP.clm.loglog.symmetric[[i]])
}

MAE.MAP.clm.cloglog.symmetric<-vector()
for (i in 1:10)
{
  MAE.MAP.clm.cloglog.symmetric[i]<-MAE(conf.mat.MAP.clm.cloglog.symmetric[[i]])
}

MAE.MAP.clm.cauchit.symmetric<-vector()
for (i in 1:10)
{
  MAE.MAP.clm.cauchit.symmetric[i]<-MAE(conf.mat.MAP.clm.cauchit.symmetric[[i]])
}

###

MAE.MAP.clm.logit.symmetric2<-vector()
for (i in 1:10)
{
  MAE.MAP.clm.logit.symmetric2[i]<-MAE(conf.mat.MAP.clm.logit.symmetric2[[i]])
}

MAE.MAP.clm.probit.symmetric2<-vector()
for (i in 1:10)
{
  MAE.MAP.clm.probit.symmetric2[i]<-MAE(conf.mat.MAP.clm.probit.symmetric2[[i]])
}

MAE.MAP.clm.loglog.symmetric2<-vector()
for (i in 1:10)
{
  MAE.MAP.clm.loglog.symmetric2[i]<-MAE(conf.mat.MAP.clm.loglog.symmetric2[[i]])
}

MAE.MAP.clm.cloglog.symmetric2<-vector()
for (i in 1:10)
{
  MAE.MAP.clm.cloglog.symmetric2[i]<-MAE(conf.mat.MAP.clm.cloglog.symmetric2[[i]])
}

MAE.MAP.clm.cauchit.symmetric2<-vector()
for (i in 1:10)
{
  MAE.MAP.clm.cauchit.symmetric2[i]<-MAE(conf.mat.MAP.clm.cauchit.symmetric2[[i]])
}

###

MAE.MAP.clm.logit.equidistant<-vector()
for (i in 1:10)
{
  MAE.MAP.clm.logit.equidistant[i]<-MAE(conf.mat.MAP.clm.logit.equidistant[[i]])
}

MAE.MAP.clm.probit.equidistant<-vector()
for (i in 1:10)
{
  MAE.MAP.clm.probit.equidistant[i]<-MAE(conf.mat.MAP.clm.probit.equidistant[[i]])
}

MAE.MAP.clm.loglog.equidistant<-vector()
for (i in 1:10)
{
  MAE.MAP.clm.loglog.equidistant[i]<-MAE(conf.mat.MAP.clm.loglog.equidistant[[i]])
}

MAE.MAP.clm.cloglog.equidistant<-vector()
for (i in 1:10)
{
  MAE.MAP.clm.cloglog.equidistant[i]<-MAE(conf.mat.MAP.clm.cloglog.equidistant[[i]])
}

MAE.MAP.clm.cauchit.equidistant<-vector()
for (i in 1:10)
{
  MAE.MAP.clm.cauchit.equidistant[i]<-MAE(conf.mat.MAP.clm.cauchit.equidistant[[i]])
}



################################################################################
### predictions with Ord.MAP criterion:
###

pred.Ord.MAP.clm.logit.flexible<-list()
for (i in 1:10)
{ pred.Ord.MAP.clm.logit.flexible[[i]]<-vector()
for (j in 1:dim(pred.test.clm.logit.flexible[[i]])[1])
{h<-min(which(cumsum(pred.test.clm.logit.flexible[[i]][j,])>=0.5))
pred.Ord.MAP.clm.logit.flexible[[i]][j]<-categories[h]
}
}

pred.Ord.MAP.clm.logit.symmetric<-list()
for (i in 1:10)
{ pred.Ord.MAP.clm.logit.symmetric[[i]]<-vector()
for (j in 1:dim(pred.test.clm.logit.symmetric[[i]])[1])
{h<-min(which(cumsum(pred.test.clm.logit.symmetric[[i]][j,])>=0.5))
pred.Ord.MAP.clm.logit.symmetric[[i]][j]<-categories[h]
}
}

pred.Ord.MAP.clm.logit.symmetric2<-list()
for (i in 1:10)
{ pred.Ord.MAP.clm.logit.symmetric2[[i]]<-vector()
for (j in 1:dim(pred.test.clm.logit.symmetric2[[i]])[1])
{h<-min(which(cumsum(pred.test.clm.logit.symmetric2[[i]][j,])>=0.5))
pred.Ord.MAP.clm.logit.symmetric2[[i]][j]<-categories[h]
}
}

pred.Ord.MAP.clm.logit.equidistant<-list()
for (i in 1:10)
{ pred.Ord.MAP.clm.logit.equidistant[[i]]<-vector()
for (j in 1:dim(pred.test.clm.logit.equidistant[[i]])[1])
{h<-min(which(cumsum(pred.test.clm.logit.equidistant[[i]][j,])>=0.5))
pred.Ord.MAP.clm.logit.equidistant[[i]][j]<-categories[h]
}
}

###

pred.Ord.MAP.clm.probit.flexible<-list()
for (i in 1:10)
{ pred.Ord.MAP.clm.probit.flexible[[i]]<-vector()
for (j in 1:dim(pred.test.clm.probit.flexible[[i]])[1])
{h<-min(which(cumsum(pred.test.clm.probit.flexible[[i]][j,])>=0.5))
pred.Ord.MAP.clm.probit.flexible[[i]][j]<-categories[h]
}
}

pred.Ord.MAP.clm.probit.symmetric<-list()
for (i in 1:10)
{ pred.Ord.MAP.clm.probit.symmetric[[i]]<-vector()
for (j in 1:dim(pred.test.clm.probit.symmetric[[i]])[1])
{h<-min(which(cumsum(pred.test.clm.probit.symmetric[[i]][j,])>=0.5))
pred.Ord.MAP.clm.probit.symmetric[[i]][j]<-categories[h]
}
}

pred.Ord.MAP.clm.probit.symmetric2<-list()
for (i in 1:10)
{ pred.Ord.MAP.clm.probit.symmetric2[[i]]<-vector()
for (j in 1:dim(pred.test.clm.probit.symmetric2[[i]])[1])
{h<-min(which(cumsum(pred.test.clm.probit.symmetric2[[i]][j,])>=0.5))
pred.Ord.MAP.clm.probit.symmetric2[[i]][j]<-categories[h]
}
}

pred.Ord.MAP.clm.probit.equidistant<-list()
for (i in 1:10)
{ pred.Ord.MAP.clm.probit.equidistant[[i]]<-vector()
for (j in 1:dim(pred.test.clm.probit.equidistant[[i]])[1])
{h<-min(which(cumsum(pred.test.clm.probit.equidistant[[i]][j,])>=0.5))
pred.Ord.MAP.clm.probit.equidistant[[i]][j]<-categories[h]
}
}

###

pred.Ord.MAP.clm.loglog.flexible<-list()
for (i in 1:10)
{ pred.Ord.MAP.clm.loglog.flexible[[i]]<-vector()
for (j in 1:dim(pred.test.clm.loglog.flexible[[i]])[1])
{h<-min(which(cumsum(pred.test.clm.loglog.flexible[[i]][j,])>=0.5))
pred.Ord.MAP.clm.loglog.flexible[[i]][j]<-categories[h]
}
}

pred.Ord.MAP.clm.loglog.symmetric<-list()
for (i in 1:10)
{ pred.Ord.MAP.clm.loglog.symmetric[[i]]<-vector()
for (j in 1:dim(pred.test.clm.loglog.symmetric[[i]])[1])
{h<-min(which(cumsum(pred.test.clm.loglog.symmetric[[i]][j,])>=0.5))
pred.Ord.MAP.clm.loglog.symmetric[[i]][j]<-categories[h]
}
}

pred.Ord.MAP.clm.loglog.symmetric2<-list()
for (i in 1:10)
{ pred.Ord.MAP.clm.loglog.symmetric2[[i]]<-vector()
for (j in 1:dim(pred.test.clm.loglog.symmetric2[[i]])[1])
{h<-min(which(cumsum(pred.test.clm.loglog.symmetric2[[i]][j,])>=0.5))
pred.Ord.MAP.clm.loglog.symmetric2[[i]][j]<-categories[h]
}
}

pred.Ord.MAP.clm.loglog.equidistant<-list()
for (i in 1:10)
{ pred.Ord.MAP.clm.loglog.equidistant[[i]]<-vector()
for (j in 1:dim(pred.test.clm.loglog.equidistant[[i]])[1])
{h<-min(which(cumsum(pred.test.clm.loglog.equidistant[[i]][j,])>=0.5))
pred.Ord.MAP.clm.loglog.equidistant[[i]][j]<-categories[h]
}
}

###

pred.Ord.MAP.clm.cloglog.flexible<-list()
for (i in 1:10)
{ pred.Ord.MAP.clm.cloglog.flexible[[i]]<-vector()
for (j in 1:dim(pred.test.clm.cloglog.flexible[[i]])[1])
{h<-min(which(cumsum(pred.test.clm.cloglog.flexible[[i]][j,])>=0.5))
pred.Ord.MAP.clm.cloglog.flexible[[i]][j]<-categories[h]
}
}

pred.Ord.MAP.clm.cloglog.symmetric<-list()
for (i in 1:10)
{ pred.Ord.MAP.clm.cloglog.symmetric[[i]]<-vector()
for (j in 1:dim(pred.test.clm.cloglog.symmetric[[i]])[1])
{h<-min(which(cumsum(pred.test.clm.cloglog.symmetric[[i]][j,])>=0.5))
pred.Ord.MAP.clm.cloglog.symmetric[[i]][j]<-categories[h]
}
}

pred.Ord.MAP.clm.cloglog.symmetric2<-list()
for (i in 1:10)
{ pred.Ord.MAP.clm.cloglog.symmetric2[[i]]<-vector()
for (j in 1:dim(pred.test.clm.cloglog.symmetric2[[i]])[1])
{h<-min(which(cumsum(pred.test.clm.cloglog.symmetric2[[i]][j,])>=0.5))
pred.Ord.MAP.clm.cloglog.symmetric2[[i]][j]<-categories[h]
}
}

pred.Ord.MAP.clm.cloglog.equidistant<-list()
for (i in 1:10)
{ pred.Ord.MAP.clm.cloglog.equidistant[[i]]<-vector()
for (j in 1:dim(pred.test.clm.cloglog.equidistant[[i]])[1])
{h<-min(which(cumsum(pred.test.clm.cloglog.equidistant[[i]][j,])>=0.5))
pred.Ord.MAP.clm.cloglog.equidistant[[i]][j]<-categories[h]
}
}

###

pred.Ord.MAP.clm.cauchit.flexible<-list()
for (i in 1:10)
{ pred.Ord.MAP.clm.cauchit.flexible[[i]]<-vector()
for (j in 1:dim(pred.test.clm.cauchit.flexible[[i]])[1])
{h<-min(which(cumsum(pred.test.clm.cauchit.flexible[[i]][j,])>=0.5))
pred.Ord.MAP.clm.cauchit.flexible[[i]][j]<-categories[h]
}
}

pred.Ord.MAP.clm.cauchit.symmetric<-list()
for (i in 1:10)
{ pred.Ord.MAP.clm.cauchit.symmetric[[i]]<-vector()
for (j in 1:dim(pred.test.clm.cauchit.symmetric[[i]])[1])
{h<-min(which(cumsum(pred.test.clm.cauchit.symmetric[[i]][j,])>=0.5))
pred.Ord.MAP.clm.cauchit.symmetric[[i]][j]<-categories[h]
}
}

pred.Ord.MAP.clm.cauchit.symmetric2<-list()
for (i in 1:10)
{ pred.Ord.MAP.clm.cauchit.symmetric2[[i]]<-vector()
for (j in 1:dim(pred.test.clm.cauchit.symmetric2[[i]])[1])
{h<-min(which(cumsum(pred.test.clm.cauchit.symmetric2[[i]][j,])>=0.5))
pred.Ord.MAP.clm.cauchit.symmetric2[[i]][j]<-categories[h]
}
}

pred.Ord.MAP.clm.cauchit.equidistant<-list()
for (i in 1:10)
{ pred.Ord.MAP.clm.cauchit.equidistant[[i]]<-vector()
for (j in 1:dim(pred.test.clm.cauchit.equidistant[[i]])[1])
{h<-min(which(cumsum(pred.test.clm.cauchit.equidistant[[i]][j,])>=0.5))
pred.Ord.MAP.clm.cauchit.equidistant[[i]][j]<-categories[h]
}
}

###
# confusion matrices

conf.mat.Ord.MAP.clm.logit.flexible<-list()
for(i in 1:10)
{
  conf.mat.Ord.MAP.clm.logit.flexible[[i]]<-mat.square(table(pred.Ord.MAP.clm.logit.flexible[[i]],test[[i]]$V5.bin.num.factor),categories)
}

conf.mat.Ord.MAP.clm.logit.symmetric<-list()
for(i in 1:10)
{
  conf.mat.Ord.MAP.clm.logit.symmetric[[i]]<-mat.square(table(pred.Ord.MAP.clm.logit.symmetric[[i]],test[[i]]$V5.bin.num.factor),categories)
}

conf.mat.Ord.MAP.clm.logit.symmetric2<-list()
for(i in 1:10)
{
  conf.mat.Ord.MAP.clm.logit.symmetric2[[i]]<-mat.square(table(pred.Ord.MAP.clm.logit.symmetric2[[i]],test[[i]]$V5.bin.num.factor),categories)
}

conf.mat.Ord.MAP.clm.logit.equidistant<-list()
for(i in 1:10)
{
  conf.mat.Ord.MAP.clm.logit.equidistant[[i]]<-mat.square(table(pred.Ord.MAP.clm.logit.equidistant[[i]],test[[i]]$V5.bin.num.factor),categories)
}

###

conf.mat.Ord.MAP.clm.probit.flexible<-list()
for(i in 1:10)
{
  conf.mat.Ord.MAP.clm.probit.flexible[[i]]<-mat.square(table(pred.Ord.MAP.clm.probit.flexible[[i]],test[[i]]$V5.bin.num.factor),categories)
}

conf.mat.Ord.MAP.clm.probit.symmetric<-list()
for(i in 1:10)
{
  conf.mat.Ord.MAP.clm.probit.symmetric[[i]]<-mat.square(table(pred.Ord.MAP.clm.probit.symmetric[[i]],test[[i]]$V5.bin.num.factor),categories)
}

conf.mat.Ord.MAP.clm.probit.symmetric2<-list()
for(i in 1:10)
{
  conf.mat.Ord.MAP.clm.probit.symmetric2[[i]]<-mat.square(table(pred.Ord.MAP.clm.probit.symmetric2[[i]],test[[i]]$V5.bin.num.factor),categories)
}

conf.mat.Ord.MAP.clm.probit.equidistant<-list()
for(i in 1:10)
{
  conf.mat.Ord.MAP.clm.probit.equidistant[[i]]<-mat.square(table(pred.Ord.MAP.clm.probit.equidistant[[i]],test[[i]]$V5.bin.num.factor),categories)
}

###

conf.mat.Ord.MAP.clm.loglog.flexible<-list()
for(i in 1:10)
{
  conf.mat.Ord.MAP.clm.loglog.flexible[[i]]<-mat.square(table(pred.Ord.MAP.clm.loglog.flexible[[i]],test[[i]]$V5.bin.num.factor),categories)
}

conf.mat.Ord.MAP.clm.loglog.symmetric<-list()
for(i in 1:10)
{
  conf.mat.Ord.MAP.clm.loglog.symmetric[[i]]<-mat.square(table(pred.Ord.MAP.clm.loglog.symmetric[[i]],test[[i]]$V5.bin.num.factor),categories)
}

conf.mat.Ord.MAP.clm.loglog.symmetric2<-list()
for(i in 1:10)
{
  conf.mat.Ord.MAP.clm.loglog.symmetric2[[i]]<-mat.square(table(pred.Ord.MAP.clm.loglog.symmetric2[[i]],test[[i]]$V5.bin.num.factor),categories)
}

conf.mat.Ord.MAP.clm.loglog.equidistant<-list()
for(i in 1:10)
{
  conf.mat.Ord.MAP.clm.loglog.equidistant[[i]]<-mat.square(table(pred.Ord.MAP.clm.loglog.equidistant[[i]],test[[i]]$V5.bin.num.factor),categories)
}

###


conf.mat.Ord.MAP.clm.cloglog.flexible<-list()
for(i in 1:10)
{
  conf.mat.Ord.MAP.clm.cloglog.flexible[[i]]<-mat.square(table(pred.Ord.MAP.clm.cloglog.flexible[[i]],test[[i]]$V5.bin.num.factor),categories)
}

conf.mat.Ord.MAP.clm.cloglog.symmetric<-list()
for(i in 1:10)
{
  conf.mat.Ord.MAP.clm.cloglog.symmetric[[i]]<-mat.square(table(pred.Ord.MAP.clm.cloglog.symmetric[[i]],test[[i]]$V5.bin.num.factor),categories)
}

conf.mat.Ord.MAP.clm.cloglog.symmetric2<-list()
for(i in 1:10)
{
  conf.mat.Ord.MAP.clm.cloglog.symmetric2[[i]]<-mat.square(table(pred.Ord.MAP.clm.cloglog.symmetric2[[i]],test[[i]]$V5.bin.num.factor),categories)
}

conf.mat.Ord.MAP.clm.cloglog.equidistant<-list()
for(i in 1:10)
{
  conf.mat.Ord.MAP.clm.cloglog.equidistant[[i]]<-mat.square(table(pred.Ord.MAP.clm.cloglog.equidistant[[i]],test[[i]]$V5.bin.num.factor),categories)
}

###

conf.mat.Ord.MAP.clm.cauchit.flexible<-list()
for(i in 1:10)
{
  conf.mat.Ord.MAP.clm.cauchit.flexible[[i]]<-mat.square(table(pred.Ord.MAP.clm.cauchit.flexible[[i]],test[[i]]$V5.bin.num.factor),categories)
}

conf.mat.Ord.MAP.clm.cauchit.symmetric<-list()
for(i in 1:10)
{
  conf.mat.Ord.MAP.clm.cauchit.symmetric[[i]]<-mat.square(table(pred.Ord.MAP.clm.cauchit.symmetric[[i]],test[[i]]$V5.bin.num.factor),categories)
}

conf.mat.Ord.MAP.clm.cauchit.symmetric2<-list()
for(i in 1:10)
{
  conf.mat.Ord.MAP.clm.cauchit.symmetric2[[i]]<-mat.square(table(pred.Ord.MAP.clm.cauchit.symmetric2[[i]],test[[i]]$V5.bin.num.factor),categories)
}

conf.mat.Ord.MAP.clm.cauchit.equidistant<-list()
for(i in 1:10)
{
  conf.mat.Ord.MAP.clm.cauchit.equidistant[[i]]<-mat.square(table(pred.Ord.MAP.clm.cauchit.equidistant[[i]],test[[i]]$V5.bin.num.factor),categories)
}

###

# MAE 

MAE.Ord.MAP.clm.logit.flexible<-vector()
for (i in 1:10)
{
  MAE.Ord.MAP.clm.logit.flexible[i]<-MAE(conf.mat.Ord.MAP.clm.logit.flexible[[i]])
}

MAE.Ord.MAP.clm.logit.symmetric<-vector()
for (i in 1:10)
{
  MAE.Ord.MAP.clm.logit.symmetric[i]<-MAE(conf.mat.Ord.MAP.clm.logit.symmetric[[i]])
}

MAE.Ord.MAP.clm.logit.symmetric2<-vector()
for (i in 1:10)
{
  MAE.Ord.MAP.clm.logit.symmetric2[i]<-MAE(conf.mat.Ord.MAP.clm.logit.symmetric2[[i]])
}

MAE.Ord.MAP.clm.logit.equidistant<-vector()
for (i in 1:10)
{
  MAE.Ord.MAP.clm.logit.equidistant[i]<-MAE(conf.mat.Ord.MAP.clm.logit.equidistant[[i]])
}

###

MAE.Ord.MAP.clm.probit.flexible<-vector()
for (i in 1:10)
{
  MAE.Ord.MAP.clm.probit.flexible[i]<-MAE(conf.mat.Ord.MAP.clm.probit.flexible[[i]])
}

MAE.Ord.MAP.clm.probit.symmetric<-vector()
for (i in 1:10)
{
  MAE.Ord.MAP.clm.probit.symmetric[i]<-MAE(conf.mat.Ord.MAP.clm.probit.symmetric[[i]])
}

MAE.Ord.MAP.clm.probit.symmetric2<-vector()
for (i in 1:10)
{
  MAE.Ord.MAP.clm.probit.symmetric2[i]<-MAE(conf.mat.Ord.MAP.clm.probit.symmetric2[[i]])
}

MAE.Ord.MAP.clm.probit.equidistant<-vector()
for (i in 1:10)
{
  MAE.Ord.MAP.clm.probit.equidistant[i]<-MAE(conf.mat.Ord.MAP.clm.probit.equidistant[[i]])
}

###


MAE.Ord.MAP.clm.loglog.flexible<-vector()
for (i in 1:10)
{
  MAE.Ord.MAP.clm.loglog.flexible[i]<-MAE(conf.mat.Ord.MAP.clm.loglog.flexible[[i]])
}

MAE.Ord.MAP.clm.loglog.symmetric<-vector()
for (i in 1:10)
{
  MAE.Ord.MAP.clm.loglog.symmetric[i]<-MAE(conf.mat.Ord.MAP.clm.loglog.symmetric[[i]])
}

MAE.Ord.MAP.clm.loglog.symmetric2<-vector()
for (i in 1:10)
{
  MAE.Ord.MAP.clm.loglog.symmetric2[i]<-MAE(conf.mat.Ord.MAP.clm.loglog.symmetric2[[i]])
}

MAE.Ord.MAP.clm.loglog.equidistant<-vector()
for (i in 1:10)
{
  MAE.Ord.MAP.clm.loglog.equidistant[i]<-MAE(conf.mat.Ord.MAP.clm.loglog.equidistant[[i]])
}

###


MAE.Ord.MAP.clm.cloglog.flexible<-vector()
for (i in 1:10)
{
  MAE.Ord.MAP.clm.cloglog.flexible[i]<-MAE(conf.mat.Ord.MAP.clm.cloglog.flexible[[i]])
}

MAE.Ord.MAP.clm.cloglog.symmetric<-vector()
for (i in 1:10)
{
  MAE.Ord.MAP.clm.cloglog.symmetric[i]<-MAE(conf.mat.Ord.MAP.clm.cloglog.symmetric[[i]])
}

MAE.Ord.MAP.clm.cloglog.symmetric2<-vector()
for (i in 1:10)
{
  MAE.Ord.MAP.clm.cloglog.symmetric2[i]<-MAE(conf.mat.Ord.MAP.clm.cloglog.symmetric2[[i]])
}

MAE.Ord.MAP.clm.cloglog.equidistant<-vector()
for (i in 1:10)
{
  MAE.Ord.MAP.clm.cloglog.equidistant[i]<-MAE(conf.mat.Ord.MAP.clm.cloglog.equidistant[[i]])
}

###


MAE.Ord.MAP.clm.cauchit.flexible<-vector()
for (i in 1:10)
{
  MAE.Ord.MAP.clm.cauchit.flexible[i]<-MAE(conf.mat.Ord.MAP.clm.cauchit.flexible[[i]])
}

MAE.Ord.MAP.clm.cauchit.symmetric<-vector()
for (i in 1:10)
{
  MAE.Ord.MAP.clm.cauchit.symmetric[i]<-MAE(conf.mat.Ord.MAP.clm.cauchit.symmetric[[i]])
}

MAE.Ord.MAP.clm.cauchit.symmetric2<-vector()
for (i in 1:10)
{
  MAE.Ord.MAP.clm.cauchit.symmetric2[i]<-MAE(conf.mat.Ord.MAP.clm.cauchit.symmetric2[[i]])
}

MAE.Ord.MAP.clm.cauchit.equidistant<-vector()
for (i in 1:10)
{
  MAE.Ord.MAP.clm.cauchit.equidistant[i]<-MAE(conf.mat.Ord.MAP.clm.cauchit.equidistant[[i]])
}


###
################################################################################
###. Comparison MAP vs Ord.MAP
###
###.  Tests MAE

p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.logit.flexible<-wilcox.test(MAE.MAP.clm.logit.flexible,MAE.Ord.MAP.clm.logit.flexible,paired = TRUE,alternative="greater")$p.value
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.logit.symmetric<-wilcox.test(MAE.MAP.clm.logit.symmetric,MAE.Ord.MAP.clm.logit.symmetric,paired = TRUE,alternative="greater")$p.value
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.logit.symmetric2<-wilcox.test(MAE.MAP.clm.logit.symmetric2,MAE.Ord.MAP.clm.logit.symmetric2,paired = TRUE,alternative="greater")$p.value
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.logit.equidistant<-wilcox.test(MAE.MAP.clm.logit.equidistant,MAE.Ord.MAP.clm.logit.equidistant,paired = TRUE,alternative="greater")$p.value

p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.logit.flexible.less<-wilcox.test(MAE.MAP.clm.logit.flexible,MAE.Ord.MAP.clm.logit.flexible,paired = TRUE,alternative="less")$p.value
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.logit.symmetric.less<-wilcox.test(MAE.MAP.clm.logit.symmetric,MAE.Ord.MAP.clm.logit.symmetric,paired = TRUE,alternative="less")$p.value
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.logit.symmetric2.less<-wilcox.test(MAE.MAP.clm.logit.symmetric2,MAE.Ord.MAP.clm.logit.symmetric2,paired = TRUE,alternative="less")$p.value
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.logit.equidistant.less<-wilcox.test(MAE.MAP.clm.logit.equidistant,MAE.Ord.MAP.clm.logit.equidistant,paired = TRUE,alternative="less")$p.value

###

p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.probit.flexible<-wilcox.test(MAE.MAP.clm.probit.flexible,MAE.Ord.MAP.clm.probit.flexible,paired = TRUE,alternative="greater")$p.value
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.probit.symmetric<-wilcox.test(MAE.MAP.clm.probit.symmetric,MAE.Ord.MAP.clm.probit.symmetric,paired = TRUE,alternative="greater")$p.value
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.probit.symmetric2<-wilcox.test(MAE.MAP.clm.probit.symmetric2,MAE.Ord.MAP.clm.probit.symmetric2,paired = TRUE,alternative="greater")$p.value
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.probit.equidistant<-wilcox.test(MAE.MAP.clm.probit.equidistant,MAE.Ord.MAP.clm.probit.equidistant,paired = TRUE,alternative="greater")$p.value

p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.probit.flexible.less<-wilcox.test(MAE.MAP.clm.probit.flexible,MAE.Ord.MAP.clm.probit.flexible,paired = TRUE,alternative="less")$p.value
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.probit.symmetric.less<-wilcox.test(MAE.MAP.clm.probit.symmetric,MAE.Ord.MAP.clm.probit.symmetric,paired = TRUE,alternative="less")$p.value
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.probit.symmetric2.less<-wilcox.test(MAE.MAP.clm.probit.symmetric2,MAE.Ord.MAP.clm.probit.symmetric2,paired = TRUE,alternative="less")$p.value
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.probit.equidistant.less<-wilcox.test(MAE.MAP.clm.probit.equidistant,MAE.Ord.MAP.clm.probit.equidistant,paired = TRUE,alternative="less")$p.value


###

p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.loglog.flexible<-wilcox.test(MAE.MAP.clm.loglog.flexible,MAE.Ord.MAP.clm.loglog.flexible,paired = TRUE,alternative="greater")$p.value
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.loglog.symmetric<-wilcox.test(MAE.MAP.clm.loglog.symmetric,MAE.Ord.MAP.clm.loglog.symmetric,paired = TRUE,alternative="greater")$p.value
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.loglog.symmetric2<-wilcox.test(MAE.MAP.clm.loglog.symmetric2,MAE.Ord.MAP.clm.loglog.symmetric2,paired = TRUE,alternative="greater")$p.value
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.loglog.equidistant<-wilcox.test(MAE.MAP.clm.loglog.equidistant,MAE.Ord.MAP.clm.loglog.equidistant,paired = TRUE,alternative="greater")$p.value

p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.loglog.flexible.less<-wilcox.test(MAE.MAP.clm.loglog.flexible,MAE.Ord.MAP.clm.loglog.flexible,paired = TRUE,alternative="less")$p.value
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.loglog.symmetric.less<-wilcox.test(MAE.MAP.clm.loglog.symmetric,MAE.Ord.MAP.clm.loglog.symmetric,paired = TRUE,alternative="less")$p.value
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.loglog.symmetric2.less<-wilcox.test(MAE.MAP.clm.loglog.symmetric2,MAE.Ord.MAP.clm.loglog.symmetric2,paired = TRUE,alternative="less")$p.value
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.loglog.equidistant.less<-wilcox.test(MAE.MAP.clm.loglog.equidistant,MAE.Ord.MAP.clm.loglog.equidistant,paired = TRUE,alternative="less")$p.value

###

p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.cloglog.flexible<-wilcox.test(MAE.MAP.clm.cloglog.flexible,MAE.Ord.MAP.clm.cloglog.flexible,paired = TRUE,alternative="greater")$p.value
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.cloglog.symmetric<-wilcox.test(MAE.MAP.clm.cloglog.symmetric,MAE.Ord.MAP.clm.cloglog.symmetric,paired = TRUE,alternative="greater")$p.value
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.cloglog.symmetric2<-wilcox.test(MAE.MAP.clm.cloglog.symmetric2,MAE.Ord.MAP.clm.cloglog.symmetric2,paired = TRUE,alternative="greater")$p.value
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.cloglog.equidistant<-wilcox.test(MAE.MAP.clm.cloglog.equidistant,MAE.Ord.MAP.clm.cloglog.equidistant,paired = TRUE,alternative="greater")$p.value

p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.cloglog.flexible.less<-wilcox.test(MAE.MAP.clm.cloglog.flexible,MAE.Ord.MAP.clm.cloglog.flexible,paired = TRUE,alternative="less")$p.value
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.cloglog.symmetric.less<-wilcox.test(MAE.MAP.clm.cloglog.symmetric,MAE.Ord.MAP.clm.cloglog.symmetric,paired = TRUE,alternative="less")$p.value
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.cloglog.symmetric2.less<-wilcox.test(MAE.MAP.clm.cloglog.symmetric2,MAE.Ord.MAP.clm.cloglog.symmetric2,paired = TRUE,alternative="less")$p.value
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.cloglog.equidistant.less<-wilcox.test(MAE.MAP.clm.cloglog.equidistant,MAE.Ord.MAP.clm.cloglog.equidistant,paired = TRUE,alternative="les")$p.value

###

p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.cauchit.flexible<-wilcox.test(MAE.MAP.clm.cauchit.flexible,MAE.Ord.MAP.clm.cauchit.flexible,paired = TRUE,alternative="greater")$p.value
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.cauchit.symmetric<-wilcox.test(MAE.MAP.clm.cauchit.symmetric,MAE.Ord.MAP.clm.cauchit.symmetric,paired = TRUE,alternative="greater")$p.value
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.cauchit.symmetric2<-wilcox.test(MAE.MAP.clm.cauchit.symmetric2,MAE.Ord.MAP.clm.cauchit.symmetric2,paired = TRUE,alternative="greater")$p.value
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.cauchit.equidistant<-wilcox.test(MAE.MAP.clm.cauchit.equidistant,MAE.Ord.MAP.clm.cauchit.equidistant,paired = TRUE,alternative="greater")$p.value

p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.cauchit.flexible.less<-wilcox.test(MAE.MAP.clm.cauchit.flexible,MAE.Ord.MAP.clm.cauchit.flexible,paired = TRUE,alternative="less")$p.value
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.cauchit.symmetric.less<-wilcox.test(MAE.MAP.clm.cauchit.symmetric,MAE.Ord.MAP.clm.cauchit.symmetric,paired = TRUE,alternative="less")$p.value
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.cauchit.symmetric2.less<-wilcox.test(MAE.MAP.clm.cauchit.symmetric2,MAE.Ord.MAP.clm.cauchit.symmetric2,paired = TRUE,alternative="less")$p.value
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.cauchit.equidistant.less<-wilcox.test(MAE.MAP.clm.cauchit.equidistant,MAE.Ord.MAP.clm.cauchit.equidistant,paired = TRUE,alternative="less")$p.value

###

p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.logit.flexible
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.logit.symmetric
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.logit.symmetric2
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.logit.equidistant

p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.logit.flexible.less
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.logit.symmetric.less
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.logit.symmetric2.less
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.logit.equidistant.less

# 

p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.probit.flexible
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.probit.symmetric
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.probit.symmetric2
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.probit.equidistant

p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.probit.flexible.less
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.probit.symmetric.less
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.probit.symmetric2.less
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.probit.equidistant.less

#

p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.loglog.flexible
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.loglog.symmetric
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.loglog.symmetric2
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.loglog.equidistant

p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.loglog.flexible.less
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.loglog.symmetric.less
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.loglog.symmetric2.less
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.loglog.equidistant.less

#

p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.cloglog.flexible
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.cloglog.symmetric
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.cloglog.symmetric2
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.cloglog.equidistant

p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.cloglog.flexible.less
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.cloglog.symmetric.less
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.cloglog.symmetric2.less
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.cloglog.equidistant.less

#

p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.cauchit.flexible
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.cauchit.symmetric
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.cauchit.symmetric2
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.cauchit.equidistant

p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.cauchit.flexible.less
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.cauchit.symmetric.less
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.cauchit.symmetric2.less
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.cauchit.equidistant.less


###############################
###############################
###   Boxplots   MAE

MAE.clm.flexible.results<-as.data.frame(cbind(MAE.MAP.clm.logit.flexible-MAE.Ord.MAP.clm.logit.flexible,
                                              #MAE.MAP.clm.probit.flexible-MAE.Ord.MAP.clm.probit.flexible,
                                              #MAE.MAP.clm.loglog.flexible-MAE.Ord.MAP.clm.loglog.flexible,
                                              #MAE.MAP.clm.cloglog.flexible-MAE.Ord.MAP.clm.cloglog.flexible,
                                              MAE.MAP.clm.cauchit.flexible-MAE.Ord.MAP.clm.cauchit.flexible))

boxplot(MAE.clm.flexible.results,
        main=paste("(d) Parkinson dataset: output V5 (6 classes). Threshold flexible"), 
        xlab=" ", ylab="MAE increment (MAP - Ord.MAP)",names=c("logistic", 
                                                               #"probit","loglog","cloglog",
                                                               "cauchit"),
        col="white")
abline(h = 0, col = 2, lwd = 2, lty=2)
# Dots
stripchart(MAE.clm.flexible.results,
           #method = "jitter",
           pch = 19,
           #col = 1:5,
           col = 1:2,
           vertical = TRUE,
           add = TRUE)

###

MAE.clm.symmetric.results<-as.data.frame(cbind(MAE.MAP.clm.logit.symmetric-MAE.Ord.MAP.clm.logit.symmetric,
                                               #MAE.MAP.clm.probit.symmetric-MAE.Ord.MAP.clm.probit.symmetric,
                                               #MAE.MAP.clm.loglog.symmetric-MAE.Ord.MAP.clm.loglog.symmetric,
                                               #MAE.MAP.clm.cloglog.symmetric-MAE.Ord.MAP.clm.cloglog.symmetric,
                                               MAE.MAP.clm.cauchit.symmetric-MAE.Ord.MAP.clm.cauchit.symmetric))

boxplot(MAE.clm.symmetric.results,
        main=paste("(d) Parkinson dataset: output V5 (6 classes). Threshold symmetric"), 
        xlab=" ", ylab="MAE increment (MAP - Ord.MAP)",names=c("logistic", 
                                                               #"probit","loglog","cloglog",
                                                               "cauchit"),
        col="white")
abline(h = 0, col = 2, lwd = 2, lty=2)
# Dots
stripchart(MAE.clm.symmetric.results,
           #method = "jitter",
           pch = 19,
           #col = 1:5,
           col = 1:2,
           vertical = TRUE,
           add = TRUE)

###

MAE.clm.symmetric2.results<-as.data.frame(cbind(MAE.MAP.clm.logit.symmetric2-MAE.Ord.MAP.clm.logit.symmetric2,
                                                #MAE.MAP.clm.probit.symmetric2-MAE.Ord.MAP.clm.probit.symmetric2,
                                                #MAE.MAP.clm.loglog.symmetric2-MAE.Ord.MAP.clm.loglog.symmetric2,
                                                #MAE.MAP.clm.cloglog.symmetric2-MAE.Ord.MAP.clm.cloglog.symmetric2,
                                                MAE.MAP.clm.cauchit.symmetric2-MAE.Ord.MAP.clm.cauchit.symmetric2))

boxplot(MAE.clm.symmetric2.results,
        main=paste("(d) Parkinson dataset: output V5 (6 classes). Threshold symmetric2"), 
        xlab=" ", ylab="MAE increment (MAP - Ord.MAP)",names=c("logistic", 
                                                               #"probit","loglog","cloglog",
                                                               "cauchit"),
        col="white")
abline(h = 0, col = 2, lwd = 2, lty=2)
# Dots
stripchart(MAE.clm.symmetric2.results,
           #method = "jitter",
           pch = 19,
           #col = 1:5,
           col = 1:2,
           vertical = TRUE,
           add = TRUE)

###

MAE.clm.equidistant.results<-as.data.frame(cbind(MAE.MAP.clm.logit.equidistant-MAE.Ord.MAP.clm.logit.equidistant,
                                                 #MAE.MAP.clm.probit.equidistant-MAE.Ord.MAP.clm.probit.equidistant,
                                                 #MAE.MAP.clm.loglog.equidistant-MAE.Ord.MAP.clm.loglog.equidistant,
                                                 #MAE.MAP.clm.cloglog.equidistant-MAE.Ord.MAP.clm.cloglog.equidistant,
                                                 MAE.MAP.clm.cauchit.equidistant-MAE.Ord.MAP.clm.cauchit.equidistant))

boxplot(MAE.clm.equidistant.results,
        main=paste("(d) Parkinson dataset: output V5 (6 classes). Threshold equidistant"), 
        xlab=" ", ylab="MAE increment (MAP - Ord.MAP)",names=c("logistic", 
                                                               #"probit","loglog","cloglog",
                                                               "cauchit"),
        col="white")
abline(h = 0, col = 2, lwd = 2, lty=2)
# Dots
stripchart(MAE.clm.equidistant.results,
           #method = "jitter",
           pch = 19,
           #col = 1:5,
           col = 1:2,
           vertical = TRUE,
           add = TRUE)


#
##
###
##
#

#### SECOND: OUTPUT VARIABLE TO PREDICT: V6.bin.num.factor

model.clm.logit.flexible<-list()
for (i in 1:10)
{model.clm.logit.flexible[[i]] <- clm(V6.bin.num.factor ~ V2 + V3 + V4 + V7 + V8 + V9 + V10 +
                                        V11 + V12 + V13 + V14 + V15 + V16 + V17 + V18 + V19 + 
                                        V20 + V21 +V22, data = training[[i]],
                                      link="logit", threshold="flexible")}


model.clm.logit.symmetric<-list()
for (i in 1:10)
{model.clm.logit.symmetric[[i]] <- clm(V6.bin.num.factor ~ V2 + V3 + V4 + V7 + V8 + V9 + V10 +
                                         V11 + V12 + V13 + V14 + V15 + V16 + V17 + V18 + V19 + 
                                         V20 + V21 +V22, data = training[[i]],
                                       link="logit", threshold="symmetric")}


model.clm.logit.symmetric2<-list()
for (i in 1:10)
{model.clm.logit.symmetric2[[i]] <- clm(V6.bin.num.factor ~ V2 + V3 + V4 + V7 + V8 + V9 + V10 +
                                          V11 + V12 + V13 + V14 + V15 + V16 + V17 + V18 + V19 + 
                                          V20 + V21 +V22, data = training[[i]],
                                        link="logit", threshold="symmetric2")}

model.clm.logit.equidistant<-list()
for (i in 1:10)
{model.clm.logit.equidistant[[i]] <- clm(V6.bin.num.factor ~ V2 + V3 + V4 + V7 + V8 + V9 + V10 +
                                           V11 + V12 + V13 + V14 + V15 + V16 + V17 + V18 + V19 + 
                                           V20 + V21 +V22, data = training[[i]],
                                         link="logit", threshold="equidistant")}


###

model.clm.probit.flexible<-list()
for (i in 1:10)
{model.clm.probit.flexible[[i]] <- clm(V6.bin.num.factor ~ V2 + V3 + V4 + V7 + V8 + V9 + V10 +
                                         V11 + V12 + V13 + V14 + V15 + V16 + V17 + V18 + V19 + 
                                         V20 + V21 +V22, data = training[[i]],
                                       link="probit", threshold="flexible")}


model.clm.probit.symmetric<-list()
for (i in 1:10)
{model.clm.probit.symmetric[[i]] <- clm(V6.bin.num.factor ~ V2 + V3 + V4 + V7 + V8 + V9 + V10 +
                                          V11 + V12 + V13 + V14 + V15 + V16 + V17 + V18 + V19 + 
                                          V20 + V21 +V22, data = training[[i]],
                                        link="probit", threshold="symmetric")}


model.clm.probit.symmetric2<-list()
for (i in 1:10)
{model.clm.probit.symmetric2[[i]] <- clm(V6.bin.num.factor ~ V2 + V3 + V4 + V7 + V8 + V9 + V10 +
                                           V11 + V12 + V13 + V14 + V15 + V16 + V17 + V18 + V19 + 
                                           V20 + V21 +V22, data = training[[i]],
                                         link="probit", threshold="symmetric2")}

model.clm.probit.equidistant<-list()
for (i in 1:10)
{model.clm.probit.equidistant[[i]] <- clm(V6.bin.num.factor ~ V2 + V3 + V4 + V7 + V8 + V9 + V10 +
                                            V11 + V12 + V13 + V14 + V15 + V16 + V17 + V18 + V19 + 
                                            V20 + V21 +V22, data = training[[i]],
                                          link="probit", threshold="equidistant")}


###

model.clm.loglog.flexible<-list()
for (i in 1:10)
{model.clm.loglog.flexible[[i]] <- clm(V6.bin.num.factor ~ V2 + V3 + V4 + V7 + V8 + V9 + V10 +
                                         V11 + V12 + V13 + V14 + V15 + V16 + V17 + V18 + V19 + 
                                         V20 + V21 +V22, data = training[[i]],
                                       link="loglog", threshold="flexible")}


model.clm.loglog.symmetric<-list()
for (i in 1:10)
{model.clm.loglog.symmetric[[i]] <- clm(V6.bin.num.factor ~ V2 + V3 + V4 + V7 + V8 + V9 + V10 +
                                          V11 + V12 + V13 + V14 + V15 + V16 + V17 + V18 + V19 + 
                                          V20 + V21 +V22, data = training[[i]],
                                        link="loglog", threshold="symmetric")}


model.clm.loglog.symmetric2<-list()
for (i in 1:10)
{model.clm.loglog.symmetric2[[i]] <- clm(V6.bin.num.factor ~ V2 + V3 + V4 + V7 + V8 + V9 + V10 +
                                           V11 + V12 + V13 + V14 + V15 + V16 + V17 + V18 + V19 + 
                                           V20 + V21 +V22, data = training[[i]],
                                         link="loglog", threshold="symmetric2")}

model.clm.loglog.equidistant<-list()
for (i in 1:10)
{model.clm.loglog.equidistant[[i]] <- clm(V6.bin.num.factor ~ V2 + V3 + V4 + V7 + V8 + V9 + V10 +
                                            V11 + V12 + V13 + V14 + V15 + V16 + V17 + V18 + V19 + 
                                            V20 + V21 +V22, data = training[[i]],
                                          link="loglog", threshold="equidistant")}


###

model.clm.cloglog.flexible<-list()
for (i in 1:10)
{model.clm.cloglog.flexible[[i]] <- clm(V6.bin.num.factor ~ V2 + V3 + V4 + V7 + V8 + V9 + V10 +
                                          V11 + V12 + V13 + V14 + V15 + V16 + V17 + V18 + V19 + 
                                          V20 + V21 +V22, data = training[[i]],
                                        link="cloglog", threshold="flexible")}


model.clm.cloglog.symmetric<-list()
for (i in 1:10)
{model.clm.cloglog.symmetric[[i]] <- clm(V6.bin.num.factor ~ V2 + V3 + V4 + V7 + V8 + V9 + V10 +
                                           V11 + V12 + V13 + V14 + V15 + V16 + V17 + V18 + V19 + 
                                           V20 + V21 +V22, data = training[[i]],
                                         link="cloglog", threshold="symmetric")}


model.clm.cloglog.symmetric2<-list()
for (i in 1:10)
{model.clm.cloglog.symmetric2[[i]] <- clm(V6.bin.num.factor ~ V2 + V3 + V4 + V7 + V8 + V9 + V10 +
                                            V11 + V12 + V13 + V14 + V15 + V16 + V17 + V18 + V19 + 
                                            V20 + V21 +V22, data = training[[i]],
                                          link="cloglog", threshold="symmetric2")}

model.clm.cloglog.equidistant<-list()
for (i in 1:10)
{model.clm.cloglog.equidistant[[i]] <- clm(V6.bin.num.factor ~ V2 + V3 + V4 + V7 + V8 + V9 + V10 +
                                             V11 + V12 + V13 + V14 + V15 + V16 + V17 + V18 + V19 + 
                                             V20 + V21 +V22, data = training[[i]],
                                           link="cloglog", threshold="equidistant")}


###

model.clm.cauchit.flexible<-list()
for (i in 1:10)
{model.clm.cauchit.flexible[[i]] <- clm(V6.bin.num.factor ~ V2 + V3 + V4 + V7 + V8 + V9 + V10 +
                                          V11 + V12 + V13 + V14 + V15 + V16 + V17 + V18 + V19 + 
                                          V20 + V21 +V22, data = training[[i]],
                                        link="cauchit", threshold="flexible")}


model.clm.cauchit.symmetric<-list()
for (i in 1:10)
{model.clm.cauchit.symmetric[[i]] <- clm(V6.bin.num.factor ~ V2 + V3 + V4 + V7 + V8 + V9 + V10 +
                                           V11 + V12 + V13 + V14 + V15 + V16 + V17 + V18 + V19 + 
                                           V20 + V21 +V22, data = training[[i]],
                                         link="cauchit", threshold="symmetric")}


model.clm.cauchit.symmetric2<-list()
for (i in 1:10)
{model.clm.cauchit.symmetric2[[i]] <- clm(V6.bin.num.factor ~ V2 + V3 + V4 + V7 + V8 + V9 + V10 +
                                            V11 + V12 + V13 + V14 + V15 + V16 + V17 + V18 + V19 + 
                                            V20 + V21 +V22, data = training[[i]],
                                          link="cauchit", threshold="symmetric2")}

model.clm.cauchit.equidistant<-list()
for (i in 1:10)
{model.clm.cauchit.equidistant[[i]] <- clm(V6.bin.num.factor ~ V2 + V3 + V4 + V7 + V8 + V9 + V10 +
                                             V11 + V12 + V13 + V14 + V15 + V16 + V17 + V18 + V19 + 
                                             V20 + V21 +V22, data = training[[i]],
                                           link="cauchit", threshold="equidistant")}



################################################################################
#
# Predictions of any of the test set
# 
#

vector.no<-c(1,5,6,23:28)

pred.test.clm.logit.flexible<-list()
for (i in 1:10)
{pred.test.clm.logit.flexible[[i]] <- predict(model.clm.logit.flexible[[i]],test[[i]][,-vector.no],type = "p")$fit
}

pred.test.clm.logit.symmetric<-list()
for (i in 1:10)
{pred.test.clm.logit.symmetric[[i]] <- predict(model.clm.logit.symmetric[[i]],test[[i]][,-vector.no],type = "p")$fit
}

pred.test.clm.logit.symmetric2<-list()
for (i in 1:10)
{pred.test.clm.logit.symmetric2[[i]] <- predict(model.clm.logit.symmetric2[[i]],test[[i]][,-vector.no],type = "p")$fit
}

pred.test.clm.logit.equidistant<-list()
for (i in 1:10)
{pred.test.clm.logit.equidistant[[i]] <- predict(model.clm.logit.equidistant[[i]],test[[i]][,-vector.no],type = "p")$fit
}

###

pred.test.clm.probit.flexible<-list()
for (i in 1:10)
{pred.test.clm.probit.flexible[[i]] <- predict(model.clm.probit.flexible[[i]],test[[i]][,-vector.no],type = "p")$fit
}

pred.test.clm.probit.symmetric<-list()
for (i in 1:10)
{pred.test.clm.probit.symmetric[[i]] <- predict(model.clm.probit.symmetric[[i]],test[[i]][,-vector.no],type = "p")$fit
}

pred.test.clm.probit.symmetric2<-list()
for (i in 1:10)
{pred.test.clm.probit.symmetric2[[i]] <- predict(model.clm.probit.symmetric2[[i]],test[[i]][,-vector.no],type = "p")$fit
}

pred.test.clm.probit.equidistant<-list()
for (i in 1:10)
{pred.test.clm.probit.equidistant[[i]] <- predict(model.clm.probit.equidistant[[i]],test[[i]][,-vector.no],type = "p")$fit
}

###

pred.test.clm.loglog.flexible<-list()
for (i in 1:10)
{pred.test.clm.loglog.flexible[[i]] <- predict(model.clm.loglog.flexible[[i]],test[[i]][,-vector.no],type = "p")$fit
}

pred.test.clm.loglog.symmetric<-list()
for (i in 1:10)
{pred.test.clm.loglog.symmetric[[i]] <- predict(model.clm.loglog.symmetric[[i]],test[[i]][,-vector.no],type = "p")$fit
}

pred.test.clm.loglog.symmetric2<-list()
for (i in 1:10)
{pred.test.clm.loglog.symmetric2[[i]] <- predict(model.clm.loglog.symmetric2[[i]],test[[i]][,-vector.no],type = "p")$fit
}

pred.test.clm.loglog.equidistant<-list()
for (i in 1:10)
{pred.test.clm.loglog.equidistant[[i]] <- predict(model.clm.loglog.equidistant[[i]],test[[i]][,-vector.no],type = "p")$fit
}

###


pred.test.clm.cloglog.flexible<-list()
for (i in 1:10)
{pred.test.clm.cloglog.flexible[[i]] <- predict(model.clm.cloglog.flexible[[i]],test[[i]][,-vector.no],type = "p")$fit
}

pred.test.clm.cloglog.symmetric<-list()
for (i in 1:10)
{pred.test.clm.cloglog.symmetric[[i]] <- predict(model.clm.cloglog.symmetric[[i]],test[[i]][,-vector.no],type = "p")$fit
}

pred.test.clm.cloglog.symmetric2<-list()
for (i in 1:10)
{pred.test.clm.cloglog.symmetric2[[i]] <- predict(model.clm.cloglog.symmetric2[[i]],test[[i]][,-vector.no],type = "p")$fit
}

pred.test.clm.cloglog.equidistant<-list()
for (i in 1:10)
{pred.test.clm.cloglog.equidistant[[i]] <- predict(model.clm.cloglog.equidistant[[i]],test[[i]][,-vector.no],type = "p")$fit
}

###


pred.test.clm.cauchit.flexible<-list()
for (i in 1:10)
{pred.test.clm.cauchit.flexible[[i]] <- predict(model.clm.cauchit.flexible[[i]],test[[i]][,-vector.no],type = "p")$fit
}

pred.test.clm.cauchit.symmetric<-list()
for (i in 1:10)
{pred.test.clm.cauchit.symmetric[[i]] <- predict(model.clm.cauchit.symmetric[[i]],test[[i]][,-vector.no],type = "p")$fit
}

pred.test.clm.cauchit.symmetric2<-list()
for (i in 1:10)
{pred.test.clm.cauchit.symmetric2[[i]] <- predict(model.clm.cauchit.symmetric2[[i]],test[[i]][,-vector.no],type = "p")$fit
}

pred.test.clm.cauchit.equidistant<-list()
for (i in 1:10)
{pred.test.clm.cauchit.equidistant[[i]] <- predict(model.clm.cauchit.equidistant[[i]],test[[i]][,-vector.no],type = "p")$fit
}

###

################################################################################
### predictions with standard MAP criterion:
###

pred.MAP.clm.logit.flexible<-list()
for (i in 1:10)
{pred.MAP.clm.logit.flexible[[i]]<-vector()
for (j in 1:dim(pred.test.clm.logit.flexible[[i]])[1])
{pred.MAP.clm.logit.flexible[[i]][j]<-categories[which.max(pred.test.clm.logit.flexible[[i]][j,])]
}
}

pred.MAP.clm.logit.symmetric<-list()
for (i in 1:10)
{pred.MAP.clm.logit.symmetric[[i]]<-vector()
for (j in 1:dim(pred.test.clm.logit.symmetric[[i]])[1])
{pred.MAP.clm.logit.symmetric[[i]][j]<-categories[which.max(pred.test.clm.logit.symmetric[[i]][j,])]
}
}

pred.MAP.clm.logit.symmetric2<-list()
for (i in 1:10)
{pred.MAP.clm.logit.symmetric2[[i]]<-vector()
for (j in 1:dim(pred.test.clm.logit.symmetric2[[i]])[1])
{pred.MAP.clm.logit.symmetric2[[i]][j]<-categories[which.max(pred.test.clm.logit.symmetric2[[i]][j,])]
}
}

pred.MAP.clm.logit.equidistant<-list()
for (i in 1:10)
{pred.MAP.clm.logit.equidistant[[i]]<-vector()
for (j in 1:dim(pred.test.clm.logit.equidistant[[i]])[1])
{pred.MAP.clm.logit.equidistant[[i]][j]<-categories[which.max(pred.test.clm.logit.equidistant[[i]][j,])]
}
}

###


pred.MAP.clm.probit.flexible<-list()
for (i in 1:10)
{pred.MAP.clm.probit.flexible[[i]]<-vector()
for (j in 1:dim(pred.test.clm.probit.flexible[[i]])[1])
{pred.MAP.clm.probit.flexible[[i]][j]<-categories[which.max(pred.test.clm.probit.flexible[[i]][j,])]
}
}

pred.MAP.clm.probit.symmetric<-list()
for (i in 1:10)
{pred.MAP.clm.probit.symmetric[[i]]<-vector()
for (j in 1:dim(pred.test.clm.probit.symmetric[[i]])[1])
{pred.MAP.clm.probit.symmetric[[i]][j]<-categories[which.max(pred.test.clm.probit.symmetric[[i]][j,])]
}
}

pred.MAP.clm.probit.symmetric2<-list()
for (i in 1:10)
{pred.MAP.clm.probit.symmetric2[[i]]<-vector()
for (j in 1:dim(pred.test.clm.probit.symmetric2[[i]])[1])
{pred.MAP.clm.probit.symmetric2[[i]][j]<-categories[which.max(pred.test.clm.probit.symmetric2[[i]][j,])]
}
}

pred.MAP.clm.probit.equidistant<-list()
for (i in 1:10)
{pred.MAP.clm.probit.equidistant[[i]]<-vector()
for (j in 1:dim(pred.test.clm.probit.equidistant[[i]])[1])
{pred.MAP.clm.probit.equidistant[[i]][j]<-categories[which.max(pred.test.clm.probit.equidistant[[i]][j,])]
}
}

###


pred.MAP.clm.loglog.flexible<-list()
for (i in 1:10)
{pred.MAP.clm.loglog.flexible[[i]]<-vector()
for (j in 1:dim(pred.test.clm.loglog.flexible[[i]])[1])
{pred.MAP.clm.loglog.flexible[[i]][j]<-categories[which.max(pred.test.clm.loglog.flexible[[i]][j,])]
}
}

pred.MAP.clm.loglog.symmetric<-list()
for (i in 1:10)
{pred.MAP.clm.loglog.symmetric[[i]]<-vector()
for (j in 1:dim(pred.test.clm.loglog.symmetric[[i]])[1])
{pred.MAP.clm.loglog.symmetric[[i]][j]<-categories[which.max(pred.test.clm.loglog.symmetric[[i]][j,])]
}
}

pred.MAP.clm.loglog.symmetric2<-list()
for (i in 1:10)
{pred.MAP.clm.loglog.symmetric2[[i]]<-vector()
for (j in 1:dim(pred.test.clm.loglog.symmetric2[[i]])[1])
{pred.MAP.clm.loglog.symmetric2[[i]][j]<-categories[which.max(pred.test.clm.loglog.symmetric2[[i]][j,])]
}
}

pred.MAP.clm.loglog.equidistant<-list()
for (i in 1:10)
{pred.MAP.clm.loglog.equidistant[[i]]<-vector()
for (j in 1:dim(pred.test.clm.loglog.equidistant[[i]])[1])
{pred.MAP.clm.loglog.equidistant[[i]][j]<-categories[which.max(pred.test.clm.loglog.equidistant[[i]][j,])]
}
}

###


pred.MAP.clm.cloglog.flexible<-list()
for (i in 1:10)
{pred.MAP.clm.cloglog.flexible[[i]]<-vector()
for (j in 1:dim(pred.test.clm.cloglog.flexible[[i]])[1])
{pred.MAP.clm.cloglog.flexible[[i]][j]<-categories[which.max(pred.test.clm.cloglog.flexible[[i]][j,])]
}
}

pred.MAP.clm.cloglog.symmetric<-list()
for (i in 1:10)
{pred.MAP.clm.cloglog.symmetric[[i]]<-vector()
for (j in 1:dim(pred.test.clm.cloglog.symmetric[[i]])[1])
{pred.MAP.clm.cloglog.symmetric[[i]][j]<-categories[which.max(pred.test.clm.cloglog.symmetric[[i]][j,])]
}
}

pred.MAP.clm.cloglog.symmetric2<-list()
for (i in 1:10)
{pred.MAP.clm.cloglog.symmetric2[[i]]<-vector()
for (j in 1:dim(pred.test.clm.cloglog.symmetric2[[i]])[1])
{pred.MAP.clm.cloglog.symmetric2[[i]][j]<-categories[which.max(pred.test.clm.cloglog.symmetric2[[i]][j,])]
}
}

pred.MAP.clm.cloglog.equidistant<-list()
for (i in 1:10)
{pred.MAP.clm.cloglog.equidistant[[i]]<-vector()
for (j in 1:dim(pred.test.clm.cloglog.equidistant[[i]])[1])
{pred.MAP.clm.cloglog.equidistant[[i]][j]<-categories[which.max(pred.test.clm.cloglog.equidistant[[i]][j,])]
}
}

###


pred.MAP.clm.cauchit.flexible<-list()
for (i in 1:10)
{pred.MAP.clm.cauchit.flexible[[i]]<-vector()
for (j in 1:dim(pred.test.clm.cauchit.flexible[[i]])[1])
{pred.MAP.clm.cauchit.flexible[[i]][j]<-categories[which.max(pred.test.clm.cauchit.flexible[[i]][j,])]
}
}

pred.MAP.clm.cauchit.symmetric<-list()
for (i in 1:10)
{pred.MAP.clm.cauchit.symmetric[[i]]<-vector()
for (j in 1:dim(pred.test.clm.cauchit.symmetric[[i]])[1])
{pred.MAP.clm.cauchit.symmetric[[i]][j]<-categories[which.max(pred.test.clm.cauchit.symmetric[[i]][j,])]
}
}

pred.MAP.clm.cauchit.symmetric2<-list()
for (i in 1:10)
{pred.MAP.clm.cauchit.symmetric2[[i]]<-vector()
for (j in 1:dim(pred.test.clm.cauchit.symmetric2[[i]])[1])
{pred.MAP.clm.cauchit.symmetric2[[i]][j]<-categories[which.max(pred.test.clm.cauchit.symmetric2[[i]][j,])]
}
}

pred.MAP.clm.cauchit.equidistant<-list()
for (i in 1:10)
{pred.MAP.clm.cauchit.equidistant[[i]]<-vector()
for (j in 1:dim(pred.test.clm.cauchit.equidistant[[i]])[1])
{pred.MAP.clm.cauchit.equidistant[[i]][j]<-categories[which.max(pred.test.clm.cauchit.equidistant[[i]][j,])]
}
}

# confusion matrices


conf.mat.MAP.clm.logit.flexible<-list()
for(i in 1:10)
{
  conf.mat.MAP.clm.logit.flexible[[i]]<-mat.square(table(pred.MAP.clm.logit.flexible[[i]],test[[i]]$V6.bin.num.factor),categories)
}


conf.mat.MAP.clm.logit.symmetric<-list()
for(i in 1:10)
{
  conf.mat.MAP.clm.logit.symmetric[[i]]<-mat.square(table(pred.MAP.clm.logit.symmetric[[i]],test[[i]]$V6.bin.num.factor),categories)
}


conf.mat.MAP.clm.logit.symmetric2<-list()
for(i in 1:10)
{
  conf.mat.MAP.clm.logit.symmetric2[[i]]<-mat.square(table(pred.MAP.clm.logit.symmetric2[[i]],test[[i]]$V6.bin.num.factor),categories)
}


conf.mat.MAP.clm.logit.equidistant<-list()
for(i in 1:10)
{
  conf.mat.MAP.clm.logit.equidistant[[i]]<-mat.square(table(pred.MAP.clm.logit.equidistant[[i]],test[[i]]$V6.bin.num.factor),categories)
}

###

conf.mat.MAP.clm.probit.flexible<-list()
for(i in 1:10)
{
  conf.mat.MAP.clm.probit.flexible[[i]]<-mat.square(table(pred.MAP.clm.probit.flexible[[i]],test[[i]]$V6.bin.num.factor),categories)
}


conf.mat.MAP.clm.probit.symmetric<-list()
for(i in 1:10)
{
  conf.mat.MAP.clm.probit.symmetric[[i]]<-mat.square(table(pred.MAP.clm.probit.symmetric[[i]],test[[i]]$V6.bin.num.factor),categories)
}


conf.mat.MAP.clm.probit.symmetric2<-list()
for(i in 1:10)
{
  conf.mat.MAP.clm.probit.symmetric2[[i]]<-mat.square(table(pred.MAP.clm.probit.symmetric2[[i]],test[[i]]$V6.bin.num.factor),categories)
}


conf.mat.MAP.clm.probit.equidistant<-list()
for(i in 1:10)
{
  conf.mat.MAP.clm.probit.equidistant[[i]]<-mat.square(table(pred.MAP.clm.probit.equidistant[[i]],test[[i]]$V6.bin.num.factor),categories)
}

###

conf.mat.MAP.clm.loglog.flexible<-list()
for(i in 1:10)
{
  conf.mat.MAP.clm.loglog.flexible[[i]]<-mat.square(table(pred.MAP.clm.loglog.flexible[[i]],test[[i]]$V6.bin.num.factor),categories)
}


conf.mat.MAP.clm.loglog.symmetric<-list()
for(i in 1:10)
{
  conf.mat.MAP.clm.loglog.symmetric[[i]]<-mat.square(table(pred.MAP.clm.loglog.symmetric[[i]],test[[i]]$V6.bin.num.factor),categories)
}


conf.mat.MAP.clm.loglog.symmetric2<-list()
for(i in 1:10)
{
  conf.mat.MAP.clm.loglog.symmetric2[[i]]<-mat.square(table(pred.MAP.clm.loglog.symmetric2[[i]],test[[i]]$V6.bin.num.factor),categories)
}


conf.mat.MAP.clm.loglog.equidistant<-list()
for(i in 1:10)
{
  conf.mat.MAP.clm.loglog.equidistant[[i]]<-mat.square(table(pred.MAP.clm.loglog.equidistant[[i]],test[[i]]$V6.bin.num.factor),categories)
}

###

conf.mat.MAP.clm.cloglog.flexible<-list()
for(i in 1:10)
{
  conf.mat.MAP.clm.cloglog.flexible[[i]]<-mat.square(table(pred.MAP.clm.cloglog.flexible[[i]],test[[i]]$V6.bin.num.factor),categories)
}


conf.mat.MAP.clm.cloglog.symmetric<-list()
for(i in 1:10)
{
  conf.mat.MAP.clm.cloglog.symmetric[[i]]<-mat.square(table(pred.MAP.clm.cloglog.symmetric[[i]],test[[i]]$V6.bin.num.factor),categories)
}


conf.mat.MAP.clm.cloglog.symmetric2<-list()
for(i in 1:10)
{
  conf.mat.MAP.clm.cloglog.symmetric2[[i]]<-mat.square(table(pred.MAP.clm.cloglog.symmetric2[[i]],test[[i]]$V6.bin.num.factor),categories)
}


conf.mat.MAP.clm.cloglog.equidistant<-list()
for(i in 1:10)
{
  conf.mat.MAP.clm.cloglog.equidistant[[i]]<-mat.square(table(pred.MAP.clm.cloglog.equidistant[[i]],test[[i]]$V6.bin.num.factor),categories)
}

###

conf.mat.MAP.clm.cauchit.flexible<-list()
for(i in 1:10)
{
  conf.mat.MAP.clm.cauchit.flexible[[i]]<-mat.square(table(pred.MAP.clm.cauchit.flexible[[i]],test[[i]]$V6.bin.num.factor),categories)
}


conf.mat.MAP.clm.cauchit.symmetric<-list()
for(i in 1:10)
{
  conf.mat.MAP.clm.cauchit.symmetric[[i]]<-mat.square(table(pred.MAP.clm.cauchit.symmetric[[i]],test[[i]]$V6.bin.num.factor),categories)
}


conf.mat.MAP.clm.cauchit.symmetric2<-list()
for(i in 1:10)
{
  conf.mat.MAP.clm.cauchit.symmetric2[[i]]<-mat.square(table(pred.MAP.clm.cauchit.symmetric2[[i]],test[[i]]$V6.bin.num.factor),categories)
}


conf.mat.MAP.clm.cauchit.equidistant<-list()
for(i in 1:10)
{
  conf.mat.MAP.clm.cauchit.equidistant[[i]]<-mat.square(table(pred.MAP.clm.cauchit.equidistant[[i]],test[[i]]$V6.bin.num.factor),categories)
}


# MAE 

MAE.MAP.clm.logit.flexible<-vector()
for (i in 1:10)
{
  MAE.MAP.clm.logit.flexible[i]<-MAE(conf.mat.MAP.clm.logit.flexible[[i]])
}

MAE.MAP.clm.probit.flexible<-vector()
for (i in 1:10)
{
  MAE.MAP.clm.probit.flexible[i]<-MAE(conf.mat.MAP.clm.probit.flexible[[i]])
}

MAE.MAP.clm.loglog.flexible<-vector()
for (i in 1:10)
{
  MAE.MAP.clm.loglog.flexible[i]<-MAE(conf.mat.MAP.clm.loglog.flexible[[i]])
}

MAE.MAP.clm.cloglog.flexible<-vector()
for (i in 1:10)
{
  MAE.MAP.clm.cloglog.flexible[i]<-MAE(conf.mat.MAP.clm.cloglog.flexible[[i]])
}

MAE.MAP.clm.cauchit.flexible<-vector()
for (i in 1:10)
{
  MAE.MAP.clm.cauchit.flexible[i]<-MAE(conf.mat.MAP.clm.cauchit.flexible[[i]])
}

###

MAE.MAP.clm.logit.symmetric<-vector()
for (i in 1:10)
{
  MAE.MAP.clm.logit.symmetric[i]<-MAE(conf.mat.MAP.clm.logit.symmetric[[i]])
}

MAE.MAP.clm.probit.symmetric<-vector()
for (i in 1:10)
{
  MAE.MAP.clm.probit.symmetric[i]<-MAE(conf.mat.MAP.clm.probit.symmetric[[i]])
}

MAE.MAP.clm.loglog.symmetric<-vector()
for (i in 1:10)
{
  MAE.MAP.clm.loglog.symmetric[i]<-MAE(conf.mat.MAP.clm.loglog.symmetric[[i]])
}

MAE.MAP.clm.cloglog.symmetric<-vector()
for (i in 1:10)
{
  MAE.MAP.clm.cloglog.symmetric[i]<-MAE(conf.mat.MAP.clm.cloglog.symmetric[[i]])
}

MAE.MAP.clm.cauchit.symmetric<-vector()
for (i in 1:10)
{
  MAE.MAP.clm.cauchit.symmetric[i]<-MAE(conf.mat.MAP.clm.cauchit.symmetric[[i]])
}

###

MAE.MAP.clm.logit.symmetric2<-vector()
for (i in 1:10)
{
  MAE.MAP.clm.logit.symmetric2[i]<-MAE(conf.mat.MAP.clm.logit.symmetric2[[i]])
}

MAE.MAP.clm.probit.symmetric2<-vector()
for (i in 1:10)
{
  MAE.MAP.clm.probit.symmetric2[i]<-MAE(conf.mat.MAP.clm.probit.symmetric2[[i]])
}

MAE.MAP.clm.loglog.symmetric2<-vector()
for (i in 1:10)
{
  MAE.MAP.clm.loglog.symmetric2[i]<-MAE(conf.mat.MAP.clm.loglog.symmetric2[[i]])
}

MAE.MAP.clm.cloglog.symmetric2<-vector()
for (i in 1:10)
{
  MAE.MAP.clm.cloglog.symmetric2[i]<-MAE(conf.mat.MAP.clm.cloglog.symmetric2[[i]])
}

MAE.MAP.clm.cauchit.symmetric2<-vector()
for (i in 1:10)
{
  MAE.MAP.clm.cauchit.symmetric2[i]<-MAE(conf.mat.MAP.clm.cauchit.symmetric2[[i]])
}

###

MAE.MAP.clm.logit.equidistant<-vector()
for (i in 1:10)
{
  MAE.MAP.clm.logit.equidistant[i]<-MAE(conf.mat.MAP.clm.logit.equidistant[[i]])
}

MAE.MAP.clm.probit.equidistant<-vector()
for (i in 1:10)
{
  MAE.MAP.clm.probit.equidistant[i]<-MAE(conf.mat.MAP.clm.probit.equidistant[[i]])
}

MAE.MAP.clm.loglog.equidistant<-vector()
for (i in 1:10)
{
  MAE.MAP.clm.loglog.equidistant[i]<-MAE(conf.mat.MAP.clm.loglog.equidistant[[i]])
}

MAE.MAP.clm.cloglog.equidistant<-vector()
for (i in 1:10)
{
  MAE.MAP.clm.cloglog.equidistant[i]<-MAE(conf.mat.MAP.clm.cloglog.equidistant[[i]])
}

MAE.MAP.clm.cauchit.equidistant<-vector()
for (i in 1:10)
{
  MAE.MAP.clm.cauchit.equidistant[i]<-MAE(conf.mat.MAP.clm.cauchit.equidistant[[i]])
}



################################################################################
### predictions with Ord.MAP criterion:
###

pred.Ord.MAP.clm.logit.flexible<-list()
for (i in 1:10)
{ pred.Ord.MAP.clm.logit.flexible[[i]]<-vector()
for (j in 1:dim(pred.test.clm.logit.flexible[[i]])[1])
{h<-min(which(cumsum(pred.test.clm.logit.flexible[[i]][j,])>=0.5))
pred.Ord.MAP.clm.logit.flexible[[i]][j]<-categories[h]
}
}

pred.Ord.MAP.clm.logit.symmetric<-list()
for (i in 1:10)
{ pred.Ord.MAP.clm.logit.symmetric[[i]]<-vector()
for (j in 1:dim(pred.test.clm.logit.symmetric[[i]])[1])
{h<-min(which(cumsum(pred.test.clm.logit.symmetric[[i]][j,])>=0.5))
pred.Ord.MAP.clm.logit.symmetric[[i]][j]<-categories[h]
}
}

pred.Ord.MAP.clm.logit.symmetric2<-list()
for (i in 1:10)
{ pred.Ord.MAP.clm.logit.symmetric2[[i]]<-vector()
for (j in 1:dim(pred.test.clm.logit.symmetric2[[i]])[1])
{h<-min(which(cumsum(pred.test.clm.logit.symmetric2[[i]][j,])>=0.5))
pred.Ord.MAP.clm.logit.symmetric2[[i]][j]<-categories[h]
}
}

pred.Ord.MAP.clm.logit.equidistant<-list()
for (i in 1:10)
{ pred.Ord.MAP.clm.logit.equidistant[[i]]<-vector()
for (j in 1:dim(pred.test.clm.logit.equidistant[[i]])[1])
{h<-min(which(cumsum(pred.test.clm.logit.equidistant[[i]][j,])>=0.5))
pred.Ord.MAP.clm.logit.equidistant[[i]][j]<-categories[h]
}
}

###

pred.Ord.MAP.clm.probit.flexible<-list()
for (i in 1:10)
{ pred.Ord.MAP.clm.probit.flexible[[i]]<-vector()
for (j in 1:dim(pred.test.clm.probit.flexible[[i]])[1])
{h<-min(which(cumsum(pred.test.clm.probit.flexible[[i]][j,])>=0.5))
pred.Ord.MAP.clm.probit.flexible[[i]][j]<-categories[h]
}
}

pred.Ord.MAP.clm.probit.symmetric<-list()
for (i in 1:10)
{ pred.Ord.MAP.clm.probit.symmetric[[i]]<-vector()
for (j in 1:dim(pred.test.clm.probit.symmetric[[i]])[1])
{h<-min(which(cumsum(pred.test.clm.probit.symmetric[[i]][j,])>=0.5))
pred.Ord.MAP.clm.probit.symmetric[[i]][j]<-categories[h]
}
}

pred.Ord.MAP.clm.probit.symmetric2<-list()
for (i in 1:10)
{ pred.Ord.MAP.clm.probit.symmetric2[[i]]<-vector()
for (j in 1:dim(pred.test.clm.probit.symmetric2[[i]])[1])
{h<-min(which(cumsum(pred.test.clm.probit.symmetric2[[i]][j,])>=0.5))
pred.Ord.MAP.clm.probit.symmetric2[[i]][j]<-categories[h]
}
}

pred.Ord.MAP.clm.probit.equidistant<-list()
for (i in 1:10)
{ pred.Ord.MAP.clm.probit.equidistant[[i]]<-vector()
for (j in 1:dim(pred.test.clm.probit.equidistant[[i]])[1])
{h<-min(which(cumsum(pred.test.clm.probit.equidistant[[i]][j,])>=0.5))
pred.Ord.MAP.clm.probit.equidistant[[i]][j]<-categories[h]
}
}

###

pred.Ord.MAP.clm.loglog.flexible<-list()
for (i in 1:10)
{ pred.Ord.MAP.clm.loglog.flexible[[i]]<-vector()
for (j in 1:dim(pred.test.clm.loglog.flexible[[i]])[1])
{h<-min(which(cumsum(pred.test.clm.loglog.flexible[[i]][j,])>=0.5))
pred.Ord.MAP.clm.loglog.flexible[[i]][j]<-categories[h]
}
}

pred.Ord.MAP.clm.loglog.symmetric<-list()
for (i in 1:10)
{ pred.Ord.MAP.clm.loglog.symmetric[[i]]<-vector()
for (j in 1:dim(pred.test.clm.loglog.symmetric[[i]])[1])
{h<-min(which(cumsum(pred.test.clm.loglog.symmetric[[i]][j,])>=0.5))
pred.Ord.MAP.clm.loglog.symmetric[[i]][j]<-categories[h]
}
}

pred.Ord.MAP.clm.loglog.symmetric2<-list()
for (i in 1:10)
{ pred.Ord.MAP.clm.loglog.symmetric2[[i]]<-vector()
for (j in 1:dim(pred.test.clm.loglog.symmetric2[[i]])[1])
{h<-min(which(cumsum(pred.test.clm.loglog.symmetric2[[i]][j,])>=0.5))
pred.Ord.MAP.clm.loglog.symmetric2[[i]][j]<-categories[h]
}
}

pred.Ord.MAP.clm.loglog.equidistant<-list()
for (i in 1:10)
{ pred.Ord.MAP.clm.loglog.equidistant[[i]]<-vector()
for (j in 1:dim(pred.test.clm.loglog.equidistant[[i]])[1])
{h<-min(which(cumsum(pred.test.clm.loglog.equidistant[[i]][j,])>=0.5))
pred.Ord.MAP.clm.loglog.equidistant[[i]][j]<-categories[h]
}
}

###

pred.Ord.MAP.clm.cloglog.flexible<-list()
for (i in 1:10)
{ pred.Ord.MAP.clm.cloglog.flexible[[i]]<-vector()
for (j in 1:dim(pred.test.clm.cloglog.flexible[[i]])[1])
{h<-min(which(cumsum(pred.test.clm.cloglog.flexible[[i]][j,])>=0.5))
pred.Ord.MAP.clm.cloglog.flexible[[i]][j]<-categories[h]
}
}

pred.Ord.MAP.clm.cloglog.symmetric<-list()
for (i in 1:10)
{ pred.Ord.MAP.clm.cloglog.symmetric[[i]]<-vector()
for (j in 1:dim(pred.test.clm.cloglog.symmetric[[i]])[1])
{h<-min(which(cumsum(pred.test.clm.cloglog.symmetric[[i]][j,])>=0.5))
pred.Ord.MAP.clm.cloglog.symmetric[[i]][j]<-categories[h]
}
}

pred.Ord.MAP.clm.cloglog.symmetric2<-list()
for (i in 1:10)
{ pred.Ord.MAP.clm.cloglog.symmetric2[[i]]<-vector()
for (j in 1:dim(pred.test.clm.cloglog.symmetric2[[i]])[1])
{h<-min(which(cumsum(pred.test.clm.cloglog.symmetric2[[i]][j,])>=0.5))
pred.Ord.MAP.clm.cloglog.symmetric2[[i]][j]<-categories[h]
}
}

pred.Ord.MAP.clm.cloglog.equidistant<-list()
for (i in 1:10)
{ pred.Ord.MAP.clm.cloglog.equidistant[[i]]<-vector()
for (j in 1:dim(pred.test.clm.cloglog.equidistant[[i]])[1])
{h<-min(which(cumsum(pred.test.clm.cloglog.equidistant[[i]][j,])>=0.5))
pred.Ord.MAP.clm.cloglog.equidistant[[i]][j]<-categories[h]
}
}

###

pred.Ord.MAP.clm.cauchit.flexible<-list()
for (i in 1:10)
{ pred.Ord.MAP.clm.cauchit.flexible[[i]]<-vector()
for (j in 1:dim(pred.test.clm.cauchit.flexible[[i]])[1])
{h<-min(which(cumsum(pred.test.clm.cauchit.flexible[[i]][j,])>=0.5))
pred.Ord.MAP.clm.cauchit.flexible[[i]][j]<-categories[h]
}
}

pred.Ord.MAP.clm.cauchit.symmetric<-list()
for (i in 1:10)
{ pred.Ord.MAP.clm.cauchit.symmetric[[i]]<-vector()
for (j in 1:dim(pred.test.clm.cauchit.symmetric[[i]])[1])
{h<-min(which(cumsum(pred.test.clm.cauchit.symmetric[[i]][j,])>=0.5))
pred.Ord.MAP.clm.cauchit.symmetric[[i]][j]<-categories[h]
}
}

pred.Ord.MAP.clm.cauchit.symmetric2<-list()
for (i in 1:10)
{ pred.Ord.MAP.clm.cauchit.symmetric2[[i]]<-vector()
for (j in 1:dim(pred.test.clm.cauchit.symmetric2[[i]])[1])
{h<-min(which(cumsum(pred.test.clm.cauchit.symmetric2[[i]][j,])>=0.5))
pred.Ord.MAP.clm.cauchit.symmetric2[[i]][j]<-categories[h]
}
}

pred.Ord.MAP.clm.cauchit.equidistant<-list()
for (i in 1:10)
{ pred.Ord.MAP.clm.cauchit.equidistant[[i]]<-vector()
for (j in 1:dim(pred.test.clm.cauchit.equidistant[[i]])[1])
{h<-min(which(cumsum(pred.test.clm.cauchit.equidistant[[i]][j,])>=0.5))
pred.Ord.MAP.clm.cauchit.equidistant[[i]][j]<-categories[h]
}
}

###
# confusion matrices

conf.mat.Ord.MAP.clm.logit.flexible<-list()
for(i in 1:10)
{
  conf.mat.Ord.MAP.clm.logit.flexible[[i]]<-mat.square(table(pred.Ord.MAP.clm.logit.flexible[[i]],test[[i]]$V6.bin.num.factor),categories)
}

conf.mat.Ord.MAP.clm.logit.symmetric<-list()
for(i in 1:10)
{
  conf.mat.Ord.MAP.clm.logit.symmetric[[i]]<-mat.square(table(pred.Ord.MAP.clm.logit.symmetric[[i]],test[[i]]$V6.bin.num.factor),categories)
}

conf.mat.Ord.MAP.clm.logit.symmetric2<-list()
for(i in 1:10)
{
  conf.mat.Ord.MAP.clm.logit.symmetric2[[i]]<-mat.square(table(pred.Ord.MAP.clm.logit.symmetric2[[i]],test[[i]]$V6.bin.num.factor),categories)
}

conf.mat.Ord.MAP.clm.logit.equidistant<-list()
for(i in 1:10)
{
  conf.mat.Ord.MAP.clm.logit.equidistant[[i]]<-mat.square(table(pred.Ord.MAP.clm.logit.equidistant[[i]],test[[i]]$V6.bin.num.factor),categories)
}

###

conf.mat.Ord.MAP.clm.probit.flexible<-list()
for(i in 1:10)
{
  conf.mat.Ord.MAP.clm.probit.flexible[[i]]<-mat.square(table(pred.Ord.MAP.clm.probit.flexible[[i]],test[[i]]$V6.bin.num.factor),categories)
}

conf.mat.Ord.MAP.clm.probit.symmetric<-list()
for(i in 1:10)
{
  conf.mat.Ord.MAP.clm.probit.symmetric[[i]]<-mat.square(table(pred.Ord.MAP.clm.probit.symmetric[[i]],test[[i]]$V6.bin.num.factor),categories)
}

conf.mat.Ord.MAP.clm.probit.symmetric2<-list()
for(i in 1:10)
{
  conf.mat.Ord.MAP.clm.probit.symmetric2[[i]]<-mat.square(table(pred.Ord.MAP.clm.probit.symmetric2[[i]],test[[i]]$V6.bin.num.factor),categories)
}

conf.mat.Ord.MAP.clm.probit.equidistant<-list()
for(i in 1:10)
{
  conf.mat.Ord.MAP.clm.probit.equidistant[[i]]<-mat.square(table(pred.Ord.MAP.clm.probit.equidistant[[i]],test[[i]]$V6.bin.num.factor),categories)
}

###

conf.mat.Ord.MAP.clm.loglog.flexible<-list()
for(i in 1:10)
{
  conf.mat.Ord.MAP.clm.loglog.flexible[[i]]<-mat.square(table(pred.Ord.MAP.clm.loglog.flexible[[i]],test[[i]]$V6.bin.num.factor),categories)
}

conf.mat.Ord.MAP.clm.loglog.symmetric<-list()
for(i in 1:10)
{
  conf.mat.Ord.MAP.clm.loglog.symmetric[[i]]<-mat.square(table(pred.Ord.MAP.clm.loglog.symmetric[[i]],test[[i]]$V6.bin.num.factor),categories)
}

conf.mat.Ord.MAP.clm.loglog.symmetric2<-list()
for(i in 1:10)
{
  conf.mat.Ord.MAP.clm.loglog.symmetric2[[i]]<-mat.square(table(pred.Ord.MAP.clm.loglog.symmetric2[[i]],test[[i]]$V6.bin.num.factor),categories)
}

conf.mat.Ord.MAP.clm.loglog.equidistant<-list()
for(i in 1:10)
{
  conf.mat.Ord.MAP.clm.loglog.equidistant[[i]]<-mat.square(table(pred.Ord.MAP.clm.loglog.equidistant[[i]],test[[i]]$V6.bin.num.factor),categories)
}

###


conf.mat.Ord.MAP.clm.cloglog.flexible<-list()
for(i in 1:10)
{
  conf.mat.Ord.MAP.clm.cloglog.flexible[[i]]<-mat.square(table(pred.Ord.MAP.clm.cloglog.flexible[[i]],test[[i]]$V6.bin.num.factor),categories)
}

conf.mat.Ord.MAP.clm.cloglog.symmetric<-list()
for(i in 1:10)
{
  conf.mat.Ord.MAP.clm.cloglog.symmetric[[i]]<-mat.square(table(pred.Ord.MAP.clm.cloglog.symmetric[[i]],test[[i]]$V6.bin.num.factor),categories)
}

conf.mat.Ord.MAP.clm.cloglog.symmetric2<-list()
for(i in 1:10)
{
  conf.mat.Ord.MAP.clm.cloglog.symmetric2[[i]]<-mat.square(table(pred.Ord.MAP.clm.cloglog.symmetric2[[i]],test[[i]]$V6.bin.num.factor),categories)
}

conf.mat.Ord.MAP.clm.cloglog.equidistant<-list()
for(i in 1:10)
{
  conf.mat.Ord.MAP.clm.cloglog.equidistant[[i]]<-mat.square(table(pred.Ord.MAP.clm.cloglog.equidistant[[i]],test[[i]]$V6.bin.num.factor),categories)
}

###

conf.mat.Ord.MAP.clm.cauchit.flexible<-list()
for(i in 1:10)
{
  conf.mat.Ord.MAP.clm.cauchit.flexible[[i]]<-mat.square(table(pred.Ord.MAP.clm.cauchit.flexible[[i]],test[[i]]$V6.bin.num.factor),categories)
}

conf.mat.Ord.MAP.clm.cauchit.symmetric<-list()
for(i in 1:10)
{
  conf.mat.Ord.MAP.clm.cauchit.symmetric[[i]]<-mat.square(table(pred.Ord.MAP.clm.cauchit.symmetric[[i]],test[[i]]$V6.bin.num.factor),categories)
}

conf.mat.Ord.MAP.clm.cauchit.symmetric2<-list()
for(i in 1:10)
{
  conf.mat.Ord.MAP.clm.cauchit.symmetric2[[i]]<-mat.square(table(pred.Ord.MAP.clm.cauchit.symmetric2[[i]],test[[i]]$V6.bin.num.factor),categories)
}

conf.mat.Ord.MAP.clm.cauchit.equidistant<-list()
for(i in 1:10)
{
  conf.mat.Ord.MAP.clm.cauchit.equidistant[[i]]<-mat.square(table(pred.Ord.MAP.clm.cauchit.equidistant[[i]],test[[i]]$V6.bin.num.factor),categories)
}

###

# MAE 

MAE.Ord.MAP.clm.logit.flexible<-vector()
for (i in 1:10)
{
  MAE.Ord.MAP.clm.logit.flexible[i]<-MAE(conf.mat.Ord.MAP.clm.logit.flexible[[i]])
}

MAE.Ord.MAP.clm.logit.symmetric<-vector()
for (i in 1:10)
{
  MAE.Ord.MAP.clm.logit.symmetric[i]<-MAE(conf.mat.Ord.MAP.clm.logit.symmetric[[i]])
}

MAE.Ord.MAP.clm.logit.symmetric2<-vector()
for (i in 1:10)
{
  MAE.Ord.MAP.clm.logit.symmetric2[i]<-MAE(conf.mat.Ord.MAP.clm.logit.symmetric2[[i]])
}

MAE.Ord.MAP.clm.logit.equidistant<-vector()
for (i in 1:10)
{
  MAE.Ord.MAP.clm.logit.equidistant[i]<-MAE(conf.mat.Ord.MAP.clm.logit.equidistant[[i]])
}

###

MAE.Ord.MAP.clm.probit.flexible<-vector()
for (i in 1:10)
{
  MAE.Ord.MAP.clm.probit.flexible[i]<-MAE(conf.mat.Ord.MAP.clm.probit.flexible[[i]])
}

MAE.Ord.MAP.clm.probit.symmetric<-vector()
for (i in 1:10)
{
  MAE.Ord.MAP.clm.probit.symmetric[i]<-MAE(conf.mat.Ord.MAP.clm.probit.symmetric[[i]])
}

MAE.Ord.MAP.clm.probit.symmetric2<-vector()
for (i in 1:10)
{
  MAE.Ord.MAP.clm.probit.symmetric2[i]<-MAE(conf.mat.Ord.MAP.clm.probit.symmetric2[[i]])
}

MAE.Ord.MAP.clm.probit.equidistant<-vector()
for (i in 1:10)
{
  MAE.Ord.MAP.clm.probit.equidistant[i]<-MAE(conf.mat.Ord.MAP.clm.probit.equidistant[[i]])
}

###


MAE.Ord.MAP.clm.loglog.flexible<-vector()
for (i in 1:10)
{
  MAE.Ord.MAP.clm.loglog.flexible[i]<-MAE(conf.mat.Ord.MAP.clm.loglog.flexible[[i]])
}

MAE.Ord.MAP.clm.loglog.symmetric<-vector()
for (i in 1:10)
{
  MAE.Ord.MAP.clm.loglog.symmetric[i]<-MAE(conf.mat.Ord.MAP.clm.loglog.symmetric[[i]])
}

MAE.Ord.MAP.clm.loglog.symmetric2<-vector()
for (i in 1:10)
{
  MAE.Ord.MAP.clm.loglog.symmetric2[i]<-MAE(conf.mat.Ord.MAP.clm.loglog.symmetric2[[i]])
}

MAE.Ord.MAP.clm.loglog.equidistant<-vector()
for (i in 1:10)
{
  MAE.Ord.MAP.clm.loglog.equidistant[i]<-MAE(conf.mat.Ord.MAP.clm.loglog.equidistant[[i]])
}

###


MAE.Ord.MAP.clm.cloglog.flexible<-vector()
for (i in 1:10)
{
  MAE.Ord.MAP.clm.cloglog.flexible[i]<-MAE(conf.mat.Ord.MAP.clm.cloglog.flexible[[i]])
}

MAE.Ord.MAP.clm.cloglog.symmetric<-vector()
for (i in 1:10)
{
  MAE.Ord.MAP.clm.cloglog.symmetric[i]<-MAE(conf.mat.Ord.MAP.clm.cloglog.symmetric[[i]])
}

MAE.Ord.MAP.clm.cloglog.symmetric2<-vector()
for (i in 1:10)
{
  MAE.Ord.MAP.clm.cloglog.symmetric2[i]<-MAE(conf.mat.Ord.MAP.clm.cloglog.symmetric2[[i]])
}

MAE.Ord.MAP.clm.cloglog.equidistant<-vector()
for (i in 1:10)
{
  MAE.Ord.MAP.clm.cloglog.equidistant[i]<-MAE(conf.mat.Ord.MAP.clm.cloglog.equidistant[[i]])
}

###


MAE.Ord.MAP.clm.cauchit.flexible<-vector()
for (i in 1:10)
{
  MAE.Ord.MAP.clm.cauchit.flexible[i]<-MAE(conf.mat.Ord.MAP.clm.cauchit.flexible[[i]])
}

MAE.Ord.MAP.clm.cauchit.symmetric<-vector()
for (i in 1:10)
{
  MAE.Ord.MAP.clm.cauchit.symmetric[i]<-MAE(conf.mat.Ord.MAP.clm.cauchit.symmetric[[i]])
}

MAE.Ord.MAP.clm.cauchit.symmetric2<-vector()
for (i in 1:10)
{
  MAE.Ord.MAP.clm.cauchit.symmetric2[i]<-MAE(conf.mat.Ord.MAP.clm.cauchit.symmetric2[[i]])
}

MAE.Ord.MAP.clm.cauchit.equidistant<-vector()
for (i in 1:10)
{
  MAE.Ord.MAP.clm.cauchit.equidistant[i]<-MAE(conf.mat.Ord.MAP.clm.cauchit.equidistant[[i]])
}


###
################################################################################
###. Comparison MAP vs Ord.MAP
###
###.  Tests MAE

p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.logit.flexible<-wilcox.test(MAE.MAP.clm.logit.flexible,MAE.Ord.MAP.clm.logit.flexible,paired = TRUE,alternative="greater")$p.value
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.logit.symmetric<-wilcox.test(MAE.MAP.clm.logit.symmetric,MAE.Ord.MAP.clm.logit.symmetric,paired = TRUE,alternative="greater")$p.value
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.logit.symmetric2<-wilcox.test(MAE.MAP.clm.logit.symmetric2,MAE.Ord.MAP.clm.logit.symmetric2,paired = TRUE,alternative="greater")$p.value
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.logit.equidistant<-wilcox.test(MAE.MAP.clm.logit.equidistant,MAE.Ord.MAP.clm.logit.equidistant,paired = TRUE,alternative="greater")$p.value

p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.logit.flexible.less<-wilcox.test(MAE.MAP.clm.logit.flexible,MAE.Ord.MAP.clm.logit.flexible,paired = TRUE,alternative="less")$p.value
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.logit.symmetric.less<-wilcox.test(MAE.MAP.clm.logit.symmetric,MAE.Ord.MAP.clm.logit.symmetric,paired = TRUE,alternative="less")$p.value
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.logit.symmetric2.less<-wilcox.test(MAE.MAP.clm.logit.symmetric2,MAE.Ord.MAP.clm.logit.symmetric2,paired = TRUE,alternative="less")$p.value
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.logit.equidistant.less<-wilcox.test(MAE.MAP.clm.logit.equidistant,MAE.Ord.MAP.clm.logit.equidistant,paired = TRUE,alternative="less")$p.value

###

p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.probit.flexible<-wilcox.test(MAE.MAP.clm.probit.flexible,MAE.Ord.MAP.clm.probit.flexible,paired = TRUE,alternative="greater")$p.value
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.probit.symmetric<-wilcox.test(MAE.MAP.clm.probit.symmetric,MAE.Ord.MAP.clm.probit.symmetric,paired = TRUE,alternative="greater")$p.value
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.probit.symmetric2<-wilcox.test(MAE.MAP.clm.probit.symmetric2,MAE.Ord.MAP.clm.probit.symmetric2,paired = TRUE,alternative="greater")$p.value
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.probit.equidistant<-wilcox.test(MAE.MAP.clm.probit.equidistant,MAE.Ord.MAP.clm.probit.equidistant,paired = TRUE,alternative="greater")$p.value

p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.probit.flexible.less<-wilcox.test(MAE.MAP.clm.probit.flexible,MAE.Ord.MAP.clm.probit.flexible,paired = TRUE,alternative="less")$p.value
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.probit.symmetric.less<-wilcox.test(MAE.MAP.clm.probit.symmetric,MAE.Ord.MAP.clm.probit.symmetric,paired = TRUE,alternative="less")$p.value
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.probit.symmetric2.less<-wilcox.test(MAE.MAP.clm.probit.symmetric2,MAE.Ord.MAP.clm.probit.symmetric2,paired = TRUE,alternative="less")$p.value
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.probit.equidistant.less<-wilcox.test(MAE.MAP.clm.probit.equidistant,MAE.Ord.MAP.clm.probit.equidistant,paired = TRUE,alternative="less")$p.value


###

p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.loglog.flexible<-wilcox.test(MAE.MAP.clm.loglog.flexible,MAE.Ord.MAP.clm.loglog.flexible,paired = TRUE,alternative="greater")$p.value
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.loglog.symmetric<-wilcox.test(MAE.MAP.clm.loglog.symmetric,MAE.Ord.MAP.clm.loglog.symmetric,paired = TRUE,alternative="greater")$p.value
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.loglog.symmetric2<-wilcox.test(MAE.MAP.clm.loglog.symmetric2,MAE.Ord.MAP.clm.loglog.symmetric2,paired = TRUE,alternative="greater")$p.value
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.loglog.equidistant<-wilcox.test(MAE.MAP.clm.loglog.equidistant,MAE.Ord.MAP.clm.loglog.equidistant,paired = TRUE,alternative="greater")$p.value

p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.loglog.flexible.less<-wilcox.test(MAE.MAP.clm.loglog.flexible,MAE.Ord.MAP.clm.loglog.flexible,paired = TRUE,alternative="less")$p.value
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.loglog.symmetric.less<-wilcox.test(MAE.MAP.clm.loglog.symmetric,MAE.Ord.MAP.clm.loglog.symmetric,paired = TRUE,alternative="less")$p.value
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.loglog.symmetric2.less<-wilcox.test(MAE.MAP.clm.loglog.symmetric2,MAE.Ord.MAP.clm.loglog.symmetric2,paired = TRUE,alternative="less")$p.value
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.loglog.equidistant.less<-wilcox.test(MAE.MAP.clm.loglog.equidistant,MAE.Ord.MAP.clm.loglog.equidistant,paired = TRUE,alternative="less")$p.value

###

p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.cloglog.flexible<-wilcox.test(MAE.MAP.clm.cloglog.flexible,MAE.Ord.MAP.clm.cloglog.flexible,paired = TRUE,alternative="greater")$p.value
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.cloglog.symmetric<-wilcox.test(MAE.MAP.clm.cloglog.symmetric,MAE.Ord.MAP.clm.cloglog.symmetric,paired = TRUE,alternative="greater")$p.value
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.cloglog.symmetric2<-wilcox.test(MAE.MAP.clm.cloglog.symmetric2,MAE.Ord.MAP.clm.cloglog.symmetric2,paired = TRUE,alternative="greater")$p.value
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.cloglog.equidistant<-wilcox.test(MAE.MAP.clm.cloglog.equidistant,MAE.Ord.MAP.clm.cloglog.equidistant,paired = TRUE,alternative="greater")$p.value

p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.cloglog.flexible.less<-wilcox.test(MAE.MAP.clm.cloglog.flexible,MAE.Ord.MAP.clm.cloglog.flexible,paired = TRUE,alternative="less")$p.value
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.cloglog.symmetric.less<-wilcox.test(MAE.MAP.clm.cloglog.symmetric,MAE.Ord.MAP.clm.cloglog.symmetric,paired = TRUE,alternative="less")$p.value
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.cloglog.symmetric2.less<-wilcox.test(MAE.MAP.clm.cloglog.symmetric2,MAE.Ord.MAP.clm.cloglog.symmetric2,paired = TRUE,alternative="less")$p.value
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.cloglog.equidistant.less<-wilcox.test(MAE.MAP.clm.cloglog.equidistant,MAE.Ord.MAP.clm.cloglog.equidistant,paired = TRUE,alternative="les")$p.value

###

p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.cauchit.flexible<-wilcox.test(MAE.MAP.clm.cauchit.flexible,MAE.Ord.MAP.clm.cauchit.flexible,paired = TRUE,alternative="greater")$p.value
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.cauchit.symmetric<-wilcox.test(MAE.MAP.clm.cauchit.symmetric,MAE.Ord.MAP.clm.cauchit.symmetric,paired = TRUE,alternative="greater")$p.value
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.cauchit.symmetric2<-wilcox.test(MAE.MAP.clm.cauchit.symmetric2,MAE.Ord.MAP.clm.cauchit.symmetric2,paired = TRUE,alternative="greater")$p.value
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.cauchit.equidistant<-wilcox.test(MAE.MAP.clm.cauchit.equidistant,MAE.Ord.MAP.clm.cauchit.equidistant,paired = TRUE,alternative="greater")$p.value

p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.cauchit.flexible.less<-wilcox.test(MAE.MAP.clm.cauchit.flexible,MAE.Ord.MAP.clm.cauchit.flexible,paired = TRUE,alternative="less")$p.value
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.cauchit.symmetric.less<-wilcox.test(MAE.MAP.clm.cauchit.symmetric,MAE.Ord.MAP.clm.cauchit.symmetric,paired = TRUE,alternative="less")$p.value
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.cauchit.symmetric2.less<-wilcox.test(MAE.MAP.clm.cauchit.symmetric2,MAE.Ord.MAP.clm.cauchit.symmetric2,paired = TRUE,alternative="less")$p.value
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.cauchit.equidistant.less<-wilcox.test(MAE.MAP.clm.cauchit.equidistant,MAE.Ord.MAP.clm.cauchit.equidistant,paired = TRUE,alternative="less")$p.value

###

p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.logit.flexible
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.logit.symmetric
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.logit.symmetric2
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.logit.equidistant

p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.logit.flexible.less
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.logit.symmetric.less
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.logit.symmetric2.less
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.logit.equidistant.less

# 

p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.probit.flexible
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.probit.symmetric
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.probit.symmetric2
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.probit.equidistant

p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.probit.flexible.less
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.probit.symmetric.less
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.probit.symmetric2.less
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.probit.equidistant.less

#

p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.loglog.flexible
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.loglog.symmetric
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.loglog.symmetric2
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.loglog.equidistant

p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.loglog.flexible.less
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.loglog.symmetric.less
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.loglog.symmetric2.less
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.loglog.equidistant.less

#

p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.cloglog.flexible
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.cloglog.symmetric
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.cloglog.symmetric2
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.cloglog.equidistant

p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.cloglog.flexible.less
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.cloglog.symmetric.less
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.cloglog.symmetric2.less
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.cloglog.equidistant.less

#

p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.cauchit.flexible
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.cauchit.symmetric
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.cauchit.symmetric2
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.cauchit.equidistant

p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.cauchit.flexible.less
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.cauchit.symmetric.less
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.cauchit.symmetric2.less
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.cauchit.equidistant.less


###############################
###############################
###   Boxplots   MAE

MAE.clm.flexible.results<-as.data.frame(cbind(MAE.MAP.clm.logit.flexible-MAE.Ord.MAP.clm.logit.flexible,
                                              #MAE.MAP.clm.probit.flexible-MAE.Ord.MAP.clm.probit.flexible,
                                              #MAE.MAP.clm.loglog.flexible-MAE.Ord.MAP.clm.loglog.flexible,
                                              #MAE.MAP.clm.cloglog.flexible-MAE.Ord.MAP.clm.cloglog.flexible,
                                              MAE.MAP.clm.cauchit.flexible-MAE.Ord.MAP.clm.cauchit.flexible))

boxplot(MAE.clm.flexible.results,
        main=paste("(d) Parkinson dataset: output V6 (6 classes). Threshold flexible"), 
        xlab=" ", ylab="MAE increment (MAP - Ord.MAP)",names=c("logistic", 
                                                               #"probit","loglog","cloglog",
                                                               "cauchit"),
        col="white")
abline(h = 0, col = 2, lwd = 2, lty=2)
# Dots
stripchart(MAE.clm.flexible.results,
           #method = "jitter",
           pch = 19,
           #col = 1:5,
           col = 1:2,
           vertical = TRUE,
           add = TRUE)

###

MAE.clm.symmetric.results<-as.data.frame(cbind(MAE.MAP.clm.logit.symmetric-MAE.Ord.MAP.clm.logit.symmetric,
                                               #MAE.MAP.clm.probit.symmetric-MAE.Ord.MAP.clm.probit.symmetric,
                                               #MAE.MAP.clm.loglog.symmetric-MAE.Ord.MAP.clm.loglog.symmetric,
                                               #MAE.MAP.clm.cloglog.symmetric-MAE.Ord.MAP.clm.cloglog.symmetric,
                                               MAE.MAP.clm.cauchit.symmetric-MAE.Ord.MAP.clm.cauchit.symmetric))

boxplot(MAE.clm.symmetric.results,
        main=paste("(d) Parkinson dataset: output V6 (6 classes). Threshold symmetric"), 
        xlab=" ", ylab="MAE increment (MAP - Ord.MAP)",names=c("logistic", 
                                                               #"probit","loglog","cloglog",
                                                               "cauchit"),
        col="white")
abline(h = 0, col = 2, lwd = 2, lty=2)
# Dots
stripchart(MAE.clm.symmetric.results,
           #method = "jitter",
           pch = 19,
           #col = 1:5,
           col = 1:2,
           vertical = TRUE,
           add = TRUE)

###

MAE.clm.symmetric2.results<-as.data.frame(cbind(MAE.MAP.clm.logit.symmetric2-MAE.Ord.MAP.clm.logit.symmetric2,
                                                #MAE.MAP.clm.probit.symmetric2-MAE.Ord.MAP.clm.probit.symmetric2,
                                                #MAE.MAP.clm.loglog.symmetric2-MAE.Ord.MAP.clm.loglog.symmetric2,
                                                #MAE.MAP.clm.cloglog.symmetric2-MAE.Ord.MAP.clm.cloglog.symmetric2,
                                                MAE.MAP.clm.cauchit.symmetric2-MAE.Ord.MAP.clm.cauchit.symmetric2))

boxplot(MAE.clm.symmetric2.results,
        main=paste("(d) Parkinson dataset: output V6 (6 classes). Threshold symmetric2"), 
        xlab=" ", ylab="MAE increment (MAP - Ord.MAP)",names=c("logistic", 
                                                               #"probit","loglog","cloglog",
                                                               "cauchit"),
        col="white")
abline(h = 0, col = 2, lwd = 2, lty=2)
# Dots
stripchart(MAE.clm.symmetric2.results,
           #method = "jitter",
           pch = 19,
           #col = 1:5,
           col = 1:2,
           vertical = TRUE,
           add = TRUE)

###

MAE.clm.equidistant.results<-as.data.frame(cbind(MAE.MAP.clm.logit.equidistant-MAE.Ord.MAP.clm.logit.equidistant,
                                                 #MAE.MAP.clm.probit.equidistant-MAE.Ord.MAP.clm.probit.equidistant,
                                                 #MAE.MAP.clm.loglog.equidistant-MAE.Ord.MAP.clm.loglog.equidistant,
                                                 #MAE.MAP.clm.cloglog.equidistant-MAE.Ord.MAP.clm.cloglog.equidistant,
                                                 MAE.MAP.clm.cauchit.equidistant-MAE.Ord.MAP.clm.cauchit.equidistant))

boxplot(MAE.clm.equidistant.results,
        main=paste("(d) Parkinson dataset: output V6 (6 classes). Threshold equidistant"), 
        xlab=" ", ylab="MAE increment (MAP - Ord.MAP)",names=c("logistic", 
                                                               #"probit","loglog","cloglog",
                                                               "cauchit"),
        col="white")
abline(h = 0, col = 2, lwd = 2, lty=2)
# Dots
stripchart(MAE.clm.equidistant.results,
           #method = "jitter",
           pch = 19,
           #col = 1:5,
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
# clm, from "ordinal" package, rfits cumulative link models (CLMs) such as the propotional odds model. 
# The model allows for various link functions and structured thresholds that restricts 
# the thresholds or cut-points to be e.g.,
# equidistant or symmetrically arranged around the central threshold(s). 
#
# clm(formula, scale, nominal, data, weights, start, subset, doFit = TRUE,
#     na.action, contrasts, model = TRUE, control=list(),
#     link = c("logit", "probit", "cloglog", "loglog", "cauchit",
#              "Aranda-Ordaz", "log-gamma"),
#     threshold = c("flexible", "symmetric", "symmetric2", "equidistant"), ...)

model.clm.logit.flexible<-list()
for (i in 1:10)
{model.clm.logit.flexible[[i]] <- clm(Class ~ . - importance - population, data = training[[i]],
                                      link="logit", threshold="flexible")}


model.clm.logit.symmetric<-list()
for (i in 1:10)
{model.clm.logit.symmetric[[i]] <- clm(Class ~ . - importance - population, data = training[[i]],
                                       link="logit", threshold="symmetric")}


model.clm.logit.symmetric2<-list()
for (i in 1:10)
{model.clm.logit.symmetric2[[i]] <- clm(Class ~ . - importance - population, data = training[[i]],
                                        link="logit", threshold="symmetric2")}

model.clm.logit.equidistant<-list()
for (i in 1:10)
{model.clm.logit.equidistant[[i]] <- clm(Class ~ . - importance - population, data = training[[i]],
                                         link="logit", threshold="equidistant")}


###

model.clm.probit.flexible<-list()
for (i in 1:10)
{model.clm.probit.flexible[[i]] <- clm(Class ~ . - importance - population, data = training[[i]],
                                       link="probit", threshold="flexible")}


model.clm.probit.symmetric<-list()
for (i in 1:10)
{model.clm.probit.symmetric[[i]] <- clm(Class ~ . - importance - population, data = training[[i]],
                                        link="probit", threshold="symmetric")}


model.clm.probit.symmetric2<-list()
for (i in 1:10)
{model.clm.probit.symmetric2[[i]] <- clm(Class ~ . - importance - population, data = training[[i]],
                                         link="probit", threshold="symmetric2")}

model.clm.probit.equidistant<-list()
for (i in 1:10)
{model.clm.probit.equidistant[[i]] <- clm(Class ~ . - importance - population, data = training[[i]],
                                          link="probit", threshold="equidistant")}


###

model.clm.loglog.flexible<-list()
for (i in 1:10)
{model.clm.loglog.flexible[[i]] <- clm(Class ~ . - importance - population, data = training[[i]],
                                       link="loglog", threshold="flexible")}


model.clm.loglog.symmetric<-list()
for (i in 1:10)
{model.clm.loglog.symmetric[[i]] <- clm(Class ~ . - importance - population, data = training[[i]],
                                        link="loglog", threshold="symmetric")}


model.clm.loglog.symmetric2<-list()
for (i in 1:10)
{model.clm.loglog.symmetric2[[i]] <- clm(Class ~ . - importance - population, data = training[[i]],
                                         link="loglog", threshold="symmetric2")}

model.clm.loglog.equidistant<-list()
for (i in 1:10)
{model.clm.loglog.equidistant[[i]] <- clm(Class ~ . - importance - population, data = training[[i]],
                                          link="loglog", threshold="equidistant")}


###

model.clm.cloglog.flexible<-list()
for (i in 1:10)
{model.clm.cloglog.flexible[[i]] <- clm(Class ~ . - importance - population, data = training[[i]],
                                        link="cloglog", threshold="flexible")}


model.clm.cloglog.symmetric<-list()
for (i in 1:10)
{model.clm.cloglog.symmetric[[i]] <- clm(Class ~ . - importance - population, data = training[[i]],
                                         link="cloglog", threshold="symmetric")}


model.clm.cloglog.symmetric2<-list()
for (i in 1:10)
{model.clm.cloglog.symmetric2[[i]] <- clm(Class ~ . - importance - population, data = training[[i]],
                                          link="cloglog", threshold="symmetric2")}

model.clm.cloglog.equidistant<-list()
for (i in 1:10)
{model.clm.cloglog.equidistant[[i]] <- clm(Class ~ . - importance - population, data = training[[i]],
                                           link="cloglog", threshold="equidistant")}


###

model.clm.cauchit.flexible<-list()
for (i in 1:10)
{model.clm.cauchit.flexible[[i]] <- clm(Class ~ . - importance - population, data = training[[i]],
                                        link="cauchit", threshold="flexible")}


model.clm.cauchit.symmetric<-list()
for (i in 1:10)
{model.clm.cauchit.symmetric[[i]] <- clm(Class ~ . - importance - population, data = training[[i]],
                                         link="cauchit", threshold="symmetric")}


model.clm.cauchit.symmetric2<-list()
for (i in 1:10)
{model.clm.cauchit.symmetric2[[i]] <- clm(Class ~ . - importance - population, data = training[[i]],
                                          link="cauchit", threshold="symmetric2")}

model.clm.cauchit.equidistant<-list()
for (i in 1:10)
{model.clm.cauchit.equidistant[[i]] <- clm(Class ~ . - importance - population, data = training[[i]],
                                           link="cauchit", threshold="equidistant")}



################################################################################
#
# Predictions of any of the test set
# 
#

pred.test.clm.logit.flexible<-list()
for (i in 1:10)
{pred.test.clm.logit.flexible[[i]] <- predict(model.clm.logit.flexible[[i]],test[[i]][, ],type = "p")$fit
}

pred.test.clm.logit.symmetric<-list()
for (i in 1:10)
{pred.test.clm.logit.symmetric[[i]] <- predict(model.clm.logit.symmetric[[i]],test[[i]][, ],type = "p")$fit
}

pred.test.clm.logit.symmetric2<-list()
for (i in 1:10)
{pred.test.clm.logit.symmetric2[[i]] <- predict(model.clm.logit.symmetric2[[i]],test[[i]][, ],type = "p")$fit
}

pred.test.clm.logit.equidistant<-list()
for (i in 1:10)
{pred.test.clm.logit.equidistant[[i]] <- predict(model.clm.logit.equidistant[[i]],test[[i]][, ],type = "p")$fit
}

###

pred.test.clm.probit.flexible<-list()
for (i in 1:10)
{pred.test.clm.probit.flexible[[i]] <- predict(model.clm.probit.flexible[[i]],test[[i]][, ],type = "p")$fit
}

pred.test.clm.probit.symmetric<-list()
for (i in 1:10)
{pred.test.clm.probit.symmetric[[i]] <- predict(model.clm.probit.symmetric[[i]],test[[i]][, ],type = "p")$fit
}

pred.test.clm.probit.symmetric2<-list()
for (i in 1:10)
{pred.test.clm.probit.symmetric2[[i]] <- predict(model.clm.probit.symmetric2[[i]],test[[i]][, ],type = "p")$fit
}

pred.test.clm.probit.equidistant<-list()
for (i in 1:10)
{pred.test.clm.probit.equidistant[[i]] <- predict(model.clm.probit.equidistant[[i]],test[[i]][, ],type = "p")$fit
}

###

pred.test.clm.loglog.flexible<-list()
for (i in 1:10)
{pred.test.clm.loglog.flexible[[i]] <- predict(model.clm.loglog.flexible[[i]],test[[i]][, ],type = "p")$fit
}

pred.test.clm.loglog.symmetric<-list()
for (i in 1:10)
{pred.test.clm.loglog.symmetric[[i]] <- predict(model.clm.loglog.symmetric[[i]],test[[i]][,],type = "p")$fit
}

pred.test.clm.loglog.symmetric2<-list()
for (i in 1:10)
{pred.test.clm.loglog.symmetric2[[i]] <- predict(model.clm.loglog.symmetric2[[i]],test[[i]][,],type = "p")$fit
}

pred.test.clm.loglog.equidistant<-list()
for (i in 1:10)
{pred.test.clm.loglog.equidistant[[i]] <- predict(model.clm.loglog.equidistant[[i]],test[[i]][, ],type = "p")$fit
}

###


pred.test.clm.cloglog.flexible<-list()
for (i in 1:10)
{pred.test.clm.cloglog.flexible[[i]] <- predict(model.clm.cloglog.flexible[[i]],test[[i]][,],type = "p")$fit
}

pred.test.clm.cloglog.symmetric<-list()
for (i in 1:10)
{pred.test.clm.cloglog.symmetric[[i]] <- predict(model.clm.cloglog.symmetric[[i]],test[[i]][,],type = "p")$fit
}

pred.test.clm.cloglog.symmetric2<-list()
for (i in 1:10)
{pred.test.clm.cloglog.symmetric2[[i]] <- predict(model.clm.cloglog.symmetric2[[i]],test[[i]][,],type = "p")$fit
}

pred.test.clm.cloglog.equidistant<-list()
for (i in 1:10)
{pred.test.clm.cloglog.equidistant[[i]] <- predict(model.clm.cloglog.equidistant[[i]],test[[i]][,],type = "p")$fit
}

###


pred.test.clm.cauchit.flexible<-list()
for (i in 1:10)
{pred.test.clm.cauchit.flexible[[i]] <- predict(model.clm.cauchit.flexible[[i]],test[[i]][,],type = "p")$fit
}

pred.test.clm.cauchit.symmetric<-list()
for (i in 1:10)
{pred.test.clm.cauchit.symmetric[[i]] <- predict(model.clm.cauchit.symmetric[[i]],test[[i]][,],type = "p")$fit
}

pred.test.clm.cauchit.symmetric2<-list()
for (i in 1:10)
{pred.test.clm.cauchit.symmetric2[[i]] <- predict(model.clm.cauchit.symmetric2[[i]],test[[i]][,],type = "p")$fit
}

pred.test.clm.cauchit.equidistant<-list()
for (i in 1:10)
{pred.test.clm.cauchit.equidistant[[i]] <- predict(model.clm.cauchit.equidistant[[i]],test[[i]][,],type = "p")$fit
}

###

################################################################################
### predictions with standard MAP criterion:
###

pred.MAP.clm.logit.flexible<-list()
for (i in 1:10)
{pred.MAP.clm.logit.flexible[[i]]<-vector()
for (j in 1:dim(pred.test.clm.logit.flexible[[i]])[1])
{pred.MAP.clm.logit.flexible[[i]][j]<-categories[which.max(pred.test.clm.logit.flexible[[i]][j,])]
}
}

pred.MAP.clm.logit.symmetric<-list()
for (i in 1:10)
{pred.MAP.clm.logit.symmetric[[i]]<-vector()
for (j in 1:dim(pred.test.clm.logit.symmetric[[i]])[1])
{pred.MAP.clm.logit.symmetric[[i]][j]<-categories[which.max(pred.test.clm.logit.symmetric[[i]][j,])]
}
}

pred.MAP.clm.logit.symmetric2<-list()
for (i in 1:10)
{pred.MAP.clm.logit.symmetric2[[i]]<-vector()
for (j in 1:dim(pred.test.clm.logit.symmetric2[[i]])[1])
{pred.MAP.clm.logit.symmetric2[[i]][j]<-categories[which.max(pred.test.clm.logit.symmetric2[[i]][j,])]
}
}

pred.MAP.clm.logit.equidistant<-list()
for (i in 1:10)
{pred.MAP.clm.logit.equidistant[[i]]<-vector()
for (j in 1:dim(pred.test.clm.logit.equidistant[[i]])[1])
{pred.MAP.clm.logit.equidistant[[i]][j]<-categories[which.max(pred.test.clm.logit.equidistant[[i]][j,])]
}
}

###


pred.MAP.clm.probit.flexible<-list()
for (i in 1:10)
{pred.MAP.clm.probit.flexible[[i]]<-vector()
for (j in 1:dim(pred.test.clm.probit.flexible[[i]])[1])
{pred.MAP.clm.probit.flexible[[i]][j]<-categories[which.max(pred.test.clm.probit.flexible[[i]][j,])]
}
}

pred.MAP.clm.probit.symmetric<-list()
for (i in 1:10)
{pred.MAP.clm.probit.symmetric[[i]]<-vector()
for (j in 1:dim(pred.test.clm.probit.symmetric[[i]])[1])
{pred.MAP.clm.probit.symmetric[[i]][j]<-categories[which.max(pred.test.clm.probit.symmetric[[i]][j,])]
}
}

pred.MAP.clm.probit.symmetric2<-list()
for (i in 1:10)
{pred.MAP.clm.probit.symmetric2[[i]]<-vector()
for (j in 1:dim(pred.test.clm.probit.symmetric2[[i]])[1])
{pred.MAP.clm.probit.symmetric2[[i]][j]<-categories[which.max(pred.test.clm.probit.symmetric2[[i]][j,])]
}
}

pred.MAP.clm.probit.equidistant<-list()
for (i in 1:10)
{pred.MAP.clm.probit.equidistant[[i]]<-vector()
for (j in 1:dim(pred.test.clm.probit.equidistant[[i]])[1])
{pred.MAP.clm.probit.equidistant[[i]][j]<-categories[which.max(pred.test.clm.probit.equidistant[[i]][j,])]
}
}

###


pred.MAP.clm.loglog.flexible<-list()
for (i in 1:10)
{pred.MAP.clm.loglog.flexible[[i]]<-vector()
for (j in 1:dim(pred.test.clm.loglog.flexible[[i]])[1])
{pred.MAP.clm.loglog.flexible[[i]][j]<-categories[which.max(pred.test.clm.loglog.flexible[[i]][j,])]
}
}

pred.MAP.clm.loglog.symmetric<-list()
for (i in 1:10)
{pred.MAP.clm.loglog.symmetric[[i]]<-vector()
for (j in 1:dim(pred.test.clm.loglog.symmetric[[i]])[1])
{pred.MAP.clm.loglog.symmetric[[i]][j]<-categories[which.max(pred.test.clm.loglog.symmetric[[i]][j,])]
}
}

pred.MAP.clm.loglog.symmetric2<-list()
for (i in 1:10)
{pred.MAP.clm.loglog.symmetric2[[i]]<-vector()
for (j in 1:dim(pred.test.clm.loglog.symmetric2[[i]])[1])
{pred.MAP.clm.loglog.symmetric2[[i]][j]<-categories[which.max(pred.test.clm.loglog.symmetric2[[i]][j,])]
}
}

pred.MAP.clm.loglog.equidistant<-list()
for (i in 1:10)
{pred.MAP.clm.loglog.equidistant[[i]]<-vector()
for (j in 1:dim(pred.test.clm.loglog.equidistant[[i]])[1])
{pred.MAP.clm.loglog.equidistant[[i]][j]<-categories[which.max(pred.test.clm.loglog.equidistant[[i]][j,])]
}
}

###


pred.MAP.clm.cloglog.flexible<-list()
for (i in 1:10)
{pred.MAP.clm.cloglog.flexible[[i]]<-vector()
for (j in 1:dim(pred.test.clm.cloglog.flexible[[i]])[1])
{pred.MAP.clm.cloglog.flexible[[i]][j]<-categories[which.max(pred.test.clm.cloglog.flexible[[i]][j,])]
}
}

pred.MAP.clm.cloglog.symmetric<-list()
for (i in 1:10)
{pred.MAP.clm.cloglog.symmetric[[i]]<-vector()
for (j in 1:dim(pred.test.clm.cloglog.symmetric[[i]])[1])
{pred.MAP.clm.cloglog.symmetric[[i]][j]<-categories[which.max(pred.test.clm.cloglog.symmetric[[i]][j,])]
}
}

pred.MAP.clm.cloglog.symmetric2<-list()
for (i in 1:10)
{pred.MAP.clm.cloglog.symmetric2[[i]]<-vector()
for (j in 1:dim(pred.test.clm.cloglog.symmetric2[[i]])[1])
{pred.MAP.clm.cloglog.symmetric2[[i]][j]<-categories[which.max(pred.test.clm.cloglog.symmetric2[[i]][j,])]
}
}

pred.MAP.clm.cloglog.equidistant<-list()
for (i in 1:10)
{pred.MAP.clm.cloglog.equidistant[[i]]<-vector()
for (j in 1:dim(pred.test.clm.cloglog.equidistant[[i]])[1])
{pred.MAP.clm.cloglog.equidistant[[i]][j]<-categories[which.max(pred.test.clm.cloglog.equidistant[[i]][j,])]
}
}

###


pred.MAP.clm.cauchit.flexible<-list()
for (i in 1:10)
{pred.MAP.clm.cauchit.flexible[[i]]<-vector()
for (j in 1:dim(pred.test.clm.cauchit.flexible[[i]])[1])
{pred.MAP.clm.cauchit.flexible[[i]][j]<-categories[which.max(pred.test.clm.cauchit.flexible[[i]][j,])]
}
}

pred.MAP.clm.cauchit.symmetric<-list()
for (i in 1:10)
{pred.MAP.clm.cauchit.symmetric[[i]]<-vector()
for (j in 1:dim(pred.test.clm.cauchit.symmetric[[i]])[1])
{pred.MAP.clm.cauchit.symmetric[[i]][j]<-categories[which.max(pred.test.clm.cauchit.symmetric[[i]][j,])]
}
}

pred.MAP.clm.cauchit.symmetric2<-list()
for (i in 1:10)
{pred.MAP.clm.cauchit.symmetric2[[i]]<-vector()
for (j in 1:dim(pred.test.clm.cauchit.symmetric2[[i]])[1])
{pred.MAP.clm.cauchit.symmetric2[[i]][j]<-categories[which.max(pred.test.clm.cauchit.symmetric2[[i]][j,])]
}
}

pred.MAP.clm.cauchit.equidistant<-list()
for (i in 1:10)
{pred.MAP.clm.cauchit.equidistant[[i]]<-vector()
for (j in 1:dim(pred.test.clm.cauchit.equidistant[[i]])[1])
{pred.MAP.clm.cauchit.equidistant[[i]][j]<-categories[which.max(pred.test.clm.cauchit.equidistant[[i]][j,])]
}
}

# confusion matrices


conf.mat.MAP.clm.logit.flexible<-list()
for(i in 1:10)
{
  conf.mat.MAP.clm.logit.flexible[[i]]<-mat.square(table(pred.MAP.clm.logit.flexible[[i]],test[[i]]$Class),categories)
}


conf.mat.MAP.clm.logit.symmetric<-list()
for(i in 1:10)
{
  conf.mat.MAP.clm.logit.symmetric[[i]]<-mat.square(table(pred.MAP.clm.logit.symmetric[[i]],test[[i]]$Class),categories)
}


conf.mat.MAP.clm.logit.symmetric2<-list()
for(i in 1:10)
{
  conf.mat.MAP.clm.logit.symmetric2[[i]]<-mat.square(table(pred.MAP.clm.logit.symmetric2[[i]],test[[i]]$Class),categories)
}


conf.mat.MAP.clm.logit.equidistant<-list()
for(i in 1:10)
{
  conf.mat.MAP.clm.logit.equidistant[[i]]<-mat.square(table(pred.MAP.clm.logit.equidistant[[i]],test[[i]]$Class),categories)
}

###

conf.mat.MAP.clm.probit.flexible<-list()
for(i in 1:10)
{
  conf.mat.MAP.clm.probit.flexible[[i]]<-mat.square(table(pred.MAP.clm.probit.flexible[[i]],test[[i]]$Class),categories)
}


conf.mat.MAP.clm.probit.symmetric<-list()
for(i in 1:10)
{
  conf.mat.MAP.clm.probit.symmetric[[i]]<-mat.square(table(pred.MAP.clm.probit.symmetric[[i]],test[[i]]$Class),categories)
}


conf.mat.MAP.clm.probit.symmetric2<-list()
for(i in 1:10)
{
  conf.mat.MAP.clm.probit.symmetric2[[i]]<-mat.square(table(pred.MAP.clm.probit.symmetric2[[i]],test[[i]]$Class),categories)
}


conf.mat.MAP.clm.probit.equidistant<-list()
for(i in 1:10)
{
  conf.mat.MAP.clm.probit.equidistant[[i]]<-mat.square(table(pred.MAP.clm.probit.equidistant[[i]],test[[i]]$Class),categories)
}

###

conf.mat.MAP.clm.loglog.flexible<-list()
for(i in 1:10)
{
  conf.mat.MAP.clm.loglog.flexible[[i]]<-mat.square(table(pred.MAP.clm.loglog.flexible[[i]],test[[i]]$Class),categories)
}


conf.mat.MAP.clm.loglog.symmetric<-list()
for(i in 1:10)
{
  conf.mat.MAP.clm.loglog.symmetric[[i]]<-mat.square(table(pred.MAP.clm.loglog.symmetric[[i]],test[[i]]$Class),categories)
}


conf.mat.MAP.clm.loglog.symmetric2<-list()
for(i in 1:10)
{
  conf.mat.MAP.clm.loglog.symmetric2[[i]]<-mat.square(table(pred.MAP.clm.loglog.symmetric2[[i]],test[[i]]$Class),categories)
}


conf.mat.MAP.clm.loglog.equidistant<-list()
for(i in 1:10)
{
  conf.mat.MAP.clm.loglog.equidistant[[i]]<-mat.square(table(pred.MAP.clm.loglog.equidistant[[i]],test[[i]]$Class),categories)
}

###

conf.mat.MAP.clm.cloglog.flexible<-list()
for(i in 1:10)
{
  conf.mat.MAP.clm.cloglog.flexible[[i]]<-mat.square(table(pred.MAP.clm.cloglog.flexible[[i]],test[[i]]$Class),categories)
}


conf.mat.MAP.clm.cloglog.symmetric<-list()
for(i in 1:10)
{
  conf.mat.MAP.clm.cloglog.symmetric[[i]]<-mat.square(table(pred.MAP.clm.cloglog.symmetric[[i]],test[[i]]$Class),categories)
}


conf.mat.MAP.clm.cloglog.symmetric2<-list()
for(i in 1:10)
{
  conf.mat.MAP.clm.cloglog.symmetric2[[i]]<-mat.square(table(pred.MAP.clm.cloglog.symmetric2[[i]],test[[i]]$Class),categories)
}


conf.mat.MAP.clm.cloglog.equidistant<-list()
for(i in 1:10)
{
  conf.mat.MAP.clm.cloglog.equidistant[[i]]<-mat.square(table(pred.MAP.clm.cloglog.equidistant[[i]],test[[i]]$Class),categories)
}

###

conf.mat.MAP.clm.cauchit.flexible<-list()
for(i in 1:10)
{
  conf.mat.MAP.clm.cauchit.flexible[[i]]<-mat.square(table(pred.MAP.clm.cauchit.flexible[[i]],test[[i]]$Class),categories)
}


conf.mat.MAP.clm.cauchit.symmetric<-list()
for(i in 1:10)
{
  conf.mat.MAP.clm.cauchit.symmetric[[i]]<-mat.square(table(pred.MAP.clm.cauchit.symmetric[[i]],test[[i]]$Class),categories)
}


conf.mat.MAP.clm.cauchit.symmetric2<-list()
for(i in 1:10)
{
  conf.mat.MAP.clm.cauchit.symmetric2[[i]]<-mat.square(table(pred.MAP.clm.cauchit.symmetric2[[i]],test[[i]]$Class),categories)
}


conf.mat.MAP.clm.cauchit.equidistant<-list()
for(i in 1:10)
{
  conf.mat.MAP.clm.cauchit.equidistant[[i]]<-mat.square(table(pred.MAP.clm.cauchit.equidistant[[i]],test[[i]]$Class),categories)
}


# MAE 

MAE.MAP.clm.logit.flexible<-vector()
for (i in 1:10)
{
  MAE.MAP.clm.logit.flexible[i]<-MAE(conf.mat.MAP.clm.logit.flexible[[i]])
}

MAE.MAP.clm.probit.flexible<-vector()
for (i in 1:10)
{
  MAE.MAP.clm.probit.flexible[i]<-MAE(conf.mat.MAP.clm.probit.flexible[[i]])
}

MAE.MAP.clm.loglog.flexible<-vector()
for (i in 1:10)
{
  MAE.MAP.clm.loglog.flexible[i]<-MAE(conf.mat.MAP.clm.loglog.flexible[[i]])
}

MAE.MAP.clm.cloglog.flexible<-vector()
for (i in 1:10)
{
  MAE.MAP.clm.cloglog.flexible[i]<-MAE(conf.mat.MAP.clm.cloglog.flexible[[i]])
}

MAE.MAP.clm.cauchit.flexible<-vector()
for (i in 1:10)
{
  MAE.MAP.clm.cauchit.flexible[i]<-MAE(conf.mat.MAP.clm.cauchit.flexible[[i]])
}

###

MAE.MAP.clm.logit.symmetric<-vector()
for (i in 1:10)
{
  MAE.MAP.clm.logit.symmetric[i]<-MAE(conf.mat.MAP.clm.logit.symmetric[[i]])
}

MAE.MAP.clm.probit.symmetric<-vector()
for (i in 1:10)
{
  MAE.MAP.clm.probit.symmetric[i]<-MAE(conf.mat.MAP.clm.probit.symmetric[[i]])
}

MAE.MAP.clm.loglog.symmetric<-vector()
for (i in 1:10)
{
  MAE.MAP.clm.loglog.symmetric[i]<-MAE(conf.mat.MAP.clm.loglog.symmetric[[i]])
}

MAE.MAP.clm.cloglog.symmetric<-vector()
for (i in 1:10)
{
  MAE.MAP.clm.cloglog.symmetric[i]<-MAE(conf.mat.MAP.clm.cloglog.symmetric[[i]])
}

MAE.MAP.clm.cauchit.symmetric<-vector()
for (i in 1:10)
{
  MAE.MAP.clm.cauchit.symmetric[i]<-MAE(conf.mat.MAP.clm.cauchit.symmetric[[i]])
}

###

MAE.MAP.clm.logit.symmetric2<-vector()
for (i in 1:10)
{
  MAE.MAP.clm.logit.symmetric2[i]<-MAE(conf.mat.MAP.clm.logit.symmetric2[[i]])
}

MAE.MAP.clm.probit.symmetric2<-vector()
for (i in 1:10)
{
  MAE.MAP.clm.probit.symmetric2[i]<-MAE(conf.mat.MAP.clm.probit.symmetric2[[i]])
}

MAE.MAP.clm.loglog.symmetric2<-vector()
for (i in 1:10)
{
  MAE.MAP.clm.loglog.symmetric2[i]<-MAE(conf.mat.MAP.clm.loglog.symmetric2[[i]])
}

MAE.MAP.clm.cloglog.symmetric2<-vector()
for (i in 1:10)
{
  MAE.MAP.clm.cloglog.symmetric2[i]<-MAE(conf.mat.MAP.clm.cloglog.symmetric2[[i]])
}

MAE.MAP.clm.cauchit.symmetric2<-vector()
for (i in 1:10)
{
  MAE.MAP.clm.cauchit.symmetric2[i]<-MAE(conf.mat.MAP.clm.cauchit.symmetric2[[i]])
}

###

MAE.MAP.clm.logit.equidistant<-vector()
for (i in 1:10)
{
  MAE.MAP.clm.logit.equidistant[i]<-MAE(conf.mat.MAP.clm.logit.equidistant[[i]])
}

MAE.MAP.clm.probit.equidistant<-vector()
for (i in 1:10)
{
  MAE.MAP.clm.probit.equidistant[i]<-MAE(conf.mat.MAP.clm.probit.equidistant[[i]])
}

MAE.MAP.clm.loglog.equidistant<-vector()
for (i in 1:10)
{
  MAE.MAP.clm.loglog.equidistant[i]<-MAE(conf.mat.MAP.clm.loglog.equidistant[[i]])
}

MAE.MAP.clm.cloglog.equidistant<-vector()
for (i in 1:10)
{
  MAE.MAP.clm.cloglog.equidistant[i]<-MAE(conf.mat.MAP.clm.cloglog.equidistant[[i]])
}

MAE.MAP.clm.cauchit.equidistant<-vector()
for (i in 1:10)
{
  MAE.MAP.clm.cauchit.equidistant[i]<-MAE(conf.mat.MAP.clm.cauchit.equidistant[[i]])
}



################################################################################
### predictions with Ord.MAP criterion:
###

pred.Ord.MAP.clm.logit.flexible<-list()
for (i in 1:10)
{ pred.Ord.MAP.clm.logit.flexible[[i]]<-vector()
for (j in 1:dim(pred.test.clm.logit.flexible[[i]])[1])
{h<-min(which(cumsum(pred.test.clm.logit.flexible[[i]][j,])>=0.5))
pred.Ord.MAP.clm.logit.flexible[[i]][j]<-categories[h]
}
}

pred.Ord.MAP.clm.logit.symmetric<-list()
for (i in 1:10)
{ pred.Ord.MAP.clm.logit.symmetric[[i]]<-vector()
for (j in 1:dim(pred.test.clm.logit.symmetric[[i]])[1])
{h<-min(which(cumsum(pred.test.clm.logit.symmetric[[i]][j,])>=0.5))
pred.Ord.MAP.clm.logit.symmetric[[i]][j]<-categories[h]
}
}

pred.Ord.MAP.clm.logit.symmetric2<-list()
for (i in 1:10)
{ pred.Ord.MAP.clm.logit.symmetric2[[i]]<-vector()
for (j in 1:dim(pred.test.clm.logit.symmetric2[[i]])[1])
{h<-min(which(cumsum(pred.test.clm.logit.symmetric2[[i]][j,])>=0.5))
pred.Ord.MAP.clm.logit.symmetric2[[i]][j]<-categories[h]
}
}

pred.Ord.MAP.clm.logit.equidistant<-list()
for (i in 1:10)
{ pred.Ord.MAP.clm.logit.equidistant[[i]]<-vector()
for (j in 1:dim(pred.test.clm.logit.equidistant[[i]])[1])
{h<-min(which(cumsum(pred.test.clm.logit.equidistant[[i]][j,])>=0.5))
pred.Ord.MAP.clm.logit.equidistant[[i]][j]<-categories[h]
}
}

###

pred.Ord.MAP.clm.probit.flexible<-list()
for (i in 1:10)
{ pred.Ord.MAP.clm.probit.flexible[[i]]<-vector()
for (j in 1:dim(pred.test.clm.probit.flexible[[i]])[1])
{h<-min(which(cumsum(pred.test.clm.probit.flexible[[i]][j,])>=0.5))
pred.Ord.MAP.clm.probit.flexible[[i]][j]<-categories[h]
}
}

pred.Ord.MAP.clm.probit.symmetric<-list()
for (i in 1:10)
{ pred.Ord.MAP.clm.probit.symmetric[[i]]<-vector()
for (j in 1:dim(pred.test.clm.probit.symmetric[[i]])[1])
{h<-min(which(cumsum(pred.test.clm.probit.symmetric[[i]][j,])>=0.5))
pred.Ord.MAP.clm.probit.symmetric[[i]][j]<-categories[h]
}
}

pred.Ord.MAP.clm.probit.symmetric2<-list()
for (i in 1:10)
{ pred.Ord.MAP.clm.probit.symmetric2[[i]]<-vector()
for (j in 1:dim(pred.test.clm.probit.symmetric2[[i]])[1])
{h<-min(which(cumsum(pred.test.clm.probit.symmetric2[[i]][j,])>=0.5))
pred.Ord.MAP.clm.probit.symmetric2[[i]][j]<-categories[h]
}
}

pred.Ord.MAP.clm.probit.equidistant<-list()
for (i in 1:10)
{ pred.Ord.MAP.clm.probit.equidistant[[i]]<-vector()
for (j in 1:dim(pred.test.clm.probit.equidistant[[i]])[1])
{h<-min(which(cumsum(pred.test.clm.probit.equidistant[[i]][j,])>=0.5))
pred.Ord.MAP.clm.probit.equidistant[[i]][j]<-categories[h]
}
}

###

pred.Ord.MAP.clm.loglog.flexible<-list()
for (i in 1:10)
{ pred.Ord.MAP.clm.loglog.flexible[[i]]<-vector()
for (j in 1:dim(pred.test.clm.loglog.flexible[[i]])[1])
{h<-min(which(cumsum(pred.test.clm.loglog.flexible[[i]][j,])>=0.5))
pred.Ord.MAP.clm.loglog.flexible[[i]][j]<-categories[h]
}
}

pred.Ord.MAP.clm.loglog.symmetric<-list()
for (i in 1:10)
{ pred.Ord.MAP.clm.loglog.symmetric[[i]]<-vector()
for (j in 1:dim(pred.test.clm.loglog.symmetric[[i]])[1])
{h<-min(which(cumsum(pred.test.clm.loglog.symmetric[[i]][j,])>=0.5))
pred.Ord.MAP.clm.loglog.symmetric[[i]][j]<-categories[h]
}
}

pred.Ord.MAP.clm.loglog.symmetric2<-list()
for (i in 1:10)
{ pred.Ord.MAP.clm.loglog.symmetric2[[i]]<-vector()
for (j in 1:dim(pred.test.clm.loglog.symmetric2[[i]])[1])
{h<-min(which(cumsum(pred.test.clm.loglog.symmetric2[[i]][j,])>=0.5))
pred.Ord.MAP.clm.loglog.symmetric2[[i]][j]<-categories[h]
}
}

pred.Ord.MAP.clm.loglog.equidistant<-list()
for (i in 1:10)
{ pred.Ord.MAP.clm.loglog.equidistant[[i]]<-vector()
for (j in 1:dim(pred.test.clm.loglog.equidistant[[i]])[1])
{h<-min(which(cumsum(pred.test.clm.loglog.equidistant[[i]][j,])>=0.5))
pred.Ord.MAP.clm.loglog.equidistant[[i]][j]<-categories[h]
}
}

###

pred.Ord.MAP.clm.cloglog.flexible<-list()
for (i in 1:10)
{ pred.Ord.MAP.clm.cloglog.flexible[[i]]<-vector()
for (j in 1:dim(pred.test.clm.cloglog.flexible[[i]])[1])
{h<-min(which(cumsum(pred.test.clm.cloglog.flexible[[i]][j,])>=0.5))
pred.Ord.MAP.clm.cloglog.flexible[[i]][j]<-categories[h]
}
}

pred.Ord.MAP.clm.cloglog.symmetric<-list()
for (i in 1:10)
{ pred.Ord.MAP.clm.cloglog.symmetric[[i]]<-vector()
for (j in 1:dim(pred.test.clm.cloglog.symmetric[[i]])[1])
{h<-min(which(cumsum(pred.test.clm.cloglog.symmetric[[i]][j,])>=0.5))
pred.Ord.MAP.clm.cloglog.symmetric[[i]][j]<-categories[h]
}
}

pred.Ord.MAP.clm.cloglog.symmetric2<-list()
for (i in 1:10)
{ pred.Ord.MAP.clm.cloglog.symmetric2[[i]]<-vector()
for (j in 1:dim(pred.test.clm.cloglog.symmetric2[[i]])[1])
{h<-min(which(cumsum(pred.test.clm.cloglog.symmetric2[[i]][j,])>=0.5))
pred.Ord.MAP.clm.cloglog.symmetric2[[i]][j]<-categories[h]
}
}

pred.Ord.MAP.clm.cloglog.equidistant<-list()
for (i in 1:10)
{ pred.Ord.MAP.clm.cloglog.equidistant[[i]]<-vector()
for (j in 1:dim(pred.test.clm.cloglog.equidistant[[i]])[1])
{h<-min(which(cumsum(pred.test.clm.cloglog.equidistant[[i]][j,])>=0.5))
pred.Ord.MAP.clm.cloglog.equidistant[[i]][j]<-categories[h]
}
}

###

pred.Ord.MAP.clm.cauchit.flexible<-list()
for (i in 1:10)
{ pred.Ord.MAP.clm.cauchit.flexible[[i]]<-vector()
for (j in 1:dim(pred.test.clm.cauchit.flexible[[i]])[1])
{h<-min(which(cumsum(pred.test.clm.cauchit.flexible[[i]][j,])>=0.5))
pred.Ord.MAP.clm.cauchit.flexible[[i]][j]<-categories[h]
}
}

pred.Ord.MAP.clm.cauchit.symmetric<-list()
for (i in 1:10)
{ pred.Ord.MAP.clm.cauchit.symmetric[[i]]<-vector()
for (j in 1:dim(pred.test.clm.cauchit.symmetric[[i]])[1])
{h<-min(which(cumsum(pred.test.clm.cauchit.symmetric[[i]][j,])>=0.5))
pred.Ord.MAP.clm.cauchit.symmetric[[i]][j]<-categories[h]
}
}

pred.Ord.MAP.clm.cauchit.symmetric2<-list()
for (i in 1:10)
{ pred.Ord.MAP.clm.cauchit.symmetric2[[i]]<-vector()
for (j in 1:dim(pred.test.clm.cauchit.symmetric2[[i]])[1])
{h<-min(which(cumsum(pred.test.clm.cauchit.symmetric2[[i]][j,])>=0.5))
pred.Ord.MAP.clm.cauchit.symmetric2[[i]][j]<-categories[h]
}
}

pred.Ord.MAP.clm.cauchit.equidistant<-list()
for (i in 1:10)
{ pred.Ord.MAP.clm.cauchit.equidistant[[i]]<-vector()
for (j in 1:dim(pred.test.clm.cauchit.equidistant[[i]])[1])
{h<-min(which(cumsum(pred.test.clm.cauchit.equidistant[[i]][j,])>=0.5))
pred.Ord.MAP.clm.cauchit.equidistant[[i]][j]<-categories[h]
}
}

###
# confusion matrices

conf.mat.Ord.MAP.clm.logit.flexible<-list()
for(i in 1:10)
{
  conf.mat.Ord.MAP.clm.logit.flexible[[i]]<-mat.square(table(pred.Ord.MAP.clm.logit.flexible[[i]],test[[i]]$Class),categories)
}

conf.mat.Ord.MAP.clm.logit.symmetric<-list()
for(i in 1:10)
{
  conf.mat.Ord.MAP.clm.logit.symmetric[[i]]<-mat.square(table(pred.Ord.MAP.clm.logit.symmetric[[i]],test[[i]]$Class),categories)
}

conf.mat.Ord.MAP.clm.logit.symmetric2<-list()
for(i in 1:10)
{
  conf.mat.Ord.MAP.clm.logit.symmetric2[[i]]<-mat.square(table(pred.Ord.MAP.clm.logit.symmetric2[[i]],test[[i]]$Class),categories)
}

conf.mat.Ord.MAP.clm.logit.equidistant<-list()
for(i in 1:10)
{
  conf.mat.Ord.MAP.clm.logit.equidistant[[i]]<-mat.square(table(pred.Ord.MAP.clm.logit.equidistant[[i]],test[[i]]$Class),categories)
}

###

conf.mat.Ord.MAP.clm.probit.flexible<-list()
for(i in 1:10)
{
  conf.mat.Ord.MAP.clm.probit.flexible[[i]]<-mat.square(table(pred.Ord.MAP.clm.probit.flexible[[i]],test[[i]]$Class),categories)
}

conf.mat.Ord.MAP.clm.probit.symmetric<-list()
for(i in 1:10)
{
  conf.mat.Ord.MAP.clm.probit.symmetric[[i]]<-mat.square(table(pred.Ord.MAP.clm.probit.symmetric[[i]],test[[i]]$Class),categories)
}

conf.mat.Ord.MAP.clm.probit.symmetric2<-list()
for(i in 1:10)
{
  conf.mat.Ord.MAP.clm.probit.symmetric2[[i]]<-mat.square(table(pred.Ord.MAP.clm.probit.symmetric2[[i]],test[[i]]$Class),categories)
}

conf.mat.Ord.MAP.clm.probit.equidistant<-list()
for(i in 1:10)
{
  conf.mat.Ord.MAP.clm.probit.equidistant[[i]]<-mat.square(table(pred.Ord.MAP.clm.probit.equidistant[[i]],test[[i]]$Class),categories)
}

###

conf.mat.Ord.MAP.clm.loglog.flexible<-list()
for(i in 1:10)
{
  conf.mat.Ord.MAP.clm.loglog.flexible[[i]]<-mat.square(table(pred.Ord.MAP.clm.loglog.flexible[[i]],test[[i]]$Class),categories)
}

conf.mat.Ord.MAP.clm.loglog.symmetric<-list()
for(i in 1:10)
{
  conf.mat.Ord.MAP.clm.loglog.symmetric[[i]]<-mat.square(table(pred.Ord.MAP.clm.loglog.symmetric[[i]],test[[i]]$Class),categories)
}

conf.mat.Ord.MAP.clm.loglog.symmetric2<-list()
for(i in 1:10)
{
  conf.mat.Ord.MAP.clm.loglog.symmetric2[[i]]<-mat.square(table(pred.Ord.MAP.clm.loglog.symmetric2[[i]],test[[i]]$Class),categories)
}

conf.mat.Ord.MAP.clm.loglog.equidistant<-list()
for(i in 1:10)
{
  conf.mat.Ord.MAP.clm.loglog.equidistant[[i]]<-mat.square(table(pred.Ord.MAP.clm.loglog.equidistant[[i]],test[[i]]$Class),categories)
}

###


conf.mat.Ord.MAP.clm.cloglog.flexible<-list()
for(i in 1:10)
{
  conf.mat.Ord.MAP.clm.cloglog.flexible[[i]]<-mat.square(table(pred.Ord.MAP.clm.cloglog.flexible[[i]],test[[i]]$Class),categories)
}

conf.mat.Ord.MAP.clm.cloglog.symmetric<-list()
for(i in 1:10)
{
  conf.mat.Ord.MAP.clm.cloglog.symmetric[[i]]<-mat.square(table(pred.Ord.MAP.clm.cloglog.symmetric[[i]],test[[i]]$Class),categories)
}

conf.mat.Ord.MAP.clm.cloglog.symmetric2<-list()
for(i in 1:10)
{
  conf.mat.Ord.MAP.clm.cloglog.symmetric2[[i]]<-mat.square(table(pred.Ord.MAP.clm.cloglog.symmetric2[[i]],test[[i]]$Class),categories)
}

conf.mat.Ord.MAP.clm.cloglog.equidistant<-list()
for(i in 1:10)
{
  conf.mat.Ord.MAP.clm.cloglog.equidistant[[i]]<-mat.square(table(pred.Ord.MAP.clm.cloglog.equidistant[[i]],test[[i]]$Class),categories)
}

###

conf.mat.Ord.MAP.clm.cauchit.flexible<-list()
for(i in 1:10)
{
  conf.mat.Ord.MAP.clm.cauchit.flexible[[i]]<-mat.square(table(pred.Ord.MAP.clm.cauchit.flexible[[i]],test[[i]]$Class),categories)
}

conf.mat.Ord.MAP.clm.cauchit.symmetric<-list()
for(i in 1:10)
{
  conf.mat.Ord.MAP.clm.cauchit.symmetric[[i]]<-mat.square(table(pred.Ord.MAP.clm.cauchit.symmetric[[i]],test[[i]]$Class),categories)
}

conf.mat.Ord.MAP.clm.cauchit.symmetric2<-list()
for(i in 1:10)
{
  conf.mat.Ord.MAP.clm.cauchit.symmetric2[[i]]<-mat.square(table(pred.Ord.MAP.clm.cauchit.symmetric2[[i]],test[[i]]$Class),categories)
}

conf.mat.Ord.MAP.clm.cauchit.equidistant<-list()
for(i in 1:10)
{
  conf.mat.Ord.MAP.clm.cauchit.equidistant[[i]]<-mat.square(table(pred.Ord.MAP.clm.cauchit.equidistant[[i]],test[[i]]$Class),categories)
}

###

# MAE 

MAE.Ord.MAP.clm.logit.flexible<-vector()
for (i in 1:10)
{
  MAE.Ord.MAP.clm.logit.flexible[i]<-MAE(conf.mat.Ord.MAP.clm.logit.flexible[[i]])
}

MAE.Ord.MAP.clm.logit.symmetric<-vector()
for (i in 1:10)
{
  MAE.Ord.MAP.clm.logit.symmetric[i]<-MAE(conf.mat.Ord.MAP.clm.logit.symmetric[[i]])
}

MAE.Ord.MAP.clm.logit.symmetric2<-vector()
for (i in 1:10)
{
  MAE.Ord.MAP.clm.logit.symmetric2[i]<-MAE(conf.mat.Ord.MAP.clm.logit.symmetric2[[i]])
}

MAE.Ord.MAP.clm.logit.equidistant<-vector()
for (i in 1:10)
{
  MAE.Ord.MAP.clm.logit.equidistant[i]<-MAE(conf.mat.Ord.MAP.clm.logit.equidistant[[i]])
}

###

MAE.Ord.MAP.clm.probit.flexible<-vector()
for (i in 1:10)
{
  MAE.Ord.MAP.clm.probit.flexible[i]<-MAE(conf.mat.Ord.MAP.clm.probit.flexible[[i]])
}

MAE.Ord.MAP.clm.probit.symmetric<-vector()
for (i in 1:10)
{
  MAE.Ord.MAP.clm.probit.symmetric[i]<-MAE(conf.mat.Ord.MAP.clm.probit.symmetric[[i]])
}

MAE.Ord.MAP.clm.probit.symmetric2<-vector()
for (i in 1:10)
{
  MAE.Ord.MAP.clm.probit.symmetric2[i]<-MAE(conf.mat.Ord.MAP.clm.probit.symmetric2[[i]])
}

MAE.Ord.MAP.clm.probit.equidistant<-vector()
for (i in 1:10)
{
  MAE.Ord.MAP.clm.probit.equidistant[i]<-MAE(conf.mat.Ord.MAP.clm.probit.equidistant[[i]])
}

###


MAE.Ord.MAP.clm.loglog.flexible<-vector()
for (i in 1:10)
{
  MAE.Ord.MAP.clm.loglog.flexible[i]<-MAE(conf.mat.Ord.MAP.clm.loglog.flexible[[i]])
}

MAE.Ord.MAP.clm.loglog.symmetric<-vector()
for (i in 1:10)
{
  MAE.Ord.MAP.clm.loglog.symmetric[i]<-MAE(conf.mat.Ord.MAP.clm.loglog.symmetric[[i]])
}

MAE.Ord.MAP.clm.loglog.symmetric2<-vector()
for (i in 1:10)
{
  MAE.Ord.MAP.clm.loglog.symmetric2[i]<-MAE(conf.mat.Ord.MAP.clm.loglog.symmetric2[[i]])
}

MAE.Ord.MAP.clm.loglog.equidistant<-vector()
for (i in 1:10)
{
  MAE.Ord.MAP.clm.loglog.equidistant[i]<-MAE(conf.mat.Ord.MAP.clm.loglog.equidistant[[i]])
}

###


MAE.Ord.MAP.clm.cloglog.flexible<-vector()
for (i in 1:10)
{
  MAE.Ord.MAP.clm.cloglog.flexible[i]<-MAE(conf.mat.Ord.MAP.clm.cloglog.flexible[[i]])
}

MAE.Ord.MAP.clm.cloglog.symmetric<-vector()
for (i in 1:10)
{
  MAE.Ord.MAP.clm.cloglog.symmetric[i]<-MAE(conf.mat.Ord.MAP.clm.cloglog.symmetric[[i]])
}

MAE.Ord.MAP.clm.cloglog.symmetric2<-vector()
for (i in 1:10)
{
  MAE.Ord.MAP.clm.cloglog.symmetric2[i]<-MAE(conf.mat.Ord.MAP.clm.cloglog.symmetric2[[i]])
}

MAE.Ord.MAP.clm.cloglog.equidistant<-vector()
for (i in 1:10)
{
  MAE.Ord.MAP.clm.cloglog.equidistant[i]<-MAE(conf.mat.Ord.MAP.clm.cloglog.equidistant[[i]])
}

###


MAE.Ord.MAP.clm.cauchit.flexible<-vector()
for (i in 1:10)
{
  MAE.Ord.MAP.clm.cauchit.flexible[i]<-MAE(conf.mat.Ord.MAP.clm.cauchit.flexible[[i]])
}

MAE.Ord.MAP.clm.cauchit.symmetric<-vector()
for (i in 1:10)
{
  MAE.Ord.MAP.clm.cauchit.symmetric[i]<-MAE(conf.mat.Ord.MAP.clm.cauchit.symmetric[[i]])
}

MAE.Ord.MAP.clm.cauchit.symmetric2<-vector()
for (i in 1:10)
{
  MAE.Ord.MAP.clm.cauchit.symmetric2[i]<-MAE(conf.mat.Ord.MAP.clm.cauchit.symmetric2[[i]])
}

MAE.Ord.MAP.clm.cauchit.equidistant<-vector()
for (i in 1:10)
{
  MAE.Ord.MAP.clm.cauchit.equidistant[i]<-MAE(conf.mat.Ord.MAP.clm.cauchit.equidistant[[i]])
}


###
################################################################################
###. Comparison MAP vs Ord.MAP
###
###.  Tests MAE

p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.logit.flexible<-wilcox.test(MAE.MAP.clm.logit.flexible,MAE.Ord.MAP.clm.logit.flexible,paired = TRUE,alternative="greater")$p.value
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.logit.symmetric<-wilcox.test(MAE.MAP.clm.logit.symmetric,MAE.Ord.MAP.clm.logit.symmetric,paired = TRUE,alternative="greater")$p.value
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.logit.symmetric2<-wilcox.test(MAE.MAP.clm.logit.symmetric2,MAE.Ord.MAP.clm.logit.symmetric2,paired = TRUE,alternative="greater")$p.value
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.logit.equidistant<-wilcox.test(MAE.MAP.clm.logit.equidistant,MAE.Ord.MAP.clm.logit.equidistant,paired = TRUE,alternative="greater")$p.value

p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.logit.flexible.less<-wilcox.test(MAE.MAP.clm.logit.flexible,MAE.Ord.MAP.clm.logit.flexible,paired = TRUE,alternative="less")$p.value
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.logit.symmetric.less<-wilcox.test(MAE.MAP.clm.logit.symmetric,MAE.Ord.MAP.clm.logit.symmetric,paired = TRUE,alternative="less")$p.value
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.logit.symmetric2.less<-wilcox.test(MAE.MAP.clm.logit.symmetric2,MAE.Ord.MAP.clm.logit.symmetric2,paired = TRUE,alternative="less")$p.value
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.logit.equidistant.less<-wilcox.test(MAE.MAP.clm.logit.equidistant,MAE.Ord.MAP.clm.logit.equidistant,paired = TRUE,alternative="less")$p.value

###

p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.probit.flexible<-wilcox.test(MAE.MAP.clm.probit.flexible,MAE.Ord.MAP.clm.probit.flexible,paired = TRUE,alternative="greater")$p.value
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.probit.symmetric<-wilcox.test(MAE.MAP.clm.probit.symmetric,MAE.Ord.MAP.clm.probit.symmetric,paired = TRUE,alternative="greater")$p.value
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.probit.symmetric2<-wilcox.test(MAE.MAP.clm.probit.symmetric2,MAE.Ord.MAP.clm.probit.symmetric2,paired = TRUE,alternative="greater")$p.value
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.probit.equidistant<-wilcox.test(MAE.MAP.clm.probit.equidistant,MAE.Ord.MAP.clm.probit.equidistant,paired = TRUE,alternative="greater")$p.value

p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.probit.flexible.less<-wilcox.test(MAE.MAP.clm.probit.flexible,MAE.Ord.MAP.clm.probit.flexible,paired = TRUE,alternative="less")$p.value
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.probit.symmetric.less<-wilcox.test(MAE.MAP.clm.probit.symmetric,MAE.Ord.MAP.clm.probit.symmetric,paired = TRUE,alternative="less")$p.value
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.probit.symmetric2.less<-wilcox.test(MAE.MAP.clm.probit.symmetric2,MAE.Ord.MAP.clm.probit.symmetric2,paired = TRUE,alternative="less")$p.value
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.probit.equidistant.less<-wilcox.test(MAE.MAP.clm.probit.equidistant,MAE.Ord.MAP.clm.probit.equidistant,paired = TRUE,alternative="less")$p.value


###

p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.loglog.flexible<-wilcox.test(MAE.MAP.clm.loglog.flexible,MAE.Ord.MAP.clm.loglog.flexible,paired = TRUE,alternative="greater")$p.value
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.loglog.symmetric<-wilcox.test(MAE.MAP.clm.loglog.symmetric,MAE.Ord.MAP.clm.loglog.symmetric,paired = TRUE,alternative="greater")$p.value
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.loglog.symmetric2<-wilcox.test(MAE.MAP.clm.loglog.symmetric2,MAE.Ord.MAP.clm.loglog.symmetric2,paired = TRUE,alternative="greater")$p.value
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.loglog.equidistant<-wilcox.test(MAE.MAP.clm.loglog.equidistant,MAE.Ord.MAP.clm.loglog.equidistant,paired = TRUE,alternative="greater")$p.value

p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.loglog.flexible.less<-wilcox.test(MAE.MAP.clm.loglog.flexible,MAE.Ord.MAP.clm.loglog.flexible,paired = TRUE,alternative="less")$p.value
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.loglog.symmetric.less<-wilcox.test(MAE.MAP.clm.loglog.symmetric,MAE.Ord.MAP.clm.loglog.symmetric,paired = TRUE,alternative="less")$p.value
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.loglog.symmetric2.less<-wilcox.test(MAE.MAP.clm.loglog.symmetric2,MAE.Ord.MAP.clm.loglog.symmetric2,paired = TRUE,alternative="less")$p.value
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.loglog.equidistant.less<-wilcox.test(MAE.MAP.clm.loglog.equidistant,MAE.Ord.MAP.clm.loglog.equidistant,paired = TRUE,alternative="less")$p.value

###

p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.cloglog.flexible<-wilcox.test(MAE.MAP.clm.cloglog.flexible,MAE.Ord.MAP.clm.cloglog.flexible,paired = TRUE,alternative="greater")$p.value
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.cloglog.symmetric<-wilcox.test(MAE.MAP.clm.cloglog.symmetric,MAE.Ord.MAP.clm.cloglog.symmetric,paired = TRUE,alternative="greater")$p.value
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.cloglog.symmetric2<-wilcox.test(MAE.MAP.clm.cloglog.symmetric2,MAE.Ord.MAP.clm.cloglog.symmetric2,paired = TRUE,alternative="greater")$p.value
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.cloglog.equidistant<-wilcox.test(MAE.MAP.clm.cloglog.equidistant,MAE.Ord.MAP.clm.cloglog.equidistant,paired = TRUE,alternative="greater")$p.value

p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.cloglog.flexible.less<-wilcox.test(MAE.MAP.clm.cloglog.flexible,MAE.Ord.MAP.clm.cloglog.flexible,paired = TRUE,alternative="less")$p.value
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.cloglog.symmetric.less<-wilcox.test(MAE.MAP.clm.cloglog.symmetric,MAE.Ord.MAP.clm.cloglog.symmetric,paired = TRUE,alternative="less")$p.value
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.cloglog.symmetric2.less<-wilcox.test(MAE.MAP.clm.cloglog.symmetric2,MAE.Ord.MAP.clm.cloglog.symmetric2,paired = TRUE,alternative="less")$p.value
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.cloglog.equidistant.less<-wilcox.test(MAE.MAP.clm.cloglog.equidistant,MAE.Ord.MAP.clm.cloglog.equidistant,paired = TRUE,alternative="les")$p.value

###

p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.cauchit.flexible<-wilcox.test(MAE.MAP.clm.cauchit.flexible,MAE.Ord.MAP.clm.cauchit.flexible,paired = TRUE,alternative="greater")$p.value
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.cauchit.symmetric<-wilcox.test(MAE.MAP.clm.cauchit.symmetric,MAE.Ord.MAP.clm.cauchit.symmetric,paired = TRUE,alternative="greater")$p.value
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.cauchit.symmetric2<-wilcox.test(MAE.MAP.clm.cauchit.symmetric2,MAE.Ord.MAP.clm.cauchit.symmetric2,paired = TRUE,alternative="greater")$p.value
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.cauchit.equidistant<-wilcox.test(MAE.MAP.clm.cauchit.equidistant,MAE.Ord.MAP.clm.cauchit.equidistant,paired = TRUE,alternative="greater")$p.value

p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.cauchit.flexible.less<-wilcox.test(MAE.MAP.clm.cauchit.flexible,MAE.Ord.MAP.clm.cauchit.flexible,paired = TRUE,alternative="less")$p.value
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.cauchit.symmetric.less<-wilcox.test(MAE.MAP.clm.cauchit.symmetric,MAE.Ord.MAP.clm.cauchit.symmetric,paired = TRUE,alternative="less")$p.value
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.cauchit.symmetric2.less<-wilcox.test(MAE.MAP.clm.cauchit.symmetric2,MAE.Ord.MAP.clm.cauchit.symmetric2,paired = TRUE,alternative="less")$p.value
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.cauchit.equidistant.less<-wilcox.test(MAE.MAP.clm.cauchit.equidistant,MAE.Ord.MAP.clm.cauchit.equidistant,paired = TRUE,alternative="less")$p.value

###

p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.logit.flexible
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.logit.symmetric
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.logit.symmetric2
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.logit.equidistant

p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.logit.flexible.less
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.logit.symmetric.less
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.logit.symmetric2.less
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.logit.equidistant.less

# 

p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.probit.flexible
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.probit.symmetric
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.probit.symmetric2
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.probit.equidistant

p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.probit.flexible.less
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.probit.symmetric.less
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.probit.symmetric2.less
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.probit.equidistant.less

#

p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.loglog.flexible
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.loglog.symmetric
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.loglog.symmetric2
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.loglog.equidistant

p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.loglog.flexible.less
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.loglog.symmetric.less
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.loglog.symmetric2.less
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.loglog.equidistant.less

#

p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.cloglog.flexible
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.cloglog.symmetric
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.cloglog.symmetric2
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.cloglog.equidistant

p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.cloglog.flexible.less
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.cloglog.symmetric.less
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.cloglog.symmetric2.less
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.cloglog.equidistant.less

#

p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.cauchit.flexible
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.cauchit.symmetric
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.cauchit.symmetric2
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.cauchit.equidistant

p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.cauchit.flexible.less
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.cauchit.symmetric.less
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.cauchit.symmetric2.less
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.cauchit.equidistant.less


###############################
###############################
###   Boxplots   MAE

MAE.clm.flexible.results<-as.data.frame(cbind(MAE.MAP.clm.logit.flexible-MAE.Ord.MAP.clm.logit.flexible,
                                              MAE.MAP.clm.probit.flexible-MAE.Ord.MAP.clm.probit.flexible,
                                              MAE.MAP.clm.loglog.flexible-MAE.Ord.MAP.clm.loglog.flexible,
                                              MAE.MAP.clm.cloglog.flexible-MAE.Ord.MAP.clm.cloglog.flexible,
                                              MAE.MAP.clm.cauchit.flexible-MAE.Ord.MAP.clm.cauchit.flexible))

boxplot(MAE.clm.flexible.results,
        main=paste("(e) CES11 dataset. Threshold flexible"), 
        xlab=" ", ylab="MAE increment (MAP - Ord.MAP)",names=c("logistic", "probit","loglog","cloglog","cauchit"),
        col="white")
abline(h = 0, col = 2, lwd = 2, lty=2)
# Dots
stripchart(MAE.clm.flexible.results,
           #method = "jitter",
           pch = 19,
           col = 1:5,
           vertical = TRUE,
           add = TRUE)

###

MAE.clm.symmetric.results<-as.data.frame(cbind(MAE.MAP.clm.logit.symmetric-MAE.Ord.MAP.clm.logit.symmetric,
                                               MAE.MAP.clm.probit.symmetric-MAE.Ord.MAP.clm.probit.symmetric,
                                               MAE.MAP.clm.loglog.symmetric-MAE.Ord.MAP.clm.loglog.symmetric,
                                               MAE.MAP.clm.cloglog.symmetric-MAE.Ord.MAP.clm.cloglog.symmetric,
                                               MAE.MAP.clm.cauchit.symmetric-MAE.Ord.MAP.clm.cauchit.symmetric))

boxplot(MAE.clm.symmetric.results,
        main=paste("(e) CES11 dataset. Threshold symmetric"), 
        xlab=" ", ylab="MAE increment (MAP - Ord.MAP)",names=c("logistic", "probit","loglog","cloglog","cauchit"),
        col="white")
abline(h = 0, col = 2, lwd = 2, lty=2)
# Dots
stripchart(MAE.clm.symmetric.results,
           #method = "jitter",
           pch = 19,
           col = 1:5,
           vertical = TRUE,
           add = TRUE)

###

MAE.clm.symmetric2.results<-as.data.frame(cbind(MAE.MAP.clm.logit.symmetric2-MAE.Ord.MAP.clm.logit.symmetric2,
                                                MAE.MAP.clm.probit.symmetric2-MAE.Ord.MAP.clm.probit.symmetric2,
                                                MAE.MAP.clm.loglog.symmetric2-MAE.Ord.MAP.clm.loglog.symmetric2,
                                                MAE.MAP.clm.cloglog.symmetric2-MAE.Ord.MAP.clm.cloglog.symmetric2,
                                                MAE.MAP.clm.cauchit.symmetric2-MAE.Ord.MAP.clm.cauchit.symmetric2))

boxplot(MAE.clm.symmetric2.results,
        main=paste("(e) CES11 dataset. Threshold symmetric2"), 
        xlab=" ", ylab="MAE increment (MAP - Ord.MAP)",names=c("logistic", "probit","loglog","cloglog","cauchit"),
        col="white")
abline(h = 0, col = 2, lwd = 2, lty=2)
# Dots
stripchart(MAE.clm.symmetric2.results,
           #method = "jitter",
           pch = 19,
           col = 1:5,
           vertical = TRUE,
           add = TRUE)

###

MAE.clm.equidistant.results<-as.data.frame(cbind(MAE.MAP.clm.logit.equidistant-MAE.Ord.MAP.clm.logit.equidistant,
                                                 MAE.MAP.clm.probit.equidistant-MAE.Ord.MAP.clm.probit.equidistant,
                                                 MAE.MAP.clm.loglog.equidistant-MAE.Ord.MAP.clm.loglog.equidistant,
                                                 MAE.MAP.clm.cloglog.equidistant-MAE.Ord.MAP.clm.cloglog.equidistant,
                                                 MAE.MAP.clm.cauchit.equidistant-MAE.Ord.MAP.clm.cauchit.equidistant))

boxplot(MAE.clm.equidistant.results,
        main=paste("(e) CES11 dataset. Threshold equidistant"), 
        xlab=" ", ylab="MAE increment (MAP - Ord.MAP)",names=c("logistic", "probit","loglog","cloglog","cauchit"),
        col="white")
abline(h = 0, col = 2, lwd = 2, lty=2)
# Dots
stripchart(MAE.clm.equidistant.results,
           #method = "jitter",
           pch = 19,
           col = 1:5,
           vertical = TRUE,
           add = TRUE)


