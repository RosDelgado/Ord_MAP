
###########################################
###########################################
#
#  EXPERIMENTAL PHASE 
#
# "ordinal" package of R 
#
# cumulative link (mixed) models also known
# as ordered regression models, proportional odds models, proportional
# hazards models for grouped survival times and ordered logit/probit/...
# models. Estimation is via maximum likelihood and mixed models are fitted
# with the Laplace approximation and adaptive Gauss-Hermite quadrature.
#
# data: wine (package "ordinal")
#
# The wine data set is adopted from Randall(1989) and from a factorial experiment 
# on factors determining the bitterness of wine. 
# Two treatment factors (temperature and contact) each have two
# levels. Temperature and contact between juice and skins can be controlled when cruching grapes
# during wine production. Nine judges each assessed wine from two bottles from each of the four
# treatment conditions, hence there are 72 observations in all.
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

p.normality.MAE.clm.logit.flexible<-shapiro.test(MAE.MAP.clm.logit.flexible-MAE.Ord.MAP.clm.logit.flexible)$p.value

if (p.normality.MAE.clm.logit.flexible >= 0.05)
{p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.logit.flexible<-t.test(MAE.MAP.clm.logit.flexible,MAE.Ord.MAP.clm.logit.flexible,paired = TRUE,alternative="greater")$p.value
} else {
  p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.logit.flexible<-wilcox.test(MAE.MAP.clm.logit.flexible,MAE.Ord.MAP.clm.logit.flexible,paired = TRUE,alternative="greater")$p.value
}


p.normality.MAE.clm.logit.symmetric<-shapiro.test(MAE.MAP.clm.logit.symmetric-MAE.Ord.MAP.clm.logit.symmetric)$p.value

if (p.normality.MAE.clm.logit.symmetric >= 0.05)
{p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.logit.symmetric<-t.test(MAE.MAP.clm.logit.symmetric,MAE.Ord.MAP.clm.logit.symmetric,paired = TRUE,alternative="greater")$p.value
} else {
  p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.logit.symmetric<-wilcox.test(MAE.MAP.clm.logit.symmetric,MAE.Ord.MAP.clm.logit.symmetric,paired = TRUE,alternative="greater")$p.value
}


p.normality.MAE.clm.logit.symmetric2<-shapiro.test(MAE.MAP.clm.logit.symmetric2-MAE.Ord.MAP.clm.logit.symmetric2)$p.value

if (p.normality.MAE.clm.logit.symmetric2 >= 0.05)
{p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.logit.symmetric2<-t.test(MAE.MAP.clm.logit.symmetric2,MAE.Ord.MAP.clm.logit.symmetric2,paired = TRUE,alternative="greater")$p.value
} else {
  p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.logit.symmetric2<-wilcox.test(MAE.MAP.clm.logit.symmetric2,MAE.Ord.MAP.clm.logit.symmetric2,paired = TRUE,alternative="greater")$p.value
}


p.normality.MAE.clm.logit.equidistant<-shapiro.test(MAE.MAP.clm.logit.equidistant-MAE.Ord.MAP.clm.logit.equidistant)$p.value

if (p.normality.MAE.clm.logit.equidistant >= 0.05)
{p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.logit.equidistant<-t.test(MAE.MAP.clm.logit.equidistant,MAE.Ord.MAP.clm.logit.equidistant,paired = TRUE,alternative="greater")$p.value
} else {
  p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.logit.equidistant<-wilcox.test(MAE.MAP.clm.logit.equidistant,MAE.Ord.MAP.clm.logit.equidistant,paired = TRUE,alternative="greater")$p.value
}

###

p.normality.MAE.clm.probit.flexible<-shapiro.test(MAE.MAP.clm.probit.flexible-MAE.Ord.MAP.clm.probit.flexible)$p.value

if (p.normality.MAE.clm.probit.flexible >= 0.05)
{p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.probit.flexible<-t.test(MAE.MAP.clm.probit.flexible,MAE.Ord.MAP.clm.probit.flexible,paired = TRUE,alternative="greater")$p.value
} else {
  p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.probit.flexible<-wilcox.test(MAE.MAP.clm.probit.flexible,MAE.Ord.MAP.clm.probit.flexible,paired = TRUE,alternative="greater")$p.value
}


p.normality.MAE.clm.probit.symmetric<-shapiro.test(MAE.MAP.clm.probit.symmetric-MAE.Ord.MAP.clm.probit.symmetric)$p.value

if (p.normality.MAE.clm.probit.symmetric >= 0.05)
{p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.probit.symmetric<-t.test(MAE.MAP.clm.probit.symmetric,MAE.Ord.MAP.clm.probit.symmetric,paired = TRUE,alternative="greater")$p.value
} else {
  p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.probit.symmetric<-wilcox.test(MAE.MAP.clm.probit.symmetric,MAE.Ord.MAP.clm.probit.symmetric,paired = TRUE,alternative="greater")$p.value
}


p.normality.MAE.clm.probit.symmetric2<-shapiro.test(MAE.MAP.clm.probit.symmetric2-MAE.Ord.MAP.clm.probit.symmetric2)$p.value

if (p.normality.MAE.clm.probit.symmetric2 >= 0.05)
{p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.probit.symmetric2<-t.test(MAE.MAP.clm.probit.symmetric2,MAE.Ord.MAP.clm.probit.symmetric2,paired = TRUE,alternative="greater")$p.value
} else {
  p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.probit.symmetric2<-wilcox.test(MAE.MAP.clm.probit.symmetric2,MAE.Ord.MAP.clm.probit.symmetric2,paired = TRUE,alternative="greater")$p.value
}


p.normality.MAE.clm.probit.equidistant<-shapiro.test(MAE.MAP.clm.probit.equidistant-MAE.Ord.MAP.clm.probit.equidistant)$p.value

if (p.normality.MAE.clm.probit.equidistant >= 0.05)
{p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.probit.equidistant<-t.test(MAE.MAP.clm.probit.equidistant,MAE.Ord.MAP.clm.probit.equidistant,paired = TRUE,alternative="greater")$p.value
} else {
  p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.probit.equidistant<-wilcox.test(MAE.MAP.clm.probit.equidistant,MAE.Ord.MAP.clm.probit.equidistant,paired = TRUE,alternative="greater")$p.value
}

###

p.normality.MAE.clm.loglog.flexible<-shapiro.test(MAE.MAP.clm.loglog.flexible-MAE.Ord.MAP.clm.loglog.flexible)$p.value

if (p.normality.MAE.clm.loglog.flexible >= 0.05)
{p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.loglog.flexible<-t.test(MAE.MAP.clm.loglog.flexible,MAE.Ord.MAP.clm.loglog.flexible,paired = TRUE,alternative="greater")$p.value
} else {
  p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.loglog.flexible<-wilcox.test(MAE.MAP.clm.loglog.flexible,MAE.Ord.MAP.clm.loglog.flexible,paired = TRUE,alternative="greater")$p.value
}


p.normality.MAE.clm.loglog.symmetric<-shapiro.test(MAE.MAP.clm.loglog.symmetric-MAE.Ord.MAP.clm.loglog.symmetric)$p.value

if (p.normality.MAE.clm.loglog.symmetric >= 0.05)
{p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.loglog.symmetric<-t.test(MAE.MAP.clm.loglog.symmetric,MAE.Ord.MAP.clm.loglog.symmetric,paired = TRUE,alternative="greater")$p.value
} else {
  p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.loglog.symmetric<-wilcox.test(MAE.MAP.clm.loglog.symmetric,MAE.Ord.MAP.clm.loglog.symmetric,paired = TRUE,alternative="greater")$p.value
}


p.normality.MAE.clm.loglog.symmetric2<-shapiro.test(MAE.MAP.clm.loglog.symmetric2-MAE.Ord.MAP.clm.loglog.symmetric2)$p.value

if (p.normality.MAE.clm.loglog.symmetric2 >= 0.05)
{p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.loglog.symmetric2<-t.test(MAE.MAP.clm.loglog.symmetric2,MAE.Ord.MAP.clm.loglog.symmetric2,paired = TRUE,alternative="greater")$p.value
} else {
  p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.loglog.symmetric2<-wilcox.test(MAE.MAP.clm.loglog.symmetric2,MAE.Ord.MAP.clm.loglog.symmetric2,paired = TRUE,alternative="greater")$p.value
}


p.normality.MAE.clm.loglog.equidistant<-shapiro.test(MAE.MAP.clm.loglog.equidistant-MAE.Ord.MAP.clm.loglog.equidistant)$p.value

if (p.normality.MAE.clm.loglog.equidistant >= 0.05)
{p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.loglog.equidistant<-t.test(MAE.MAP.clm.loglog.equidistant,MAE.Ord.MAP.clm.loglog.equidistant,paired = TRUE,alternative="greater")$p.value
} else {
  p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.loglog.equidistant<-wilcox.test(MAE.MAP.clm.loglog.equidistant,MAE.Ord.MAP.clm.loglog.equidistant,paired = TRUE,alternative="greater")$p.value
}

###


p.normality.MAE.clm.cloglog.flexible<-shapiro.test(MAE.MAP.clm.cloglog.flexible-MAE.Ord.MAP.clm.cloglog.flexible)$p.value

if (p.normality.MAE.clm.cloglog.flexible >= 0.05)
{p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.cloglog.flexible<-t.test(MAE.MAP.clm.cloglog.flexible,MAE.Ord.MAP.clm.cloglog.flexible,paired = TRUE,alternative="greater")$p.value
} else {
  p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.cloglog.flexible<-wilcox.test(MAE.MAP.clm.cloglog.flexible,MAE.Ord.MAP.clm.cloglog.flexible,paired = TRUE,alternative="greater")$p.value
}


p.normality.MAE.clm.cloglog.symmetric<-shapiro.test(MAE.MAP.clm.cloglog.symmetric-MAE.Ord.MAP.clm.cloglog.symmetric)$p.value

if (p.normality.MAE.clm.cloglog.symmetric >= 0.05)
{p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.cloglog.symmetric<-t.test(MAE.MAP.clm.cloglog.symmetric,MAE.Ord.MAP.clm.cloglog.symmetric,paired = TRUE,alternative="greater")$p.value
} else {
  p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.cloglog.symmetric<-wilcox.test(MAE.MAP.clm.cloglog.symmetric,MAE.Ord.MAP.clm.cloglog.symmetric,paired = TRUE,alternative="greater")$p.value
}


p.normality.MAE.clm.cloglog.symmetric2<-shapiro.test(MAE.MAP.clm.cloglog.symmetric2-MAE.Ord.MAP.clm.cloglog.symmetric2)$p.value

if (p.normality.MAE.clm.cloglog.symmetric2 >= 0.05)
{p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.cloglog.symmetric2<-t.test(MAE.MAP.clm.cloglog.symmetric2,MAE.Ord.MAP.clm.cloglog.symmetric2,paired = TRUE,alternative="greater")$p.value
} else {
  p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.cloglog.symmetric2<-wilcox.test(MAE.MAP.clm.cloglog.symmetric2,MAE.Ord.MAP.clm.cloglog.symmetric2,paired = TRUE,alternative="greater")$p.value
}


p.normality.MAE.clm.cloglog.equidistant<-shapiro.test(MAE.MAP.clm.cloglog.equidistant-MAE.Ord.MAP.clm.cloglog.equidistant)$p.value

if (p.normality.MAE.clm.cloglog.equidistant >= 0.05)
{p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.cloglog.equidistant<-t.test(MAE.MAP.clm.cloglog.equidistant,MAE.Ord.MAP.clm.cloglog.equidistant,paired = TRUE,alternative="greater")$p.value
} else {
  p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.cloglog.equidistant<-wilcox.test(MAE.MAP.clm.cloglog.equidistant,MAE.Ord.MAP.clm.cloglog.equidistant,paired = TRUE,alternative="greater")$p.value
}

###

p.normality.MAE.clm.cauchit.flexible<-shapiro.test(MAE.MAP.clm.cauchit.flexible-MAE.Ord.MAP.clm.cauchit.flexible)$p.value

if (p.normality.MAE.clm.cauchit.flexible >= 0.05)
{p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.cauchit.flexible<-t.test(MAE.MAP.clm.cauchit.flexible,MAE.Ord.MAP.clm.cauchit.flexible,paired = TRUE,alternative="greater")$p.value
} else {
  p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.cauchit.flexible<-wilcox.test(MAE.MAP.clm.cauchit.flexible,MAE.Ord.MAP.clm.cauchit.flexible,paired = TRUE,alternative="greater")$p.value
}


p.normality.MAE.clm.cauchit.symmetric<-shapiro.test(MAE.MAP.clm.cauchit.symmetric-MAE.Ord.MAP.clm.cauchit.symmetric)$p.value

if (p.normality.MAE.clm.cauchit.symmetric >= 0.05)
{p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.cauchit.symmetric<-t.test(MAE.MAP.clm.cauchit.symmetric,MAE.Ord.MAP.clm.cauchit.symmetric,paired = TRUE,alternative="greater")$p.value
} else {
  p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.cauchit.symmetric<-wilcox.test(MAE.MAP.clm.cauchit.symmetric,MAE.Ord.MAP.clm.cauchit.symmetric,paired = TRUE,alternative="greater")$p.value
}


p.normality.MAE.clm.cauchit.symmetric2<-shapiro.test(MAE.MAP.clm.cauchit.symmetric2-MAE.Ord.MAP.clm.cauchit.symmetric2)$p.value

if (p.normality.MAE.clm.cauchit.symmetric2 >= 0.05)
{p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.cauchit.symmetric2<-t.test(MAE.MAP.clm.cauchit.symmetric2,MAE.Ord.MAP.clm.cauchit.symmetric2,paired = TRUE,alternative="greater")$p.value
} else {
  p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.cauchit.symmetric2<-wilcox.test(MAE.MAP.clm.cauchit.symmetric2,MAE.Ord.MAP.clm.cauchit.symmetric2,paired = TRUE,alternative="greater")$p.value
}


p.normality.MAE.clm.cauchit.equidistant<-shapiro.test(MAE.MAP.clm.cauchit.equidistant-MAE.Ord.MAP.clm.cauchit.equidistant)$p.value

if (p.normality.MAE.clm.cauchit.equidistant >= 0.05)
{p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.cauchit.equidistant<-t.test(MAE.MAP.clm.cauchit.equidistant,MAE.Ord.MAP.clm.cauchit.equidistant,paired = TRUE,alternative="greater")$p.value
} else {
  p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.cauchit.equidistant<-wilcox.test(MAE.MAP.clm.cauchit.equidistant,MAE.Ord.MAP.clm.cauchit.equidistant,paired = TRUE,alternative="greater")$p.value
}

###

p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.logit.flexible
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.logit.symmetric
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.logit.symmetric2
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.logit.equidistant

p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.probit.flexible
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.probit.symmetric
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.probit.symmetric2
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.probit.equidistant

p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.loglog.flexible
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.loglog.symmetric
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.loglog.symmetric2
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.loglog.equidistant

p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.cloglog.flexible
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.cloglog.symmetric
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.cloglog.symmetric2
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.cloglog.equidistant

p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.cauchit.flexible
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.cauchit.symmetric
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.cauchit.symmetric2
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.clm.cauchit.equidistant



###############################
###############################
###   Boxplots   MAE

MAE.clm.flexible.results<-as.data.frame(cbind(MAE.MAP.clm.logit.flexible-MAE.Ord.MAP.clm.logit.flexible,
                           MAE.MAP.clm.probit.flexible-MAE.Ord.MAP.clm.probit.flexible,
                           MAE.MAP.clm.loglog.flexible-MAE.Ord.MAP.clm.loglog.flexible,
                           MAE.MAP.clm.cloglog.flexible-MAE.Ord.MAP.clm.cloglog.flexible,
                           MAE.MAP.clm.cauchit.flexible-MAE.Ord.MAP.clm.cauchit.flexible))

boxplot(MAE.clm.flexible.results,
        main=paste("threshold flexible"), 
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
        main=paste("threshold symmetric"), 
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
        main=paste("threshold symmetric2"), 
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
        main=paste("threshold equidistant"), 
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



#########.  END
