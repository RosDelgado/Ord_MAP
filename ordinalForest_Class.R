
###########################################
###########################################
#
#  EXPERIMENTAL PHASE 
#
# The ordinal forest (OF) method allows ordinal regression with high-dimensional 
# and low-dimensional data. After having constructed an OF prediction rule using 
# a training dataset, it can be used to predict the values of the ordinal target 
# variable for new observations. Moreover, by means of the (permutation-based) 
# variable importance measure of OF, it is also possible to rank the covariates 
# with respect to their importances in the prediction of the values of the ordinal target variable.
# OF is presented in Hornung (2020).

# ordinal random forest using ordfor function from the ordinalForest package. 
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
# wave inversions and/or ST elevation or depression of > 0.05 mV); 0 = normal)
# thalach. numeric. Maximum heart rate achieved
# exang. factor. Exercise induced angina (1 = yes; 0 = no)
# oldpeak. numeric. ST depression induced by exercise relative to rest
# Class. factor. Ordinal target variable - severity of coronary artery disease (determined using
# angiograms) (1 = no disease; 2 = degree 1; 3 = degree 2; 4 = degree 3; 5 = degree 4)

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

install.packages("ordinalForest")

library(ordinalForest)
# Load example dataset:
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
v.b<-seq(from=5, to=14, by=1)    # values for the hyperparameter nbest (default=10)
 

model.ord.forest<-list()
for (k in 1:length(v.b))
{ model.ord.forest[[k]]<-list()
  for (i in 1:10)
  {model.ord.forest[[k]][[i]] <- ordfor(depvar="Class", data=training[[i]], nsets=v.a, nbest=v.b[k],ntreeperdiv=100, 
                                   ntreefinal=5000, npermtrial=500,perffunction = "probability")
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

p.normality.MAE.ord.forest<-vector()
for (k in 1:length(v.b))
{p.normality.MAE.ord.forest[k]<-shapiro.test(MAE.MAP.ord.forest[[k]]-MAE.Ord.MAP.ord.forest[[k]])$p.value
}


p.val.test.MAE.MAP.greater.MAE.ORD.MAP.ord.forest<-vector()

for (k in 1:length(v.b))
{if (p.normality.MAE.ord.forest[k] >= 0.05)
{p.val.test.MAE.MAP.greater.MAE.ORD.MAP.ord.forest[k]<-t.test(MAE.MAP.ord.forest[[k]],MAE.Ord.MAP.ord.forest[[k]],paired = TRUE,alternative="greater")$p.value
} else {
  p.val.test.MAE.MAP.greater.MAE.ORD.MAP.ord.forest[k]<-wilcox.test(MAE.MAP.ord.forest[[k]],MAE.Ord.MAP.ord.forest[[k]],paired = TRUE,alternative="greater")$p.value
}
}


p.val.test.MAE.MAP.greater.MAE.ORD.MAP.ord.forest




###############################
###############################
###   Boxplots   MAE


MAE.results.ord.forest<-as.data.frame(cbind(MAE.MAP.ord.forest[[1]]-MAE.Ord.MAP.ord.forest[[1]],
                                            MAE.MAP.ord.forest[[2]]-MAE.Ord.MAP.ord.forest[[2]],
                                            MAE.MAP.ord.forest[[3]]-MAE.Ord.MAP.ord.forest[[3]],
                                            MAE.MAP.ord.forest[[4]]-MAE.Ord.MAP.ord.forest[[4]],
                                            MAE.MAP.ord.forest[[5]]-MAE.Ord.MAP.ord.forest[[5]],
                                            MAE.MAP.ord.forest[[6]]-MAE.Ord.MAP.ord.forest[[6]],
                                            MAE.MAP.ord.forest[[7]]-MAE.Ord.MAP.ord.forest[[7]],
                                            MAE.MAP.ord.forest[[8]]-MAE.Ord.MAP.ord.forest[[8]],
                                            MAE.MAP.ord.forest[[9]]-MAE.Ord.MAP.ord.forest[[9]],
                                            MAE.MAP.ord.forest[[10]]-MAE.Ord.MAP.ord.forest[[10]]
                                            ))

boxplot(MAE.results.ord.forest,
        main=paste(" "), 
        xlab="nbest hyperparmeter", ylab="MAE increment (MAP - Ord.MAP)",names=c("5", "6", "7", "8", "9", "10", "11", "12", "13", "14"),
        col="white")
abline(h = 0, col = 2, lwd = 2, lty=2)
# Dots
stripchart(MAE.results.ord.forest,
           #method = "jitter",
           pch = 19,
           col = 1:5,
           vertical = TRUE,
           add = TRUE)



#########.  END
