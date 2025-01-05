########################################
############   train function, from caret library
############
############ This function sets up a grid of tuning parameters for a number of 
############ classification and regression routines, fits each model and 
############ calculates a resampling based performance measure.
############ Uses "trainControl" argument from caret
###########
########################################

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

################################################################################

source("mat_square.R")

library(arules)   # for "discretize".
library(doParallel)
registerDoParallel(cores=6)

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

#arules::discretize(parkinson$V5, method = "frequency", breaks=5)

summary(parkinson$V5)
cut.points<-c(0,13,18,24,29,1000)

V5.bin <- arules::discretize(parkinson$V5, method = "fixed", breaks=cut.points, infinity=TRUE)
table(V5.bin)
parkinson<-as.data.frame(cbind(parkinson,V5.bin))


levels.V5.int<-names(table(parkinson$V5.bin))
levels(parkinson$V5.bin)<-c("<13","[13,18)","[18,24)","[24,29)",">=29")
levels.V5.ordinal.encod<-sort(unique(as.numeric(parkinson$V5.bin)))

par(mfrow = c(1, 1))
plot(parkinson$V5.bin,parkinson$V5, xlab="Intervals for motor UPDRS", ylab="motor UPDRS")

parkinson$V5.bin.num<-as.numeric(parkinson$V5.bin)
table(parkinson$V5.bin)
table(parkinson$V5.bin.num)

parkinson$V5.bin.num.factor<-as.factor(parkinson$V5.bin.num)   # classes in number but factor type 

################################################################################
####### Error functions for summaryFunction argument, trainControl function
#######

standard.MAE.ord.parkinson<-function(data,lev=levels.V5.ordinal.encod,model=NULL)
  # data=dataframe with columns "obs" and "pred" of character/factor type
  # lev=character string with outcome factor levels
{Conf.mat<- mat.square(table(as.numeric(data$pred),as.numeric(data$obs)),lev)
  value<-MAE(Conf.mat)
  c(MAE.ord=value)
}


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
################################################################################
########## caret::train. Resampling method: cross-validation
################################################################################
library(caret)

mtry<-sqrt(ncol(df)-4)
ntree<-3

fitControl.Accuracy <- trainControl(
  method = "cv",
  number = 3,
  search="random")



fitControl.MAE.parkinson <- trainControl(
                           method = "cv",
                           number = 3,
                           ## Evaluate performance using 
                           ## the following function
                           summaryFunction = standard.MAE.ord.parkinson)




tuned.rf.caret.Accuracy.parkinson<-list()
tuned.rf.caret.MAE.parkinson<-list()


for (i in 1:10)
{ set.seed(12345)
  tuned.rf.caret.Accuracy.parkinson[[i]] <- caret::train(training[[i]][ ,features], 
                                                 training[[i]][ ,25], 
                                                 method="rf", 
                                                 metric="Accuracy",
                                                 tuneLength=10,
                                                 trControl=fitControl.Accuracy)
  
  print(i)}




##############

for (i in 1:10)
{ set.seed(12345)
  tuned.rf.caret.MAE.parkinson[[i]] <- caret::train(training[[i]][ ,features], 
                                               training[[i]][ ,25], 
                                               method="rf",
                                          metric="MAE.ord",
                                          maximize=FALSE,
                                               tuneLength=10,
                                               trControl=fitControl.MAE.parkinson)
  
  print(i)}


################################################################################
#
# Predictions of any of the test set
# 
#

pred.test.rf.caret.Accuracy.parkinson<-list()
for (i in 1:10)
{pred.test.rf.caret.Accuracy.parkinson[[i]] <- predict(tuned.rf.caret.Accuracy.parkinson[[i]],test[[i]][,features],type="prob")
}


pred.test.rf.caret.MAE.parkinson<-list()
for (i in 1:10)
{pred.test.rf.caret.MAE.parkinson[[i]] <- predict(tuned.rf.caret.MAE.parkinson[[i]],test[[i]][,features],type="prob")
}

###

################################################################################
### predictions with standard MAP criterion:
###

pred.MAP.rf.caret.Accuracy.parkinson<-list()
for (i in 1:10)
{pred.MAP.rf.caret.Accuracy.parkinson[[i]]<-vector()
for (j in 1:dim(pred.test.rf.caret.Accuracy.parkinson[[i]])[1])
{pred.MAP.rf.caret.Accuracy.parkinson[[i]][j]<-categories[which.max(pred.test.rf.caret.Accuracy.parkinson[[i]][j,])]
}
}

pred.MAP.rf.caret.MAE.parkinson<-list()
for (i in 1:10)
{pred.MAP.rf.caret.MAE.parkinson[[i]]<-vector()
for (j in 1:dim(pred.test.rf.caret.MAE.parkinson[[i]])[1])
{pred.MAP.rf.caret.MAE.parkinson[[i]][j]<-categories[which.max(pred.test.rf.caret.MAE.parkinson[[i]][j,])]
}
}


# confusion matrices


conf.mat.MAP.rf.caret.Accuracy.parkinson<-list()
for(i in 1:10)
{
  conf.mat.MAP.rf.caret.Accuracy.parkinson[[i]]<-mat.square(table(pred.MAP.rf.caret.Accuracy.parkinson[[i]],test[[i]][[25]]),categories)
}


conf.mat.MAP.rf.caret.MAE.parkinson<-list()
for(i in 1:10)
{
  conf.mat.MAP.rf.caret.MAE.parkinson[[i]]<-mat.square(table(pred.MAP.rf.caret.MAE.parkinson[[i]],test[[i]][[25]]),categories)
}



# MAE 

MAE.MAP.rf.caret.Accuracy.parkinson<-vector()
for (i in 1:10)
{
  MAE.MAP.rf.caret.Accuracy.parkinson[i]<-MAE(conf.mat.MAP.rf.caret.Accuracy.parkinson[[i]])
}


MAE.MAP.rf.caret.MAE.parkinson<-vector()
for (i in 1:10)
{
  MAE.MAP.rf.caret.MAE.parkinson[i]<-MAE(conf.mat.MAP.rf.caret.MAE.parkinson[[i]])
}


###

################################################################################
### predictions with Ord.MAP criterion:
###

pred.Ord.MAP.rf.caret.Accuracy.parkinson<-list()
for (i in 1:10)
{pred.Ord.MAP.rf.caret.Accuracy.parkinson[[i]]<-vector()
for (j in 1:dim(pred.test.rf.caret.Accuracy.parkinson[[i]])[1])
{h<-min(which(cumsum(as.vector(pred.test.rf.caret.Accuracy.parkinson[[i]][j,]))>=0.5))
pred.Ord.MAP.rf.caret.Accuracy.parkinson[[i]][j]<-categories[h]
}
}


pred.Ord.MAP.rf.caret.MAE.parkinson<-list()
for (i in 1:10)
{pred.Ord.MAP.rf.caret.MAE.parkinson[[i]]<-vector()
for (j in 1:dim(pred.test.rf.caret.MAE.parkinson[[i]])[1])
{h<-min(which(cumsum(as.vector(pred.test.rf.caret.MAE.parkinson[[i]][j,]))>=0.5))
pred.Ord.MAP.rf.caret.MAE.parkinson[[i]][j]<-categories[h]
}
}


###
# confusion matrices


conf.mat.Ord.MAP.rf.caret.Accuracy.parkinson<-list()
for(i in 1:10)
{
  conf.mat.Ord.MAP.rf.caret.Accuracy.parkinson[[i]]<-mat.square(table(pred.Ord.MAP.rf.caret.Accuracy.parkinson[[i]],test[[i]][[25]]),categories)
}



conf.mat.Ord.MAP.rf.caret.MAE.parkinson<-list()
for(i in 1:10)
{
  conf.mat.Ord.MAP.rf.caret.MAE.parkinson[[i]]<-mat.square(table(pred.Ord.MAP.rf.caret.MAE.parkinson[[i]],test[[i]][[25]]),categories)
}


# MAE 


MAE.Ord.MAP.rf.caret.Accuracy.parkinson<-vector()
for (i in 1:10)
{
  MAE.Ord.MAP.rf.caret.Accuracy.parkinson[i]<-MAE(conf.mat.Ord.MAP.rf.caret.Accuracy.parkinson[[i]])
}


MAE.Ord.MAP.rf.caret.MAE.parkinson<-vector()
for (i in 1:10)
{
  MAE.Ord.MAP.rf.caret.MAE.parkinson[i]<-MAE(conf.mat.Ord.MAP.rf.caret.MAE.parkinson[[i]])
}



###
################################################################################
###. Comparison MAP vs Ord.MAP
###
###.  Tests MAE

p.normality.MAE.rf.caret.Accuracy.parkinson<-shapiro.test(MAE.MAP.rf.caret.Accuracy.parkinson-MAE.Ord.MAP.rf.caret.Accuracy.parkinson)$p.value

if (p.normality.MAE.rf.caret.Accuracy.parkinson >= 0.05)
{p.val.test.MAE.MAP.greater.MAE.ORD.MAP.rf.caret.Accuracy.parkinson<-t.test(MAE.MAP.rf.caret.Accuracy.parkinson,MAE.Ord.MAP.rf.caret.Accuracy.parkinson,paired = TRUE,alternative="greater")$p.value
} else {
  p.val.test.MAE.MAP.greater.MAE.ORD.MAP.rf.caret.Accuracy.parkinson<-wilcox.test(MAE.MAP.rf.caret.Accuracy.parkinson,MAE.Ord.MAP.rf.caret.Accuracy.parkinson,paired = TRUE,alternative="greater")$p.value
}


p.normality.MAE.rf.caret.MAE.parkinson<-shapiro.test(MAE.MAP.rf.caret.MAE.parkinson-MAE.Ord.MAP.rf.caret.MAE.parkinson)$p.value

if (p.normality.MAE.rf.caret.MAE.parkinson >= 0.05)
{p.val.test.MAE.MAP.greater.MAE.ORD.MAP.rf.caret.MAE.parkinson<-t.test(MAE.MAP.rf.caret.MAE.parkinson,MAE.Ord.MAP.rf.caret.MAE.parkinson,paired = TRUE,alternative="greater")$p.value
} else {
  p.val.test.MAE.MAP.greater.MAE.ORD.MAP.rf.caret.MAE.parkinson<-wilcox.test(MAE.MAP.rf.caret.MAE.parkinson,MAE.Ord.MAP.rf.caret.MAE.parkinson,paired = TRUE,alternative="greater")$p.value
}



p.val.test.MAE.MAP.greater.MAE.ORD.MAP.rf.caret.Accuracy.parkinson
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.rf.caret.MAE.parkinson


###############################
###############################
###   Boxplots   MAE

MAE.rf.caret.parkinson.results<-as.data.frame(cbind(MAE.MAP.rf.caret.Accuracy.parkinson-MAE.Ord.MAP.rf.caret.Accuracy.parkinson,
                                                    MAE.MAP.rf.caret.MAE.parkinson-MAE.Ord.MAP.rf.caret.MAE.parkinson))

boxplot(MAE.rf.caret.parkinson.results,
        main=paste(" "), 
        xlab="Summary Function", ylab="MAE increment (MAP - Ord.MAP)",names=c("Accuracy", "MAE"),
        col="white")
        
abline(h = 0, col = 2, lwd = 2, lty=2)
# Dots
stripchart(MAE.rf.caret.parkinson.results,
           #method = "jitter",
           pch = 19,
           col = 1:5,
           vertical = TRUE,
           add = TRUE)


##################################################################################
##################################################################################
##################################################################################


#### SECOND: OUTPUT VARIABLE TO PREDICT: V6: total_UPDRS (ranges from 0 to 176)

#arules::discretize(parkinson$V6, method = "frequency", breaks=5)

cut.points.v6<-c(0,20,25,30,40,1000)

V6.bin <- arules::discretize(parkinson$V6, method = "fixed", breaks=cut.points.v6, infinity=TRUE)
table(V6.bin)
parkinson<-as.data.frame(cbind(parkinson,V6.bin))


levels.V6.int<-names(table(parkinson$V6.bin))
levels(parkinson$V6.bin)<-c("<20","[20,25)","[25,30)","[30,40)",">=40")
levels.V6.ordinal.encod<-sort(unique(as.numeric(parkinson$V6.bin)))

par(mfrow = c(1, 1))
plot(parkinson$V6.bin,parkinson$V6, xlab="Intervals for total UPDRS", ylab="total UPDRS")

parkinson$V6.bin.num<-as.numeric(parkinson$V6.bin)
table(parkinson$V6.bin)
table(parkinson$V6.bin.num)

parkinson$V6.bin.num.factor<-as.factor(parkinson$V6.bin.num)   # classes in number but factor type 

################################################################################
####### Error functions for summaryFunction argument, trainControl function
#######

standard.MAE.ord.parkinson<-function(data,lev=levels.V5.ordinal.encod,model=NULL)
  # data=dataframe with columns "obs" and "pred" of character/factor type
  # lev=character string with outcome factor levels
{Conf.mat<- mat.square(table(as.numeric(data$pred),as.numeric(data$obs)),lev)
value<-MAE(Conf.mat)
c(MAE.ord=value)
}


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
################################################################################
########## caret::train. Resampling method: cross-validation
################################################################################
library(caret)

mtry<-sqrt(ncol(df)-4)
ntree<-3

fitControl.Accuracy <- trainControl(
  method = "cv",
  number = 3,
  search="random")



fitControl.MAE.parkinson <- trainControl(
  method = "cv",
  number = 3,
  ## Evaluate performance using 
  ## the following function
  summaryFunction = standard.MAE.ord.parkinson)




tuned.rf.caret.Accuracy.parkinson.v6<-list()
tuned.rf.caret.MAE.parkinson.v6<-list()


for (i in 1:10)
{ set.seed(12345)
  tuned.rf.caret.Accuracy.parkinson.v6[[i]] <- caret::train(training[[i]][ ,features], 
                                                         training[[i]][ ,28], 
                                                         method="rf", 
                                                         metric="Accuracy",
                                                         tuneLength=10,
                                                         trControl=fitControl.Accuracy)
  
  print(i)}


##############

for (i in 1:10)
{ set.seed(12345)
  tuned.rf.caret.MAE.parkinson.v6[[i]] <- caret::train(training[[i]][ ,features], 
                                                    training[[i]][ ,28], 
                                                    method="rf",
                                                    metric="MAE.ord",
                                                    maximize=FALSE,
                                                    tuneLength=10,
                                                    trControl=fitControl.MAE.parkinson)
  
  print(i)}



################################################################################
#
# Predictions of any of the test set
# 
#

pred.test.rf.caret.Accuracy.parkinson.v6<-list()
for (i in 1:10)
{pred.test.rf.caret.Accuracy.parkinson.v6[[i]] <- predict(tuned.rf.caret.Accuracy.parkinson.v6[[i]],test[[i]][,features],type="prob")
}


pred.test.rf.caret.MAE.parkinson.v6<-list()
for (i in 1:10)
{pred.test.rf.caret.MAE.parkinson.v6[[i]] <- predict(tuned.rf.caret.MAE.parkinson.v6[[i]],test[[i]][,features],type="prob")
}

###

################################################################################
### predictions with standard MAP criterion:
###

pred.MAP.rf.caret.Accuracy.parkinson.v6<-list()
for (i in 1:10)
{pred.MAP.rf.caret.Accuracy.parkinson.v6[[i]]<-vector()
for (j in 1:dim(pred.test.rf.caret.Accuracy.parkinson.v6[[i]])[1])
{pred.MAP.rf.caret.Accuracy.parkinson.v6[[i]][j]<-categories[which.max(pred.test.rf.caret.Accuracy.parkinson.v6[[i]][j,])]
}
}

pred.MAP.rf.caret.MAE.parkinson.v6<-list()
for (i in 1:10)
{pred.MAP.rf.caret.MAE.parkinson.v6[[i]]<-vector()
for (j in 1:dim(pred.test.rf.caret.MAE.parkinson.v6[[i]])[1])
{pred.MAP.rf.caret.MAE.parkinson.v6[[i]][j]<-categories[which.max(pred.test.rf.caret.MAE.parkinson.v6[[i]][j,])]
}
}


# confusion matrices


conf.mat.MAP.rf.caret.Accuracy.parkinson.v6<-list()
for(i in 1:10)
{
  conf.mat.MAP.rf.caret.Accuracy.parkinson.v6[[i]]<-mat.square(table(pred.MAP.rf.caret.Accuracy.parkinson.v6[[i]],test[[i]][[28]]),categories)
}


conf.mat.MAP.rf.caret.MAE.parkinson.v6<-list()
for(i in 1:10)
{
  conf.mat.MAP.rf.caret.MAE.parkinson.v6[[i]]<-mat.square(table(pred.MAP.rf.caret.MAE.parkinson.v6[[i]],test[[i]][[28]]),categories)
}



# MAE 

MAE.MAP.rf.caret.Accuracy.parkinson.v6<-vector()
for (i in 1:10)
{
  MAE.MAP.rf.caret.Accuracy.parkinson.v6[i]<-MAE(conf.mat.MAP.rf.caret.Accuracy.parkinson.v6[[i]])
}


MAE.MAP.rf.caret.MAE.parkinson.v6<-vector()
for (i in 1:10)
{
  MAE.MAP.rf.caret.MAE.parkinson.v6[i]<-MAE(conf.mat.MAP.rf.caret.MAE.parkinson.v6[[i]])
}


###

################################################################################
### predictions with Ord.MAP criterion:
###

pred.Ord.MAP.rf.caret.Accuracy.parkinson.v6<-list()
for (i in 1:10)
{pred.Ord.MAP.rf.caret.Accuracy.parkinson.v6[[i]]<-vector()
for (j in 1:dim(pred.test.rf.caret.Accuracy.parkinson.v6[[i]])[1])
{h<-min(which(cumsum(as.vector(pred.test.rf.caret.Accuracy.parkinson.v6[[i]][j,]))>=0.5))
pred.Ord.MAP.rf.caret.Accuracy.parkinson.v6[[i]][j]<-categories[h]
}
}


pred.Ord.MAP.rf.caret.MAE.parkinson.v6<-list()
for (i in 1:10)
{pred.Ord.MAP.rf.caret.MAE.parkinson.v6[[i]]<-vector()
for (j in 1:dim(pred.test.rf.caret.MAE.parkinson.v6[[i]])[1])
{h<-min(which(cumsum(as.vector(pred.test.rf.caret.MAE.parkinson.v6[[i]][j,]))>=0.5))
pred.Ord.MAP.rf.caret.MAE.parkinson.v6[[i]][j]<-categories[h]
}
}


###
# confusion matrices


conf.mat.Ord.MAP.rf.caret.Accuracy.parkinson.v6<-list()
for(i in 1:10)
{
  conf.mat.Ord.MAP.rf.caret.Accuracy.parkinson.v6[[i]]<-mat.square(table(pred.Ord.MAP.rf.caret.Accuracy.parkinson.v6[[i]],test[[i]][[28]]),categories)
}



conf.mat.Ord.MAP.rf.caret.MAE.parkinson.v6<-list()
for(i in 1:10)
{
  conf.mat.Ord.MAP.rf.caret.MAE.parkinson.v6[[i]]<-mat.square(table(pred.Ord.MAP.rf.caret.MAE.parkinson.v6[[i]],test[[i]][[28]]),categories)
}


# MAE 


MAE.Ord.MAP.rf.caret.Accuracy.parkinson.v6<-vector()
for (i in 1:10)
{
  MAE.Ord.MAP.rf.caret.Accuracy.parkinson.v6[i]<-MAE(conf.mat.Ord.MAP.rf.caret.Accuracy.parkinson.v6[[i]])
}


MAE.Ord.MAP.rf.caret.MAE.parkinson.v6<-vector()
for (i in 1:10)
{
  MAE.Ord.MAP.rf.caret.MAE.parkinson.v6[i]<-MAE(conf.mat.Ord.MAP.rf.caret.MAE.parkinson.v6[[i]])
}



###
################################################################################
###. Comparison MAP vs Ord.MAP
###
###.  Tests MAE

p.normality.MAE.rf.caret.Accuracy.parkinson.v6<-shapiro.test(MAE.MAP.rf.caret.Accuracy.parkinson.v6-MAE.Ord.MAP.rf.caret.Accuracy.parkinson.v6)$p.value

if (p.normality.MAE.rf.caret.Accuracy.parkinson.v6 >= 0.05)
{p.val.test.MAE.MAP.greater.MAE.ORD.MAP.rf.caret.Accuracy.parkinson.v6<-t.test(MAE.MAP.rf.caret.Accuracy.parkinson.v6,MAE.Ord.MAP.rf.caret.Accuracy.parkinson.v6,paired = TRUE,alternative="greater")$p.value
} else {
  p.val.test.MAE.MAP.greater.MAE.ORD.MAP.rf.caret.Accuracy.parkinson.v6<-wilcox.test(MAE.MAP.rf.caret.Accuracy.parkinson.v6,MAE.Ord.MAP.rf.caret.Accuracy.parkinson.v6,paired = TRUE,alternative="greater")$p.value
}


p.normality.MAE.rf.caret.MAE.parkinson.v6<-shapiro.test(MAE.MAP.rf.caret.MAE.parkinson.v6-MAE.Ord.MAP.rf.caret.MAE.parkinson.v6)$p.value

if (p.normality.MAE.rf.caret.MAE.parkinson.v6 >= 0.05)
{p.val.test.MAE.MAP.greater.MAE.ORD.MAP.rf.caret.MAE.parkinson.v6<-t.test(MAE.MAP.rf.caret.MAE.parkinson.v6,MAE.Ord.MAP.rf.caret.MAE.parkinson.v6,paired = TRUE,alternative="greater")$p.value
} else {
  p.val.test.MAE.MAP.greater.MAE.ORD.MAP.rf.caret.MAE.parkinson.v6<-wilcox.test(MAE.MAP.rf.caret.MAE.parkinson.v6,MAE.Ord.MAP.rf.caret.MAE.parkinson.v6,paired = TRUE,alternative="greater")$p.value
}



p.val.test.MAE.MAP.greater.MAE.ORD.MAP.rf.caret.Accuracy.parkinson.v6
p.val.test.MAE.MAP.greater.MAE.ORD.MAP.rf.caret.MAE.parkinson.v6


###############################
###############################
###   Boxplots   MAE

MAE.rf.caret.parkinson.results.v6<-as.data.frame(cbind(MAE.MAP.rf.caret.Accuracy.parkinson.v6-MAE.Ord.MAP.rf.caret.Accuracy.parkinson.v6,
                                                    MAE.MAP.rf.caret.MAE.parkinson.v6-MAE.Ord.MAP.rf.caret.MAE.parkinson.v6))

boxplot(MAE.rf.caret.parkinson.results.v6,
        main=paste(" "), 
        xlab="Summary Function", ylab="MAE increment (MAP - Ord.MAP)",names=c("Accuracy", "MAE"),
        col="white")

abline(h = 0, col = 2, lwd = 2, lty=2)
# Dots
stripchart(MAE.rf.caret.parkinson.results.v6,
           #method = "jitter",
           pch = 19,
           col = 1:5,
           vertical = TRUE,
           add = TRUE)


###########. END




