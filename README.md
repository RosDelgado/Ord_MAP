# Ord_MAP

R Scripts for the experimentation phase of the manuscript:

"Ord-MAP criterion: extending binary MAP for Ordinal Classification"
by Rosario Delgado (submitted to publication)
This repository was developed by Rosario Delgado (Universitat Autònoma de Barcelona, 2025).


For experimental evaluation of the Ord-MAP criterion for ordinal classification (versus the traditional MAP criterion), in Section 5 of the paper diverse datasets and classifiers are considered. Even when classifiers are tailored for ordinal data, Ord-MAP achieves statistically significant improvements in prediction, with respect to the Mean Absolute Error (MAE) metric. 

We consider four real-world datasets:

# World Values Surveys (WVS) dataset. 
This dataset is sourced from the "carData" R package (Dataset to accompany J. Fox and S. Weisberg, An R Companion to Applied Regression, Third Edition, Sage (2019). https://doi.org/10.32614/CRAN.package.carData). Target variable Poverty with three categories: Too Little, About Right, Too Much. 

Script: ordinal_logistic_regression_poverty.R

Classifier: ordinal logistic regression using "polr" function from the "MASS" R package. 

# Wine dataset. 
This dataset is available in the "ordinal" R package (https://doi.org/10.32614/CRAN.package.ordinal). The target variable is Rating, ordinal with five values, from 1 to 5. 

Script: ordinal_package_wine.R

Classifier: cumulative link (mixed) models also known as ordered regression models, proportional odds models, proportional hazards models for grouped survival times and ordered logit/probit/... models, using the "clm" function from "ordinal" R package. 

# Hearth dataset.
This dataset is included in the "ordinalForest" R package (https://doi.org/10.32614/CRAN.package.ordinalForest). The target variable Class has five ordered categories, from 1 to 5. 

Script: ordinalForest_Class.R

Classifier: ordinal random forest using the "ordfor" function from the "ordinalForest" R package. 

# Parkinson dataset. 
This dataset is available from the UC Irvine Machine Learning Repository(https://archive.ics.uci.edu/dataset/189/parkinsons+telemonitoring). Two target variables are considered, which are continuous and have been discretized: 
• v5:motor UPDRS (clinician’s motor Unified Parkinson’s Disease Rating Scale). This scale is a widely used score to track disease progression. Categories 1 to 5: < 13, [13, 18), [18, 24), [24, 29), ≥ 29.
• v6:total UPDRS (clinician’s total Unified Parkinson’s Disease Rating Scale). 
Categories 1 to 5: < 20, [20, 25), [25, 30), [30, 40), ≥ 40.

Script: train_caret_rf_PARKINSON.R

Classifier: we use the "train" function, from the "caret" R library to set up a grid of tuning parameters for a number of 
classification and regression routines, fits each model and calculates a resampling based performance measure.
Uses "trainControl" argument from "caret".

# Specific R scripts
The four scripts use the script mat_square.R (introduced here), which converts any matrix in a square matrix with desired row/column labels, by adding zeros if needed.


# Author
Rosario Delgado (Universitat Autònoma de Barcelona, Spain, 2025).


