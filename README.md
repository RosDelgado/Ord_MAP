# Ord_MAP

R Scripts for the experimentation phase of the manuscript:

"Ord-MAP criterion: extending binary MAP for Ordinal Classification"
by Rosario Delgado (submitted to publication)
This repository was developed by Rosario Delgado (Universitat Autònoma de Barcelona, 2025).


For experimental evaluation of the Ord-MAP criterion for ordinal classification (versus the traditional MAP criterion), in Section 5 of the paper diverse datasets and classifiers are considered. Even when classifiers are tailored for ordinal data, Ord-MAP achieves statistically significant improvements in prediction, with respect to the Mean Absolute Error (MAE) metric. 



We consider three real-world datasets:

Facial age dataset (https://www.kaggle.com/datasets/frabbisw/facial-age)
Abalone dataset (https://archive.ics.uci.edu/dataset/1/abalone)
Parkinson dataset (https://archive.ics.uci.edu/dataset/189/parkinson+telemonitoring)
It uses the content in https://github.com/giuliabinotto/ IntervalScaleClassification, which correspond to Section 4 fot the same paper, where scripts facilitate the computation of two ordinal metrics, Mean Absolute Error (MAE) and Total Cost (TC), alongside their interval scale counterparts introduced in the paper, with a specific section designed to address scenarios in which the rightmost interval is unbounded.

Description
From_png_to_dataframe.R
This script allows to load face files (.png) obtained from https://www.kaggle.com/datasets/frabbisw/facial-age and transform them into a dataframe, with 9,673 rows corresponding to face pictures, and 32x32+1=1025 columns, the last one being "age", while the others are the features V1,..., V1024. The dataframe is save as "faces.grey.32.Rda".

train_caret_rf_FACES.R
This script loads "faces.grey.32.Rda" and develops the experimental phase explained in Section 5 of the paper, corresponding use the interval-scale metrics to tuning hyper-parameter mtry for random forest using the caret library.

tune_control_e1071_knn_FACES.R
This script is similar to the previous one, but corresponds to tuning hyper-parameter k for k-nearest neighbors (knn) using e1071 library.

train_caret_rf_ABALONE.R
Similar to "train_caret_rf_FACES.R" but for the Abalone dataset.

tune_control_e1071_knn_ABALONE.R
Similar to "tune_control_e1071_knn_FACES.R" but for the Abalone dataset.

train_caret_rf_PARKINSON.R
Similar to "train_caret_rf_FACES.R" but for the Parkinson dataset.

tune_control_e1071_knn_PARKINSON.R
Similar to "tune_control_e1071_knn_FACES.R" but for the Parkinson dataset.

Requirements
General R libraries
The following libraries are needed:

magick, stringr, mdatools, png and utils (used by "From_png_to_dataframe.R")

arules (used by "train_caret_rf.R" and "tune_control_e1071.R")

caret (used by "train_caret_rf.R")

e1071 and class (used by "tune_control_e1071.R")

Specific R scripts
mat_square.R (introduced here): converts any matrix in a square matrix with desired row/column labels, by adding zeros if needed.

From https://github.com/giuliabinotto/ IntervalScaleClassification

MAE.R: computes MAE and normalized SMAE metrics.

MAEintervals.R: computes MAE.int and SMAE.int metrics.

Authors
Giulia Binotto & Rosario Delgado (Universitat Autònoma de Barcelona, Spain, 2024).

About
No description, website, or topics provided.
Resources
 Readme
License
 GPL-3.0 license
 Activity
Stars
 0 stars
Watchers
 1 watching
Forks
 0 forks
Releases
No releases published
Create a new release
Packages
No packages published
Publish your first package
Languages
R
100.0%
Footer
© 2025 GitHub, Inc.
Footer navigation
Terms
Privacy
Security
Status
Docs
Contact
Manage cookiesg the cases proposed in George et al. (2016).
The Simulations.R script includes the code to generate the plots of section 4.2.
Authors
Giulia Binotto & Rosario Delgado (Universitat Autònoma de Barcelona, Spain).
