#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jan 15 07:30:04 2023

@author: lucas
"""

from mol_data_prep import *
from sklearn.feature_selection import mutual_info_classif


compounds = pd.read_csv("~/virtual_screening_pipeline/compounds_test.csv")


###### Steps        

#### Load molecule dataset into object, the csv should contain the molecules in the first column and the variable of interest in the second

dataset = mol_dataset(compounds)

#### Defining the problem defines which method wil be used for synthetic sample generation
dataset.problem = 'classification'

### Set NA action for the y values, can be set to mean, zero or median
### Every missing value will be substituted by zero, the mean or median value

dataset.na_action_y('zero')
dataset.na_action_x = 'mean'

### transform into y into log(y)

dataset.log_y()

## Binarize the data for classification
## takes the cutoff as arguments

dataset.binarize(5)

### Calculates features for the training set using mordred
dataset.calculate_features('mordred', 'training')

#### Split sets, takes as argument the size os the test set
#### if clusters=True is used, a vector of the cluster to which each molecule belong should also be passed as argument
#### if clusters=True is used the test set will be taken equally from each cluster

dataset.split_sets(0.2)


## Can be set to normalization or standardization, applying min-max or z-scores respectivelly
## The scaler will be saved to normalize further sets

dataset.normalization_type = 'standardize'

## Apply normalization/standardization
dataset.normalize_training()

####### Iterativelly removes correlated features whose correlations are higher than the absolute value of the threshold

dataset.remove_correlated_training(threshold=0.9)

### Runs features selection. Takes as arguments a method and the number of features to be selected
dataset.feature_selection(mutual_info_classif, 10)
dataset.filter_selected_features_training()


####### load molecules into external dataset and apply the same pipeline
dataset.load_external(compounds)
dataset.calculate_features('mordred', 'external')
dataset.normalize_external()
dataset.remove_correlated_external()
dataset.filter_selected_features_external()

