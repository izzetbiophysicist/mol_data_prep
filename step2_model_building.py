#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jan 15 07:30:04 2023

@author: lucas
"""

from sklearn.feature_selection import mutual_info_classif
import os
import pickle
import pandas as pd

os.chdir('/home/lucas/virtual_screening_pipeline/') 

from mol_data_prep_6 import *

#### Sklearn
from sklearn.model_selection import GridSearchCV
from sklearn.ensemble import RandomForestRegressor
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import accuracy_score
from sklearn.metrics import f1_score
from sklearn.metrics import precision_score
from sklearn.metrics import recall_score
from sklearn.metrics import roc_curve, auc

import copy


def run_dataset(compounds):
    model_log = pd.DataFrame(columns=['features','synthetic', 'p_train','p_test','accuracy', 'f1_score', 'precision_score', 'recall_score', 'AUC', 'youden'])

    
    ###### Steps        
    
    #### Load molecule dataset into object, the csv should contain the molecules in the first column and the variable of interest in the second
    
    dataset = mol_dataset(compounds)
    
    #### Defining the problem defines which method wil be used for synthetic sample generation
    dataset.problem = 'classification'
    
    ### Set NA action for the y values, can be set to mean, zero or median
    ### Every missing value will be substituted by zero, the mean or median value
    
    dataset.na_action_y('zero')
    dataset.na_action_x = 'mean'
    
    #### Clean the smiles
    dataset.clean_smiles()
    
    ##### deal with duplicates removing all molecules that have identical tanimoto coefficient
    ##### can also be used to remove molecules have identical smiles using 'string' instead of 'tanimoto'
    ##### the y value of the duplicates is compressed into a single value using the 'mean' or the 'median' of the y values
    dataset.remove_dups('tanimoto', 'median')
    
    
    ### transform into y into log(y)
    
    #dataset.log_y()
    
    ## Binarize the data for classification
    ## takes the cutoff as arguments
    
    dataset.binarize(1000)
    
    
    #### Split sets, takes as argument the size os the test set
    #### if clusters=True is used, a vector of the cluster to which each molecule belong should also be passed as argument
    #### if clusters=True is used the test set will be taken equally from each cluster
    dataset.calculate_features('mordred', 'train_test')
    
    to_cluster = dataset.training_features
    
    to_cluster_solved = dataset.solve_non_numeric(to_cluster, na_action='mean')
    
    
    scaler = StandardScaler()
    scaler.fit(to_cluster_solved)
    to_cluster_norm = scaler.transform(to_cluster_solved)
    
    to_cluster_PCA = pca_dataset(to_cluster_norm, 3)
    
    
    kmeans = KMeans(init="random", n_clusters=kmeans_screen(to_cluster_PCA[0], 50, 20)+1, n_init=10, max_iter=300, random_state=42)
    kclust = kmeans.fit(to_cluster_PCA[0])
    
    
    
    ### Calculates features for the training set using mordred
    
    
    ### split_sets
    dataset.split_sets(0.2, clusters=kclust.labels_, use_cluster=True)
    dataset.training_set = dataset.solve_non_numeric(dataset.training_set[0], na_action='mean'), dataset.training_set[1]
    
    sample_strategies = ['none','smote']
    
    for synth_sample in sample_strategies:
        if synth_sample == 'smote':
            dataset.synthetic()
            
        
        ## Can be set to normalization or standardization, applying min-max or z-scores respectivelly
        ## The scaler will be saved to normalize further sets
        
        dataset.normalization_type = 'standardize'
        
        ## Apply normalization/standardization
        dataset.normalize_training()
        
        ####### Iterativelly removes correlated features whose correlations are higher than the absolute value of the threshold
        
        
        dataset.remove_correlated_training(threshold=0.9)
    
        
        
        ###### Save dataset into another object for iterative feature selection
        
        ### Runs features selection. Takes as arguments a method and the number of features to be selected
        

        feature_list = [10, 50, 100, 200]
        #feature_list = [10, 50]
        #feature_list = list(range(10, np.shape(dataset.training_set[0])[1]+1, 90))
        #feature_list.append(np.shape(dataset.training_set[0])[1])
        
        for n_features in feature_list:
            
            
            ### Copy dataset and select features
            
            dataset_tolearn = copy.deepcopy(dataset)
            dataset_tolearn.feature_selection(mutual_info_classif, k=n_features)
            dataset_tolearn.filter_selected_features_training()
        
        
            ### Save dataset with pickle for future stuff
        
            with open('dataset_'+synth_sample+'_'+str(n_features)+'_features', 'wb') as handle:
                pickle.dump(dataset, handle, protocol=pickle.HIGHEST_PROTOCOL)
        
            #### Random forest
        
        # =============================================================================
            n_estimators = [int(x) for x in np.linspace(start = 200, stop = 2000, num = 5)]
            max_features = ['auto', 'sqrt']
            max_depth = [int(x) for x in np.linspace(10, 110, num = 5)]
            max_depth.append(None)
            min_samples_split = [2, 5, 12,20]
            min_samples_leaf = [1, 2, 4, 8]
            bootstrap = [True, False]
            
            
            grid = {'n_estimators': n_estimators,
                            'max_features': max_features,
                            'max_depth': max_depth,
                            'min_samples_split': min_samples_split,
                            'min_samples_leaf': min_samples_leaf,
                            'bootstrap': bootstrap}
            
            #### Random Forest model
            rf = RandomForestClassifier()
            rf_search = GridSearchCV(estimator = rf, param_grid = grid, verbose=2, cv=3, n_jobs = -1)   
            rf_search.fit(dataset_tolearn.training_set[0], dataset_tolearn.training_set[1])    
            pred = rf_search.predict_proba(dataset_tolearn.test_set[0])
            
        
            #### Save model with pickle
            
            with open('RF_model_'+str(n_features)+'_features_'+synth_sample, 'wb') as handle:
                pickle.dump(rf_search, handle, protocol=pickle.HIGHEST_PROTOCOL)
                
        
            ### ROC and Youden
        
            fpr, tpr, threshold = roc_curve(dataset_tolearn.test_set[1], pred[:,1])
            roc_auc = auc(fpr, tpr)
            optimal_idx = np.argmax(tpr - fpr)
            optimal_threshold = threshold[optimal_idx]
            
            binary = [1 if i > optimal_threshold else 0 for i in pred[:,1]]
            
            accuracy = accuracy_score(dataset_tolearn.test_set[1],binary)
            f1 = f1_score(dataset_tolearn.test_set[1],binary)
            precision = precision_score(dataset_tolearn.test_set[1],binary)
            recall = recall_score(dataset_tolearn.test_set[1],binary)
            
            #### Other data to save
            ## Percentage of actives in training and test set
            p_test = sum(dataset_tolearn.test_set[1])/len(dataset_tolearn.test_set[1])
            p_train = sum(dataset_tolearn.training_set[1])/len(dataset_tolearn.training_set[1])
            
            
            new_data_log = pd.DataFrame({'features':n_features,'synthetic':synth_sample, 'p_train':p_train,'p_test':p_test,'accuracy':accuracy, 'f1_score':f1, 'precision_score':precision, 'recall_score':recall, 'AUC':roc_auc, 'youden':optimal_threshold}, index=[0])
            model_log = pd.concat([model_log, new_data_log], ignore_index=True)
            model_log.to_csv('tmp_model_log.csv')
    
    
    model_log.to_csv('model_log.csv')
    
### List all directories
all_entries = os.listdir('/home/lucas/virtual_screening_pipeline/datasets/')

processed_dirs = []

for curr_dir in all_entries:
    os.chdir('/home/lucas/virtual_screening_pipeline/datasets/'+curr_dir) 
    
    compounds = pd.read_csv("dataset.csv")
    
    run_dataset(compounds)
    
    processed_dirs.append(curr_dir)
    
    # Write the list of processed directories to a file
    with open("/home/lucas/virtual_screening_pipeline/meta_progression.txt", "w") as file:
        for dir_name in processed_dirs:
            file.write(dir_name + "\n")
    
    
