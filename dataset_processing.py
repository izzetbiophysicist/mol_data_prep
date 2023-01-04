# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import pandas as pd
import numpy as np
import rdkit as rd
import argparse
from mordred import Calculator
from mordred import descriptors
from random import sample
from sklearn.model_selection import train_test_split
from numpy.random import randint

compounds = pd.read_csv("~/virtual_screening_pipeline/compounds2.csv")

class dataset:
    def __init__(self, compounds):
        self.molecules = compounds[compounds.columns[0]]
        self.molecules = self.molecules.tolist()
        self.y = compounds[compounds.columns[1]]
        self.y = self.y.tolist() 
        
        self.external = []
        
        self.training = []
        
        self.validation = []
        
        self.training_features = []
                
        
    ##########################
    ### Define what to do to NA, NaN and Inf values
    ##########################
    
    def na_action(self, action):
        if action == 'zero':
            for cp in range(len(self.y)):
                if self.y[cp] == 'Inf' or self.y[cp] == 'NA' or self.y[cp] == 'NaN' or self.y[cp] == '' or self.y[cp] == ' ':
                    self.y[cp] = 0
        
        if action == 'mean':
            
            ### verify which are float and take mean
            self.y_num = [a for a in self.y if type(a) == float]
            mean_y = np.mean(self.y_num)
            
            for cp in range(len(self.y)):
                if self.y[cp] == 'Inf' or self.y[cp] == 'NA' or self.y[cp] == 'NaN' or self.y[cp] == '' or self.y[cp] == ' ':
                    self.y[cp] = mean_y
                    
        if action == 'median':
            
            ### verify which are float and take mean
            self.y_num = [a for a in self.y if type(a) == float]
            median_y = np.median(self.y_num)
            
            for cp in range(len(self.y)):
                if self.y[cp] == 'Inf' or self.y[cp] == 'NA' or self.y[cp] == 'NaN' or self.y[cp] == '' or self.y[cp] == ' ':
                    self.y[cp] = median_y
          
    ##########################
    ### Binarize data for classification
    ##########################
    
    def binarize(self, cutoff):
        self.binary = ['active' if a <= cutoff else 'inactive' for a in self.y]

    ##########################
    ### apply log to y
    ##########################

    def log_y(self):
        self.y = [np.log(a) for a in self.y]   
    
    #########################
    ### separate sets
    #########################
    
    
    def split_sets(self, use_cluster, training_size, external_size):
       #### use_cluster - use clustering for creating the test set, if True, previous clustering analysis must be carried out
       #### training_size - size of the training set
       #### external_size - size of the external dataset, can be zero. In this case only training and validation are used
       def sample_set(molecules, y, size):
           ##### Sample external
           set_sample = sample(range(len(molecules)), size)
           
           mol_sample = [molecules[a] for a in set_sample]
           y_sample = [y[a] for a in set_sample]
           
           
           rest_mol = [molecules[a] for a in range(len(molecules)) if a not in set_sample]
           rest_y = [y[a] for a in range(len(y)) if a not in set_sample]
           
           return mol_sample, y_sample, rest_mol, rest_y
       
       if training_size + external_size > 1:
            print('wrong group sizes')
       else:    
            
           #### Separate sets without clustering
           
           if use_cluster == False:
               if external_size > 0:
                    ### sample external set
                    ext_n = round(external_size*len(self.molecules))
                    train_n = round(training_size*len(self.molecules))
                    
                    external = sample_set(self.molecules, self.y, ext_n)
                    training = sample_set(external[2], external[3], train_n)
                    validation = sample_set(training[2], training[3], len(self.molecules)-(ext_n+train_n))
                
                    ### returns external
                    self.external = external[0], external[1]
                    self.training = training[0], training[1]
                    self.validation = validation[0], validation[1]

                   
               if external_size <= 0:
                    ### sample external set
                    external = sample_set(self.molecules, self.y, external_size)
                    training = sample_set(external[2], external[3], training_size)
                    validation = sample_set(training[2], training[3], 1-training_size)
                
                    ### returns external
                    self.training = training[0], training[1]
                    self.validation = validation[0], validation[1]
                  
               
         
            
    #########################
    ### calculate descriptors
    #########################    

    def calculate_features_training(self, feature_type):
            
        def calculate_mordred_training(self):
            
        ### Calculate descriptors for training set
            
            mol_array = [rd.Chem.MolFromSmiles(m, sanitize=True) for m in self.training[0]]
            calc = Calculator(descriptors, ignore_3D=True)
            return calc.pandas(mol_array)
        
        #### check feature type
        if feature_type in ['mordred','amino']:
            if feature_type == 'mordred':
                self.training_features = calculate_mordred_training(self)
        
        
bla = dataset(compounds)
bla.na_action('zero')
bla.binarize(30)
bla.log_y()

bla.split_sets(False, 0.7, 0.1)
bla.training[0]
bla.calculate_features_training('mordred')

bla.training_features
