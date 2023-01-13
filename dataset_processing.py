# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""


### computational chemistry stuff
import rdkit as rd
from mordred import Calculator
from mordred import descriptors

### Pandas, random and numpy
import pandas as pd
import numpy as np
from random import sample

from numpy.random import randint

from math import isnan

### Plot
from matplotlib import pyplot
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import pyplot as plt


### SKlearn
from sklearn.preprocessing import MinMaxScaler
from sklearn.preprocessing import StandardScaler
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA
from sklearn.metrics import silhouette_score
from sklearn.model_selection import ParameterGrid
from sklearn.feature_selection import SelectKBest 
from sklearn.feature_selection import chi2
from sklearn.feature_selection import f_classif
from sklearn.feature_selection import f_regression
from sklearn.feature_selection import mutual_info_regression
from sklearn.feature_selection import mutual_info_classif
from sklearn.model_selection import train_test_split
from sklearn.linear_model import LinearRegression
from imblearn.over_sampling import SMOTE



compounds = pd.read_csv("~/virtual_screening_pipeline/compounds2.csv")






class dataset:
    def __init__(self, compounds):
        self.molecules = compounds[compounds.columns[0]]
        self.molecules = self.molecules.tolist()
        
        ### set to be tested
        self.external = []
        self.external_features = []
        
        self.y = compounds[compounds.columns[1]]
        self.y = self.y.tolist() 
        
        self.normalization_type = 'None'

        self.na_action_x = []
        self.problem = 'None'
        
        self.model_list = []
        
        self.isbinary = 'False'

        
        
    ##########################
    ### Define what to do to NA, NaN and Inf values
    ##########################
    
    def na_action_y(self, action):
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
        self.isbinary = 'True'

    ##########################
    ### apply log to y
    ##########################

    def log_y(self):
        self.y = [np.log(a) for a in self.y]   
        
        
        
    ########################
    ### Normalization
    ########################
    def solve_non_numeric(self,feature_dataframe, na_action):
        #### Deal with non-numeric values
        if na_action == 'mean':
        
            norm = feature_dataframe.to_numpy()
            for j in range(np.shape(norm)[1]):
                error = [i for i in range(np.shape(norm)[0]) if type(norm[i,j]) != float and type(norm[i,j]) != int]            
                
                ### action if all column if empty
                if len(error) == np.shape(norm)[0]:
                    for to_replace in error:
                        norm[to_replace][j] = 0
                
                ## else, calculate
                else:
                    
                    to_av = [norm[i,j] for i in range(np.shape(norm)[0]) if i not in error]
                    col_av = np.mean(to_av)
                    
                    for to_replace in error:
                        norm[to_replace][j] = col_av
                        
        if na_action == 'median':
            norm = feature_dataframe.to_numpy()
            for j in range(np.shape(norm)[1]):
                error = [i for i in range(np.shape(norm)[0]) if type(norm[i,j]) != float and type(norm[i,j]) != int]            
                
                ### action if all column if empty
                if len(error) == np.shape(norm)[0]:
                    for to_replace in error:
                        norm[to_replace][j] = 0
                
                ## else, calculate
                else:
                    
                    to_av = [norm[i,j] for i in range(np.shape(norm)[0]) if i not in error]
                    col_av = np.median(to_av)
                    
                    for to_replace in error:
                        norm[to_replace][j] = col_av
                    
        if na_action == 'zero':
            norm = feature_dataframe.to_numpy()
            for j in range(np.shape(norm)[1]):
                error = [i for i in range(np.shape(norm)[0]) if type(norm[i,j]) != float and type(norm[i,j]) != int]            
                
                for to_replace in error:
                    norm[to_replace][j] = 0
                    
        return norm
    
    def normalize_set(self, set_to_norm, na_action):
        
        
        if self.normalization_type == 'normalize':
            set_to_norm = self.solve_non_numeric(set_to_norm, na_action)
            
            scaler = MinMaxScaler()
            scaler.fit(set_to_norm)
            
            set_to_norm = scaler.transform(set_to_norm)
            return set_to_norm, scaler
        
        
        if self.normalization_type == 'standardize':
            set_to_norm = self.solve_non_numeric(set_to_norm, na_action)
            
            scaler = StandardScaler()
            scaler.fit(set_to_norm)
            
            set_to_norm = scaler.transform(set_to_norm)
            return set_to_norm, scaler
        
    def normalize_training(self):
        normalized = self.normalize_set(self.training_set[0], self.na_action_x)
        self.training_set = normalized[0], self.training_set[1]
        self.scaler = normalized[1]
        
        
        test_solve = self.solve_non_numeric(self.test_set[0], self.na_action_x)
        self.test_set = self.scaler.transform(test_solve), self.test_set[1]
        
        
    
    #########################
    ### separate sets
    #########################

    
    def split_sets(self, test_size, clusters=[], use_synthetic=False, use_cluster=False, smote_k=3):
       #### use_cluster - use clustering for creating the test set, if True, previous clustering analysis must be carried out
       #### training_size - size of the training set
       #### external_size - size of the external dataset, can be zero. In this case only training and validation are used
        def sample_set(features, y, size):
            ##### Sample external
            set_sample = sample(range(len(features)), size)
            
            mol_sample = features.iloc[set_sample,]
            y_sample = [y[a] for a in set_sample]
            
            rest = [i for i in range(len(features)) if i not in set_sample]
            rest_mol = features.iloc[rest,]
            rest_y = [y[a] for a in range(len(y)) if a not in set_sample]
            
            return mol_sample, y_sample, rest_mol, rest_y
        
            #### Separate sets without clustering
           
        if use_cluster == False:
            if test_size > 0:
                 ### sample external set
                 ext_n = round(test_size*len(self.training_features))
                 
                 
                 test = sample_set(self.training_features, self.y, ext_n)
                 
                 ### returns external
                 self.training_features = test[2], test[3]
                 self.test_features = test[0], test[1]
        

        
            
        if use_cluster == True:
            if len(clusters) != len(self.training_features):
                print('cluster vector does not match the dataset')
            
            if test_size > 0:
                
                 ### test N
                 ext_n = round(test_size*len(self.training_features))
                 
                 cluster_set = set(clusters)
                 
                 if ext_n < len(cluster_set):
                     print('test samples =< clusters')
                 
                 else:
                     ### loop over set
                                          
                     ###### This part takes 1 sample from each cluster until the ext_n reaches zero
                     
                     set_sample=[]
                     counter = 0
                     #### Equally distributed among clusters
                
                     while counter != ext_n:
                        for current_cl in cluster_set:
                            cl_mols = [i for i in range(len(clusters)) if clusters[i] == current_cl]
                            sampled = sample(cl_mols, 1)[0]
                            set_sample.append(sampled)
                            cl_mols = [cl_mols[i] for i in range(len(cl_mols)) if cl_mols[i] != sampled]
                            counter+=1
                            
                            if counter == ext_n:
                                break
                     
                     y_test = [self.y[i] for i in range(len(self.y)) if i in set_sample]      
                     self.test_set = self.training_features.loc[set_sample], y_test
                     y_training = [self.y[i] for i in range(len(self.y)) if i not in set_sample]
                     self.training_set = self.training_features.drop(set_sample), y_training
                     
        if use_synthetic == True and self.problem == 'classification':
            oversample = SMOTE(k_neighbors=smote_k)
            overs_training = oversample.fit_resample(self.training_set[0], self.training_set[1])
            self.training_set = overs_training[0], overs_training[1]
        
        if use_synthetic == True and self.problem == 'regression':
            print('Synthetic Sampling not implemented for regression yet')
               
    #### calculate descriptors with mordred
    def calculate_mordred(self, mol):
        
        
        mol_array = [rd.Chem.MolFromSmiles(m, sanitize=True) for m in mol]
        calc = Calculator(descriptors, ignore_3D=True)
        return calc.pandas(mol_array)    

    
    def calculate_features(self, feature_type, set_to_calc):
        if set_to_calc == 'training':      
        
            #### check feature type
            if feature_type in ['mordred','amino']:
                if feature_type == 'mordred':
                    self.training_features = self.calculate_mordred(self.molecules)
                    
        if set_to_calc == 'external':       
        
            #### check feature type
            if feature_type in ['mordred','amino']:
                if feature_type == 'mordred':
                    self.external_features = self.calculate_mordred(self.external)


    
    ############
    #### SAVE NORMALIZATION DATA FOR TEST AND SCREENING SET NORMALIZATION
    #############
    def load_external(self, external_compounds):
        self.external = external_compounds[external_compounds.columns[0]]
        self.external = self.external.tolist()
        
    def normalize_external(self):
        
        external_solve = self.solve_non_numeric(self.external_features, self.na_action_x)
        self.external_set = self.scaler.transform(external_solve)
    
      
    #########
    ### Feature selection
    #########
    def feature_selection(self, method, k):
        selector = SelectKBest(method, k)
        features = selector.fit_transform(self.training_set[0], self.training_set[1])
        
        feature_index = [i for i in range(len(selector.get_support())) if selector.get_support()[i] == True]
        
        column_names = list(pd.DataFrame(self.training_features).columns.values)
        feature_names = [column_names[i] for i in feature_index] 
        
        
        self.feature_selection = features, feature_index, feature_names
        
    def filter_selected_features_training(self):
        self.training_set = self.training_set[0][:,self.feature_selection[1]], self.training_set[1]
        self.test_set = self.test_set[0][:,self.feature_selection[1]], self.test_set[1]
        
    def filter_selected_features_external(self):
        self.external_set = self.external_set[:,self.feature_selection[1]]
        
    
            
    def detect_correlated(self, data, threshold):
        
        #### Does not support dataframes
        #### Feature dataset is converted to numpy array during normalization
        
        cormat=np.corrcoef(np.transpose(data))
        for i in range(np.shape(cormat)[0]):
            for j in range(np.shape(cormat)[1]):
                if isnan(cormat[i,j]):
                    cormat[i,j] = 0
                
                if i == j:
                    cormat[i,j] = 0
        
        redundant = []
        
        
        i=0
        j=0
        
        while i <= np.shape(cormat)[0]:
            while j <= np.shape(cormat)[1]:
                if np.absolute(cormat[i,j]) > threshold:
                    redundant.append(j)
                    cormat=np.delete(cormat, j, 1)
                    
                j+=1
            i+=1
        
        non_redundant = [i for i in range(np.shape(data)[1]) if i not in redundant]
        
        return non_redundant
    
    def remove_correlated_training(self, threshold):
        
        data = self.training_set[0]
        
        self.non_redundant = self.detect_correlated(self.training_set[0], threshold) 
        
        self.training_set = data[:,self.non_redundant], self.training_set[1]
        self.test_set = data[:,self.non_redundant], self.test_set[1]

    def remove_correlated_external(self):
        
        self.external_set = self.external_set[:,self.non_redundant]
        
    
    ################
    #### create new model   
    ################
    
            



###### Steps        

bla=[]
#### Load dataset into object
compounds = compounds.loc[1:9]
compounds=compounds.append(compounds.loc[1])
compounds=compounds.append(compounds.loc[1])


bla = dataset(compounds)
bla.problem = 'regression'

### Set NA action
bla.na_action_y('zero')
bla.na_action_x = 'mean'

### transform into log
bla.log_y()

## Binarize
bla.binarize(30)



### Split sets and calculate features for training set
bla.calculate_features('mordred', 'training')


#test=bla.normalize_set(bla.training_features,'mean')[0]
clusters=[0,0,0,0,0,0,1,1,1,2,2,2]
bla.split_sets(0.5, clusters, use_cluster=True, use_synthetic=True)


## set normalization type and normalize training set
bla.normalization_type = 'normalize'

######## normalize training set and define scaler
bla.normalize_training()

####### Remove correlated
### 

bla.remove_correlated_training(threshold=0.9)

### Selection
bla.feature_selection(f_regression, 5)
bla.filter_selected_features_training()









####### load and proccess external data
bla.load_external(compounds)
bla.calculate_features('mordred', 'external')
bla.normalize_external()
bla.remove_correlated_external()
bla.filter_selected_features_external()

###### Train models

### Linear regression
### Random Forest
### Neural net
### Gradient boost 











data = bla.training_set[0]
threshold=0.9







bla.pca_set(3)

### Cluster
bla.kmeans_scan_k(5,1)
bla.kmeans(3)







### Load external
bla.load_external(compounds)
bla.calculate_features('mordred','external')
bla.normalize_external('mean')

bla.split_sets(True, 0.1, [0,0,0,0,0,0,1,1,1,1,1,1])

bla.plot_pca()

bla.split_sets(False, 0.5)

#kmeans_model = KMeans()     # instantiating KMeans model
#parameters = [2, 3, 4]
#parameter_grid = ParameterGrid({'n_clusters': parameters})
#silhouette_scores = []

#for p in parameter_grid:
#    print(p)
#    kmeans_model.set_params(**p)    # set current hyper parameter
#    kmeans_model.fit(bla.pca_molecules)   
#    ss = silhouette_score(bla.pca_molecules, kmeans_model.labels_)   # calculate silhouette_score
#    silhouette_scores += [ss]       # store all the scores

#def plot_clusters(self):
    

    
def pca_set(self, n_components):
    
    set_to_calc = bla.normalize_set(bla.training_features, bla.na_action_x)[0]
    
    set_to_calc2 = np.transpose(set_to_calc)
    
    pca = PCA(n_components=n_components, svd_solver='arpack')
    self.pca_molecules = pca.fit(set_to_calc2)
    
    
