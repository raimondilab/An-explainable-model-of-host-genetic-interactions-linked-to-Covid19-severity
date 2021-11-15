#### 1) Stratified K-fold split of sample cohort into train and test sets

'''In this task we embedded a strategy for variant screening in, stratified k-fold cross-validation (using the StratifiedKFold function from the scikit-learn library) scheme to generate 5 random splits, into a training and testing test, of the original dataset. Each fold was constituted by an 80 % training set which was also employed for variant screening and a 20 % testing set. The variants in the test set were curated from the variants screened in the training set.'''

# import relevant libraries 
#!/data/SW/anaconda3/envs/myenv/bin/python

# import relevant libraries 
import os,sys,operator
import gzip, glob
import pickle, math
import numpy as np
import pandas as pd
from scipy.stats import chisquare, fisher_exact
from statsmodels.stats.contingency_tables import Table2x2
from sklearn.model_selection import StratifiedKFold, train_test_split


# Import phenotype dataset
fs="\t"
# path="/data/DB/Covid19_WES/NEW_SET2/"  # path to zip vcf WES files in remote server repository.
path="/Users/Hp/Desktop/2000_cohort_analysis/" # path to csv phenotype file in local computer

df = pd.read_csv(path+"pheno_sample_improved_grading_2000_patients.csv", delimiter=',')

# define the dependent & independent variables 
y = df.delta_grading
X = df.drop('delta_grading', 1)

# define stratified k-split
skf = StratifiedKFold(n_splits=5)

# Stratified K-fold cross-validation 
df['kfold'] = -1
df = df.sample(frac=1).reset_index(drop=True)

# set the fold training and testing set sampling units 
fold_no = 0
for train_index, test_index in skf.split(df, y):
    train = df.loc[train_index,:]
    test = df.loc[test_index,:]
    print('Fold',str(fold_no),
          'Class Ratio:',
          sum(test['delta_grading'])/len(test['delta_grading']))
    fold_no += 1
    train_filename = 'train_split_' + str(fold_no) + '.csv'
    test_filename = 'test_split_' + str(fold_no) + '.csv'
    train.to_csv(path + 'pheno_sample_patients_train_2000_patients' + train_filename, index=False)
    test.to_csv(path + 'pheno_sample_patients_train_2000_patients' + test_filename, index=False)
    