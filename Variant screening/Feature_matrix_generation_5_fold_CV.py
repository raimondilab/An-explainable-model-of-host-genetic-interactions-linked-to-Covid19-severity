
#### 3) Feature matrix generation 
'''For each split, we generated a feature matrix for the training set by assigning the allele counts of each screened variant for each sample of the training: i.e. 0 
for genotype 0/0, 1 for genotypes 1/0 or 0/1, 2 for genotype 1/1. The feature matrix for the test set was defined by considering only variants identified as significant
after screening the training set of the corresponding split and by assigning the allele count of each sample of the test set. We also included as additional features age,
which was normalized, and gender, which was binarized by setting males to 0 and females to 1. Severe patients from group “3+4+5” were given the classification label “1”, 
the asymptomatic patients from group 0 were given the label “0”.''' 

# Note we re-run this script 10-times to generate the 5-folds of 80 % screened variant training set and 5-folds of 20 % unscreened test sets variants identified 
# from the training sets of each folds. 

#!/data/SW/anaconda3/envs/myenv/bin/python

import os,sys,operator
import gzip, glob
import pickle, math
import pandas as pd
import gzip
import csv
import json


fs="\t"
#path="/data/DB/Covid19_WES/NEW_SET/feature_count_matrices&codes_1381_variants_grade43vs0/"

path="/data/DB/Covid19_WES/NEW_SET2/variant_stat_grading_improved_version_5-fold_CV_2000_cohort_latest/"
# assign phenotype dataset names
list_of_names = ['pheno_sample_patients_train_split_1','pheno_sample_patients_test_split_1',
                 'pheno_sample_patients_train_split_2','pheno_sample_patients_test_split_2', 
                'pheno_sample_patients_train_split_3','pheno_sample_patients_test_split_3',
                 'pheno_sample_patients_train_split_4','pheno_sample_patients_test_split_4', 
                 'pheno_sample_patients_train_split_5','pheno_sample_patients_test_split_5']
 
# create empty list
li = []
phenotype_list = []
filter_train = []
# append datasets into the list
for i in range(len(list_of_names)):
    df = pd.read_csv(path+list_of_names[i]+".csv", delimiter=',', quotechar='"', index_col='sample_ID')
    phenotype_list.append(df)

# define phenotype train and testing list 
phenotype_train_list = [phenotype_list[0], phenotype_list[2], phenotype_list[4], 
                      phenotype_list[6], phenotype_list[8]]

phenotype_test_list = [phenotype_list[1], phenotype_list[3], phenotype_list[5], 
                     phenotype_list[7], phenotype_list[9]]

#for df1 in phenotype_train_list:

df1=pd.read_csv(path+"pheno_sample_patients_train_split_1.csv", dtype=object, comment="#")
grade54 = df1[(df1['grading']=="5")]['sample'].tolist() + df1[(df1['grading']=="4")]['sample'].tolist() 
grade3 = df1[(df1['grading']=="3")]['sample'].tolist()
grade543 = grade54 + grade3
grade0 = df1[(df1['grading']=="0")]['sample'].tolist()


#print (grade4)
#print (grade0)
path="/data/DB/Covid19_WES/NEW_SET2/"
grade543_col=[s+"_GT" for s in grade543]
grade0_col=[s+"_GT" for s in grade0]

column_filter=['Chr','Start', 'End', 'Ref', 'Alt', 'Gene.refGene', 'ExonicFunc.refGene', 'AAChange.refGene', 'CADD_phred']
for l in open("../GT_columns.txt", "rt"):
    col=l.strip("\n")
    if col in grade543_col or col in grade0_col:
        column_filter.append(col)
fmask = os.path.join(path, "*.csv.gz")
###Listing all the files in the specified path ending with *.featureCounts-GENCODE.txt
 ###Reading dataframe
for infile in glob.glob(path+"/*.csv.gz"):
    df=pd.read_csv(infile, dtype=object, comment="#",  usecols=column_filter)
    ###Ge-tting the column names corresponding to each individual groups (e.g. in this case grade0 and grade4)
    li.append(df)
    df = pd.concat(li, axis=0, ignore_index=True)
#print(df)
df = pd.DataFrame(df)
# import the significant filtered variants for grade 43 vs0
path="/data/DB/Covid19_WES/NEW_SET2/variant_stat_grading_improved_version_5-fold_CV_2000_cohort_latest/"


df1=pd.read_csv(path+"filter_variant_stat_grading_improved_trainset_fold_1_CV.csv", dtype=object, comment="#")
df_new = pd.DataFrame(df1, columns=['Chr', 'Start', 'End', 'Ref', 'Alt', 'Gene.refGene'])

# import the training and test set feature count matrices without covariates phenotypic features 
#df_new.head()
df = pd.merge(df, df_new, on=['Chr', 'Start', 'End', 'Ref', 'Alt', 'Gene.refGene'], how='inner')
#select columns that ends with genotype ('_GT') from the dervied variant dataset 
df = df.loc[:, df.columns.str.endswith('_GT')]

# remove _GT from columns and transpose columns to rows 
df.columns = df.columns.str.rstrip('_GT')
df = pd.concat([df_new['Gene.refGene'], df], axis=1)
df = df.T 

df.columns = df.iloc[0]
df = df.iloc[1:]
#print the feature count matrix
#Creating a dictionary 
#mapping = {'0/0': 0, '0/1': 1, './.':0, 'somatic':0, '1|1':2, '0|1':1, '1/1':2, '0|0':0}
# values where key matches 
df[df == '0/0'] = 0
df[df == '0/1'] = 1
df[df == './.'] = 0
df[df == '.|.'] = 0
df[df == 'somatic'] = 0
df[df == '1|1'] = 2
df[df == '0|1'] = 1
df[df == '1/0'] = 1

df[df == '0|0'] = 0
df[df == '1/1'] = 2

path="/data/DB/Covid19_WES/NEW_SET2/variant_stat_grading_improved_version_5-fold_CV_2000_cohort_latest/"

df.to_csv(path + "feature_count_fold_1_train.csv", index = True)
