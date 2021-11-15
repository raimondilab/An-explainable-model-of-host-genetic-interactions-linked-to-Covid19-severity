### 4) 
''' We used this script to generate the feature matrix from 3000 cohort follow-up dataset. Here we considered only the handful identified 16 full supported variants that consistently maintained non-zero weights 
from decision tree like models (Random Forest and XGBoost) classifier in the stratified 5-fold 2000 cohort analyses. We generated this feature matrix to validate the saved models from 2000 cohort analysis'''

# Note: this script was run once to generate the feature matrix to validate the saved models from 2000 cohort study. 

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
li = []
path="/data/DB/Covid19_WES/NEW_SET3/variant_stat_extractions_new_data/"

df1=pd.read_csv(path+"phenotype_adjustedbyage_based_grading_5_3000_cohort_only.csv", dtype=object, comment="#")
grade54 = df1[(df1['grading_5']=="5")]['sample'].tolist() + df1[(df1['grading_5']=="4")]['sample'].tolist() 
grade3 = df1[(df1['grading_5']=="3")]['sample'].tolist()
grade543 = grade54 + grade3
grade0 = df1[(df1['grading_5']=="0")]['sample'].tolist()


#print (grade4)
#print (grade0)
path="/data/DB/Covid19_WES/NEW_SET3/variant_stat_extractions_new_data/"
#grade543_col=[s+"_GT" for s in grade543]
#grade0_col=[s+"_GT" for s in grade0]

column_filter=['Chr','Start', 'End', 'Ref', 'Alt', 'Gene.refGene', 'ExonicFunc.refGene', 'AAChange.refGene', 'CADD_phred']
for l in open("../column_GT_new_grading_5.txt", "rt"):
    col=l.strip("\n")
    if col in grade543 or col in grade0:
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
path="/data/DB/Covid19_WES/NEW_SET3/variant_stat_extractions_new_data/"


df1=pd.read_csv(path+"filter_16_best_fully_supported_variants.csv", dtype=object, comment="#")
df_new = pd.DataFrame(df1, columns=['Chr', 'Start', 'End', 'Ref', 'Alt', 'Gene.refGene'])

# import the training and test set feature count matrices without covariates phenotypic features 
#df_new.head()
df = pd.merge(df, df_new, on=['Chr', 'Start', 'End', 'Ref', 'Alt', 'Gene.refGene'], how='inner')
#select columns that ends with genotype ('_GT') from the dervied variant dataset 
df = df.loc[:, df.columns.str.endswith('_hg38')]

# remove _GT from columns and transpose columns to rows 
#df.columns = df.columns.str.rstrip('_hg38')
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

path="/data/DB/Covid19_WES/NEW_SET3/variant_stat_extractions_new_data/"

df.to_csv(path + "feature_count_16_fully_supported_variants_latest.csv", index = True)
