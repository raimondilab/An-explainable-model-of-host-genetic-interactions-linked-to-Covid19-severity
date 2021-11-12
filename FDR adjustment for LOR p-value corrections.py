#### (2.1) : We used the fisher's exact test to make adjustment (corrections) on the p-value column of the extracted variants from (2) for each of the stratified 5-folds. 

# Note this was re-run five-times to make corrections in each stratified 5-foldCV extracted variants p-values.


#!/data/SW/anaconda3/envs/myenv/bin/python

# Import relevant libraries 
import os,sys,operator
import gzip, glob
import pickle, math
import pandas as pd
from scipy.stats import chisquare, fisher_exact
from statsmodels.stats.contingency_tables import Table2x2
from statsmodels.stats.multitest import fdrcorrection

fs="\t"
path="/data/DB/Covid19_WES/NEW_SET2/"

df1=pd.read_csv(path+"pheno_sample_patients_train_split_1.csv", dtype=object, comment="#")
grade54 = df1[(df1['grading']=="5")]['sample'].tolist() + df1[(df1['grading']=="4")]['sample'].tolist()
grade3 = df1[(df1['grading']=="3")]['sample'].tolist()

grade543 = grade54 + grade3
grade0 = df1[(df1['grading']=="0")]['sample'].tolist()


df=pd.read_csv("variant_stat_grading_improved_trainset_fold_1.tsv", dtype=object, sep="\t")
df["adjP"]=fdrcorrection(df["LOR Pvalue"].astype(float))[1]

#df.to_csv(path + "variant_stat_grading_improved_trainset_fold_1_adj.tsv", index=False, sep="\t")
