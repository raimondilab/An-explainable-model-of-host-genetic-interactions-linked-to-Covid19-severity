### 3.2) 
'''We used this script to extract relevant information (Chr, Ref, Alt, genotypes) from zipped VCF files of the 24 chromosomes (chromosome 1 - Y) 3000 cohort 
follow-up dataset that can further be useful in feature matrix generation. The extracted information were merged with their relevant zipped CSV information files for all the 24 chromosomes.''' 


#!/data/SW/anaconda3/envs/myenv/bin/python
import os,sys,operator
import gzip, glob
import pickle, math
import re
from functools import partial
import warnings
import vcf
import pandas as pd
from scipy.stats import chisquare, fisher_exact
from statsmodels.stats.contingency_tables import Table2x2
from pandas.api.types import CategoricalDtype
import numpy as np
import pathlib

# import scikit-allel
import allel
import sys
import pandas as pd
from gzip import open as gzopen
from vcf_to_dataframe import vcf_to_dataframe
import numpy as np

fs="\t"
path="/data/DB/Covid19_WES/NEW_SET3/"

# Select genotype present in phenotype information 
genotype = []
for l in open("../column_GT_new_grading_5.txt", "rt"):
    as_list = l.split(", ")
    genotype.append(as_list[0].replace("\n", ""))

for infile in glob.glob(path+"joint_chrY.annovar.vcf.gz"):
    df = vcf_to_dataframe(infile, keep_samples = genotype,keep_format_data=False)
    
    df.rename(columns = {'chrom':'Chrom'}, inplace = True)

#merge converted vcf to dataframe file and corresponding csv files 
column_filter=['Chr','Start', 'End', 'Ref', 'Alt', 'Gene.refGene', 'ExonicFunc.refGene', 'AAChange.refGene', 'CADD_phred']
for infile in glob.glob(path+"joint_chrY.annovar.csv.gz"):
    #print (infile)
    #continue
    ###Reading dataframe
    df_1=pd.read_csv(infile, dtype=object, comment="#", usecols=column_filter)
    # merge dataframes 
    df_new = pd.concat([df_1, df], axis=1)
    
    # delete irrelevant columnns 
   
    df_new = df_new.drop(['Chrom', 'pos', 'ref', 'alt', 'info', 'format' ], axis=1)
    pd.set_option('display.max_rows', None)
    #print(df_new)
    path="/data/DB/Covid19_WES/NEW_SET3/variant_stat_extractions_new_data/"
    #df_new.to_csv(path+ "joint_chrY.annovar_merged.csv")
