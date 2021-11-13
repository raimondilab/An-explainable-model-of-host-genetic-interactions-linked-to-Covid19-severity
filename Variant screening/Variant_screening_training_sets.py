
#### 2) Variant Screening on the training set

'''We employed a Log-Odds Ratio (LOR) statistics to perform case-control association and to screen variants associated with either severe or asymptomatic patients
in each of the training sets for each of the five folds generated'''

# Note: this script was re-used 5-times to extract filtered training set variants for each stratified 5-fold CV.

#!/data/SW/anaconda3/envs/myenv/bin/python

# import relevant libraries 
import os,sys,operator
import gzip, glob
import pickle, math
import numpy as np
import pandas as pd
from scipy.stats import chisquare, fisher_exact
from statsmodels.stats.contingency_tables import Table2x2

fs="\t"
path="/data/DB/Covid19_WES/NEW_SET2/"

# assign phenotype dataset names
list_of_names = ['pheno_sample_patients_train_split_1','pheno_sample_patients_test_split_1',
                 'pheno_sample_patients_train_split_2','pheno_sample_patients_test_split_2', 
                'pheno_sample_patients_train_split_3','pheno_sample_patients_test_split_3',
                 'pheno_sample_patients_train_split_4','pheno_sample_patients_test_split_4', 
                 'pheno_sample_patients_train_split_5','pheno_sample_patients_test_split_5']
 
# create empty list
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

grade543_col=[s+"_GT" for s in grade543]
grade0_col=[s+"_GT" for s in grade0]

column_filter=['Chr','Start', 'End', 'Ref', 'Alt', 'Gene.refGene', 'ExonicFunc.refGene', 'AAChange.refGene', 'CADD_phred']
for l in open("../GT_columns.txt", "rt"):
    col=l.strip("\n")
    if col in grade543_col or col in grade0_col:
        column_filter.append(col)
#print (column_filter)  

###Odds ration calculation implemented based on the following course:
###http://courses.washington.edu/b516/lectures_2009/Odds_Ratios.pdf
print ('Chr'+fs+'Start'+fs+'End'+fs+'Ref'+fs+'Alt'+fs+'Gene.refGene'+fs+"ExonicFunc.refGene"+fs+"AAChange.refGene"+fs+"CADD_phred"+fs+"AF543"+fs+"AF0"+fs+"LogAFRatio"+fs+"MUT543"+fs+"MUT0"+fs+"WT543"+fs+"WT0"+fs+"Odds-Ratio"+fs+"Log Odds-Ratio"+fs+"LOR Pvalue"+fs+"LOR CI") 
###Listing all the files in the specified path ending with *.featureCounts-GENCODE.txt
path="/data/DB/Covid19_WES/NEW_SET2/"
for infile in glob.glob(path+"*.csv.gz"):
    #print (infile)
    #continue
    ###Reading dataframe
    df=pd.read_csv(infile, dtype=object, comment="#", usecols=column_filter)
    ###Getting the column names corresponding to each individual groups (e.g. in this case grade0 and grade4)
    #print (grade0_col)
    #print (grade43_col)
    ###Iterating over each dataframe's row
    for index, row in df.iterrows():
        AN0=0###Total allele number
        AC0=0###Allele counts
        AF0=0###Allele frequency
        AN543=0###Total allele number
        AC543=0###Allele counts
        AF543=0###Allele frequency
        WT0=0
        WT543=0
        MUT0=0
        MUT543=0
        for col in grade0_col:
            #print (col, row[col])
            ####First calculating the nr of WT and mutated patients
            WT0+=row[col].count("0")
            WT0+=row[col].count("somatic")
            MUT0+=row[col].count("1")
            #if row[col].find("0/0") != -1:
            #    WT0.append(col)
            #if row[col].find("1") != -1:
            #    MUT0.append(col)
            AN0+=row[col].count("0")
            AN0+=row[col].count("1")        ###To get AN I'm counting everything but ./., meaning the gene allele is not there
            AN0+=row[col].count("somatic")
            AC0+=row[col].count("1")###To get AC, I'm counting all the instances where we find the mutated allele, i.e. 1 (at this stage, I'm making no distinction between heterozygotes or homozygites)
            if AN0 == 0:
                AF0=0
            else:
                AF0=float(AC0)/float(AN0)
        for col in grade543_col:
            WT543+=row[col].count("0")
            WT543+=row[col].count("somatic")
            MUT543+=row[col].count("1")
            #if row[col].find("0/0") != -1:
            #    WT4.append(col)
            #if row[col].find("1") != -1:
            #    MUT4.append(col)
            AN543+=row[col].count("0")
            AN543+=row[col].count("1")
            AN543+=row[col].count("somatic")
            AC543+=row[col].count("1")
            if AN543 == 0:
                AF543=0
            else:
                AF543=float(AC543)/float(AN543)
            #print (col, row[col], row[col].count("0"))
        #print (AF4, AF0)
        ####Try both chi2 and Fisher's exact test to evaluate the significance of enrichment
        #chi2, chi2P=chisquare([len(set(MUT4)),len(set(MUT0))],f_exp=[len(set(WT4)),len(set(WT0))])
        #se=math.sqrt((1/MUT4)+(1/MUT0)+(1/WT4)+(1/WT0))
        #OR, FisherP=fisher_exact([[MUT4,MUT0],[WT4,WT0]])
        ctable=Table2x2(np.asarray([[MUT543,MUT0],[WT543,WT0]]))
        lor_CI=ctable.log_oddsratio_confint()
        lor_P=ctable.log_oddsratio_pvalue()
        lor=ctable.log_oddsratio
        OR=ctable.oddsratio
        lAFor=0
        if AF543 == 0:
            AF543=1/(673*2) ### This is a workaround to use the log function even when AFs are 0, i.e. using Lowest possible allele frequency in this dataset
        if AF0 == 0:
            AF0=1/(673*2) 
            
        lAFor=math.log(AF543/AF0)/math.log(2) ###Calculating the Log Odds Ratio
        
        
        print(row['Chr']+fs+row['Start']+fs+row['End']+fs+row['Ref']+fs+row['Alt']+fs+row['Gene.refGene']+fs+str(row["ExonicFunc.refGene"])+fs+str(row["AAChange.refGene"])+fs+str(row["CADD_phred"])+fs+str(AF543)+fs+str(AF0)+fs+str(lAFor)+fs+str(MUT543)+fs+str(MUT0)+fs+str(WT543)+fs+str(WT0)+fs+str(OR)+fs+str(lor)+fs+str(lor_P)+fs+str(lor_CI))


##outdf.to_csv("variant_stat.tsv", index=False)
