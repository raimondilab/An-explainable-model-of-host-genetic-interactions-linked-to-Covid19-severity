#### 2.2.1) filter significant variants from 2.1 using the cut-off = 6.0, p-value<0.05 and |LOR| ≥ 1.0 

# Note this was re-run five-times to extract significant fitered variants from each of the stratified 5-foldCV. 


#!/data/SW/anaconda3/envs/myenv/bin/python

# import relevant libraries 
import os,sys,operator
import gzip, glob
import pickle, math
import pandas as pd
from scipy.stats import chisquare, fisher_exact

fs="\t"
path="/data/DB/Covid19_WES/NEW_SET2/"
li = []
lor_cutoff=6.0
for infile in glob.glob(path+"/*.tsv.gz"): # zipped tsv FDR correction file for each stratified 5-foldCV training set.
    df = csv.reader(infile, delimiter="\t")
    l = li.append(df)
    #if (float(l.split(fs)[10]) > lor_cutoff or float(l.split(fs)[10]) < -lor_cutoff): # We tried several cut-offs
    if l.split(fs)[18] != "LOR Pvalue" and float(l.split(fs)[18]) < 0.05 and abs(float(l.split(fs)[17])) > 1 and l.split(fs)[6] != "synonymous SNV":
    #if l.split(fs)[17] != "LOR Pvalue" and float(l.split(fs)[17]) < 0.001 and abs(float(l.split(fs)[16])) > 0.5 and l.split(fs)[6] != "synonymous SNV":
    #if l.split(fs)[21].find("adjP") == -1 and float(l.split(fs)[21]) < 0.1:
        #print (l.split(fs)[5].split(";")[0]+fs+l.split(fs)[8]+fs+l.split(fs)[9]+fs+l.split(fs)[10])
        print (l.strip("\n"))
