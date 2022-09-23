#!/data/SW/anaconda3/envs/myenv/bin/python


import os, sys, math
import pandas as pd
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import*
from matplotlib import pyplot as plt
from mpl_toolkits.axes_grid import make_axes_locatable
from mpl_toolkits import *
from matplotlib import cm
#import adjustText




df=pd.read_csv("merged_disease_variants_updated_sig.csv", dtype=object)


#print (df['Pathway name'])
#print (df['Pathway name'], -np.log10(df['Entities FDR_x'].astype(float)),  -np.log10(df['Entities FDR_y'].astype(float)))
#cutoff=15
#df_sort=df.sort_values(by = 'FDR_x') 
#df["FDR"] = df['FDR'].str.replace(',', '.')

cat="respiratory or thoracic disease"
#cat="immune system disease"
#cat="pancreas disease"
#cat="phenotype"

#df1=df.loc[(df['traiCategory'] == cat) & ((df['XGB_nonzero'].astype(float) > 0) | (df['RF_nonzero'].astype(float) > 0) )]
df1=df.loc[(df['traiCategory'] == cat) & ((df['XGB_nonzero'].astype(float) > 0) )]


#p=-np.log10(df1['pval'].astype(float)).iloc[::-1]
#print (p)
p=-np.log10(df1['pval'].astype(float))
odds=df1['OddsRatio'].astype(float)
lors=df1['Log Odds-Ratio'].astype(float)
nonzero_w=np.array(df1['XGB_nonzero'].astype(float))+np.array(df1['RF_nonzero'].astype(float))
#labels=df1['Gene.refGene'].astype(str)+"\n"+df1['traitReported'].astype(str)
labels=df1['Gene.refGene'].astype(str)
df1['traitReported']=df1['traitReported'].str.replace("\(more controls excluded\)",'')
df1['traitReported']=df1['traitReported'].str.replace("\(mode\)",'')
traits=df1['traitReported'].astype(str)
#path=[]
#adjp=[]
for a1,a2,a3,a4,a5 in zip(p, odds, lors,nonzero_w, labels):
  print(a1,a2,a3,a4,a5.replace("\n",","))


#ind1 = np.arange(len(path))  # the x locations for the groups
#print (len(ind1), len(ind2))
width = 0.7      # the width of the bars
fig = plt.figure(figsize=(15,15))

ax0 = fig.add_axes([0.1, 0.1,0.8,0.8])

print (min(p), max(p))
ax0.scatter(p, odds, s=nonzero_w*100, c=lors, cmap="coolwarm")
for x,y,w,l,t in zip(p,odds,nonzero_w,labels, traits):
  if w>=5: 
    ax0.text(x, y, s=l, fontsize=math.log(w)*15)
    ax0.text(x+0.03, y-0.03, s=t, fontsize=math.log(w)*5, c="gray")
    
  elif x >= 3 and y >=2.5: 
    ax0.text(x, y, s=l, fontsize=math.log(w)*15)
    ax0.text(x+0.03, y-0.03, s=t, fontsize=math.log(w)*5, c="gray")
#major_ticksy = np.arange(0, len(path), 1)
#minor_ticksy=[]
#for t in major_ticksy:
#  minor_ticksy.append(t+0.5)                                            
#ax0.set_yticks(major_ticksy)


#ax0.set_yticklabels(ylab, horizontalalignment = "left")
#ax0.set_yticklabels(ylab)

ax0.tick_params(axis='x', which='major', labelsize=20)
ax0.tick_params(axis='y', which='major', labelsize=20)
#ax0.yaxis.set_ticks_position('left')
#ax0.xaxis.set_ticks_position('top')

plt.xlabel('-log10(PheWAS_Pvalue)', fontsize=30)
plt.ylabel('PheWAS_OddsRatio', fontsize=30)
#print (adjp_1_nonan)
#print (adjp_2_nonan)

#max_range=0
#if max(adjp_1_nonan) > max(adjp_2_nonan):
#	max_range=max(adjp_1_nonan)
#else:
#	max_range=max(adjp_2_nonan)

#ax0.set_xlim([0,max_range])

#ax0.tight_layout()

plt.draw()
plt.savefig('scatter_%s.eps' % (cat))
plt.savefig('scatter_%s.png' % (cat))


