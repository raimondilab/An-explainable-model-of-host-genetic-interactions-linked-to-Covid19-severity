#### 7.1) Pathway analysis: This script was used to plot bar-chart of pathway enrichment clusters (non-zero variants from stratified 5-foldCVs decision tree like models) 2000 cohort. 

#!/data/SW/anaconda3/envs/myenv/bin/python


import os, sys
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




df=pd.read_csv("variant_stat_xgb_w1_genes_pathways_modules.txt", dtype=object, sep="\t")


#print (df['Pathway name'])
#print (df['Pathway name'], -np.log10(df['Entities FDR_x'].astype(float)),  -np.log10(df['Entities FDR_y'].astype(float)))
cutoff=15
#df_sort=df.sort_values(by = 'FDR_x') 
df["FDR"] = df['FDR'].str.replace(',', '.')

df1=df.loc[df['Module'] == str(sys.argv[1])]

path1=df1['GeneSet'][:cutoff].iloc[::-1]
path1=[p[:-3] for p in path1]

adjp_1=-np.log10(df1['FDR'].astype(float))[:cutoff].iloc[::-1]

path=[]
adjp=[]
for p,ap in zip(path1, adjp_1):
  #p=path1[ii].upper()
  #ap=adjp_1[ii]
  if p not in path:
    path.append(p)
    adjp.append(ap)

print (path)


#print (df1["FDR"])

#print (adjp_1)
#for ii in range(len(path)):
#	print (path[ii], adjp_1[ii], adjp_2[ii])

ind1 = np.arange(len(path))  # the x locations for the groups
#print (len(ind1), len(ind2))
width = 0.7      # the width of the bars
fig = plt.figure(figsize=(20,15))

ax0 = fig.add_axes([0.55, 0.1,0.4,0.8])
ax0.barh(ind1, adjp, width, color="gray", edgecolor="gray", linewidth=0.1)

major_ticksy = np.arange(0, len(path), 1)
minor_ticksy=[]
for t in major_ticksy:
  minor_ticksy.append(t+0.5)                                            
ax0.set_yticks(major_ticksy)

ylab=[]
for p in path:
  print (p)
  if p.find('HDR through Homologous Recombination (HRR) or Single Strand Annealing (SSA)') != -1:
    ylab.append(p.replace('HDR through Homologous Recombination (HRR) or Single Strand Annealing (SSA)', 'HDR through HRR or SSA'))
  else:
    ylab.append(p)

#ax0.set_yticklabels(ylab, horizontalalignment = "left")
ax0.set_yticklabels(ylab)

ax0.tick_params(axis='x', which='major', labelsize=27.5, pad=15)
ax0.tick_params(axis='y', which='major', labelsize=40, direction="in")
ax0.spines['top'].set_visible(True)
ax0.spines['right'].set_visible(False)
ax0.spines['bottom'].set_visible(False)
ax0.spines['top'].set_linewidth(3)
ax0.spines['left'].set_visible(True)
ax0.spines['left'].set_linewidth(3)
ax0.yaxis.set_ticks_position('left')
ax0.xaxis.set_ticks_position('top')

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
plt.savefig('module_%s.eps' % (sys.argv[1]))
plt.savefig('module_%s.png' % (sys.argv[1]))
plt.savefig('module_%s.pdf' % (sys.argv[1]), format="pdf", bbox_inches="tight")


