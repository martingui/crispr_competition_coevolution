#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 19 10:10:42 2020

@author: guillemet
"""

import pandas as pd
import os
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import math
import copy
from matplotlib.ticker import MaxNLocator
from sklearn.linear_model import LinearRegression
os.chdir('/home/guillemet/Documents/crispr/work/Martin/scripts')
colors= sns.color_palette()
#########frequency of mutation in protospacer


#### R

dfpam=pd.read_csv("../../Aarhus/data/PAM_positions.csv")
dfmuts=pd.DataFrame()
for i in range(1,9):
    dfmut=pd.read_csv('../data/FreeBayes/R_seq/R'+str(i)+'_data.csv', delimiter='\t')
    dfmut.TIME=dfmut.TIME.apply(lambda x: 'T'+x[1] if 'R' in x else x )
    dfmut['rep']='R'+str(i)
    dfmuts=pd.concat([dfmuts,dfmut])
dfmuts['escape_pos']=-1
for i in range(len(dfmuts)):
    pos=dfmuts.POS.iloc[i]
    filt_df=dfpam[((dfpam.escape_start+1<=pos) & (dfpam.escape_end>=pos)) | ((dfpam.escape_start+1>=pos) & (dfpam.escape_end<=pos))]
    if len(filt_df)==1:
        start=filt_df.escape_start.iloc[0]
        end=filt_df.escape_end.iloc[0]
        if filt_df.strand.iloc[0] == '+':
            dfmuts.escape_pos.iloc[i] = end - pos + 1
        elif filt_df.strand.iloc[0] == '-':
            dfmuts.escape_pos.iloc[i] = pos - start -1
    elif len(filt_df)>1:
        dfmuts.escape_pos.iloc[i] =-0.1
        
ndfmuts=dfmuts[dfmuts.escape_pos>=0]
ndfmuts=ndfmuts.loc[ndfmuts.FREQ>0.1]
sns.displot(data=ndfmuts, x='escape_pos', bins=[x+0.5 for x in range(0,38)])

dfmutsr=copy.deepcopy(dfmuts)

tndfmuts=ndfmuts.loc[ndfmuts.escape_pos>25]
tdfmuts=dfmuts.loc[dfmuts.escape_pos>25]


#### W

dfpam=pd.read_csv("../../Aarhus/data/PAM_positions.csv")
dfmuts=pd.DataFrame()
for i in range(1,9):
    dfmut=pd.read_csv('../data/FreeBayes/W_seq/W'+str(i)+'_data.csv', delimiter='\t')
    dfmut.TIME=dfmut.TIME.apply(lambda x: 'T'+x[1] if 'R' in x else x )
    dfmut['rep']='W'+str(i)
    dfmuts=pd.concat([dfmuts,dfmut])
dfmuts['escape_pos']=-1
for i in range(len(dfmuts)):
    pos=dfmuts.POS.iloc[i]
    filt_df=dfpam[((dfpam.escape_start+1<=pos) & (dfpam.escape_end>=pos)) | ((dfpam.escape_start+1>=pos) & (dfpam.escape_end<=pos))]
    if len(filt_df)==1:
        start=filt_df.escape_start.iloc[0]
        end=filt_df.escape_end.iloc[0]
        if filt_df.strand.iloc[0] == '+':
            dfmuts.escape_pos.iloc[i] = end - pos + 1
        elif filt_df.strand.iloc[0] == '-':
            dfmuts.escape_pos.iloc[i] = pos - start -1
    elif len(filt_df)>1:
        dfmuts.escape_pos.iloc[i] =-0.1
        
ndfmuts=dfmuts[dfmuts.escape_pos>=0]
ndfmuts=ndfmuts.loc[ndfmuts.FREQ>0.1]
sns.displot(data=ndfmuts, x='escape_pos', bins=[x+0.5 for x in range(0,38)])
            
dfmutsw=copy.deepcopy(dfmuts)
[31149,1209,971,25461,7037,27013,34608,31065,31725,21039,24343,3233,29998,23084,16236,30386]


#### R0

dfpam=pd.read_csv("../../Aarhus/data/PAM_positions.csv")
dfmuts=pd.read_csv('../data/FreeBayes/Other/T0_R_data.csv', delimiter='\t')
dfmut.TIME=dfmut.TIME.apply(lambda x: 'T'+x[1] if 'R' in x else x )
dfmuts['escape_pos']=-1
for i in range(len(dfmuts)):
    pos=dfmuts.POS.iloc[i]
    filt_df=dfpam[((dfpam.escape_start+1<=pos) & (dfpam.escape_end>=pos)) | ((dfpam.escape_start+1>=pos) & (dfpam.escape_end<=pos))]
    if len(filt_df)==1:
        start=filt_df.escape_start.iloc[0]
        end=filt_df.escape_end.iloc[0]
        if filt_df.strand.iloc[0] == '+':
            dfmuts.escape_pos.iloc[i] = end - pos + 1
        elif filt_df.strand.iloc[0] == '-':
            dfmuts.escape_pos.iloc[i] = pos - start -1
    elif len(filt_df)>1:
        dfmuts.escape_pos.iloc[i] =-0.1
dfmuts0r=copy.deepcopy(dfmuts)


#### W0

dfpam=pd.read_csv("../../Aarhus/data/PAM_positions.csv")
dfmuts=pd.read_csv('../data/FreeBayes/Other/T0_W_data.csv', delimiter='\t')
dfmut.TIME=dfmut.TIME.apply(lambda x: 'T'+x[1] if 'R' in x else x )
dfmuts['escape_pos']=-1
for i in range(len(dfmuts)):
    pos=dfmuts.POS.iloc[i]
    filt_df=dfpam[((dfpam.escape_start+1<=pos) & (dfpam.escape_end>=pos)) | ((dfpam.escape_start+1>=pos) & (dfpam.escape_end<=pos))]
    if len(filt_df)==1:
        start=filt_df.escape_start.iloc[0]
        end=filt_df.escape_end.iloc[0]
        if filt_df.strand.iloc[0] == '+':
            dfmuts.escape_pos.iloc[i] = end - pos + 1
        elif filt_df.strand.iloc[0] == '-':
            dfmuts.escape_pos.iloc[i] = pos - start -1
    elif len(filt_df)>1:
        dfmuts.escape_pos.iloc[i] =-0.1
        
dfmuts0w=copy.deepcopy(dfmuts)

##### Filter: not in t0, > 5 reads, no duplicates
    

## Plot phage mutations along genome
dfmutsall=pd.concat([dfmutsw, dfmutsr])

dfmutsr=dfmutsr.loc[dfmutsr.TYPE!='complex']
dfmutsw=dfmutsw.loc[dfmutsw.TYPE!='complex']
x=1

for dfs in [[dfmutsw, dfmuts0w, 'W', colors[1]], [dfmutsr, dfmuts0r, 'R', colors[3]]]:
    dfmutsall=dfs[0]
    dfmutsall=dfmutsall.reset_index(drop=True)
    dfmuts0=dfs[1]
    dfmuts0=dfmuts0.reset_index(drop=True)
    newindex=[]
    noindex=[]
    for i in dfmutsall.index:
        for iii in dfmuts0.index:
            if (dfmutsall.ALT.loc[i] == dfmuts0.ALT.loc[iii]) and (dfmutsall.REF.loc[i] == dfmuts0.REF.loc[iii]) and (dfmutsall.POS.loc[i] == dfmuts0.POS.loc[iii]):
                noindex+=[i]
            else:
                if i not in newindex:
                    newindex+=[i]
                    
    newindex=[x for x in newindex if x not in noindex]
    
    dfnewmuts=dfmutsall.loc[newindex]
    dfnewmuts.POS=dfnewmuts.POS/1000
    dfnewmuts=dfnewmuts.drop_duplicates(subset=['ALT','REF','POS'])
    dfnewmuts=dfnewmuts.loc[dfnewmuts.AO>5]
    
    fig, ax=plt.subplots()
    sns.set_style('ticks')
    sns.histplot(data=dfnewmuts, x='POS', kde=False, bins=range(0,37,2), color=dfs[3], ax=ax, alpha=0.35)
    sns.histplot(data=dfnewmuts.loc[dfnewmuts.escape_pos>=-0.1], x='POS', kde=False, bins=range(0,37,2), color=dfs[3], ax=ax)
    sns.set_style('ticks')
    plt.plot(windf.center_x, windf.y, 'k--', linewidth=1)
    plt.ylabel('Number of mutations per 2-kb', fontsize=12)
    plt.xlabel('Position in the phage genome (kb)', fontsize=12)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    sns.despine()
    plt.savefig('../steps/phage_mutation/along_genome_'+dfs[2]+'.png', dpi=300, bbox_inches='tight')
    plt.savefig('/home/guillemet/Documents/crispr/final_figures/S5_'+str(x)+'.png', dpi=300, bbox_inches='tight')
    
    plt.figure()
    sns.set_style('ticks')
    sns.displot(data=dfnewmuts, x='escape_pos', bins=[x+0.5 for x in range(0,38)], color=dfs[3])
    plt.xlabel('Position in the protospacer', fontsize=12)
    plt.ylabel('Number of mutations', fontsize=12)
    plt.plot([0.5,4.5],[-0.1,-0.1], 'k')
    plt.plot([6.5,27.5],[-0.1,-0.1], 'k')
    plt.text(2.5, -0.5, 'PAM', ha='center', va='center', fontsize=12)
    plt.text(17, -0.5, 'Seed', ha='center', va='center', fontsize=12)
    sns.despine()
    plt.ylim(-0.8,15)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.savefig('../steps/phage_mutation/along_protospacer_'+dfs[2]+'.png', dpi=300, bbox_inches='tight')
    plt.savefig('/home/guillemet/Documents/crispr/final_figures/S6_'+str(x)+'.png', dpi=300, bbox_inches='tight')
    x=2

    

