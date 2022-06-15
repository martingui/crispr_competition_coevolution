#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 11 10:12:24 2021

@author: guillemet
"""
import pandas as pd
import os
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import math
import copy
os.chdir('/home/guillemet/Documents/crispr/work/Martin/scripts')


colors= sns.color_palette()
sns.set_style("white")

bacter_all=pd.read_csv('../data/nBacteria_genos_filled.csv', sep=',', header=0) 
bacter_all=bacter_all.iloc[:,1:]
phage_bacterR=pd.read_csv("../steps/df/phage_bacterR.csv")
phage_bacterW=pd.read_csv("../steps/df/phage_bacterW.csv")
phage_bacterW.loc[phage_bacterW.f_variant>1, 'f_variant']=1
phage_bacterC=pd.read_csv("../steps/df/phage_bacterC.csv")

phage_bacter=pd.concat([phage_bacterR, phage_bacterW])




phagesW1=pd.read_csv('../data/phages/W1.tsv', sep='\t', header=0) 
phagesW1['rep']='W1'
phagesW1['cond']='W'
phagesW2=pd.read_csv('../data/phages/W2.tsv', sep='\t', header=0) 
phagesW2['rep']='W2'
phagesW2['cond']='W'
phagesW3=pd.read_csv('../data/phages/W3.tsv', sep='\t', header=0) 
phagesW3['rep']='W3'
phagesW3['cond']='W'
phagesW4=pd.read_csv('../data/phages/W4.tsv', sep='\t', header=0) 
phagesW4['rep']='W4'
phagesW4['cond']='W'
phagesW5=pd.read_csv('../data/phages/W5.tsv', sep='\t', header=0) 
phagesW5['rep']='W5'
phagesW5['cond']='W'
phagesW6=pd.read_csv('../data/phages/W6.tsv', sep='\t', header=0) 
phagesW6['rep']='W6'
phagesW6['cond']='W'
phagesW7=pd.read_csv('../data/phages/W7.tsv', sep='\t', header=0) 
phagesW7['rep']='W7'
phagesW7['cond']='W'
phagesW8=pd.read_csv('../data/phages/W8.tsv', sep='\t', header=0) 
phagesW8['rep']='W8'
phagesW8['cond']='W'
phagesR1=pd.read_csv('../data/phages/R1.tsv', sep='\t', header=0) 
phagesR1['rep']='R1'
phagesR1['cond']='R'
phagesR2=pd.read_csv('../data/phages/R2.tsv', sep='\t', header=0) 
phagesR2['rep']='R2'
phagesR2['cond']='R'
phagesR3=pd.read_csv('../data/phages/R3.tsv', sep='\t', header=0) 
phagesR3['rep']='R3'
phagesR3['cond']='R'
phagesR4=pd.read_csv('../data/phages/R4.tsv', sep='\t', header=0) 
phagesR4['rep']='R4'
phagesR4['cond']='R'
phagesR5=pd.read_csv('../data/phages/R5.tsv', sep='\t', header=0) 
phagesR5['rep']='R5'
phagesR5['cond']='R'
phagesR6=pd.read_csv('../data/phages/R6.tsv', sep='\t', header=0) 
phagesR6['rep']='R6'
phagesR6['cond']='R'
phagesR7=pd.read_csv('../data/phages/R7.tsv', sep='\t', header=0) 
phagesR7['rep']='R7'
phagesR7['cond']='R'
phagesR8=pd.read_csv('../data/phages/R8.tsv', sep='\t', header=0) 
phagesR8['rep']='R8'
phagesR8['cond']='R'
phages_all=pd.concat([phagesW1, phagesW2, phagesW3, phagesW4, phagesW5, phagesW6, phagesW7, phagesW7, phagesW8
                      , phagesR1, phagesR2, phagesR3, phagesR4, phagesR5, phagesR6, phagesR7, phagesR8])

def adaptation(r, df, phag):
    df=df[df.freq>0]
    df=df[df.Rep==r]
    adapts=[]
    for time in range(0,5):
        pha=phag[phag.time==time]
        pha=pha[pha.rep==r]
        adapt=0
        tdf=df[df.Time==time]
        for gen in set(tdf.Genotype):
            g=gen.split('-')[1:-1]
            if len(g) == 0:
                pass
            elif len(g) == 1:
                sp=int(g[0].split('_')[1])
                adapt+= tdf.loc[tdf.Genotype==gen, 'freq'].iloc[0] * pha[pha.spacer==sp].f_variant.iloc[0]
            else:
                toadd=tdf.loc[tdf.Genotype==gen, 'freq'].iloc[0]
                for sp in g:
                    sp=int(sp.split('_')[1])
                    toadd=toadd*pha[pha.spacer==sp].f_variant.iloc[0]
                adapt+=toadd
        adapts+=[adapt]
    return(adapts)




#R adaptation  
value=[]
times=[]
groups=[]
for i in ['R1','R2','R3','R4','R5','R6','R7','R8']:
    value+=adaptation(i, bacter_all, phage_bacter)
    times+=[0,1,2,3,4]
    groups+=[i]*5
dfr=pd.DataFrame({'value':value, 'times':times, 'groups':groups}) 

#W adaptation  
value=[]
times=[]
groups=[]
for i in ['W1','W2','W3','W4','W5','W6','W7','W8']:
    value+=adaptation(i, bacter_all, phage_bacter)
    times+=[0,1,2,3,4]
    groups+=[i]*5
dfw=pd.DataFrame({'value':value, 'times':times, 'groups':groups})  
dfall1=pd.concat([dfr,dfw])

#little for spatial adapt
dfr_here=dfr
dfw_here=dfw



def spatial_adaptation(r, df, phages):
    df=df[df.freq>0]
    df=df[df.Rep==r]
    adapts=[]
    for time in range(0,5):
        pha=phages
        pha=pha[pha.rep!=r]
        pha=pha[pha.cond==r[0]]
        adapt=0
        tdf=df[df.Time==time]
        for gen in set(tdf.Genotype):
            g=gen.split('-')[1:-1]
            if len(g) == 0:
                pass
            elif len(g) == 1:
                sp=int(g[0].split('_')[1])
                toadd=tdf.loc[tdf.Genotype==gen, 'freq'].iloc[0] 
                tomup=0
                for proto in set(pha.site):
                    if abs(proto-sp)<18:
                        tomup+=pha[pha.site==proto]["f_"+str(time)].sum()/7
                        if tomup>1:
                            print(sp)
                            print(tomup)
                adapt+= toadd * tomup
            else:
                toadd=tdf.loc[tdf.Genotype==gen, 'freq'].iloc[0]
                for sp in g:
                    tomup=0
                    sp=int(sp.split('_')[1])
                    for proto in set(pha.site):
                        if abs(proto-sp)<18:
                            tomup+=pha[pha.site==proto]["f_"+str(time)].sum()/7
                    #print(toadd, tomup)              
                    toadd=toadd*tomup
                adapt+=toadd

        adapts+=[adapt]
    return(adapts)
    

a=spatial_adaptation('R2', bacter_all, phages_all)


#R adaptation  
value=[]
times=[]
groups=[]
for i in ['R1','R2','R3','R4','R5','R6','R7','R8']:
    value+=spatial_adaptation(i, bacter_all, phages_all)
    times+=[0,1,2,3,4]
    groups+=[i]*5
dfr_away=pd.DataFrame({'value':value, 'times':times, 'groups':groups}) 


#W adaptation  
value=[]
times=[]
groups=[]
for i in ['W1','W2','W3','W4','W5','W6','W7','W8']:
    value+=spatial_adaptation(i, bacter_all, phages_all)
    times+=[0,1,2,3,4]
    groups+=[i]*5
dfw_away=pd.DataFrame({'value':value, 'times':times, 'groups':groups})  


dfr_here["where"]='here'
dfw_here["where"]='here'
dfr_away["where"]='away'
dfw_away["where"]='away'

dfw_spatial=pd.concat([dfw_here,dfw_away])
dfr_spatial=pd.concat([dfr_here,dfr_away])

plt.figure()
sns.lineplot(x=dfw_spatial.times, y=dfw_spatial.value, hue=dfw_spatial["where"])
plt.ylim([0,1])
plt.ylabel("Phage adaptation")
plt.xlabel("Days")
plt.legend(loc='upper left')
plt.title("Spatial adaptation in M populations")
plt.xticks(np.arange(0, 5, step=1))
#plt.savefig('../steps/Adaptation/spatial_adaptation_M.pdf')

plt.figure()
sns.lineplot(x=dfr_spatial.times, y=dfr_spatial.value, hue=dfr_spatial["where"])
plt.ylim([0,1])
plt.ylabel("Phage adaptation")
plt.xlabel("Days")
plt.legend(loc='upper left')
plt.title("Spatial adaptation in D populations")
plt.xticks(np.arange(0, 5, step=1))
#plt.savefig('../steps/Adaptation/spatial_adaptation_D.pdf')

df_spatialr=copy.deepcopy(dfr_spatial.loc[dfr_spatial['where']=='here'])
df_spatialr=df_spatialr.drop(columns=['where'])
df_spatialr.columns=['here','time','group']
df_spatialr['away']=dfr_spatial.loc[dfr_spatial['where']=='away']['value']
df_spatialr['diff']=df_spatialr['here']-df_spatialr['away']
df_spatialr['cond']='D'

df_spatialw=copy.deepcopy(dfw_spatial.loc[dfw_spatial['where']=='here'])
df_spatialw=df_spatialw.drop(columns=['where'])
df_spatialw.columns=['here','time','group']
df_spatialw['away']=dfw_spatial.loc[dfw_spatial['where']=='away']['value']
df_spatialw['diff']=df_spatialw['here']-df_spatialw['away']
df_spatialw['cond']='M'




df_spatial_all=pd.concat([df_spatialr, df_spatialw])
plt.figure(figsize=[5.4,5], dpi=300)
sns.set_style('ticks')
#sns.catplot(data=df_spatial_all, x='time', y='diff', hue='cond',hue_order=['M','D'], palette=['darkorange','darkred'], \
#            linestyles=["--", "--"], kind="point", dodge=True)
sns.lineplot(data=df_spatial_all.loc[(df_spatial_all.cond=='M') | (df_spatial_all.cond=='D')], x='time', y='diff', hue='cond',hue_order=['M','D'], palette=[colors[1],colors[3]], ci=95, n_boot=1000, seed=10)
plt.ylabel("Phage local adaptation", fontsize=12)
plt.xlabel("Time (days)", fontsize=12)
plt.legend(loc='center')
plt.ylim(-0.35, 0.7)
ax=plt.gca()
ax.get_legend().remove()
'''handles, labels = plt.gca().get_legend_handles_labels()
plt.legend(handles=handles[1:], labels=labels[1:])'''
plt.xticks(np.arange(0, 5, step=1))
plt.plot([0,4], [0,0], 'k:')
sns.despine()
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
fig=plt.gcf()
plt.subplots_adjust(left=0.151)
plt.savefig('../steps/Adaptation/spatial_adaptation_diff_ppt.pdf')
plt.savefig('/home/guillemet/Documents/crispr/final_figures/6_1.pdf', bbox_inches='tight')



