#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 24 14:29:27 2020

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


nbacter_all=pd.read_csv('../data/nBacteria_genos_filled.csv', sep=',', header=0) 
nbacter_all=nbacter_all.iloc[:,1:]

phage_bacterR=pd.read_csv("../steps/df/phage_bacterR.csv")
phage_bacterW=pd.read_csv("../steps/df/phage_bacterW.csv")
phage_bacterC=pd.read_csv("../steps/df/phage_bacterC.csv")

phage_bacter=pd.concat([phage_bacterR, phage_bacterW])

colors= sns.color_palette()


def subf(x):
    try:
        return(x.split('-')[1].split('_')[1])
    except:
        return('0')


sns.set_style("white")

bim=['A7','A8','A9','A10','B7','B8','B9',"C2","C3","C4","C6",'C7','C8','D1','D5','D6']
spacer=[31149,1209,971,25461,7037,27013,34608,31065,31725,21039,24343,3233,29998,23084,16236,30386]
spacers=[31149,1209,971,25461,7037,27013,34608,31065,31725,21039,24343,3233,29998,23084,16236,30386]
eof=[2.5e-5,1.7e-5,2.5e-5,4.8e-5,1.3e-5,1.7e-6,3.8e-6,7.6E-5,1E-6,7.3e-5,3.1e-6,3.6e-4,3.4e-5,4.9e-6,1.9e-5,1.0e-5]
eof_sd=[1.7e-5, 1.2e-5, 0.5e-5, 2.1e-5, 0.6e-5, 2.1e-6, 2.6e-6, 2.5e-5, 0, 1.5e-5, 1.8e-6, 1.9e-4, 1.2e-5, 2.8e-6, 1.9e-5, 1.4e-5]
pe=[0.36458333, 0.33333333, 0.10416667, 0.59722222, 0.27430556, 0.10069444, 0.39583333, 0.24305556,0.02083333, 0.69097222, 0.09722222, 0.90972222, 0.73958333, 0.11805556, 0.90277778, 0.62847222]
myorder=[31725, 1209, 27013, 30386, 3233, 25461, 24343, 29998, 31149, 23084, 21039, 34608, 7037, 16236, 31065, 971]
myordercolor=["#D6404EC8","#E65948C8","#F47346C8","#F99655C8","#FDB567C8","#FDD07DC8","#FEE695C8","#FEF6B1C8","#F8FCB4C8","#EBF7A0C8","#D3ED9BC8","#B4E0A2C8","#91D3A4C8","#6DC4A4C8","#50A9AFC8","#348BBBC8"]
mr_sd=[9.30E-08,6.80E-08,3.10E-08,2.20E-07,2.70E-08,4.60E-09,5.20E-08,6.50E-08,-1.20E-09,2.50E-07,1.00E-08,8.60E-07,2.90E-07,2.10E-08,6.20E-07,1.60E-07]
mr=[2.10E-07,1.90E-07,4.80E-08,4.30E-07,1.50E-07,4.80E-08,2.50E-07,1.30E-07,9.40E-09,6.00E-07,4.60E-08,1.10E-06,7.10E-07,5.70E-08,1.20E-06,5.10E-07]


df_color=pd.DataFrame({'sp':spacer, 'ordersp':myorder, 'eof':eof, 'color':myordercolor, 'pe':pe, 'mr':mr, 'mr_sd':mr_sd, 'ncolor':[-1]*16})
for i in df_color.index:
    sp=df_color.loc[i,'sp']
    ncol=df_color.loc[df_color.ordersp==sp, 'color']
    df_color.loc[i,'ncolor']=ncol.iloc[0]

df_color=df_color.sort_values('mr').reset_index(drop=True)

df=pd.DataFrame({'bim':bim, 'spacer':spacer, 'eof':eof, 'pe':pe, 'eof_sd':eof_sd, 'mr':mr, 'mr_sd':mr_sd})
#df['r']=0
df['w1']=0
df['w2']=0
df['w3']=0
df['w4']=0
df['sdw1']=0
df['sdw2']=0
df['sdw3']=0
df['sdw4']=0
df['fw1']=0
df['fw2']=0
df['fw3']=0
df['fw4']=0

t1w=phage_bacterW[(phage_bacterW.spacer.isin(spacer)) & (phage_bacterW.time==3)]
t1r=phage_bacterR[(phage_bacterR.spacer.isin(spacer)) & (phage_bacterR.time==3)]


spa=[]
means=[]
sds=[]
for time in range(5):
    for sp in spacer:
        t1w=phage_bacterW[(phage_bacterW.spacer.isin(spacer)) & (phage_bacterW.time==time)]
        mean=t1w.loc[t1w.spacer==sp,'f_variant'].mean()
        sd=t1w.loc[t1w.spacer==sp,'f_variant'].std()
        spa+=[sp]
        means+=[mean]
        sds+=[sd]
        df.loc[df.spacer==sp, 'w'+str(time)]=mean
        df.loc[df.spacer==sp, 'sdw'+str(time)]=sd


## sort by mean frequency in control

control=phage_bacterC.loc[phage_bacterC.spacer!=0]
means_control=control.groupby('spacer').freq.mean().sort_values().index

df.spacer = df.spacer.astype("category")
df.spacer.cat.set_categories(means_control, inplace=True)
df = df.sort_values(["spacer"])

i=0
fig, axs = plt.subplots(ncols=16, figsize=(10,2))
for spacer in df.spacer:
    sns.lineplot(x=[0,1,2,3,4], y=df.loc[df.spacer==spacer].loc[:,['w0','w1','w2','w3','w4']].iloc[0].tolist(), color=df_color.loc[df_color.sp==spacer,'ncolor'].iloc[0], ax=axs[i])
    ax=plt.gca()
    axs[i].set(ylim=(0,1.02))
    axs[i].xaxis.set_major_locator(MaxNLocator(integer=True))
    sns.despine()
    if i!=0:
        axs[i].get_yaxis().set_ticks([])
    i+=1
plt.savefig('../steps/spacers/1l/16_acquisitions_mean_control_W.png', dpi=300)


from scipy.stats import spearmanr
from scipy.stats import pearsonr



means_control=control.groupby('spacer').freq.mean().sort_values()
df['mean_variant']=df.loc[:,['w1','w2','w3','w4']].apply(np.mean, axis=1)
df.index=df.spacer
df['mean_freq']=means_control[df.index]
plt.figure()
sns.scatterplot(data=df, x='mean_freq', y='mean_variant', color=colors[1])
ax=plt.gca()
ax.set(xscale='log')
sns.despine()
plt.xlim(0.009,1.05)
plt.ylim(-0.025,1.05)
plt.xlabel('Mean frequency of host in the control', fontsize=12)
plt.ylabel('Mean frequency of phage mutation', fontsize=12)
plt.xticks(ticks=[0.01,0.1,1], labels=[0.01,0.1,1])
ax.set_xticks([0.01,0.1,1])
plt.tick_params(axis='x', bottom=True)
plt.tick_params(left=True)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)

print(spearmanr(a=df.mr, b=df.mean_variant))
print(spearmanr(a=df.mean_freq, b=df.mean_variant))
print(pearsonr(x=np.log(df.mr), y=df.mean_variant))
print(pearsonr(x=np.log(df.mean_freq), y=df.mean_variant))
reg0=LinearRegression().fit(np.log(np.array(df.mean_freq)).reshape((-1,1)), np.array(df.mean_variant))
plt.plot([0.01,1],[reg0.intercept_+math.log(0.01)*reg0.coef_,reg0.intercept_+math.log(1)*reg0.coef_], 'k:')
pear=pearsonr(x=df.mean_freq, y=df.mean_variant)
pear1=pear[1]
#plt.text(0.03,0.8,'Pearson\'s r='+str(round(pear[0],2)), size=12)
#plt.text(0.03,0.69,'p-value='+f"{pear1:.1E}", size=12)
#plt.savefig('../steps/spacers/means/W_mean_freq_log_newfig.png', dpi=300)
plt.savefig('/home/guillemet/Documents/crispr/final_figures/S9_1.pdf', bbox_inches='tight')



means_control=control.groupby('spacer').freq.mean().sort_values()
df['mean_variant']=df.loc[:,['w1','w2','w3','w4']].apply(np.mean, axis=1)
df.index=df.spacer
df['mean_freq']=means_control[df.index]


plt.figure()
sns.scatterplot(data=df, x='mr', y='mean_variant', color=colors[1])
pear=pearsonr(x=df.mr, y=df.mean_variant)
df=df.loc[df.index!=31725]
reg0=LinearRegression().fit(np.log(np.array(df.mr)).reshape((-1,1)), np.array(df.mean_variant))
plt.plot([1e-8,2e-6],[reg0.intercept_+math.log(1e-8)*reg0.coef_,reg0.intercept_+math.log(2e-6)*reg0.coef_], 'k:')


ax=plt.gca()
ax.set(xscale='log')
sns.despine()
plt.xlim(0.8e-8,1.5e-6)
plt.ylim(-0.025,1.05)
plt.xlabel('Protospacer mutation rate', fontsize=12)
plt.ylabel('Mean frequency of phage mutation', fontsize=12)
#plt.xticks(ticks=[0.01,0.1,1], labels=[0.01,0.1,1])
#ax.set_xticks([0.01,0.1,1])
plt.tick_params(axis='x', bottom=True)
plt.tick_params(left=True)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
#plt.text(1e-7,0.8,'Pearson\'s r='+str(round(pear[0],2)), size=12)
#plt.text(1e-7,0.69,'p-value='+str(round(pear[1], 2)), size=12)
plt.savefig('../steps/spacers/means/W_mr_freq_newfig.pdf')
plt.savefig('/home/guillemet/Documents/crispr/final_figures/S9_3.pdf', bbox_inches='tight')

from scipy.stats import spearmanr
from scipy.stats import pearsonr
print(pearsonr(x=df.mr, y=df.mean_variant))
print(pearsonr(x=np.log(df.mean_freq), y=df.mean_variant))







################ R pop
spacer=[31149,1209,971,25461,7037,27013,34608,31065,31725,21039,24343,3233,29998,23084,16236,30386]


df=pd.DataFrame({'bim':bim, 'spacer':spacer, 'eof':eof, 'eof_sd':eof_sd, 'mr':mr, 'mr_sd':mr_sd})
#df['r']=0
df['w1']=0
df['w2']=0
df['w3']=0
df['w4']=0
df['w0']=0
df['sdw1']=0
df['sdw2']=0
df['sdw3']=0
df['sdw4']=0
df['sdw0']=0
df['fw1']=0
df['fw2']=0
df['fw3']=0
df['fw4']=0
df['fw0']=0

t1w=phage_bacterW[(phage_bacterW.spacer.isin(spacer)) & (phage_bacterW.time==3)]



spa=[]
means=[]
sds=[]
for time in range(5):
    for sp in spacer:
        t1w=phage_bacterR[(phage_bacterR.spacer.isin(spacer)) & (phage_bacterR.time==time)]
        mean=t1w.loc[t1w.spacer==sp,'f_variant'].mean()
        sd=t1w.loc[t1w.spacer==sp,'f_variant'].std()
        spa+=[sp]
        means+=[mean]
        sds+=[sd]
        df.loc[df.spacer==sp, 'w'+str(time)]=mean
        df.loc[df.spacer==sp, 'sdw'+str(time)]=sd


#sort by eof DEFAULT

#sort by position
#df=df.sort_values(['spacer']).reset_index(drop=True)


df=df.sort_values('mr')
df=df.reset_index(drop=True)


## sort by mean frequency in control

control=phage_bacterC.loc[phage_bacterC.spacer!=0]
means_control=control.groupby('spacer').freq.mean().sort_values().index

df.spacer = df.spacer.astype("category")
df.spacer.cat.set_categories(means_control, inplace=True)
df = df.sort_values(["spacer"])

i=0
fig, axs = plt.subplots(ncols=16, figsize=(10,2))
for spacer in df.spacer:
    sns.lineplot(x=[0,1,2,3,4], y=df.loc[df.spacer==spacer].loc[:,['w0','w1','w2','w3','w4']].iloc[0].tolist(), color=df_color.loc[df_color.sp==spacer,'ncolor'].iloc[0], ax=axs[i])
    ax=plt.gca()
    axs[i].set(ylim=(0,1.02))
    axs[i].xaxis.set_major_locator(MaxNLocator(integer=True))
    sns.despine()
    if i!=0:
        axs[i].get_yaxis().set_ticks([])
    i+=1
plt.savefig('../steps/spacers/1l/16_acquisitions_mean_control_R.png', dpi=300)

means_control=control.groupby('spacer').freq.mean().sort_values()
df['mean_variant']=df.loc[:,['w1','w2','w3','w4']].apply(np.mean, axis=1)
df.index=df.spacer
df['mean_freq']=means_control[df.index]
plt.figure()
sns.scatterplot(data=df, x='mean_freq', y='mean_variant', color=colors[3])
ax=plt.gca()
ax.set(xscale='log')
sns.despine()
plt.xlim(0.009,1.05)
plt.ylim(-0.025,1.05)
plt.xlabel('Mean frequency of host in the control', fontsize=12)
plt.ylabel('Mean frequency of phage mutation', fontsize=12)
plt.xticks(ticks=[0.01,0.1,1], labels=[0.01,0.1,1])
ax.set_xticks([0.01,0.1,1])
plt.tick_params(axis='x', bottom=True)
plt.tick_params(left=True)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
reg0=LinearRegression().fit(np.log(np.array(df.mean_freq)).reshape((-1,1)), np.array(df.mean_variant))
plt.plot([0.01,1],[reg0.intercept_+math.log(0.01)*reg0.coef_,reg0.intercept_+math.log(1)*reg0.coef_], 'k:')
pear=pearsonr(x=df.mean_freq, y=df.mean_variant)
pear1=pear[1]
#plt.text(0.11,0.8,'Pearson\'s r='+str(round(pear[0],2)), size=12)
#plt.text(0.11,0.69,'p-value='+f"{pear1:.1E}", size=12)
plt.savefig('../steps/spacers/means/R_mean_freq_log_newfig.pdf')
plt.savefig('/home/guillemet/Documents/crispr/final_figures/S9_2.pdf', bbox_inches='tight')

means_control=control.groupby('spacer').freq.mean().sort_values()
df['mean_variant']=df.loc[:,['w1','w2','w3','w4']].apply(np.mean, axis=1)
df.index=df.spacer
df['mean_freq']=means_control[df.index]
plt.figure()
sns.scatterplot(data=df, x='mr', y='mean_variant', color=colors[3])
pear=pearsonr(x=df.mr, y=df.mean_variant)
df=df.loc[df.index!=31725]
reg0=LinearRegression().fit(np.log(np.array(df.mr)).reshape((-1,1)), np.array(df.mean_variant))
plt.plot([1e-8,2e-6],[reg0.intercept_+math.log(1e-8)*reg0.coef_,reg0.intercept_+math.log(2e-6)*reg0.coef_], 'k:')
ax=plt.gca()
ax.set(xscale='log')
sns.despine()
plt.xlim(0.8e-8,1.5e-6)
plt.ylim(-0.025,1.05)
plt.xlabel('Protospacer mutation rate', fontsize=12)
plt.ylabel('Mean frequency of phage mutation', fontsize=12)
#plt.xticks(ticks=[0.01,0.1,1], labels=[0.01,0.1,1])
#ax.set_xticks([0.01,0.1,1])
plt.tick_params(axis='x', bottom=True)
plt.tick_params(left=True)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
#plt.text(1e-7,0.8,'Pearson\'s r='+str(round(pear[0],2)), size=12)
#plt.text(1e-7,0.69,'p-value='+str(round(pear[1], 2)), size=12)
plt.savefig('../steps/spacers/means/R_mr_freq_newfig.pdf')
plt.savefig('/home/guillemet/Documents/crispr/final_figures/S9_4.pdf', bbox_inches='tight')

from scipy.stats import spearmanr
from scipy.stats import pearsonr
print(spearmanr(a=df.mr, b=df.mean_variant))
print(spearmanr(a=df.mean_freq, b=df.mean_variant))
print(pearsonr(x=np.log(df.mr), y=df.mean_variant))
print(pearsonr(x=np.log(df.mean_freq), y=df.mean_variant))






        
        
bacter_all=pd.read_csv('../data/nBacteria_genos_filled.csv', sep=',', header=0) 
bacter_all=bacter_all.iloc[:,1:]      
        
bacter_all=bacter_all.drop_duplicates()
newdf=pd.DataFrame(columns=['Genotype','freq','Rep','Time'])
bims=['0','971', '16236', '31065', '7037', '21039', '34608', '29998', '23084', '25461', '24343', '31149', '3233', '30386', '27013', '1209', '31725']
for rep in set(bacter_all.Rep):
    for time in set(bacter_all.Time):
        for bim in bims:
            subdf=bacter_all.loc[(bacter_all.Time==time) & (bacter_all.Rep==rep)]
            subdf.Genotype=subdf.Genotype.apply(lambda x: subf(x))
            newf=subdf.loc[subdf.Genotype==bim].freq.sum()
            newdf.loc[len(newdf)]=[bim,newf,rep,time]
newdf.columns=['spacer','f_spacer','rep','time']
