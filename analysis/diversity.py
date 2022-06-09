#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 11 10:11:36 2021

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
sns.set_style('ticks')
colors= sns.color_palette()
sns.set_style("white")
os.chdir('/home/guillemet/Documents/crispr/work/Martin/scripts')
bacter_all=pd.read_csv('../data/nBacteria_genos_filled.csv', sep=',', header=0) 
bacter_all=bacter_all.iloc[:,1:]
phage_bacterR=pd.read_csv("../steps/df/phage_bacterR.csv")
phage_bacterW=pd.read_csv("../steps/df/phage_bacterW.csv")
phage_bacterC=pd.read_csv("../steps/df/phage_bacterC.csv")

phage_bacter=pd.concat([phage_bacterR, phage_bacterW])



def diversity(bacter_all=bacter_all):
    #bacter_all=bacter_all.drop(columns=['Spacer'])
    bacter_all=bacter_all.drop_duplicates()
    
    hsreps=[]
    for time in range(5):
        for rep in ['C1','C2','C3','C4','C5','C6','C7']:
            hs=1
            for gen in set(bacter_all[(bacter_all.Rep==rep) & (bacter_all.Time==time)].Genotype):
                hs = hs - (bacter_all[(bacter_all.Rep==rep) & (bacter_all.Genotype==gen) & (bacter_all.Time==time)].freq.iloc[0])**2
            hsreps+=[hs]
    eff_c=[1/(1-x) for x in hsreps]
    print(eff_c)
    
    
    
    hsreps=[]
    for time in range(5):
        for rep in ['R1','R2','R3','R4','R5','R6','R7','R8']:
            hs=1
            for gen in set(bacter_all[(bacter_all.Rep==rep) & (bacter_all.Time==time)].Genotype):
                hs = hs - (bacter_all[(bacter_all.Rep==rep) & (bacter_all.Genotype==gen) & (bacter_all.Time==time)].freq.iloc[0])**2
            hsreps+=[hs]
    eff_r=[1/(1-x) for x in hsreps]
    print(eff_r)
        
    hsreps=[]    
    for time in range(5):
        for rep in ['W1','W2','W3','W4','W5','W6','W7','W8']:
            hs=1
            for gen in set(bacter_all[(bacter_all.Rep==rep) & (bacter_all.Time==time)].Genotype):
                hs = hs - (bacter_all[(bacter_all.Rep==rep) & (bacter_all.Genotype==gen) & (bacter_all.Time==time)].freq.iloc[0])**2
            hsreps+=[hs]
    eff_w=[1/(1-x) for x in hsreps]
    print(eff_w)

        
    print(len(eff_w+eff_r+eff_c))    
    print(len(8*[0]+8*[1]+8*[2]+8*[3]+8*[4]+8*[0]+8*[1]+8*[2]+8*[3]+8*[4]+7*[0]+7*[1]+7*[2]+7*[3]+7*[4]))
    df_fst=pd.DataFrame({'fst':eff_w+eff_r+eff_c, 'time':8*[0]+8*[1]+8*[2]+8*[3]+8*[4]+8*[0]+8*[1]+8*[2]+8*[3]+8*[4]+7*[0]+7*[1]+7*[2]+7*[3]+7*[4], 'cond':['W']*40+['R']*40+['C']*35, 'rep':5*['R1','R2','R3','R4','R5','R6','R7','R8']+5*['W1','W2','W3','W4','W5','W6','W7','W8']+5*['C1','C2','C3','C4','C5','C6','C7']})
    df_fst_plot=copy.deepcopy(df_fst)
    df_fst_plot.cond=df_fst_plot.cond.apply(lambda x : 'M' if x=='W' else ('D' if x=='R' else x))
    print(df_fst_plot)
    plt.figure()
    sns.set_style('ticks')
    sns.lineplot(data=df_fst_plot.loc[df_fst_plot.cond!='X'], x='time', y='fst', hue='cond', estimator=None, units='rep', hue_order=['C','M','D'], palette=[colors[0],colors[1],colors[3]], linewidth=0.8)
    sns.scatterplot(data=df_fst_plot.loc[df_fst_plot.cond!='X'], x='time', y='fst', hue='cond',style='cond' ,markers=['o','s','P'], estimator=None, units='rep', hue_order=['C','M','D'], palette=[colors[0],colors[1],colors[3]])
    plt.ylabel('Bacteria diversity', fontsize=12)
    handles, labels = plt.gca().get_legend_handles_labels()
    labels=['A','B','C','A','B','C']
    plt.legend(handles=handles[3:], labels=labels[3:])
    ax=plt.gca()
    ax.xaxis.set_major_locator(MaxNLocator(integer=True))
    plt.xlabel('Time (days)', fontsize=12)
    plt.ylim(0,18)
    sns.despine()
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.savefig('../steps/Diversity/bacteria_diversity_CMD_ppt.png', dpi=300, bbox_inches='tight')
    plt.savefig('/home/guillemet/Documents/crispr/final_figures/3.png', dpi=300, bbox_inches='tight')

    return(df_fst_plot)

df=diversity()

def subf(x):
    try:
        return(x.split('-')[1].split('_')[1])
    except:
        return('0')

def bim_diversity(bacter_all=bacter_all):
    #bacter_all=bacter_all.drop(columns=['Spacer'])
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
                
    bacter_all=newdf       
    hsreps=[]
    for time in range(5):
        for rep in ['C1','C2','C3','C4','C5','C6','C7']:
            hs=1
            for gen in set(bacter_all[(bacter_all.Rep==rep) & (bacter_all.Time==time)].Genotype):
                hs = hs - (bacter_all[(bacter_all.Rep==rep) & (bacter_all.Genotype==gen) & (bacter_all.Time==time)].freq.iloc[0])**2
            hsreps+=[hs]
    eff_c=[1/(1-x) for x in hsreps]
    print(eff_c)
    
    
    
    hsreps=[]
    for time in range(5):
        for rep in ['R1','R2','R3','R4','R5','R6','R7','R8']:
            hs=1
            for gen in set(bacter_all[(bacter_all.Rep==rep) & (bacter_all.Time==time)].Genotype):
                hs = hs - (bacter_all[(bacter_all.Rep==rep) & (bacter_all.Genotype==gen) & (bacter_all.Time==time)].freq.iloc[0])**2
            hsreps+=[hs]
    eff_r=[1/(1-x) for x in hsreps]
    print(eff_r)
        
    hsreps=[]    
    for time in range(5):
        for rep in ['W1','W2','W3','W4','W5','W6','W7','W8']:
            hs=1
            for gen in set(bacter_all[(bacter_all.Rep==rep) & (bacter_all.Time==time)].Genotype):
                hs = hs - (bacter_all[(bacter_all.Rep==rep) & (bacter_all.Genotype==gen) & (bacter_all.Time==time)].freq.iloc[0])**2
            hsreps+=[hs]
    eff_w=[1/(1-x) for x in hsreps]
    print(eff_w)

        

    df_fst=pd.DataFrame({'fst':eff_w+eff_r+eff_c, 'time':8*[0]+8*[1]+8*[2]+8*[3]+8*[4]+8*[0]+8*[1]+8*[2]+8*[3]+8*[4]+7*[0]+7*[1]+7*[2]+7*[3]+7*[4], 'cond':['W']*40+['R']*40+['C']*35, 'rep':5*['R1','R2','R3','R4','R5','R6','R7','R8']+5*['W1','W2','W3','W4','W5','W6','W7','W8']+5*['C1','C2','C3','C4','C5','C6','C7']})

    df_fst_plot=copy.deepcopy(df_fst)
    df_fst_plot.cond=df_fst_plot.cond.apply(lambda x : 'M' if x=='W' else ('D' if x=='R' else x))

    plt.figure()

    sns.set_style('ticks')
    sns.lineplot(data=df_fst_plot.loc[df_fst_plot.cond!='X'], x='time', y='fst', hue='cond', estimator=None, units='rep', hue_order=['C','M','D'], palette=[colors[0],colors[1],colors[3]], linewidth=0.8)
    sns.scatterplot(data=df_fst_plot.loc[df_fst_plot.cond!='X'], x='time', y='fst', hue='cond',style='cond' ,markers=['o','s','P'], estimator=None, units='rep', hue_order=['C','M','D'], palette=[colors[0],colors[1],colors[3]])

    plt.ylabel('First spacer diversity', fontsize=12)
    handles, labels = plt.gca().get_legend_handles_labels()
    handles, labels = plt.gca().get_legend_handles_labels()
    labels=['A','B','C','A','B','C']
    plt.legend(handles=handles[3:], labels=labels[3:])
    ax=plt.gca()
    plt.ylim(0,18)
    ax.xaxis.set_major_locator(MaxNLocator(integer=True))
    plt.xlabel('Time (days)', fontsize=12)
    sns.despine()
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.savefig('../steps/Diversity/bacteria_diversity_original_CMD_ppt.png', dpi=300, bbox_inches='tight')
    plt.savefig('/home/guillemet/Documents/crispr/final_figures/S3.png', dpi=300, bbox_inches='tight')

    return(df_fst_plot)

bim_diversity()


def ci(p):
    """
    Return 2-sided symmetric confidence interval specified
    by p.
    """
    u_pval = (1+p)/2.
    l_pval = (1-u_pval)
    l_indx = int(np.floor(n*l_pval))
    u_indx = int(np.floor(n*u_pval))
    return(simulations[l_indx],simulations[u_indx])

def bootstrap(data, n=1000, func=np.mean):
    """
    Generate `n` bootstrap samples, evaluating `func`
    at each resampling. `bootstrap` returns a function,
    which can be called to obtain confidence intervals
    of interest.
    """
    simulations = list()
    sample_size = len(data)
    xbar_init = np.mean(data)
    for c in range(n):
        itersample = np.random.choice(data, size=sample_size, replace=True)
        simulations.append(func(itersample))
    simulations.sort()

    return(ci)

for time in set(df.time):
    for cond in set(df.cond):
        boot=bootstrap(df.loc[(df.time==time)&(df.cond==cond),'fst'], n=5000)
        cintervals = boot(.95)
        print(time, cond)
        print((cintervals[1] + cintervals[0])/2,(cintervals[1] - cintervals[0])/2)
        print(cintervals)
        
import statsmodels.api as sm
from statsmodels.formula.api import ols

for t in range(5):
    print(t)
    model = ols('fst ~ cond', data=df[df.time==t]).fit()
    anova_table = sm.stats.anova_lm(model, typ=2)
    print(anova_table)

