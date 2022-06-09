#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 11 10:56:47 2021

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


#Make bacter_all

bacter_all=pd.read_csv('../data/Bacteria_genos.csv')

def rep_change(x):
    if 'B' in x:
        return('C'+x[-1])
    elif 'PR' in x:
        return(x[1:])
    elif 'Pw' in x:
        return('W'+x[-1])

def fill_bacter_all(bacter_all):
    '''
    Parameters
    ----------
    bacter_all_file : bacteria genotypes file

    Returns
    -------
    The bacteria genotypes file with containing lines with frequencies of 0 when the genotype is not found
    at this time but found at other times in the replicate
    '''
    bacter_all=bacter_all.drop(columns=['Cond'])
    for i in bacter_all.index:
        sp=bacter_all.loc[i,'Spacer']
        geno=bacter_all.loc[i,'Genotype']
        for rep in set(bacter_all.Rep):
            for time in range(5):
                if len(bacter_all[(bacter_all['Spacer']==sp) & (bacter_all['Genotype']==geno) & (bacter_all['Rep']==rep) & (bacter_all['Time']==time)])==0:
                    bacter_all.loc[len(bacter_all)]=[geno,0,rep,time,sp]
    bacter_all=bacter_all.sort_values(['Rep','Spacer','Time'])
    return(bacter_all)

bacter_all=fill_bacter_all(bacter_all)
bacter_all.to_csv('../data/Bacteria_genos_filled.csv')


import copy
nbacter_all=copy.deepcopy(bacter_all)
nbacter_all=nbacter_all.drop(columns=['Spacer'])
nbacter_all=nbacter_all.drop_duplicates()
for gen in set(nbacter_all.Genotype):
    for time in range(4):
        for rep in set(nbacter_all.Rep):
            nbacter_all.loc[(nbacter_all.Genotype==gen) & (nbacter_all.Time==time) & (nbacter_all.Rep==rep), 'freq_growth']= \
                nbacter_all.loc[(nbacter_all.Genotype==gen) & (nbacter_all.Time==time+1) & (nbacter_all.Rep==rep), 'freq'].iloc[0] - \
                nbacter_all.loc[(nbacter_all.Genotype==gen) & (nbacter_all.Time==time) & (nbacter_all.Rep==rep), 'freq'].iloc[0]
nbacter_all.to_csv('../data/nBacteria_genos_filled.csv')




#make phage_bacter




def make_spacer_variant_table(rep,bacter_all_file='../data/Bacteria_genos_filled.csv'):
    '''

    Parameters
    ----------
    phage_file : R1 to R8.csv / W1 to W8.csv
    bacter_all_file : bacteria genotypes file
    Rep : which replicate (should match the name of phage_file)

    Returns
    -------
    A table with: time
        spacer
        frequency of the spacer
        number of genotype with this spacer
        frequency of the phage variant bypassing this spacer
        number of different mutations bypassing this spacer

    '''
    phage_file='../data/phages/'+rep+'.tsv'
    if 'C' not in rep:
        phage=pd.read_csv(phage_file, sep='\t', header=0)
    bacter_all=pd.read_csv(bacter_all_file, sep=',',header=0)
    if 'C' in rep:
        bacter_all.Rep=bacter_all.Rep.apply(rep_change)
    bacter_all=bacter_all[bacter_all['Rep']==rep]
    bacter_all=bacter_all.reset_index(drop=True)
    #bacter_all=fill_bacter_all(bacter_all) 
    phage_bacter=pd.DataFrame(columns=['rep','time','spacer','f_spacer','n_spacer','f_variant','n_variant',])
    #phage_bacter.loc[len(phage_bacter)]=[1,2,3,4,5]
    p=19
    # Associate variant to bacteria genotype

    for time in [0,1,2,3,4]:
        bacter=bacter_all[bacter_all['Time']==time]
        for spacer in bacter['Spacer']:
            f=0
            n=0
            if 'C' not in rep:
                for proto_index in phage.index:
                    if spacer>phage.loc[proto_index, 'site']-p and spacer<phage.loc[proto_index, 'site']+p and phage.loc[proto_index,'f_'+str(time)]!=0:
                        n+=1
                        f+=phage.loc[proto_index, 'f_'+str(time)]
            #list_f=bacter.loc[(bacter['Spacer']==spacer) & (bacter['Rep']==rep) & (bacter['Time']==time) & (bacter['freq']!=0)]['freq']
            list_f=bacter.loc[(bacter['Spacer']==spacer) & (bacter['Rep']==rep) & (bacter['Time']==time)]['freq']
            if spacer not in phage_bacter[phage_bacter.time==time]['spacer']:
                phage_bacter.loc[len(phage_bacter)]=[rep,time,spacer,sum(list_f),len(list_f),f,n]
    phage_bacter=phage_bacter.sort_values(['spacer','time'])
    phage_bacter=phage_bacter.drop_duplicates()
    phage_bacter=phage_bacter.reset_index(drop=True)
        
    
    #growth computation
    phage_bacter['spacer_growth']=0
    phage_bacter['variant_growth']=0
    for i in phage_bacter.index:
        if phage_bacter.time[i]!=4:
            #print(phage_bacter.loc[i,])
            #print(phage_bacter.loc[i+1,])
            
            phage_bacter.loc[i,'spacer_growth']=phage_bacter.loc[i+1, 'f_spacer']-phage_bacter.loc[i, 'f_spacer']
            phage_bacter.loc[i,'variant_growth']=phage_bacter.loc[i+1, 'f_variant']-phage_bacter.loc[i, 'f_variant']
    # 0 variants for DGCC
    for row in phage_bacter.index:
        if phage_bacter.loc[row,'spacer']==0:
            phage_bacter.loc[row,'n_spacer']=0
            phage_bacter.loc[row,'n_variant']=0
            phage_bacter.loc[row,'f_variant']=0
            
    #Normalise for sum freq = 1
    return(phage_bacter)



phage_bacterW=pd.concat([make_spacer_variant_table('W1'),
                        make_spacer_variant_table('W2'),
                        make_spacer_variant_table('W3'),
                        make_spacer_variant_table('W4'),
                        make_spacer_variant_table('W5'),
                        make_spacer_variant_table('W6'),
                        make_spacer_variant_table('W7'),
                        make_spacer_variant_table('W8')], ignore_index=True)

phage_bacterR=pd.concat([make_spacer_variant_table('R1'),
                        make_spacer_variant_table('R2'),
                        make_spacer_variant_table('R3'),
                        make_spacer_variant_table('R4'),
                        make_spacer_variant_table('R5'),
                        make_spacer_variant_table('R6'),
                        make_spacer_variant_table('R7'),
                        make_spacer_variant_table('R8')], ignore_index=True)



phage_bacterC=pd.read_csv('../data/Bacteria_genos.csv')
phage_bacterC=phage_bacterC[phage_bacterC.Cond=='C']
phage_bacterC=fill_bacter_all(phage_bacterC)
phage_bacterC.columns=['genotype', 'freq','rep','time','spacer']
phage_bacterC_mean=pd.DataFrame(columns=['genotype', 'freq','time','spacer'])
for geno in set(phage_bacterC.genotype):
    for t in set(phage_bacterC.time):
        meanf=phage_bacterC[(phage_bacterC.genotype==geno) & (phage_bacterC.time==t)].freq.mean()
        phage_bacterC_mean.loc[len(phage_bacterC_mean)]=[geno,meanf,t,phage_bacterC[phage_bacterC.genotype==geno].spacer.iloc[0]]
    
phage_bacter=pd.concat([phage_bacterR, phage_bacterW])
