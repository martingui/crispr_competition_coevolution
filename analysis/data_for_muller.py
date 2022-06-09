import pandas as pd
import os
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import math
import copy
os.chdir('/home/guillemet/Documents/crispr/work/Martin/scripts')

bacter_all=pd.read_csv('../data/Bacteria_genos_filled.csv', sep=',', header=0, index_col=0)



def data_for_muller():
    for rep in ['R1','R2','R3','R4','R5','R6','R7','R8','W1','W2','W3','W4','W5','W6','W7','W8','C1','C2','C3','C4','C5','C6','C7']:
        pops=pd.read_csv("../data/Bacteria_pop.csv")
        df=bacter_all[bacter_all.Rep==rep]
        df=df[df.freq>0]
        df=df.drop(columns=['Spacer','Rep'])
        df=df.drop_duplicates()
        df=df.reset_index(drop=True)
        for i in df.index:
            a=df.loc[i,'Genotype'].split('-')
            if len(a)==2:
                df.loc[i,'Genotype']='DGCC'
            if len(a)==3:
                df.loc[i,'Genotype']=str(a[1].split('_')[1])
            if len(a)==4:
                df.loc[i,'Genotype']=   str(a[1].split('_')[1])+"-"+str(a[2].split('_')[1])
            if len(a)==5:
                df.loc[i,'Genotype']=str(a[1].split('_')[1])+"-"+str(a[2].split('_')[1])+"-"+str(a[3].split('_')[1])
        
        
        for gen in set(df.Genotype):
            for t in range(5):
                if gen not in list(df[df.Time==t].Genotype):
                    df.loc[len(df)]=[gen,0,t]  
        for gen in set(df.Genotype):
                for t in range(5):
                    if t!=0:
                        if df[(df.Genotype==gen)&(df.Time==t)].freq.iloc[0]!=0 and df[(df.Genotype==gen)&(df.Time==t-1)].freq.iloc[0] == 0:
                            i=df[(df.Genotype==gen)&(df.Time==t-1)].index[0]
                    if t!=4:
                        if df[(df.Genotype==gen)&(df.Time==t)].freq.iloc[0]!=0 and df[(df.Genotype==gen)&(df.Time==t+1)].freq.iloc[0] == 0:
                            i=df[(df.Genotype==gen)&(df.Time==t+1)].index[0]
                            df.loc[i,'freq'] = 0
                            if t==3:
                                df.loc[i,'freq'] = 1e-8
                            
        
        ### ADD TOTAL POP AND MULTIPLY
        df['pop']=0
        for i in df.index:
            df.loc[i,'pop']=math.log10(pops.loc[(pops.time==df.Time[i]) & (pops.rep==rep)]['count'].iloc[0])
        df.columns=["Identity","Population","Generation","tot"]
        df.Population=df.Population*df.tot
        df=df.drop(columns='tot')
        df=df.sort_values('Generation')
        
        ### INTERPOLATION
        ndf=pd.DataFrame(columns=df.columns)
        for gen in set(df.Identity):
            interp_f = scipy.interpolate.PchipInterpolator(df[df.Identity==gen].Generation, df[df.Identity==gen].Population)
            nx=np.linspace(0,4,101)
            for t in nx:
                ndf.loc[len(ndf)]=[gen, interp_f(t), t]
        ndf=ndf.sort_values(["Identity","Generation"])
        ndf.to_csv("../steps/muller/csv/apop_"+rep+".csv")
        edges=pd.DataFrame(columns=["Parent","Identity"])
        for i in set(df.Identity):
            for j in set(df.Identity):
                if i!=j:
                    cond=True
                    if len(i.split('-')) == len(j.split('-'))+1:
                        maxx=i
                        minn=j
                    elif len(j.split('-')) == len(i.split('-'))+1:
                        maxx=j
                        minn=i
                    else:
                        cond=False
                    if cond:
                        if maxx.split('-')[:-1] == minn.split('-'):
                            edges.loc[len(edges)]=[minn,maxx]
    
            if len(i.split('-'))==1 and i!='DGCC':
                edges.loc[len(edges)]=['DGCC',i]
        for i in set(df.Identity):
            if i!='DGCC' and i not in list(edges.Identity):
                parent=0
                for spa in range(1,len(i.split('-'))):
                    if '-'.join(i.split('-')[:spa]) in set(df.Identity):
                        parent='-'.join(i.split('-')[:spa])
                if parent == 0:
                    parent='DGCC'
                edges.loc[len(edges)]=[parent,i]
                        
        edges=edges.drop_duplicates()
        edges.to_csv("../steps/muller/csv/edges_"+rep+".csv")
        
data_for_muller()
                        