import pandas as pd
import os
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import math
import copy
from matplotlib.ticker import MaxNLocator
from sklearn.linear_model import LinearRegression

colors= sns.color_palette()


bacter_all=pd.read_csv('../data/nBacteria_genos_filled.csv', sep=',', header=0) 
bacter_all=bacter_all.iloc[:,1:]

phage_bacterR=pd.read_csv("../steps/df/phage_bacterR.csv")
phage_bacterW=pd.read_csv("../steps/df/phage_bacterW.csv")
phage_bacterW.loc[phage_bacterW.f_variant>1, 'f_variant']=1
phage_bacterC=pd.read_csv("../steps/df/phage_bacterC.csv")

phage_bacter=pd.concat([phage_bacterR, phage_bacterW])

sns.set_style("white")

def fst_phages_jk(jacknife=False, toremove=-1, Qst=False):
    import copy
    if toremove>0:
        r=7
    elif toremove==-1:
        r=8
    phages_read=pd.DataFrame()
    listrep=list(range(1,9))
    listrep.remove(toremove)
    listrep=['W'+str(x) for x in listrep]+['R'+str(x) for x in listrep]
    for i in listrep:
        temp=pd.read_csv('../data/FreeBayes/'+i[0]+"_seq/"+i+"_data.csv", delimiter="\t")
        temp['rep']=i
        phages_read=pd.concat([phages_read, temp])
    phages_read.TIME=phages_read.TIME.apply(lambda x: x if 'R' not in x else 'T'+x[1])
    
    p=19
    print(r)
    for cond in ['W','R']:
        diff_df=pd.DataFrame(columns = ['time','locus','allele','nbar','nc','pbar','s2','a','b'])
        baseRdf=phages_read[phages_read.rep.str.contains(cond)]
        baseRdf=baseRdf.reset_index(drop=True)
        Rdf=baseRdf
        
        if Qst==True:
            all_Rdf=pd.DataFrame()
            for repn in range(1,9):
                for time in range(1,5):
                    rep=cond+str(repn)
                    Rdf=baseRdf.loc[(baseRdf.rep==rep) & (baseRdf.TIME=='T'+str(time))]
                    dict_freq={}
                    dict_dp={}
                    dict_count={}
                    gens=[]
                    tokeep=[]
                    for i in range(len(Rdf.index)):
                        pos=Rdf.iloc[i]['POS']
                        for gen in set(phage_bacter[phage_bacter.rep.str.contains(rep)].spacer):
                            if abs(pos-gen)<p:
                                tokeep+=[pos]
                                gens+=[gen]
                                if gen not in list(dict_freq.keys()):
                                    dict_freq[gen]=Rdf.iloc[i]['FREQ']
                                    dict_dp[gen]=Rdf.iloc[i]['DP']
                                    dict_count[gen]=1
                                else:
                                    dict_freq[gen]+=Rdf.iloc[i]['FREQ']
                                    dict_dp[gen]+=Rdf.iloc[i]['DP']
                                    dict_count[gen]+=1
                    nRdf=pd.DataFrame(columns=['unnamed','AO','DP','TYPE','TIME','ALT','REF','POS','FREQ','rep'])                
                    for gen in dict_freq.keys():
                        n=dict_count[gen]
                        nRdf.loc[len(nRdf)]=['null',round(((dict_freq[gen]/n)*dict_dp[gen])/n),dict_dp[gen]/n,'mut', 'T'+str(time),'mut','wt', gen, dict_freq[gen]/n, rep]
                    all_Rdf=pd.concat([all_Rdf, nRdf])
                    Rdf=all_Rdf
                
        #### add the locus for every rep with frequency 0 if not already there, add a DP=mean DP where this locus is found
        allpos=set(Rdf.POS)
        dicallpos={}
        for pos in allpos:
                dicallpos[pos]=Rdf.loc[Rdf.POS==pos, 'DP'].mean()
        for time in range(5):
            pdf=Rdf[(Rdf.rep.str.contains(cond)) & (Rdf.TIME=='T'+str(time))]
            pdf=pdf.reset_index(drop=True)
            reprange=list(range(1,9))
            reprange.remove(toremove)
            for rep in reprange:
                for pos in allpos:
                    if pos not in list(pdf[pdf.rep==cond+str(rep)].POS):
                        pdf.loc[len(pdf)]=[pos, dicallpos[pos], dicallpos[pos], 'wt', 'T'+str(time), 'wt', 'wt', pos, 1, cond+str(rep)]
            #Add the wt allele for each mutation
            for site in set(pdf.POS):
                site_df=pdf[pdf.POS==site]
                site_df=site_df.reset_index(drop=True)
                for rep in set(site_df.rep):
                    repsite_df=site_df[site_df.rep==rep]
                    if repsite_df.TYPE.iloc[0]=='wt':
                        pass
                    else:
                        sumAO=repsite_df.AO.sum()
                        sumfreq=repsite_df.FREQ.sum()
                        site_df.loc[len(site_df)]=repsite_df.iloc[0]
                        site_df.loc[len(site_df)-1,'AO'] = site_df.loc[len(site_df)-1].DP - sumAO
                        site_df.loc[len(site_df)-1,'FREQ'] = 1 - sumfreq
                        site_df.loc[len(site_df)-1,'ALT'] = 'wt'
                        site_df.loc[len(site_df)-1,'TYPE'] = 'wt'
                
                #
                nbar=site_df[site_df.TYPE=='wt'].DP.mean()
                nc=(r*nbar - ( sum([z**2 for z in site_df[site_df.TYPE=='wt'].DP]) / (r*nbar) ) ) / (r-1)
                for allele in set(site_df.ALT):
                    for rep in set(site_df.rep):
                        if len(site_df[(site_df.rep==rep) & (site_df.ALT==allele)])<1:
                            site_df.loc[len(site_df)]=[site,0,site_df.loc[site_df.ALT=='wt','DP'].mean(), 'not', 'T'+str(time), allele,'wt',pos,0,rep]
                    pbar=site_df[site_df.ALT==allele].FREQ.sum()/r
                    s2= sum(site_df[site_df.ALT==allele].DP * (site_df[site_df.ALT==allele].FREQ - pbar)**2)  / ((r-1) * nbar)
                    diff_df.loc[len(diff_df)]=[time,site, allele, nbar,nc,pbar,s2,-1,-1]
                    
        diff_df['a'] = (diff_df.nbar / diff_df.nc) * \
            ( diff_df.s2 - (1/(diff_df.nbar - 1)) * \
             (diff_df.pbar * (1 - diff_df.pbar) - ((r-1)/r) * diff_df.s2) )
                
        diff_df['b'] = ( diff_df.nbar / (diff_df.nbar -1)) * \
            (diff_df.pbar * (1 - diff_df.pbar) - ((r-1)/r) * diff_df.s2)
          
        if cond == 'R':
            diff_dfr=copy.deepcopy(diff_df)
        elif cond=='W':
            diff_dfw=copy.deepcopy(diff_df)     
    return(diff_dfw, diff_dfr)

df_fst_phage_all=pd.DataFrame()
for rem in range(1,9):
    diff_dfw, diff_dfr=fst_phages_jk(jacknife=True, toremove=rem, Qst=True)    
    fst_phager=[]
    for time in range(5):
        fst_phager+=[diff_dfr[diff_dfr.time==time].a.sum() / (diff_dfr[diff_dfr.time==time].a.sum() + diff_dfr[diff_dfr.time==time].b.sum())]
    fst_phager[0]=0    
    
    fst_phagew=[]
    for time in range(5):
        fst_phagew+=[diff_dfw[diff_dfw.time==time].a.sum() / (diff_dfw[diff_dfw.time==time].a.sum() + diff_dfw[diff_dfw.time==time].b.sum())]
    fst_phagew[0]=0     
    
    df_fst_phage=pd.DataFrame({'fst':fst_phagew+fst_phager, 'time':list(range(5))+list(range(5)), 'cond':['W']*5+['R']*5})    
    df_fst_phage_all=pd.concat([df_fst_phage_all,df_fst_phage])
    
    
#df_fst_phage_all.to_csv('../steps/Fst/jk_fst_phage_rw.csv')


df_fst_phage_all=pd.read_csv('../steps/Fst/jk_fst_phage_rw.csv')
plt.figure()
sns.set_style('ticks')
sns.lineplot(data=df_fst_phage_all, x='time', y='fst', hue='cond',palette=[colors[0],colors[1],colors[3]], hue_order=['C','W','R'])
ax=plt.gca()
ax.xaxis.set_major_locator(MaxNLocator(integer=True))
ax.get_legend().remove()
plt.xlabel('Time (days)',fontsize=12)
plt.ylabel('Fst',fontsize=12)
plt.ylim(-0.02,0.6)
sns.despine()
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.savefig('../steps/Fst/jk_fst_phage_ppt.pdf', bbox_inches='tight')
plt.savefig('/home/guillemet/Documents/crispr/final_figures/S6_1.pdf', bbox_inches='tight')


df_fst_phage_all=pd.read_csv('../steps/Fst/jk_qst_phage_rw.csv')
plt.figure()
sns.set_style('ticks')
sns.lineplot(data=df_fst_phage_all, x='time', y='fst', hue='cond',palette=[colors[0],colors[1],colors[3]], hue_order=['C','W','R'])
ax=plt.gca()
ax.xaxis.set_major_locator(MaxNLocator(integer=True))
ax.get_legend().remove()
plt.xlabel('Time (days)',fontsize=12)
plt.ylabel('Qst',fontsize=12)
plt.ylim(-0.02,0.6)
sns.despine()
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.savefig('../steps/Fst/jk_qst_phage_ppt.pdf', bbox_inches='tight')
plt.savefig('/home/guillemet/Documents/crispr/final_figures/S6_2.pdf', bbox_inches='tight')


def fst_bacteria_jk(bact=bacter_all, toremove=-1):
    print('b')
    bacter_df=copy.deepcopy(bact)
    bacter_df=bacter_df.drop_duplicates()
    bacter_df=bacter_df.loc[bacter_df.Rep.str.contains(str(toremove)) == False]
    bacter_df=bacter_df.reset_index(drop=True)
    
    for cond in ['R','W','C']:
        for gen in set(bacter_df[bacter_df.Rep.str.contains(cond)].Genotype):
            for time in range(5):
                for rep in set(bacter_df[bacter_df.Rep.str.contains(cond)].Rep):
                    if len(bacter_df[(bacter_df.Rep==rep) & (bacter_df.Time==time) & (bacter_df.Genotype==gen)])<1:
                        bacter_df.loc[len(bacter_df)]=[gen, 0, rep, time, 0]
    
    bacter_df['DP']=-1
    data=pd.DataFrame()
    
    for time in range(1,5):
        for rep in ['PR1','PR2','PR3','PR4','PR5','PR6','PR7','PR8','Pw1','Pw2','Pw3','Pw4','Pw5','Pw6','Pw7','Pw8','B1','B2','B3','B4','B5','B6','B7']:
            data=pd.read_csv("../data/samples/"+rep+'T'+str(time)+"/inter4_replaced_sort_count",header=None, delimiter='\s+')
            data.columns=['Count','Genotype']
            if 'PR' in rep:
                bacter_df.loc[(bacter_df.Rep=='R'+rep[2]) & (bacter_df.Time==time),'DP']=data.Count.sum()
            elif 'Pw' in rep:
                bacter_df.loc[(bacter_df.Rep=='W'+rep[2]) & (bacter_df.Time==time),'DP']=data.Count.sum()
            elif 'B' in rep:
                bacter_df.loc[(bacter_df.Rep=='C'+rep[1]) & (bacter_df.Time==time),'DP']=data.Count.sum()
    
    
    for cond in ['R','W','C']:
        if cond == 'R' or cond == 'W':
            if toremove>0:
                r=7
            elif toremove==-1:
                r=8
        elif cond == 'C':
            if toremove>0:
                r=6
            elif toremove==-1:
                r=7
        diff_df=pd.DataFrame(columns = ['time','genotype','nbar','nc','pbar','s2','a','b'])
        for time in range(1,5):
            pdf=bacter_df[(bacter_df.Rep.str.contains(cond)) & (bacter_df.Time==time)]
            pdf=pdf.reset_index(drop=True)
            site_df=pdf
            for gen in set(site_df.Genotype):   
                nbar=site_df[site_df.Genotype==gen].DP.mean()
                nc=(r*nbar - ( sum([z**2 for z in site_df[site_df.Genotype==gen].DP]) / (r*nbar) ) ) / (r-1)
                pbar=site_df[site_df.Genotype==gen].freq.sum()/r
                s2= sum(site_df[site_df.Genotype==gen].DP * (site_df[site_df.Genotype==gen].freq - pbar)**2)  / ((r-1) * nbar)
                diff_df.loc[len(diff_df)]=[time,gen, nbar,nc,pbar,s2,-1,-1]
                    
        diff_df['a'] = (diff_df.nbar / diff_df.nc) * \
            ( diff_df.s2 - (1/(diff_df.nbar - 1)) * \
             (diff_df.pbar * (1 - diff_df.pbar) - ((r-1)/r) * diff_df.s2) )
                
        diff_df['b'] = ( diff_df.nbar / (diff_df.nbar -1)) * \
            (diff_df.pbar * (1 - diff_df.pbar) - ((r-1)/r) * diff_df.s2)
          
    
        if cond == 'R':
            diff_dfr=copy.deepcopy(diff_df)
        elif cond=='W':
            diff_dfw=copy.deepcopy(diff_df)
        elif cond=='C':
            diff_dfc=copy.deepcopy(diff_df)
            
    return(diff_dfw, diff_dfr, diff_dfc)

df_fst_bacter_all=pd.DataFrame()
for rem in range(1,9):
    print(rem)
    diff_dfw, diff_dfr, diff_dfc= fst_bacteria_jk(bact=bacter_all, toremove=rem)
    fst_bacterr=[]
    for time in range(5):
        fst_bacterr+=[diff_dfr[diff_dfr.time==time].a.sum() / (diff_dfr[diff_dfr.time==time].a.sum() + diff_dfr[diff_dfr.time==time].b.sum())]
    fst_bacterr[0]=0    
    
    fst_bacterw=[]
    for time in range(5):
        fst_bacterw+=[diff_dfw[diff_dfw.time==time].a.sum() / (diff_dfw[diff_dfw.time==time].a.sum() + diff_dfw[diff_dfw.time==time].b.sum())]
    fst_bacterw[0]=0  
    
    if not rem==8:
        fst_bacterc=[]
        for time in range(5):
            fst_bacterc+=[diff_dfc[diff_dfc.time==time].a.sum() / (diff_dfc[diff_dfc.time==time].a.sum() + diff_dfc[diff_dfc.time==time].b.sum())]
        fst_bacterc[0]=0     
    if not rem == 8:   
        df_fst_bacter=pd.DataFrame({'fst':fst_bacterr+fst_bacterw+fst_bacterc, 'time':list(range(5))+list(range(5))+list(range(5)), 'cond':['R']*5+['W']*5+['C']*5})    
        df_fst_bacter_all=pd.concat([df_fst_bacter_all,df_fst_bacter])
    else:
        df_fst_bacter=pd.DataFrame({'fst':fst_bacterr+fst_bacterw, 'time':list(range(5))+list(range(5)), 'cond':['R']*5+['W']*5})    
        df_fst_bacter_all=pd.concat([df_fst_bacter_all,df_fst_bacter])

#df_fst_bacter_all.to_csv('../steps/Fst/jk_bact_rwc.csv')
df_fst_bacter_all=pd.read_csv('../steps/Fst/jk_bact_rwc.csv')
plt.figure()
sns.set_style('ticks')
sns.lineplot(data=df_fst_bacter_all.loc[df_fst_bacter_all.cond.str.contains('L')==False], x='time', y='fst', hue='cond',palette=[colors[0],colors[1],colors[3]], hue_order=['C','W','R'])
ax=plt.gca()
ax.xaxis.set_major_locator(MaxNLocator(integer=True))
ax.get_legend().remove()
plt.xlabel('Time (days)',fontsize=12)
plt.ylabel('Fst',fontsize=12)
plt.ylim(-0.03,0.75)
sns.despine()
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.savefig('../steps/Fst/ajk_bact_ppt_nod.pdf', bbox_inches='tight')
plt.savefig('/home/guillemet/Documents/crispr/final_figures/S4_1.pdf', bbox_inches='tight')

