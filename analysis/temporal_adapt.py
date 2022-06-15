import pandas as pd
import os
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import math
import copy
from matplotlib.ticker import MaxNLocator
os.chdir('/home/guillemet/Documents/crispr/work/Martin/scripts')

colors= sns.color_palette()
sns.set_style("white")

nbacter_all=pd.read_csv('../data/nBacteria_genos_filled.csv', sep=',', header=0) 
nbacter_all=nbacter_all.iloc[:,1:]

phage_bacterR=pd.read_csv("../steps/df/phage_bacterR.csv")
phage_bacterW=pd.read_csv("../steps/df/phage_bacterW.csv")
phage_bacterW.loc[phage_bacterW.f_variant>1, 'f_variant']=1
phage_bacterC=pd.read_csv("../steps/df/phage_bacterC.csv")

phage_bacter=pd.concat([phage_bacterR, phage_bacterW])

                
def tau_adaptation(r, df, phag):
    df=df[df.freq>0]
    df=df[df.Rep==rep]
    adapts=[]
    for tp in range(5):
        for tb in range(5):
            pha=phag[phag.time==tp]
            pha=pha[pha.rep==r]
            adapt=0
            tdf=df[df.Time==tb]
            for gen in set(tdf.Genotype):
                g=gen.split('-')[1:-1]
                if len(g) == 0:
                    adapt += tdf.loc[tdf.Genotype==gen, 'freq'].iloc[0]
                elif len(g) == 1:
                    sp=int(g[0].split('_')[1])
                    adapt+= tdf.loc[tdf.Genotype==gen, 'freq'].iloc[0] * pha[pha.spacer==sp].f_variant.iloc[0]
                else:
                    toadd=tdf.loc[tdf.Genotype==gen, 'freq'].iloc[0]
                    for sp in g:
                        sp=int(sp.split('_')[1])
                        toadd=toadd*pha[pha.spacer==sp].f_variant.iloc[0]    ##no linkage when on
                    sp=int(g[-1].split('_')[1])
                    #toadd=toadd*pha[pha.spacer==sp].f_variant.iloc[0]     ##full linkage when on
                    adapt+=toadd
            adapts+=[adapt]
    df_out=pd.DataFrame({'cond':25*[r[0]], 'rep':25*[r], 'tbact':[0,1,2,3,4]*5, 'tphage':[0]*5+[1]*5+[2]*5+[3]*5+[4]*5, 'value':adapts})
    return(df_out)


def tau_adaptation_no3(r, df, phag):
    df=df[df.freq>0]
    df=df[df.Rep==rep]
    adapts=[]
    for tp in range(5):
        for tb in range(5):
            pha=phag[phag.time==tp]
            pha=pha[pha.rep==r]
            adapt=0
            tdf=df[df.Time==tb]
            for gen in set(tdf.Genotype):
                g=gen.split('-')[1:-1]
                if len(g) == 0:
                    adapt += tdf.loc[tdf.Genotype==gen, 'freq'].iloc[0]
                elif len(g) == 1:
                    sp=int(g[0].split('_')[1])
                    if pha[pha.spacer==sp].f_variant.iloc[0]>0.5 and tdf.loc[tdf.Genotype==gen, 'freq'].iloc[0]>0.1 and tb>=3:
                        pass
                    else:
                        adapt+= tdf.loc[tdf.Genotype==gen, 'freq'].iloc[0] * pha[pha.spacer==sp].f_variant.iloc[0]
                else:
                    toadd=tdf.loc[tdf.Genotype==gen, 'freq'].iloc[0]
                    for sp in g:
                        sp=int(sp.split('_')[1]) ##no linkage when on
                        if sp==g[-1] and pha[pha.spacer==sp].f_variant.iloc[0]>0.5 and tdf.loc[tdf.Genotype==gen, 'freq'].iloc[0]>0.1 and tb>=3:
                            toadd=0   
                        else:   
                            toadd=toadd*pha[pha.spacer==sp].f_variant.iloc[0]
#                    sp=int(g[-1].split('_')[1]) ## full linkage when on
#                    if pha[pha.spacer==sp].f_variant.iloc[0]>0.5 and tdf.loc[tdf.Genotype==gen, 'freq'].iloc[0]>0.1 and tb>=3:
#                        toadd=0
#                    else:
#                        toadd=toadd*pha[pha.spacer==sp].f_variant.iloc[0]
                    adapt+=toadd
            adapts+=[adapt]
    df_out=pd.DataFrame({'cond':25*[r[0]], 'rep':25*[r], 'tbact':[0,1,2,3,4]*5, 'tphage':[0]*5+[1]*5+[2]*5+[3]*5+[4]*5, 'value':adapts})
    return(df_out)
    
all_df=pd.DataFrame()
for rep in ['W1','W2','W3','W4','W5','W6','W7','W8']:
    #all_df=pd.concat([all_df,tau_adaptation_no3(rep,nbacter_all, phage_bacterW)])
    all_df=pd.concat([all_df,tau_adaptation(rep,nbacter_all, phage_bacterW)])
for rep in ['R1','R2','R3','R4','R5','R6','R7','R8']:
    #all_df=pd.concat([all_df,tau_adaptation_no3(rep,nbacter_all, phage_bacterR)])
    all_df=pd.concat([all_df,tau_adaptation(rep,nbacter_all, phage_bacterR)])
    
all_df_meanr=pd.DataFrame({'tbact':[0,1,2,3,4]*5, 'tphage':[0]*5+[1]*5+[2]*5+[3]*5+[4]*5, 'mvalue':25*[-1]})
all_df_meanw=pd.DataFrame({'tbact':[0,1,2,3,4]*5, 'tphage':[0]*5+[1]*5+[2]*5+[3]*5+[4]*5, 'mvalue':25*[-1]})

for p in range(5):
    for b in range(5):
        all_df_meanr.loc[(all_df_meanr.tbact==b) & (all_df_meanr.tphage==p), 'mvalue']=all_df.loc[(all_df.cond=='R') & (all_df.tbact==b) & (all_df.tphage==p), 'value'].mean()
for p in range(5):
    for b in range(5):
        all_df_meanw.loc[(all_df_meanw.tbact==b) & (all_df_meanw.tphage==p), 'mvalue']=all_df.loc[(all_df.cond=='W') & (all_df.tbact==b) & (all_df.tphage==p), 'value'].mean()

mat_meanr=pd.DataFrame(index=range(5), columns=range(5))
mat_meanw=pd.DataFrame(index=range(5), columns=range(5))
for i in range(5):
    for j in range(5):
        mat_meanr.loc[j,i]=all_df_meanr.loc[(all_df_meanr.tphage==i) & (all_df_meanr.tbact==j), 'mvalue'].iloc[0]
        mat_meanw.loc[j,i]=all_df_meanw.loc[(all_df_meanw.tphage==i) & (all_df_meanw.tbact==j), 'mvalue'].iloc[0]

mat_meanr=mat_meanr.astype(float)

mat_meanw=mat_meanw.astype(float)

plt.figure(figsize=(6.3,5))
sns.heatmap(mat_meanr, vmin=0, vmax=0.65)
ax=plt.gca()
ax.invert_yaxis()
plt.ylabel('Bacteria time (Days)', fontsize=12)
plt.xlabel('Phages time (Days)', fontsize=12)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.savefig('../steps/Adaptation/temporal/heatmap_ppt_D_no_linkage.pdf', bbox_inches='tight')
plt.savefig('/home/guillemet/Documents/crispr/final_figures/S8_2.pdf', bbox_inches='tight')
plt.figure(figsize=(6.3,5))
sns.heatmap(mat_meanw, vmin=0, vmax=0.65)
ax=plt.gca()
ax.invert_yaxis()
plt.ylabel('Bacteria time (Days)', fontsize=12)
plt.xlabel('Phages time (Days)', fontsize=12)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.savefig('../steps/Adaptation/temporal/heatmap_ppt_M_no_linkage.pdf', bbox_inches='tight')
plt.savefig('/home/guillemet/Documents/crispr/final_figures/S8_1.pdf', bbox_inches='tight')

all_df=pd.DataFrame()
for rep in ['W1','W2','W3','W4','W5','W6','W7','W8']:
    all_df=pd.concat([all_df,tau_adaptation_no3(rep,nbacter_all, phage_bacterW)])
    #all_df=pd.concat([all_df,tau_adaptation(rep,nbacter_all, phage_bacterW)])
for rep in ['R1','R2','R3','R4','R5','R6','R7','R8']:
    all_df=pd.concat([all_df,tau_adaptation_no3(rep,nbacter_all, phage_bacterR)])
    #all_df=pd.concat([all_df,tau_adaptation(rep,nbacter_all, phage_bacterR)])
    
all_df_meanr=pd.DataFrame({'tbact':[0,1,2,3,4]*5, 'tphage':[0]*5+[1]*5+[2]*5+[3]*5+[4]*5, 'mvalue':25*[-1]})
all_df_meanw=pd.DataFrame({'tbact':[0,1,2,3,4]*5, 'tphage':[0]*5+[1]*5+[2]*5+[3]*5+[4]*5, 'mvalue':25*[-1]})

for p in range(5):
    for b in range(5):
        all_df_meanr.loc[(all_df_meanr.tbact==b) & (all_df_meanr.tphage==p), 'mvalue']=all_df.loc[(all_df.cond=='R') & (all_df.tbact==b) & (all_df.tphage==p), 'value'].mean()
for p in range(5):
    for b in range(5):
        all_df_meanw.loc[(all_df_meanw.tbact==b) & (all_df_meanw.tphage==p), 'mvalue']=all_df.loc[(all_df.cond=='W') & (all_df.tbact==b) & (all_df.tphage==p), 'value'].mean()

mat_meanr=pd.DataFrame(index=range(5), columns=range(5))
mat_meanw=pd.DataFrame(index=range(5), columns=range(5))
for i in range(5):
    for j in range(5):
        mat_meanr.loc[j,i]=all_df_meanr.loc[(all_df_meanr.tphage==i) & (all_df_meanr.tbact==j), 'mvalue'].iloc[0]
        mat_meanw.loc[j,i]=all_df_meanw.loc[(all_df_meanw.tphage==i) & (all_df_meanw.tbact==j), 'mvalue'].iloc[0]

mat_meanr=mat_meanr.astype(float)

mat_meanw=mat_meanw.astype(float)

plt.figure(figsize=(6.3,5))
sns.heatmap(mat_meanr, vmin=0, vmax=0.65)
ax=plt.gca()
ax.invert_yaxis()
plt.ylabel('Bacteria time (Days)', fontsize=12)
plt.xlabel('Phages time (Days)', fontsize=12)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.savefig('../steps/Adaptation/temporal/heatmap_nocr3_ppt_D_no_linkage.pdf', bbox_inches='tight')
plt.savefig('/home/guillemet/Documents/crispr/final_figures/S8_4.pdf', bbox_inches='tight')

plt.figure(figsize=(6.3,5))
sns.heatmap(mat_meanw, vmin=0, vmax=0.65)
ax=plt.gca()
ax.invert_yaxis()
plt.ylabel('Bacteria time (Days)', fontsize=12)
plt.xlabel('Phages time (Days)', fontsize=12)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.savefig('../steps/Adaptation/temporal/heatmap_nocr3_ppt_M_no_linkage.pdf', bbox_inches='tight')
plt.savefig('/home/guillemet/Documents/crispr/final_figures/S8_3.pdf', bbox_inches='tight')











all_df['tau']=all_df.tbact - all_df.tphage
all_df_w = all_df.loc[all_df.rep.str.contains('W')]
all_df_r = all_df.loc[all_df.rep.str.contains('R')]

#tau plot with both treatments

plt.figure(figsize=[5.4,5], dpi=300)
sns.set_style('ticks')
sns.lineplot(data=all_df_w, x='tau', y='value', color=colors[1])
plt.ylabel('Phage temporal adaptation', fontsize=12)
plt.xlabel('Time delay between bacteria and phage (days)', fontsize=12)
plt.ylim(0,0.55)
sns.lineplot(data=all_df_r, x='tau', y='value', color=colors[3])
from matplotlib.collections import PolyCollection
ax=plt.gca()
for child in ax.findobj(PolyCollection):
    child.set_alpha(0.15)

sns.despine()
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.savefig('../steps/Adaptation/temporal/aadaptation_papier_nolinkage.pdf')
plt.savefig('/home/guillemet/Documents/crispr/final_figures/6_2.pdf', bbox_inches='tight')



plt.figure(figsize=[5.4,5], dpi=300)
plt.plot([1,6], [1,2])
plt.ylabel('Phage fitness', fontsize=12)
plt.xlabel('Time delay between bacteria and phage (days)', fontsize=12)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
sns.despine()
plt.savefig('/home/guillemet/Documents/crispr/final_figures/aaaa6_2.pdf', bbox_inches='tight')



