#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 13 12:04:24 2021

@author: guillemet
"""

def fst_jk(toremove, bacter_all=bacter_all):
    #bacter_all=bacter_all.drop(columns=['Spacer'])
    bacter_all=bacter_all.drop_duplicates()
    bacter_all=bacter_all.loc[bacter_all.Rep.str.contains(str(toremove))==False]
    hssr=[]
    for time in range(5):
        hsreps=[]
        for rep in [x for x in set(bacter_all.Rep) if 'R' in x]:
            hs=1
            for gen in set(bacter_all[(bacter_all.Rep==rep) & (bacter_all.Time==time)].Genotype):
                hs = hs - (bacter_all[(bacter_all.Rep==rep) & (bacter_all.Genotype==gen) & (bacter_all.Time==time)].freq.iloc[0])**2
            hsreps+=[hs]
        hssr+=[sum(hsreps)/7]
    print('hssr')    
    htsr=[]
    for time in range(5):
        hsgen=[]
        for gen in set(bacter_all[(bacter_all.Rep.str.contains('R')) & (bacter_all.Time==time)].Genotype):
            hs=0
            for rep in [x for x in set(bacter_all.Rep) if 'R' in x]:
                try: 
                    hs = hs + (bacter_all[(bacter_all.Rep==rep) & (bacter_all.Genotype==gen) & (bacter_all.Time==time)].freq.iloc[0])
                except:
                    pass
            hsgen+=[(hs/7)**2]
        htsr+=[1-sum(hsgen)]
    print('htsr')    
    Dsr=[]
    fstsr=[]
    for time in range(5):
        fstsr+=[(7/6)*(htsr[time]-hssr[time]) / (htsr[time])]
        Dsr+=[(7/6)*(htsr[time]-hssr[time]) / (1-hssr[time])]
        
        
    hssw=[]
    for time in range(5):
        hsreps=[]
        for rep in [x for x in set(bacter_all.Rep) if 'W' in x]:
            hs=1
            for gen in set(bacter_all[(bacter_all.Rep==rep) & (bacter_all.Time==time)].Genotype):
                hs = hs - (bacter_all[(bacter_all.Rep==rep) & (bacter_all.Genotype==gen) & (bacter_all.Time==time)].freq.iloc[0])**2
            hsreps+=[hs]
        hssw+=[sum(hsreps)/7]
    print('hssw')    
    htsw=[]
    for time in range(5):
        hsgen=[]
        for gen in set(bacter_all[(bacter_all.Rep.str.contains('R')) & (bacter_all.Time==time)].Genotype):
            hs=0
            for rep in [x for x in set(bacter_all.Rep) if 'W' in x]:
                try: 
                    hs = hs + (bacter_all[(bacter_all.Rep==rep) & (bacter_all.Genotype==gen) & (bacter_all.Time==time)].freq.iloc[0])
                except:
                    pass
            hsgen+=[(hs/7)**2]
        htsw+=[1-sum(hsgen)]
    print('htsw')    
    Dsw=[]
    fstsw=[]
    for time in range(5):
        fstsw+=[(7/6)*(htsw[time]-hssw[time]) / (htsw[time])]
        Dsw+=[(7/6)*(htsw[time]-hssw[time]) / (1-hssw[time])]
        
    hssc=[]
    for time in range(5):
        hsreps=[]
        for rep in [x for x in set(bacter_all.Rep) if 'C' in x]:
            hs=1
            for gen in set(bacter_all[(bacter_all.Rep==rep) & (bacter_all.Time==time)].Genotype):
                hs = hs - (bacter_all[(bacter_all.Rep==rep) & (bacter_all.Genotype==gen) & (bacter_all.Time==time)].freq.iloc[0])**2
            hsreps+=[hs]
        hssc+=[sum(hsreps)/6]
    print('hssc')    
    htsc=[]
    for time in range(5):
        hsgen=[]
        for gen in set(bacter_all[(bacter_all.Rep.str.contains('C')) & (bacter_all.Time==time)].Genotype):
            hs=0
            for rep in [x for x in set(bacter_all.Rep) if 'C' in x]:
                try: 
                    hs = hs + (bacter_all[(bacter_all.Rep==rep) & (bacter_all.Genotype==gen) & (bacter_all.Time==time)].freq.iloc[0])
                except:
                    pass
            hsgen+=[(hs/6)**2]
        htsc+=[1-sum(hsgen)]
    print('htsc')    
    Dsc=[]
    fstsc=[]
    for time in range(5):
        fstsc+=[(6/5)*(htsc[time]-hssc[time]) / (htsc[time])]
        Dsc+=[(6/5)*(htsc[time]-hssc[time]) / (1-hssc[time])]

    if toremove != 8:
        df_D=pd.DataFrame({'D':Dsr+Dsw+Dsc, 'time':[0,1,2,3,4,0,1,2,3,4,0,1,2,3,4], 'cond':['R']*5+['W']*5+['C']*5})
    elif toremove == 8:
        df_D=pd.DataFrame({'D':Dsr+Dsw, 'time':[0,1,2,3,4,0,1,2,3,4], 'cond':['R']*5+['W']*5})
    

    return(df_D)

df_D_all=fst_jk(toremove=1)
for i in range(2,9):
    print(i)
    df_D_all= pd.concat([df_D_all, fst_jk(toremove=i)])
    

#df_D_all.to_csv('../steps/Fst/D_bacteria.csv')
df_D_all = pd.read_csv('../steps/Fst/D_bacteria.csv')
df_D_plot=copy.deepcopy(df_D_all)
df_D_plot.cond=df_D_plot.cond.apply(lambda x : 'M' if x=='W' else ('D' if x=='R' else x))
plt.figure()
sns.set_style('ticks')
sns.lineplot(data=df_D_plot, x='time', y='D', hue='cond', hue_order=['C','M','D'], palette=[colors[0], colors[1], colors[3]])
plt.xticks([0,1,2,3,4])
plt.xlabel('Time (days)',fontsize=12)
plt.ylabel('D',fontsize=12)
plt.ylim(-0.03,1)
handles, labels = plt.gca().get_legend_handles_labels()
plt.legend(handles=handles[1:], labels=labels[1:])
ax=plt.gca()
ax.get_legend().remove()
sns.despine()
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.savefig('../steps/Fst/D_bacteria_heterozygosities_CMD,jk.png', dpi=300, bbox_inches='tight')
plt.savefig('/home/guillemet/Documents/crispr/final_figures/S4_2.png', dpi=300, bbox_inches='tight')



