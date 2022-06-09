import pandas as pd
import os
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import math
import copy
from matplotlib.ticker import MaxNLocator

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

os.chdir('/home/guillemet/Documents/crispr/work/Martin/scripts')
colors= sns.color_palette()
sns.set_style("white")
phage_bacterC=pd.read_csv("../steps/df/phage_bacterC.csv")
df=phage_bacterC


def compute_fitness(t0, t1):
    return(np.log10(t1 * (1-t0)/(t0 * (1-t1))))

spacers=[]
reps=[]
t0=[]
t1=[]
for spacer in set(df.spacer):
    for rep in set(df.rep):
        dfg = df.loc[(df.spacer==spacer) & (df.rep==rep)]
        spacers+=[spacer]
        reps+=[rep]
        t0+=[dfg.loc[dfg.time==0,'freq'].iloc[0]]
        t1+=[dfg.loc[dfg.time==1,'freq'].iloc[0]]
        
ndf=pd.DataFrame({'spacer':spacers, 'rep':reps, 't0':t0, 't1':t1})
ndf['fitness'] = compute_fitness(ndf.t0, ndf.t1)
mean = ndf.loc[ndf.spacer==0, 'fitness'].mean()
ndf['fitness']=ndf['fitness']-mean
toplot=ndf.groupby('spacer').fitness.mean()
plt.figure()
sns.histplot(x=toplot, binwidth=0.1, binrange=(-1.05, 0.35))
ax=plt.gca()
ax.yaxis.set_major_locator(MaxNLocator(integer=True))
plt.plot([0,0],[1.05,2.5], 'k')
plt.text(0, 2.7, 'Susceptible strain', ha='center')
plt.ylabel('Number of strain')
plt.xlabel('Fitness')
plt.savefig('../steps/fitness/fitness_day1.png', dpi=300)

dic_color={}
for i in range(16):
    dic_color[df_color.loc[i, 'sp']] = df_color.loc[i, 'ncolor']
dic_color[0]='grey'


def getxy(i):
    i=15-i
    Z=0.007
    if i in [0,1,2,3,4,7,10,12,13,15]:
        y=0.5
    if i in[5,8,11,14]:
        y=1.5
    if i ==9:
        y=2.5
    if i==6:
        y=3.5
    if i==0:
        x=-1-Z
    if i==1:
        x=-0.9-Z
    if i==2:
        x=-0.8-Z
    if i==3:
        x=-0.7-Z
    if i in [4,5]:
        x=-0.6-Z
    if i in [6,7,8,9]:
        x=-0.5-Z
    if i ==10:
        x=-0.4-Z
    if i in [11,12]:
        x=-0.3-Z
    if i in [13,14]:
        x=-0.2-Z
    if i == 15:
        x=0.3-Z
    return(x,y)


    

plt.figure()
df=pd.DataFrame({'spacer':toplot.index, 'fitness':toplot.tolist()})
ax=plt.gca()
df=df.sort_values('fitness')
sns.histplot(data=df, ax=ax, stat="count", multiple="stack",
             x="fitness", kde=False,
             palette=dic_color, hue="spacer",
             element="bars", legend=False, binwidth=0.1, binrange=(-1.05, 0.35))
plt.xlabel('Fitness relative to the wild type',fontsize=12)
plt.ylabel('Number of strain',fontsize=12)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
ax.set_xticks([-1,-0.8,-0.6,-0.4,-0.2,0,0.2,0.4])
plt.tick_params(axis='x', bottom=True)
plt.tick_params(left=True)
plt.plot([0,0],[1.05,2.5], 'k')
plt.text(0, 2.7, 'Wild-Type strain', ha='center')
ax.yaxis.set_major_locator(MaxNLocator(integer=True))
sns.despine()

for i in range(len(myorder)-1,-1,-1):
    print(i)
    plt.text(getxy(i)[0], getxy(i)[1],myorder[i], rotation=270, ha='center',va='center')
plt.savefig('../steps/fitness/fitness_day1_color.png', dpi=300, bbox_inches='tight')
plt.savefig('/home/guillemet/Documents/crispr/final_figures/S1_labels.png', dpi=300, bbox_inches='tight')
plt.savefig('/home/guillemet/Documents/crispr/final_figures/S1.png', dpi=300, bbox_inches='tight')




