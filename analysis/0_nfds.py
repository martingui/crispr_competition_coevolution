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

sns.set_style("white")
def growth_reg_03_genotype(phage_bacter,title):
    '''
    
    Parameters
    ----------
    phage_bacter : phage_bacterR or phage_bacterW
    title : 'R' or 'W'

    print a plot of spacer frequency between t,t+1 change as f(spacer frequency )
    '''
    nbacter_all=phage_bacter
    nbacter_all=nbacter_all.loc[nbacter_all.Rep.str.contains(title)]
                    
    fig=plt.figure()
    nbacter_all=nbacter_all[(nbacter_all.freq!=0) & (nbacter_all.freq_growth!=0)]
    df03=nbacter_all[nbacter_all.Time!=4]
    
    if title=='C':
        couleurs=[colors[0]]
    elif title=='W':
        couleurs=[colors[1]]
    elif title=='R':
        couleurs=[colors[3]]
    
    for i in range(3):
        couleurs.append(tuple([x+(1-x)*0.35 for x in  couleurs[-1]]))
        
    couleurs.reverse()
    print(couleurs)
    
    plt.figure(figsize=(4.2,4))
    plt.fill([-1,-1,2],[1,-2,-2], color=(0.97, 0.97, 0.97), zorder=0)
    plt.fill([-1,2,2],[2,-1,2], color=(0.97, 0.97, 0.97), zorder=0)
    sns.scatterplot(df03.freq, df03.freq_growth, hue=df03.Time, palette=couleurs, hue_order=[0,1,2,3]) 
    plt.ylim(-1,1)
    plt.xlim(0,1)

    
    reg0=LinearRegression().fit(np.array(df03.freq).reshape((-1,1)),np.array(df03.freq_growth))
    plt.plot([0,1],[reg0.intercept_,reg0.intercept_+reg0.coef_], 'k:')
    print(reg0.coef_)
    plt.xlabel('Genotype frequency at day t',fontsize=12)
    plt.ylabel(r'$\Delta$'+' genotype frequency between t and t+1',fontsize=12)
    legend = plt.gca().legend()
    plt.legend(title='Time (days)',fontsize=12)
    #legend.texts[0].set_text("Day")
    
    #
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.tick_params(axis='x', bottom=True)
    plt.tick_params(left=True)
    plt.savefig('../steps/Regression-NFDS/genotype_'+title+'_ppt.pdf', bbox_inches='tight')
    if title=='C':
        num=1
    elif title=="W":
        num=2
    elif title=='R':
        num=3
    plt.savefig('/home/guillemet/Documents/crispr/final_figures/5_'+str(num)+'.pdf', bbox_inches='tight')
    plt.show()
    return

growth_reg_03_genotype(nbacter_all,'W')
growth_reg_03_genotype(nbacter_all,'R')
growth_reg_03_genotype(nbacter_all,'C')

handles, labels = plt.gca().get_legend_handles_labels()
plt.legend(handles=handles[1:], labels=labels[1:])



import statsmodels.formula.api as smf
import statsmodels.api as sm
reg_testdf=copy.deepcopy(nbacter_all)
reg_testdf=reg_testdf.loc[reg_testdf.Time!=4]

reg_testdf['cond']=reg_testdf['Rep']
reg_testdf['cond']=reg_testdf['cond'].apply(lambda x: str(x[0]))
reg_testdf=reg_testdf.loc[reg_testdf.cond!='W']
lm1 = smf.ols(formula='freq_growth ~ freq * cond', data=reg_testdf).fit()
aov=sm.stats.anova_lm(lm1, type=2)
print(aov)


