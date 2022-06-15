import pandas as pd
import os
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import math
import copy
from matplotlib.ticker import MaxNLocator
from sklearn.linear_model import LinearRegression
import Bio
import re
import random
from Bio import SeqIO
os.chdir('/home/guillemet/Documents/crispr/work/Martin/scripts')


p=18
record = SeqIO.read("../data/TheOneWeUse_T0.fasta", "fasta")

windows3=[]
pattern=re.compile('[A-Z][G][G][A-Z][G]')
for m in pattern.finditer(str(record.seq)):
    windows3+=[[m.start()+5-2*p, m.start()+5]]
    
    
pattern=re.compile('[A-Z][G][G][A-Z][G]')
for m in pattern.finditer(str(record.seq.reverse_complement())):
    windows3+=[[len(record.seq)-(m.start()), len(record.seq)-m.start()+2*p]]
    
    



dfmutr0=pd.read_csv('../data/FreeBayes/Other/T0_R_data.csv',  delimiter='\t')
dfmutr0=copy.deepcopy(dfmutr0)


dfmutw0=pd.read_csv('../data/FreeBayes/Other/T0_W_data.csv',  delimiter='\t')
dfmutw0=copy.deepcopy(dfmutw0)
#dfmutr0=dfmutr0.loc[dfmutr0.FREQ>0.2]

dfpam=pd.read_csv("../../Aarhus/data/PAM_positions.csv")
dfmuts=pd.DataFrame()
for i in range(1,9):
    dfmut=pd.read_csv('../data/FreeBayes/R_seq/R'+str(i)+'_data.csv', delimiter='\t')
    dfmut.TIME=dfmut.TIME.apply(lambda x: 'T'+x[1] if 'R' in x else x )
    dfmut['rep']='R'+str(i)
    dfmuts=pd.concat([dfmuts,dfmut])
dfmuts['escape_pos']=-1
newdf=pd.DataFrame(columns=dfmuts.columns)
for i in range(len(dfmuts)):
    pos=dfmuts.POS.iloc[i]
    filt_df=dfpam[((dfpam.escape_start<pos) & (dfpam.escape_end>pos)) | ((dfpam.escape_start>pos) & (dfpam.escape_end<pos))]
    if len(filt_df)==1:
        start=filt_df.escape_start.iloc[0]
        end=filt_df.escape_end.iloc[0]
        if filt_df.strand.iloc[0] == '+':
            dfmuts.escape_pos.iloc[i] = end - pos + 1
        elif filt_df.strand.iloc[0] == '-':
            dfmuts.escape_pos.iloc[i] = pos - start + 1
        else:
            dfmuts.escape_pos.iloc[i] = -1000
    elif len(filt_df)>1:
        start=filt_df.escape_start.iloc[0]
        end=filt_df.escape_end.iloc[0]
        if filt_df.strand.iloc[0] == '+':
            dfmuts.escape_pos.iloc[i] = end - pos + 1
        elif filt_df.strand.iloc[0] == '-':
            dfmuts.escape_pos.iloc[i] = pos - start + 1
        else:
            dfmuts.escape_pos.iloc[i] = -1000
        newdf.loc[len(newdf)]=dfmuts.iloc[i]
        start=filt_df.escape_start.iloc[1]
        end=filt_df.escape_end.iloc[1]
        if filt_df.strand.iloc[1] == '+':
            newdf.escape_pos.loc[len(newdf)-1] = end - pos + 1
        elif filt_df.strand.iloc[1] == '-':
            newdf.escape_pos.loc[len(newdf)-1] = pos - start + 1
        else:
            newdf.escape_pos.loc[len(newdf)-1] = -1000


dfmuts=pd.concat([dfmuts, newdf])
dfmuts=dfmuts.loc[dfmuts.POS!=34230]

def sort_escape(x):
    if x<0:
        return('non-escape')
    elif x>0:
        return('escape')
    
def is_cr3(x):
    for win in windows3:
        if x>win[0] and x<win[1]:
            return True
    return False
    
dfmuts.escape_pos=dfmuts.escape_pos.apply(sort_escape)
dfmuts=dfmuts.drop_duplicates()
dfmutsDall=copy.deepcopy(dfmuts)
#dfmuts=dfmuts.loc[dfmuts.POS.isin(dfmutr0.POS)==False]
dfmuts=dfmuts.loc[dfmuts.TYPE.isin(['snp','mnp'])]
dfmuts=dfmuts.loc[dfmuts.TIME=='T4']
dfmuts['cr3']=dfmuts.POS.apply(is_cr3)
dfmutsD=dfmuts
df=dfmuts.loc[dfmuts.escape_pos=='non-escape']



df=df.sort_values('FREQ')
df['x']=range(len(df),0,-1)
df.x=df.x.apply(lambda x: x/len(df))

dfx=df.FREQ.apply(lambda x: math.log10(x))
dfy=df.x.apply(lambda x: math.log10(x))
reg0=LinearRegression().fit(np.array(dfx).reshape((-1,1)),np.array(dfy))
print(reg0.coef_)
c=reg0.coef_
i=reg0.intercept_


df1=dfmuts.loc[dfmuts.escape_pos=='escape']

df1=df1.sort_values('FREQ')
df1['x']=range(len(df1),0,-1)
df1.x=df1.x.apply(lambda x: x/len(df1))

df1x=df1.FREQ.apply(lambda x: math.log10(x))
df1y=df1.x.apply(lambda x: math.log10(x))
reg0=LinearRegression().fit(np.array(df1x).reshape((-1,1)),np.array(df1y))
print(reg0.coef_)
c1=reg0.coef_
i1=reg0.intercept_



hist, edges = np.histogram(df.FREQ, [x/10 for x in range(0,11)])
freq = hist/float(hist.sum())
plt.figure()
plt.bar([(x/10)+0.005 for x in range(0,10)], freq, width=0.09, align="edge", ec="k" )

plt.ylim(0,0.7)
plt.text(0.4,0.6,str(hist.sum())+" mutations")
plt.title('non-escape D')
plt.savefig('../steps/mutations/D_non_escape.pdf')


oldhist, edges = np.histogram(df.FREQ, [x/10 for x in range(0,11)])
histcr3, edges = np.histogram(df.loc[df.cr3==True,'FREQ'], [x/10 for x in range(0,11)])
freqcr3 = histcr3/float(oldhist.sum())
hist, edges = np.histogram(df.loc[df.cr3==False,'FREQ'], [x/10 for x in range(0,11)])
freq = hist/float(oldhist.sum())
plt.figure()
plt.bar([(x/10)+0.005 for x in range(0,10)], freq, width=0.09, align="edge", ec="k" )
plt.bar([(x/10)+0.005 for x in range(0,10)], freqcr3, width=0.09, bottom=freq,align="edge", ec="r" )
plt.ylim(0,0.7)
hist, edges = np.histogram(df.FREQ, [x/10 for x in range(0,11)])
plt.text(0.4,0.6,str(oldhist.sum())+" mutations")
plt.title('non-escape D')
plt.savefig('../steps/mutations/D_non_escape_cr3.pdf')


hist, edges = np.histogram(df1.FREQ, [x/10 for x in range(0,11)])
freq = hist/float(hist.sum())
plt.figure()
plt.bar([(x/10)+0.005 for x in range(0,10)], freq, width=0.09, align="edge", ec="k" )
plt.ylim(0,0.7)
plt.text(0.4,0.6,str(hist.sum())+" mutations")
plt.title('escape D')
plt.savefig('../steps/mutations/D_escape.pdf')


###

###





dfpam=pd.read_csv("../../Aarhus/data/PAM_positions.csv")
dfmuts=pd.DataFrame()
for i in range(1,9):
    dfmut=pd.read_csv('../data/FreeBayes/W_seq/W'+str(i)+'_data.csv', delimiter='\t')
    dfmut.TIME=dfmut.TIME.apply(lambda x: 'T'+x[1] if 'R' in x else x )
    dfmut['rep']='W'+str(i)
    dfmuts=pd.concat([dfmuts,dfmut])
dfmuts['escape_pos']=-1
newdf=pd.DataFrame(columns=dfmuts.columns)
for i in range(len(dfmuts)):
    pos=dfmuts.POS.iloc[i]
    filt_df=dfpam[((dfpam.escape_start<pos) & (dfpam.escape_end>pos)) | ((dfpam.escape_start>pos) & (dfpam.escape_end<pos))]
    if len(filt_df)==1:
        start=filt_df.escape_start.iloc[0]
        end=filt_df.escape_end.iloc[0]
        if filt_df.strand.iloc[0] == '+':
            dfmuts.escape_pos.iloc[i] = end - pos + 1
        elif filt_df.strand.iloc[0] == '-':
            dfmuts.escape_pos.iloc[i] = pos - start + 1
        else:
            dfmuts.escape_pos.iloc[i] = -1000
    elif len(filt_df)>1:
        start=filt_df.escape_start.iloc[0]
        end=filt_df.escape_end.iloc[0]
        if filt_df.strand.iloc[0] == '+':
            dfmuts.escape_pos.iloc[i] = end - pos + 1
        elif filt_df.strand.iloc[0] == '-':
            dfmuts.escape_pos.iloc[i] = pos - start + 1
        else:
            dfmuts.escape_pos.iloc[i] = -1000
        newdf.loc[len(newdf)]=dfmuts.iloc[i]
        start=filt_df.escape_start.iloc[1]
        end=filt_df.escape_end.iloc[1]
        if filt_df.strand.iloc[1] == '+':
            newdf.escape_pos.loc[len(newdf)-1] = end - pos + 1
        elif filt_df.strand.iloc[1] == '-':
            newdf.escape_pos.loc[len(newdf)-1] = pos - start + 1
        else:
            newdf.escape_pos.loc[len(newdf)-1] = -1000


dfmuts=pd.concat([dfmuts, newdf])
dfmuts=dfmuts.loc[dfmuts.POS!=34230]

def sort_escape(x):
    if x<0:
        return('non-escape')
    elif x>0:
        return('escape')
    
dfmuts.escape_pos=dfmuts.escape_pos.apply(sort_escape)
dfmuts=dfmuts.drop_duplicates()
dfmutsMall=copy.deepcopy(dfmuts)
dfmuts=dfmuts.loc[dfmuts.POS.isin(dfmutr0.POS)==False]
dfmuts=dfmuts.loc[dfmuts.TYPE.isin(['snp','mnp'])]
dfmuts=dfmuts.loc[dfmuts.TIME=='T4']
dfmuts['cr3']=dfmuts.POS.apply(is_cr3)
dfmutsM=dfmuts
df=dfmuts.loc[dfmuts.escape_pos=='non-escape']


df=df.sort_values('FREQ')
df['x']=range(len(df),0,-1)
df.x=df.x.apply(lambda x: x/len(df))

dfx=df.FREQ.apply(lambda x: math.log10(x))
dfy=df.x.apply(lambda x: math.log10(x))
reg0=LinearRegression().fit(np.array(dfx).reshape((-1,1)),np.array(dfy))
print(reg0.coef_)
c=reg0.coef_
i=reg0.intercept_


df1=dfmuts.loc[dfmuts.escape_pos=='escape']

df1=df1.sort_values('FREQ')
df1['x']=range(len(df1),0,-1)
df1.x=df1.x.apply(lambda x: x/len(df1))

df1x=df1.FREQ.apply(lambda x: math.log10(x))
df1y=df1.x.apply(lambda x: math.log10(x))
reg0=LinearRegression().fit(np.array(df1x).reshape((-1,1)),np.array(df1y))
print(reg0.coef_)
c1=reg0.coef_
i1=reg0.intercept_



hist, edges = np.histogram(df.FREQ, [x/10 for x in range(0,11)])
freq = hist/float(hist.sum())
plt.figure()
plt.bar([(x/10)+0.005 for x in range(0,10)], freq, width=0.09, align="edge", ec="k" )
plt.ylim(0,0.7)
plt.text(0.4,0.6,str(hist.sum())+" mutations")
plt.title('non-escape M')
plt.savefig('../steps/mutations/M_non_escape.pdf')


hist, edges = np.histogram(df1.FREQ, [x/10 for x in range(0,11)])
freq = hist/float(hist.sum())
plt.figure()
plt.bar([(x/10)+0.005 for x in range(0,10)], freq, width=0.09, align="edge", ec="k" )
plt.ylim(0,0.7)
plt.text(0.4,0.6,str(hist.sum())+" mutations")
plt.title('escape M')
plt.savefig('../steps/mutations/M_escape.pdf')


oldhist, edges = np.histogram(df.FREQ, [x/10 for x in range(0,11)])
histcr3, edges = np.histogram(df.loc[df.cr3==True,'FREQ'], [x/10 for x in range(0,11)])
freqcr3 = histcr3/float(oldhist.sum())
hist, edges = np.histogram(df.loc[df.cr3==False,'FREQ'], [x/10 for x in range(0,11)])
freq = hist/float(oldhist.sum())
plt.figure()
plt.bar([(x/10)+0.005 for x in range(0,10)], freq, width=0.09, align="edge", ec="k" )
plt.bar([(x/10)+0.005 for x in range(0,10)], freqcr3, width=0.09, bottom=freq,align="edge", ec="r" )
plt.ylim(0,0.7)
hist, edges = np.histogram(df.FREQ, [x/10 for x in range(0,11)])
plt.text(0.4,0.6,str(oldhist.sum())+" mutations")
plt.title('non-escape M')
plt.savefig('../steps/mutations/M_non_escape_cr3.pdf')




#### number of phages muta >.05

dfmutsMall=pd.concat([dfmutw0, dfmutsMall])
dfmutsMall=dfmutsMall.drop_duplicates()
dfmutsMall=dfmutsMall.loc[dfmutsMall.FREQ>0.025]
dfmutsDall=pd.concat([dfmutr0, dfmutsDall])
dfmutsDall=dfmutsDall.drop_duplicates()
dfmutsDall=dfmutsDall.loc[dfmutsDall.FREQ>0.025]

cpt=[0]*5
cpt_esc=[0]*5
for i in range(5):
    cpt[i]=len(dfmutsMall.loc[dfmutsMall.TIME=='T'+str(i)])
    cpt_esc[i]=len(dfmutsMall.loc[(dfmutsMall.TIME=='T'+str(i)) & (dfmutsMall.escape_pos=='escape')])

time=[0]*8
rep=list(range(1,9))
val=[len(dfmutsMall.loc[(dfmutsMall.TIME=='-W')])]*8
for i in range(1,5):
    for r in range(1,9):
        val+=[len(dfmutsMall.loc[(dfmutsMall.TIME=='T'+str(i)) & (dfmutsMall.rep=='W'+str(r))])]
        rep+=[r]
        time+=[i]
dfmall=pd.DataFrame({'time':time, 'rep':rep, 'val':val})

time=[0]*8
rep=list(range(1,9))
val=[len(dfmutsMall.loc[(dfmutsMall.TIME=='-W') & (dfmutsMall.escape_pos=='escape')])]*8
for i in range(1,5):
    for r in range(1,9):
        val+=[len(dfmutsMall.loc[(dfmutsMall.TIME=='T'+str(i)) & (dfmutsMall.rep=='W'+str(r)) & (dfmutsMall.escape_pos=='escape')])]
        rep+=[r]
        time+=[i]
dfmallesc=pd.DataFrame({'time':time, 'rep':rep, 'val':val})

colors= sns.color_palette()

plt.figure()
plt.plot(range(5),cpt, color=colors[1], ls='--')
plt.plot(range(5),cpt_esc, color=colors[1])
sns.despine()
ax=plt.gca()
ax.xaxis.set_major_locator(MaxNLocator(integer=True))
plt.xlabel('Time (days)', fontsize=12)
plt.ylabel('Number of phage mutations', fontsize=12)

plt.figure()
sns.lineplot(data=dfmall, x='time', y='val', color=colors[1], alpha=0.6, err_kws={'alpha':0.05})
sns.lineplot(data=dfmallesc, x='time', y='val', color=colors[1])
ax=plt.gca()
ax.lines[0].set_linestyle("--")
ax.xaxis.set_major_locator(MaxNLocator(integer=True))
plt.ylim(0,50)
plt.xlabel('Time (days)', fontsize=12)
plt.ylabel('Number of phage mutations', fontsize=12)
sns.despine()
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.savefig('../steps/phage_mutation/number_of_mutations_M.pdf', bbox_inches='tight')
plt.savefig('/home/guillemet/Documents/crispr/final_figures/S7_1.pdf', bbox_inches='tight')


cpt=[0]*5
cpt_esc=[0]*5
for i in range(5):
    cpt[i]=len(dfmutsDall.loc[dfmutsDall.TIME=='T'+str(i)])
    cpt_esc[i]=len(dfmutsDall.loc[(dfmutsDall.TIME=='T'+str(i)) & (dfmutsDall.escape_pos=='escape')])

time=[0]*8
rep=list(range(1,9))
val=[len(dfmutsDall.loc[(dfmutsDall.TIME=='-M')])]*8
for i in range(1,5):
    for r in range(1,9):
        val+=[len(dfmutsDall.loc[(dfmutsDall.TIME=='T'+str(i)) & (dfmutsDall.rep=='R'+str(r))])]
        rep+=[r]
        time+=[i]
dfpall=pd.DataFrame({'time':time, 'rep':rep, 'val':val})

time=[0]*8
rep=list(range(1,9))
val=[16]*8
for i in range(1,5):
    for r in range(1,9):
        val+=[len(dfmutsDall.loc[(dfmutsDall.TIME=='T'+str(i)) & (dfmutsDall.rep=='R'+str(r)) & (dfmutsDall.escape_pos=='escape')])]
        rep+=[r]
        time+=[i]
dfpallesc=pd.DataFrame({'time':time, 'rep':rep, 'val':val})

plt.figure()
plt.plot(range(5),cpt, color=colors[3], ls='--')
plt.plot(range(5),cpt_esc, color=colors[3])
sns.despine()
ax=plt.gca()
ax.xaxis.set_major_locator(MaxNLocator(integer=True))
plt.xlabel('Time (days)', fontsize=12)
plt.ylabel('Number of phage mutations', fontsize=12)


plt.figure()
sns.lineplot(data=dfpall, x='time', y='val', color=colors[3], alpha=0.6, err_kws={'alpha':0.05})
sns.lineplot(data=dfpallesc, x='time', y='val', color=colors[3])
ax=plt.gca()
ax.lines[0].set_linestyle("--")
ax.xaxis.set_major_locator(MaxNLocator(integer=True))
plt.ylim(0,50)
plt.xlabel('Time (days)', fontsize=12)
plt.ylabel('Number of phage mutations', fontsize=12)
sns.despine()
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.savefig('../steps/phage_mutation/number_of_mutations_P.pdf', bbox_inches='tight')
plt.savefig('/home/guillemet/Documents/crispr/final_figures/S7_2.pdf', bbox_inches='tight')



