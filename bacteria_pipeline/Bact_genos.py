import numpy as np
import pandas as pd
import re
import os
import copy

#set current directory to save file of this script
os.chdir(os.path.dirname(os.path.abspath("__file__")))
os.chdir('/home/guillemet/Documents/crispr/work/Martin/scripts/Bacteria')

"""
Script that read in blasted bacterial spacers. The PAM dictionary is used to find the 
middle position of the targeted escape region in the phage. The bacterial genotype is then 
constituted of genomic position(s), that is the middle nucleotide in an escape region. 

Input: inter4_replaced_replaced_sort_count for all samples (W1-W8 and R1-R8, for all T1-T4)

output: Bacteria_genos.csv 		A dataframe with:
										Genotype: Dict name of spacer
										freq: Frequency of that genotype
										Rep: Replicat (W1-W8 or R1-R8)
										Cond: WT for W samples and Mut for R samples
										Time: 1-4
										Spacer: Middle position of escape region
		no_match_R.txt 			A txt with percentage of non-matched genotypes for all R samples
		no_match_W.txt			A txt with percentage of non-matched genotypes for all W samples

"""

## Path to sample folder with the bacteria blast

threshold=0.01

t0=["Mix-B0-A","Mix-B0-B","Mix-B0-C"]
def make_data0():
    path = "../../data/samples/"
    data=pd.DataFrame()
    for rep in t0:
        data=pd.concat([data,pd.read_csv(path+rep+"/inter4_replaced_sort_count",header=None, delimiter='\s+')], ignore_index=True)
    data.columns=['count','Genotype']
    data0=pd.DataFrame(columns=['count','Genotype'])
    for spacer in set(data.Genotype):
        data0.loc[len(data0)]=[sum(data[data.Genotype==spacer]['count']), spacer]
    sum_all=sum(data0['count'])
    data0=data0[data0['count']>threshold*sum_all]
    sum_sub=sum(data0['count'])           
    data0['freq']=data0['count']/sum_sub     
    data0['Time']=0  
    data0['Spacer']=data0.Genotype.apply(lambda x: x.split('-')[1])
    data0['Spacer']=data0['Spacer'].apply(lambda x: x.split('_')[1] if x!='fin' else 0)
    #data0.to_csv('../../data/data0.csv', sep=',')

def BactGenotypes(samples, cond, data0, PAM_dict):
    non_match_geno = {}
    data_list = []
    
    for sample in samples:
        path = "../../data/samples/"
        path = "".join([path,sample,"/inter4_replaced_sort_count"])
        
        with open(path) as file:
            genos = [line.lstrip(' ').rstrip("\n").split(" ") for line in file]
        
        for geno in genos:
            geno[0]=int(geno[0])
        
        bacteria_genos = list()
        non_matches = list()
        bacteria_genos_count = list()
        non_matches_count = list()
        for geno in genos: 
            l = [n.start() for n in re.finditer('-', geno[1])]
            if any(value > 15 for value in [l[x] - l[x-1] for x in range(len(l))][1:]):
                non_matches_count.append(int(geno[0]))
                non_matches.append(geno)
            else: 
                bacteria_genos_count.append(int(geno[0]))
                bacteria_genos.append(geno)
                
        non_match_geno[sample.strip("P")] = (sum(non_matches_count)/(sum(bacteria_genos_count)+sum(non_matches_count)))*100
        df = pd.DataFrame(bacteria_genos, columns=["Count","Genotype"])
        df = df.assign(freq = df['Count'] / sum(df['Count']))
        df = df.loc[df['freq'] > 0.01] ## Applying filter of 0.01 to genotype frequencies
        df.freq = df['Count'] / sum(df['Count'])
        
        df = df.drop(['Count'], axis=1)
        
        df['Time'] = [int(sample[-1:])] * len(df)
        df['Cond'] = [cond] * len(df)
        df['Rep'] = [sample.strip("P")[:-2].upper()] * len(df)
        
        
        if int(sample[-1:]) == 1:
            data0['Cond'] = [cond] * len(data0)
            data0['Rep'] = [sample.strip("P")[:-2].upper()] * len(data0)
            
            df = pd.concat([data0,df]) ## Appending T0 genotypes to dataframe if time point T1
        
        df = df.reset_index()
            
        pos_list = list()

        for i in range(len(df)):
            pos = re.findall(r'\_(.*?)\-', df[i:i+1]["Genotype"][i])
            if len(pos) == 0:
                pos_list.append(['0'])
            else:
                pos_list.append(pos)
        
        new_df = pd.DataFrame(pos_list, index=[df.Genotype,df.freq,df.Rep,df.Cond,df.Time]).stack()

        new_df = new_df.reset_index()
        
        new_df = new_df.drop(["level_5"], axis = 1)

        new_df = new_df.rename(columns= {0 : "Pos"})
        
        ## Make new column with the middle position in the spacer region
        spacers = []

        for i in range(len(new_df.Pos)):
            if len(PAM_dict.POS[PAM_dict.POS == new_df.Pos[i]].index.tolist()) == 0:
                spacers.append(0)
            else: 
                index = PAM_dict.POS[PAM_dict.POS == new_df.Pos[i]].index.tolist()[0]
                if PAM_dict.strand[index] == "-":
                    spacers.append(int(new_df.Pos[i])+17)
                elif PAM_dict.strand[index] == "+":
                    spacers.append(int(new_df.Pos[i])-14)
                    
        
        new_df = new_df.drop(["Pos"], axis = 1)
        new_df["Spacer"] = spacers
        
        data_list.append(new_df)
    
    big_data = pd.concat(data_list)
    big_data = big_data.reset_index(drop=True)
    
    return big_data, non_match_geno


## Read in PAM_dict, to look up the spacer name and find the middle position of its escape 
# region in the phage genome.
PAM_dict = pd.read_csv("../../data/PAM_positions.csv")
pos_list = list()
for seq_name in PAM_dict["seq_name"]:
    pos_list.append(re.findall(r'PAM_(.*?)_NC', seq_name)[0])
PAM_dict['POS'] = pos_list


## Read in T0 data to append to all samples. 
# data0.csv is made based on MIX-BO-A, MIX-BO-B, MIX-BO-C. They all include the same 17 genotypes
# when filtering for freq > 0.01
# data0.csv hold the average frequency for the 17 genotypes, in these 3 samples.
data0 = pd.read_csv("../../data/data0.csv")


R_samples = ["PR1T1","PR1T2","PR1T3","PR1T4",
"PR2T1","PR2T2","PR2T3","PR2T4",
"PR3T1","PR3T2","PR3T3","PR3T4",
"PR4T1","PR4T2","PR4T3","PR4T4",
"PR5T1","PR5T2","PR5T3","PR5T4",
"PR6T1","PR6T2","PR6T3","PR6T4",
"PR7T1","PR7T2","PR7T3","PR7T4",
"PR8T1","PR8T2","PR8T3","PR8T4"]

W_samples = ["Pw1T1","Pw1T2","Pw1T3","Pw1T4",
"Pw2T1","Pw2T2","Pw2T3","Pw2T4",
"Pw3T1","Pw3T2","Pw3T3","Pw3T4",
"Pw4T1","Pw4T2","Pw4T3","Pw4T4",
"Pw5T1","Pw5T2","Pw5T3","Pw5T4",
"Pw6T1","Pw6T2","Pw6T3","Pw6T4",
"Pw7T1","Pw7T2","Pw7T3","Pw7T4",
"Pw8T1","Pw8T2","Pw8T3","Pw8T4"]

C_samples = ['B'+str(x)+'T'+str(y) for x in list(range(1,8)) for y in list(range(1,5))]

R_data, no_match_R = BactGenotypes(R_samples, "Mut", data0, PAM_dict)
W_data, no_match_W = BactGenotypes(W_samples, "WT", data0, PAM_dict)
C_data, no_match_C = BactGenotypes(C_samples, "C", data0, PAM_dict)

all_data = pd.concat([W_data, R_data, C_data])

### THIS CHANGES THE NAME OF THE GENOTYPE, THE SPACER INSIDE ARE NAMED WITH THE MIDDLE 
### POSITION ON THE PHAGE GENOME AND NOT THE START OR END
bacter_all=copy.deepcopy(all_data)
bacter_all=bacter_all.reset_index(drop=True)
for i in bacter_all.index:
    gen=bacter_all.loc[i,'Genotype'].split('-')[1:-1]
    l=[]
    for pam in gen:
        pos=int(pam.split('_')[1])
        if pos+17 in list(bacter_all.Spacer):
            l+=['PAM_'+str(pos+17)]
        elif pos-14 in list(bacter_all.Spacer):
            l+=['PAM_'+str(pos-14)]
        else:
            pass
    nl='-'.join(l)
    nl='start-'+nl+'-fin'
    if nl == 'start--fin':
        nl='start-fin'
    bacter_all.loc[i,'Genotype']=nl
    
    if 'B' in bacter_all.loc[i,'Rep']:
        bacter_all.loc[i,'Rep']='C'+bacter_all.loc[i,'Rep'][1]




all_data.to_csv('../../data/Bacteria_genos.csv', index=False)

f = open("../../data/no_match_R.txt", "w")
f.write(str(no_match_R))
f.close()

f = open("../../data/no_match_W.txt", "w")
f.write(str(no_match_W))
f.close()

R_data=R_data.drop(columns='Spacer')
R_data=R_data.drop_duplicates()
sums=[]
for i in set(R_data.Rep):
    for j in set(R_data.Time):
        print(sum(R_data[(R_data.Rep==i) &(R_data.Time==j)]['freq']))




