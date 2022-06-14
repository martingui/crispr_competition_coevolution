#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 19 14:01:01 2020

@author: guillemet
"""

import os
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
os.chdir('/home/guillemet/Documents/crispr/work/Martin/scripts')
colors= sns.color_palette()
Phages = pd.read_csv("../data/Phages.csv")
Bacteria = pd.read_csv("../data/Bacteria.csv")
Phages.loc[Phages.Time==0, 'Phages']=Phages.loc[Phages.Time==0, 'Phages']/10
Bacteria.loc[Bacteria.Time==0, 'Bacteries']=Bacteria.loc[Bacteria.Time==0, 'Bacteries']/10

Phages.Type=Phages.Type.apply(lambda x : 'M' if x=='WT' else ('D' if x=='Mut' else x))


Bacteria.Type=Bacteria.Type.apply(lambda x : 'M' if x=='WT' else ('D' if x=='Mut' else x))


plt.figure()
sns.set_style('ticks')
sns.lineplot(x='Time', y='Bacteries', hue='Type', units='Echantillons',data=Bacteria.loc[Bacteria.Type!='x'], estimator=None, hue_order=['C','M','D'], \
             palette=[colors[0],colors[1],colors[3]], \
             linewidth=0.8)
sns.scatterplot(x='Time', y='Bacteries', hue='Type',style='Type' ,markers=['o','s','P'], units='Echantillons',data=Bacteria.loc[Bacteria.Type!='x'], estimator=None, hue_order=['C','M','D'], \
             palette=[colors[0],colors[1],colors[3]])
plt.yscale("log")
plt.ylabel("Bacteria density (cfu/mL)", fontsize=12)
plt.xlabel("Time (days)", fontsize=12)
plt.ylim(10,1.5e10)
handles, labels = plt.gca().get_legend_handles_labels()
handles, labels = plt.gca().get_legend_handles_labels()
labels=['A','B','C','A','B','C']
plt.legend(handles=handles[3:], labels=labels[3:])

ax=plt.gca()
ax.xaxis.set_major_locator(MaxNLocator(integer=True))
sns.despine()
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.savefig('../steps/Demography/bacteria_demo_ppt.png', dpi=300, bbox_inches='tight')
plt.savefig('/home/guillemet/Documents/crispr/final_figures/2_1.png', dpi=300, bbox_inches='tight')


plt.figure()
sns.set_style('ticks')
sns.set_style('ticks')
sns.lineplot(x='Time', y='Phages', hue='Type', units='Echantillons',data=Phages[Phages.Type!='C'], estimator=None, hue_order=['M','D'], \
             palette=[colors[1],colors[3]], \
             linewidth=0.8)
sns.scatterplot(x='Time', y='Phages', hue='Type', style='Type',markers=['s','P'],units='Echantillons',data=Phages[Phages.Type!='C'], estimator=None, hue_order=['M','D'], \
             palette=[colors[1],colors[3]])    
plt.yscale("log")
plt.ylabel("Phage density (pfu/mL)", fontsize=12)
plt.xlabel("Time (days)", fontsize=12)
plt.ylim(10,1.5e10)
handles, labels = plt.gca().get_legend_handles_labels()
labels=['B','C','B','C']
plt.legend(handles=handles[2:], labels=labels[2:])
ax=plt.gca()
# ax.get_legend().remove()
ax.xaxis.set_major_locator(MaxNLocator(integer=True))
sns.despine()
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.savefig('../steps/Demography/phages_demo_ppt_nod.png', dpi=300, bbox_inches='tight')
plt.savefig('/home/guillemet/Documents/crispr/final_figures/2_2.png', dpi=300, bbox_inches='tight')