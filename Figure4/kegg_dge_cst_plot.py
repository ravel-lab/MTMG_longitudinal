#!/usr/bin/env python3

#importing packages to be used
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import sys
import matplotlib
import random
import matplotlib.patches as mpatches

pd.options.mode.chained_assignment = None
matplotlib.rc('font', serif='Helvetica')
matplotlib.rc('axes',linewidth=2)

data = pd.read_csv("MT_MG_KEGG_DGE.csv",sep=",",index_col=0)

#filtering out KOs which did not appear in the MT data
data['MT_FC'] = data['MT_FC'].fillna(-999)
data = data[data['MT_FC'] != -999]


stacked_fig, stacked_axs = plt.subplots(3,1, figsize=(4.25,12), facecolor='w', edgecolor='k')
stacked_fig.subplots_adjust(left=0.18,right=0.92,bottom=0.18,top=0.95, wspace=0.3,hspace=0.3)

stacked_axs = stacked_axs.ravel()

plot_count = 0

for compar in ['I_III','I_IV','III_IV']:

    data_compar = data[data['Comparison']==compar]
    
    data_compar_mgsig = data_compar[data_compar['MG_p'] <= 0.05]
    data_compar_mginsig = data_compar[data_compar['MG_p'] > 0.05]

    data_compar_mgsig['logMG_p'] = np.log10(data_compar_mgsig['MG_p'])*-1.0
    data_compar_mgsig['logMT_p'] = np.log10(data_compar_mgsig['MT_p'])*-1.0
    data_compar_mginsig['logMG_p'] = np.log10(data_compar_mginsig['MG_p'])*-1.0
    data_compar_mginsig['logMT_p'] = np.log10(data_compar_mginsig['MT_p'])*-1.0

    stacked_axs[plot_count].plot(data_compar_mgsig['MT_FC'],data_compar_mgsig['logMT_p'],'o',c='#327fa8',ms=1.75)
    stacked_axs[plot_count].plot(data_compar_mginsig['MT_FC'],data_compar_mginsig['logMT_p'],'d',mfc='none',c='#e8b03f',ms=1.75)

    stacked_axs[plot_count].set_xlim(-10,10)
    stacked_axs[plot_count].set_ylim(0,25)
    stacked_axs[plot_count].set_xticks([-10,-8,-6,-4,-2,0,2,4,6,8,10])
    stacked_axs[plot_count].set_yticks([0,5,10,15,20,25])
    stacked_axs[plot_count].set_xticklabels(['10','8','6','4','2','0','2','4','6','8','10'],fontsize=14)
    stacked_axs[plot_count].set_yticklabels(['0','5','10','15','20','25'],fontsize=14)
    stacked_axs[plot_count].set_ylabel('-log10 p-value',fontsize=16)
    stacked_axs[plot_count].set_xlabel('Fold Change',fontsize=16)
    stacked_axs[plot_count].axvline(x=0,color="k")
    stacked_axs[plot_count].set_title('%s                     %s' %(compar.split("_")[0],compar.split("_")[1]),fontsize=18,y=0.85)
    stacked_axs[plot_count].invert_xaxis()

    plot_count+=1

stacked_fig.savefig('KEGG_DGE.pdf')