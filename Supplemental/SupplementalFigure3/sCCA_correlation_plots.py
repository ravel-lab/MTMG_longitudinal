#!/usr/bin/python

#importing packages to be used
import  matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import pandas as pd
import numpy as np
import sys
import matplotlib 
import seaborn as sns
import scipy.cluster.hierarchy as hac
import operator
from scipy.spatial.distance import pdist, squareform
from scipy.cluster.hierarchy import fcluster
from matplotlib import rc, rcParams

rc('text', usetex=True)

lc_data = pd.read_csv("Lactobacillus_crispatus_axesCov.csv",sep=",")
li_data = pd.read_csv("Lactobacillus_iners_axesCov.csv",sep=",")
gv_data = pd.read_csv("Gardnerella_vaginalis_axesCov.csv",sep=",")

stacked_fig, stacked_axs = plt.subplots(3,1, figsize=(4,8), facecolor='w', edgecolor='k')
stacked_fig.subplots_adjust(left=0.2,right=0.95,bottom=0.05,top=0.95, wspace=0.4,hspace=0.4)
stacked_axs = stacked_axs.ravel()

stacked_axs[0].scatter(lc_data['x_pred'],lc_data['y_pred'],c='#ff0000',s=14)
stacked_axs[0].set_xlabel('Taxa Canonical Variable',fontsize=10)
stacked_axs[0].set_ylabel('VOG Canonical Variable',fontsize=10)
stacked_axs[0].set_title(r'\textit{L. crispatus} gene expr.',fontdict={'fontsize':12})
stacked_axs[0].text(-0.1,1.1,r"a",fontsize=18,ha='right',transform=stacked_axs[0].transAxes)

stacked_axs[1].scatter(li_data['x_pred'],li_data['y_pred'],c='#ff8c00',s=14)
stacked_axs[1].set_xlabel('Taxa Canonical Variable',fontsize=10)
stacked_axs[1].set_ylabel('VOG Canonical Variable',fontsize=10)
stacked_axs[1].set_title(r'\textit{L. iners} gene expr.',fontdict={'fontsize':12})
stacked_axs[1].text(-0.1,1.1,r"b",fontsize=18,ha='right',transform=stacked_axs[1].transAxes)

stacked_axs[2].scatter(gv_data['x_pred'],gv_data['y_pred'],c='#20b2aa',s=14)
stacked_axs[2].set_xlabel('Taxa Canonical Variable',fontsize=10)
stacked_axs[2].set_ylabel('VOG Canonical Variable',fontsize=10)
stacked_axs[2].set_title(r'\textit{Gardnerella} gene expr.',fontdict={'fontsize':12})
stacked_axs[2].text(-0.1,1.1,r"c",fontsize=18,ha='right',transform=stacked_axs[2].transAxes)

stacked_fig.savefig('sparseCCA_correlation_plot.pdf')