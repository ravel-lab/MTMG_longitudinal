#!/usr/bin/env python3

#importing packages to be used
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import sys
import matplotlib.colors as clr
import matplotlib.patches as mpatches

#file with PCA coords of the samples
pca_coords = pd.read_csv(sys.argv[1],sep=",")

pca_coords_simp = pca_coords.groupby('Replicate').mean()
color_dict = {'chl':'#808080','Lcr6':'#eb4034','Lcr12':'#ff8178','Lin6':'#ffb300','Lin12':'#fad682'}
#####make one figure of average PCA coordinates

stacked_fig, stacked_axs = plt.subplots(1,1, figsize=(6,8), facecolor='w', edgecolor='k')
stacked_fig.subplots_adjust(left=0.10,right=0.9,bottom=0.1,top=0.9, wspace=0.1,hspace=0.1)

for category in color_dict:

	cat_data = pca_coords_simp[pca_coords_simp['Trajectory']==category]

    scatter = stacked_axs.scatter(x=cat_data['PC1'],y=cat_data['PC2'],c=color_dict[category],s=16)


stacked_axs.yaxis.grid(color='#ededed',lw=0.75)
stacked_axs.xaxis.grid(color='#ededed',lw=0.75)


stacked_axs.set_axisbelow(True)

stacked_fig.colorbar(scatter,shrink=0.5,orientation='horizontal')

stacked_fig.savefig("PSTIN_rnaSeq_PCA.pdf")
