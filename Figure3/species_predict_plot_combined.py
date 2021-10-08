#!/usr/bin/env python3

#importing packages to be used
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import pandas as pd
import numpy as np
import sys
import matplotlib
import seaborn as sns
import random
import matplotlib.patches as mpatches
pd.options.mode.chained_assignment = None

taxa_color_scheme = {'Lactobacillus_crispatus':'#ff0000','Gardnerella_vaginalis':'#20b2aa','g_Lactobacillus':'#eef06c',
                     'Lactobacillus_iners':'#ff8c00','Lactobacillus_gasseri':'#7fff00','Lactobacillus_jensenii':'#333333',
                     'Enterococcus_faecalis':'#bbffff','Raoultella_planticola':'#ffc2ff','g_Peptoniphilus':'#CCCC00',
                     'Sneathia_sanguinegens':'#c7aa8f','Atopobium_vaginae':'#0000cd','g_Atopobium':'#0000cd',
                     'Lacotbacillus_helveticus':'#00ccff','Mageeibacillus_indolicus':'#3cb371','g_Anaerococcus':'#87cefa','g_Gardnerella':'#20b2aa',
                     'Megasphaera_sp_type_1':'#0000ff','Streptococcus_agalactiae':'#ff738a','g_Megasphaera':'#008B45',
                     'Streptococcus_oralis':'#ff66ff','Prevotella_bivia':'#bfbfbf','Aerococcus_christensenii':'#bebebe',
                     'Anaerococcus_tetradius':'#87cefa','Gemella':'#daa520','Prevotella_genogroup_1':'#b0b0b0',
                     'Lactobacillus_vaginalis':'#ffffff','other':'#808080','Bifidobacterium_longum':'#c1ffc1',
                     'Bifidobacterium_breve':'#c1ffc1','Eggerthella':'#DE7710','Mycoplasma_hominis':'#10DE4E',
                     'Porphyromonas_bennonis':'#DE4310','Eubacterium_saphenum':'#8F10DE','Fusobacterium_nucleatum':'#CD853F',
                     'Fusobacterium_gonidiaformans':'#CD853F','Streptococcus_anginosus':'#ffc0cb','Peptostreptococcus_anaerobius':'#DEDB10',
                     'Arcanobacterium_phocae':'#8c10de','Bacteroides_uniformis':'#de1058','Ureaplasma_parvum':'#9999ff',
                     'Peptoniphilus_harei':'#CCCC00','Mobiluncus_mulieris':'#f08080','Megasphaera_sp._type_2':'#008B45',
                     'BVAB1':'#b31900','BVAB2':'#3bb16f','g_Escherichia.Shigella':'#12456b','Peptoniphilus_lacrimalis':'#d2b48c',
                     'Veillonella_montpellierensis':'#ff8c69','Prevotella_genogroup_3':'#b36200','Parvimonas_micra':'#cdcd00',
                     'Corynebacterium_accolens':'#ffff00','Finegoldia_magna':'#800080','Prevotella_genogroup_2':'#b0b0b0',
                     'Staphylococcus_epidermidis':'#ffffff','Prevotella_timonensis':'#f8fca9','g_Streptococcus':'#ffc0cb','g_Bifidobacterium':'#c1ffc1','g_Enterococcus':'#bbffff'
                     ,'g_Staphylococcus':'#800080','g_Finegoldia':'#800080','g_Prevotella':'#bfbfbf','g_Sneathia':'#d1e8eb','g_Aerococcus':'#e8e1ba'
                     ,'g_Leptotrichia':'#7ca386','g_Veillonella':'#499996','g_Dialister':'#997dbd','g_Corynebacterium_1':'#ffff00'
                     ,'g_Varibaculum':'#094717','g_Delftia':'#090c47','g_Corynebacterium_1':'#ffff00','Prevotella_amnii':'#ac8fc7','Sneathia_amnii':'#6963ff'}

def get_color(taxa):

    chars = '0123456789ABCDEF'

    if taxa in taxa_color_scheme:
        
        taxa_color = taxa_color_scheme[taxa]

    else:
        taxa_color_scheme[taxa] = '#'+''.join(random.sample(chars,6))
        taxa_color = taxa_color_scheme[taxa]

    return taxa_color


matplotlib.rc('font', serif='Helvetica Neue')
matplotlib.rc('axes',linewidth=2)
data = pd.read_csv("all_species_expr_predict.csv",sep=",")

species_list = list(data['Species'].unique())
stacked_fig, stacked_axs = plt.subplots(1,1, figsize=(8,8), facecolor='w', edgecolor='k')
stacked_fig.subplots_adjust(left=0.10,right=0.9,bottom=0.1,top=0.9, wspace=0.1,hspace=0.1)

for species in species_list:

    taxa_color = get_color(species)

    species_data = data[data['Species']==species]
    species_data = species_data[species_data['LogPriorExpr'].notna()]
    species_data = species_data[species_data['LogChange'].notna()]
    
    stacked_axs.scatter(x=species_data['LogPriorExpr'],y=species_data['LogChange'],c="#adadad",s=4,marker="o")

stacked_axs.set_xlim([-8,8])
stacked_axs.set_ylim([-8,8])
stacked_axs.set_xlabel('Log 10 Prior Rel. Expression',fontsize=20)
stacked_axs.set_ylabel('Log 10 Change in Rel. Abundandce',fontsize=20)
stacked_axs.set_xticks([-8,-6,-4,-2,0,2,4,6,8])
stacked_axs.set_xticklabels(['-8','-6','-4','-2','0','2','4','6','8'],fontsize=16)
stacked_axs.set_yticks([-8,-6,-4,-2,0,2,4,6,8])
stacked_axs.set_yticklabels(['-8','-6','-4','-2','0','2','4','6','8'],fontsize=16)

#stacked_axs.axhline(y=0,color="k",linewidth=2)
x_vals = np.array(stacked_axs.get_xlim())
intercept = -0.51429
slope = 0.6253872
y_vals = intercept + slope * x_vals
stacked_axs.plot(x_vals, y_vals, '--',lw=3,c="k")

stacked_fig.savefig('all_predict_plot.pdf')