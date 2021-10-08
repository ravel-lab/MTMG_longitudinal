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
matplotlib.rc('font', serif='Helvetica Neue')
matplotlib.rc('axes',linewidth=2)

taxa_color_scheme = {'Lactobacillus_crispatus':'#ff0000','Gardnerella_vaginalis':'#20b2aa','g_Lactobacillus':'#eef06c',
                     'Lactobacillus_iners':'#ff8c00','Lactobacillus_gasseri':'#7fff00','Lactobacillus_jensenii':'#333333',
                     'Enterococcus_faecalis':'#bbffff','Raoultella_planticola':'#ffc2ff','g_Peptoniphilus':'#CCCC00',
                     'Sneathia_sanguinegens':'#c7aa8f','Atopobium_vaginae':'#0000cd','g_Atopobium':'#0000cd',
                     'Lacotbacillus_helveticus':'#00ccff','Mageeibacillus_indolicus':'#3cb371','g_Anaerococcus':'#87cefa','g_Gardnerella':'#20b2aa',
                     'Megasphaera_genomosp.':'#0000ff','Streptococcus_agalactiae':'#ff738a','Megasphaera_genomosp.':'#008B45',
                     'Streptococcus_oralis':'#ff66ff','Prevotella_bivia':'#bfbfbf','Aerococcus_christensenii':'#bebebe',
                     'Anaerococcus_tetradius':'#87cefa','Gemella':'#daa520','Prevotella_genogroup_1':'#b0b0b0','Prevotella_disiens':'#9bbade',
                     'Lactobacillus_vaginalis':'#ffffff','other':'#808080','Bifidobacterium_longum':'#c1ffc1','Prevotella_buccalis':'#debdff',
                     'Bifidobacterium_breve':'#c1ffc1','Eggerthella':'#DE7710','Mycoplasma_hominis':'#10DE4E',
                     'Porphyromonas_uenonis':'#DE4310','Eubacterium_saphenum':'#8F10DE','Fusobacterium_nucleatum':'#CD853F',
                     'Fusobacterium_gonidiaformans':'#CD853F','Streptococcus_anginosus':'#ffc0cb','Peptostreptococcus_anaerobius':'#DEDB10',
                     'Arcanobacterium_phocae':'#8c10de','Bacteroides_uniformis':'#de1058','Ureaplasma_parvum':'#9999ff','Prevotella_corporis':'#ad8bab',
                     'Peptoniphilus_harei':'#CCCC00','Mobiluncus_mulieris':'#f08080','Megasphaera_sp._type_2':'#008B45',
                     'BVAB1':'#b31900','BVAB2':'#3bb16f','g_Escherichia.Shigella':'#12456b','Peptoniphilus_lacrimalis':'#803849',
                     'Veillonella_montpellierensis':'#ff8c69','Prevotella_genogroup_3':'#b36200','Parvimonas_micra':'#cdcd00',
                     'Corynebacterium_accolens':'#ffff00','Finegoldia_magna':'#800080','Prevotella_genogroup_2':'#b0b0b0','Propionibacterium_sp.':'#5f7362',
                     'Staphylococcus_epidermidis':'#ffffff','Prevotella_timonensis':'#f8fca9','g_Streptococcus':'#ffc0cb','g_Bifidobacterium':'#c1ffc1','g_Enterococcus':'#bbffff'
                     ,'g_Staphylococcus':'#800080','g_Finegoldia':'#800080','Prevotella_sp.':'#ffd9b3','g_Sneathia':'#d1e8eb','g_Aerococcus':'#e8e1ba'
                     ,'g_Leptotrichia':'#7ca386','g_Veillonella':'#499996','g_Dialister':'#997dbd','g_Corynebacterium_1':'#ffff00'
                     ,'g_Varibaculum':'#094717','g_Delftia':'#090c47','g_Corynebacterium_1':'#ffff00','Prevotella_amnii':'#ac8fc7','Sneathia_amnii':'#6963ff'}


species_names = ["$\it{S. anginosus}$","$\it{F. magna}$","$\it{M. mulieris}$","$\it{G. vaginalis}$","$\it{A. vaginae}$","$\it{P. corporis}$","$\it{P. bivia}$","$\it{P. disiens}$","$\it{L. crispatus}$",
                 "$\it{B. breve}$","$\it{Propionibacterium}$","$\it{L. jensenii}$","$\it{Ca. L. vaginae}$","$\it{M. indolicus}$","$\it{L. gasseri}$","$\it{P. lacrimalis}$","$\it{L. iners}$","$\it{Megasphaera}$",
                "$\it{Prevotella}$","$\it{P. buccalis}$","$\it{P. timonensis}$","$\it{P. amnii}$","$\it{P. uenonis}$","$\it{S. sanguinegens}$","$\it{S. amnii}$"]
def get_color(taxa):

    chars = '0123456789ABCDEF'

    if taxa in taxa_color_scheme:
        
        taxa_color = taxa_color_scheme[taxa]

    else:
        taxa_color_scheme[taxa] = '#'+''.join(random.sample(chars,6))
        taxa_color = taxa_color_scheme[taxa]

    return taxa_color

data = pd.read_csv("species_coef.csv",sep=",")
data['UL'] = data.apply(lambda row: row['Xint'] + row['Std.Error'], axis=1)
data['LL'] = data.apply(lambda row: row['Xint'] - row['Std.Error'], axis=1)

species_list = list(data['Coef'].unique())
stacked_fig, stacked_axs = plt.subplots(1,1, figsize=(10,6), facecolor='w', edgecolor='k')
stacked_fig.subplots_adjust(left=0.1,right=0.9,bottom=0.3,top=0.9, wspace=0.1,hspace=0.1)

for species in species_list:

    species_data = data[data['Coef']==species]
    taxa_color = get_color(species)

    stacked_axs.plot([species_data['x_loc'].values,species_data['x_loc'].values],[species_data['UL'].values,species_data['LL'].values],c='#adadad',lw=2)
    stacked_axs.scatter(y=species_data['Xint'],x=species_data['x_loc'],c=taxa_color,s=35,zorder=5)

stacked_axs.set_ylim(-0.55,3.15)
stacked_axs.set_xlim(0.5,25.5)
stacked_axs.set_yticks([-.5,0,0.5,1,1.5,2,2.5,3])
stacked_axs.set_yticklabels([-.5,0,0.5,1,1.5,2,2.5,3],fontsize=18)
stacked_axs.set_xticks([1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25])
stacked_axs.set_xticklabels(species_names,fontsize=17,rotation=50,ha='right')
stacked_axs.set_ylabel("Log 10 relative expression\n to maintain abundance",fontsize=20)
stacked_axs.yaxis.grid(color='gray')
stacked_axs.set_axisbelow(True)

stacked_fig.savefig("Species_pred_inters.pdf")
