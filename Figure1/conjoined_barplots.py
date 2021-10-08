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
from skbio.diversity.alpha import shannon


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
matplotlib.rc('axes',linewidth=1)

#reading the data in 
mg_data = pd.read_csv(sys.argv[1],sep="\t")
mt_data = pd.read_csv(sys.argv[2],sep="\t")

xtick_list = list(mg_data['loc'])
xtick_list.extend(x for x in list(mt_data['loc']) if x not in xtick_list)

keep_columns = ['Lactobacillus_crispatus','Lactobacillus_iners','Gardnerella_vaginalis','Lactobacillus_gasseri','Lactobacillus_jensenii','BVAB1','Atopobium_vaginae','Sneathia_sanguinegens','Sneathia_amnii','Prevotella_timonensis','Prevotella_bivia','Mageeibacillus_indolicus','Bifidobacterium_breve','Streptococcus_agalactiae','Streptococcus_anginosus','Prevotella_amnii']

bar_width = 1

plot_count=0


legend_entries = {}

stacked_fig, stacked_axs = plt.subplots(2,1, figsize=(14,10), facecolor='w', edgecolor='k')
stacked_fig.subplots_adjust(left=0.1,right=0.8,bottom=0.1,top=0.9, wspace=0.05,hspace=0.1)
stacked_axs = stacked_axs.ravel()

mg_data_micro = mg_data[keep_columns]

mg_data_micro['other'] = mg_data_micro.apply(lambda row: 1.0 - row.sum(), axis=1)

mg_data_sub = pd.concat([mg_data[mg_data.columns[0:9]],mg_data_micro],axis=1,sort=False)

mg_data_sub['bottom_count'] = pd.Series([0.0 for x in range(len(mg_data_sub.index))], index=mg_data_sub.index)

for taxa in range(9,len(keep_columns)+10):

    taxa = mg_data_sub.columns[taxa]

    print(taxa)

    taxa_color = get_color(taxa)

    if taxa not in legend_entries:

        legend_entries[taxa] = taxa_color_scheme[taxa]
    
    stacked_axs[0].bar(mg_data_sub['loc'],mg_data_sub[taxa],width=bar_width,bottom=mg_data_sub['bottom_count'],color=taxa_color,clip_on=False)
    
    mg_data_sub['bottom_count'] = mg_data_sub['bottom_count'] + mg_data_sub[taxa]

stacked_axs[0].set_ylim([0.0,1.0])
stacked_axs[0].tick_params(width=2)
stacked_axs[0].set_xticks([], minor=False)
stacked_axs[0].set_yticks([0,0.2,0.4,0.6,0.8,1.0],minor=False)
stacked_axs[0].set_yticklabels(['0','0.2','0.4','0.6','0.8','1.0'],minor=False,fontsize=20)
stacked_axs[0].set_xticklabels([],minor=False)
stacked_axs[0].set_xlabel('',fontsize=24)
stacked_axs[0].set_title('Metagenomic data',fontsize=20,loc='left')
stacked_axs[0].set_ylabel('Relative abundance',fontsize=20)

stacked_axs[0].yaxis.grid(color='gray')
stacked_axs[0].set_axisbelow(True)

#stacked_axs[0].vlines(line_loc,linewidth=1.5,ymin=-0.1,ymax=1,clip_on=False)

####MT data bars
mt_data_micro = mt_data[keep_columns]

mt_data_micro['other'] = mt_data_micro.apply(lambda row: 1.0 - row.sum(), axis=1)

mt_data_sub = pd.concat([mt_data[mt_data.columns[0:8]],mt_data_micro],axis=1,sort=False)

mt_data_sub['bottom_count'] = pd.Series([0.0 for x in range(len(mt_data_sub.index))], index=mt_data_sub.index)

for taxa in range(8,len(keep_columns)+9):

    taxa = mt_data_sub.columns[taxa]

    print(taxa)

    taxa_color = get_color(taxa)

    if taxa not in legend_entries:

        legend_entries[taxa] = taxa_color_scheme[taxa]
    
    stacked_axs[1].bar(mt_data_sub['loc'],mt_data_sub[taxa],width=bar_width,bottom=mt_data_sub['bottom_count'],color=taxa_color,clip_on=False)
    
    mt_data_sub['bottom_count'] = mt_data_sub['bottom_count'] + mt_data_sub[taxa]

stacked_axs[1].set_ylim([0.0,1.0])
stacked_axs[1].tick_params(width=2)
stacked_axs[1].set_xticks([], minor=False)
stacked_axs[1].set_yticks([0,0.2,0.4,0.6,0.8,1.0],minor=False)
stacked_axs[1].set_yticklabels(['0','0.2','0.4','0.6','0.8','1.0'],minor=False,fontsize=20)
stacked_axs[1].set_xticklabels([],minor=False)
stacked_axs[1].set_xlabel('Metatranscriptomic data',fontsize=20)
stacked_axs[1].xaxis.set_label_coords(.12, -0.05)
stacked_axs[1].set_ylabel('Relative expression',fontsize=20)
stacked_axs[1].xaxis.tick_top()

stacked_axs[1].yaxis.grid(color='gray')
stacked_axs[1].set_axisbelow(True)

taxa_text_legend = {'Lactobacillus_iners':'L. iners','Lactobacillus_crispatus':'L. crispatus','Gardnerella_vaginalis':'G. vaginalis','Lactobacillus_jensenii':'L. jensenii',
                    'Atopobium_vaginae':"A. vaginae","Lactobacillus_gasseri":'L. gasseri','BVAB1':'BVAB1','Streptococcus_agalactiae':'S. agalactiae','Sneathia_sanguinegens':'S. sanguinegens',
                    'Sneathia_amnii':'S. amnii','Prevotella_timonensis':'P. timonensis', 'Prevotella_bivia':'P. bivia','Mageeibacillus_indolicus':'M. indolicus','Streptococcus_anginosus':'S. aginosus',
                    'Prevotella_amnii':'P. amnii','Bifidobacterium_breve':'B. breve','other':'other'}

patch_list = []
for taxa in legend_entries:
    taxa_text = taxa_text_legend[taxa]

    if taxa_text == "BVAB1" or taxa_text == "other":
        data_key = mpatches.Patch(color=legend_entries[taxa],label="%s" %(taxa_text))
    
    else:    
        data_key = mpatches.Patch(color=legend_entries[taxa],label="$\it{%s}$" %(taxa_text))
    
    patch_list.append(data_key)

stacked_fig.legend(handles=patch_list,loc=7 ,ncol=1,fontsize=14)
stacked_fig.savefig('mg_mt_barplots.pdf')









