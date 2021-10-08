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
from numpy.polynomial.polynomial import polyfit
from numpy import inf

pd.options.mode.chained_assignment = None

#taxa color dictionary
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

species_color_set = {'Lcrispatus':'#ff0000','Liners':'#ff8c00','Lgasseri':'#7fff00','Ljensenii':'#FFFF00','NonLacto':'#0000cd'}


species_list = ['Lactobacillus_iners','Lactobacillus_crispatus','Gardnerella_vaginalis','Lactobacillus_gasseri','Sneathia_amnii','Lactobacillus_jensenii','Prevotella_amnii','Prevotella_timonensis','BVAB1',
'Prevotella_bivia','Sneathia_sanguinegens','Prevotella_buccalis','Bifidobacterium_breve','Megasphaera_genomosp.','Prevotella_sp.','Propionibacterium_sp.',
'Mageeibacillus_indolicus','Atopobium_vaginae','Porphyromonas_uenonis','Mobiluncus_mulieris','Peptoniphilus_lacrimalis','Prevotella_disiens','Prevotella_corporis','Streptococcus_anginosus','Finegoldia_magna']

species_text = [r'\textit{Lactobacillus iners}',r'\textit{Lactobacillus crispatus}',r'\textit{Gardnerella vaginalis}',r'\textit{Lactobacillus gasseri}',r'\textit{Sneathia amnii}',
                r'\textit{Lactobacillus jensenii}',r'\textit{Prevotella amnii}',r'\textit{Prevotella timonensis}',r'Ca. \textit{Lachnocurva vaginae}',r'\textit{Prevotella bivia}',
                r'\textit{Sneathia sanguinegens}',r'\textit{Prevotella buccalis}',r'\textit{Bifidobacterium breve}',r'\textit{Megasphaera genomosp.}',r'\textit{Prevotella sp.}',
                r'\textit{Propionibacterium sp.}',r'\textit{Mageeibacillus indolicus}',r'\textit{Atopobium vaginae}',r'\textit{Porphyromonas uenonis}',r'\textit{Mobiluncus mulieris}',
                r'\textit{Peptoniphilus lacrimalis}',r'\textit{Prevotella disiens}',r'\textit{Prevotella corporis}',r'\textit{Streptococcus anginosus}',r'\textit{Finegoldia magna}']

#species_list = ['Lactobacillus_iners']
#setting font and axes globally
matplotlib.rc('font', serif='Helvetica Neue')
matplotlib.rc('axes',linewidth=2)
matplotlib.rc('text', usetex=True)

#reading in 16S data e.g. hmp_u01_combined_061318_mf.xlsx
MG_data = pd.read_csv(sys.argv[1],sep="\t")
MG_data = MG_data[MG_data['MT_avail'] == 1]

#dropping rows without 16S data
MT_data = pd.read_csv(sys.argv[2],sep="\t")
MT_data = MT_data[MT_data['MG_avail'] == 1]

relative_expression = pd.DataFrame(columns=["UID",'SID','Count','Species','Ratio','LogAbund'])

plot_count = 0 
stacked_fig, stacked_axs = plt.subplots(5,5,figsize=(20,16), facecolor='w', edgecolor='k')
stacked_fig.subplots_adjust(left=0.05,right=0.95,bottom=0.05,top=0.95, wspace=0.2,hspace=0.2)
stacked_axs = stacked_axs.ravel()


for species in species_list:
    print(species)

    species_MG_data = MG_data[['UID','SID','Stability','Community_type','Count','%s' %(species)]]
    species_MG_data.rename(columns={'%s'%(species):"MG"},inplace=True)
    species_MT_data = MT_data[['UID','%s' %(species)]]
    species_MT_data.rename(columns={'%s'%(species):"MT"},inplace=True)
    species_MG_data = species_MG_data.replace(0,np.nan)
    species_MG_data = species_MG_data.fillna(0.00000001)

    species_combo = pd.merge(left=species_MG_data,right=species_MT_data,on="UID",how="inner")
    species_combo['Species'] = species
    species_combo['Ratio'] = np.log10(species_combo['MT']/species_combo['MG'])
    species_combo['LogAbund'] = np.log10(species_combo['MG'])
    species_combo = species_combo.drop(['MG','MT'],axis=1)
    relative_expression = relative_expression.append(species_combo)

    species_combo['color'] = species_combo['Community_type'].map(species_color_set)
    
    #species_combo = species_combo[species_combo['Community_type'] == 'NonLacto']
    
    #species_combo_stable = species_combo[species_combo['Stability'] == "Stable"]
    #species_combo_unstable = species_combo[species_combo['Stability'] == "Unstable"]
    species_combo['Ratio'] = species_combo['Ratio'].replace(-inf,-4)
    species_combo['Ratio'] = species_combo['Ratio'].replace(0,-4)
    b, m = polyfit(species_combo['LogAbund'], species_combo['Ratio'], 1)
    x_int = -b/m
    stacked_axs[plot_count].scatter(x=species_combo['LogAbund'],y=species_combo['Ratio'],c=taxa_color_scheme[species],s=10,marker="o")
    stacked_axs[plot_count].plot(species_combo['LogAbund'], b + m * species_combo['LogAbund'], '-',c='k',lw=2)
    #stacked_axs.scatter(x=species_combo_stable['LogAbund'],y=species_combo_stable['Ratio'],c='#4287f5',s=6,marker="o")
    #stacked_axs.scatter(x=species_combo_unstable['LogAbund'],y=species_combo_unstable['Ratio'],c='#56d177',s=6,marker="o")

    #stacked_axs.set_title('%s' %(species))
    stacked_axs[plot_count].set_xlim([-7,0])
    stacked_axs[plot_count].set_xticks([-7,-6,-5,-4,-3,-2,-1,0])
    stacked_axs[plot_count].set_yticks([-6,-4,-2,0,2,4,6])
    stacked_axs[plot_count].set_xticklabels([-7,-6,-5,-4,-3,-2,-1,0],fontsize=8)
    stacked_axs[plot_count].set_yticklabels([-6,-4,-2,0,2,4,6],fontsize=8)
    stacked_axs[plot_count].set_ylim([-6,6])
    stacked_axs[plot_count].set_xlabel('Log 10 Rel. Abundance',fontsize=10)
    stacked_axs[plot_count].set_ylabel('Log 10 Expr/Abund',fontsize=10)
    stacked_axs[plot_count].axhline(y=0,color="#adadad",linewidth=1,linestyle='-')
    stacked_axs[plot_count].axvline(x=x_int,color="#adadad",linewidth=1)
    stacked_axs[plot_count].text(.05,0.05,species_text[plot_count],fontsize=12,ha='left',transform=stacked_axs[plot_count].transAxes)

    plot_count+=1

stacked_fig.savefig('All_taxa_MGMT_plot.pdf')
plt.close()

#relative_expression.to_csv("species_relative_exp.csv",sep=",",index=None)