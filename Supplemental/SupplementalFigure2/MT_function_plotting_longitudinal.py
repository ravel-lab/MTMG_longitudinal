#!/usr/bin/python3

import sys
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
import matplotlib
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
                     'Arcanobacterium_phocae':'#8c10de','Bacteroides_uniformis':'#de1058','Ureaplasma_parvum':'#9999ff','Prevotella_corporis':'#d63c87',
                     'Peptoniphilus_harei':'#CCCC00','Mobiluncus_mulieris':'#f08080','Megasphaera_sp._type_2':'#008B45',
                     'BVAB1':'#b31900','BVAB2':'#3bb16f','g_Escherichia.Shigella':'#12456b','Peptoniphilus_lacrimalis':'#803849',
                     'Veillonella_montpellierensis':'#ff8c69','Prevotella_genogroup_3':'#b36200','Parvimonas_micra':'#cdcd00',
                     'Corynebacterium_accolens':'#ffff00','Finegoldia_magna':'#800080','Prevotella_genogroup_2':'#b0b0b0','Propionibacterium_sp.':'#5f7362',
                     'Staphylococcus_epidermidis':'#ffffff','Prevotella_timonensis':'#f8fca9','g_Streptococcus':'#ffc0cb','g_Bifidobacterium':'#c1ffc1','g_Enterococcus':'#bbffff'
                     ,'g_Staphylococcus':'#800080','g_Finegoldia':'#800080','Prevotella_sp.':'#ffd9b3','g_Sneathia':'#d1e8eb','g_Aerococcus':'#e8e1ba'
                     ,'g_Leptotrichia':'#7ca386','g_Veillonella':'#499996','g_Dialister':'#997dbd','g_Corynebacterium_1':'#ffff00'
                     ,'g_Varibaculum':'#094717','g_Delftia':'#090c47','g_Corynebacterium_1':'#ffff00','Prevotella_amnii':'#ac8fc7','Sneathia_amnii':'#6963ff'}


species_names_key = {'Lactobacillus_crispatus':"$\it{L. crispatus}$",'Lactobacillus_iners':"$\it{L. iners}$",'Lactobacillus_jensenii':"$\it{L. jensenii}$",'Lactobacillus_gasseri':"$\it{L. gasseri}$",'Gardnerella_vaginalis':"$\it{G. vaginalis}$",
                     'BVAB1':"$\it{Ca. L. vaginae}$",'Atopobium_vaginae':"$\it{A. vaginae}$",'Mageeibacillus_indolicus':"$\it{M. indolicus}$",'Prevotella_amnii':"$\it{P. amnii}$",
                     'Prevotella_bivia':"$\it{P. bivia}$",'Prevotella_buccalis':"$\it{P. buccalis}$",'Prevotella_corporis':"$\it{P. corporis}$",'Prevotella_disiens':"$\it{P. disiens}$",'Prevotella_sp.':"$\it{Prevotella}$",
                     'Prevotella_timonensis':"$\it{P. timonensis}$",'Sneathia_amnii':"$\it{S. amnii}$",'Sneathia_sanguinegens':"$\it{S. sanguinegens}$",'Finegoldia_magna':"$\it{F. magna}$",
                     'Megasphaera_genomosp.':"$\it{Megasphaera}$",'Mobiluncus_mulieris':"$\it{M. mulieris}$",'Peptoniphilus_lacrimalis':"$\it{P. lacrimalis}$",'Porphyromonas_uenonis':"$\it{P. uenonis}$",
                     'Propionibacterium_sp.':"$\it{Propionibacterium}$",'Streptococcus_anginosus':"$\it{S. anginosus}$",'Bifidobacterium_breve':"$\it{B. breve}$",'other':'Other'}

species_names = list(species_names_key.keys())
def get_color(taxa):

    chars = '0123456789ABCDEF'

    if taxa in taxa_color_scheme:
        
        taxa_color = taxa_color_scheme[taxa]

    else:
        taxa_color_scheme[taxa] = '#'+''.join(random.sample(chars,6))
        taxa_color = taxa_color_scheme[taxa]

    return taxa_color

CST_color_dict = {"I":'#d44646',"II":'#7ac472',"III":'#d48326',"IV":'#262fd4',"V":'#525254'}

function_abund_data = pd.read_csv(sys.argv[1],sep=",")
#function_abund_data = function_abund_data[function_abund_data['CST'] != "II"]
#function_abund_data = function_abund_data[function_abund_data['CST'] != "V"]
function_abund_data = function_abund_data.drop(['Cytolysin','Nactylgalactosaminidase'],axis=1)
functions = list(function_abund_data.columns[5:])
subjects = list(set(function_abund_data['SID']))

count=0 

stacked_fig, stacked_axs = plt.subplots(2,4, figsize=(16,8), facecolor='w', edgecolor='k')
stacked_fig.subplots_adjust(left=0.1,right=0.9,bottom=0.1,top=0.9,hspace = 0.4, wspace=0.4)
stacked_axs = stacked_axs.ravel()

for function in functions:

    for subject in subjects:

        subject_data = function_abund_data[function_abund_data['SID']==subject]

        stacked_axs[count].plot(subject_data['Timepoint'],subject_data[function],linewidth=1,c='#000000',linestyle='--')
        for i in range(len(subject_data.index)):
    
            stacked_axs[count].scatter(subject_data.iloc[i]['Timepoint'],subject_data.iloc[i][function],c=subject_data.iloc[i]['Color'])

    stacked_axs[count].set_yscale("log")
    if function == 'Glycogen_debranch':
        stacked_axs[count].set_ylim([0.0001,0.01])
        stacked_axs[count].set_yticks([0.0001,0.001,0.01])
        stacked_axs[count].set_yticklabels(['100','1000','10,000'],fontsize=12)

    elif function == 'Sulfatase':
        stacked_axs[count].set_ylim([0.000001,0.01])
        stacked_axs[count].set_yticks([0.000001,0.00001,0.0001,0.001,0.01])
        stacked_axs[count].set_yticklabels(['1','10','100','1000','10,000'],fontsize=12)
    else:
        stacked_axs[count].set_ylim([0.000000001,0.01])
        stacked_axs[count].set_yticks([0.00000001,0.0000001,0.000001,0.00001,0.0001,0.001,0.01])
        stacked_axs[count].set_yticklabels(['0.01','0.1','1','10','100','1000','10,000'],fontsize=12)

    stacked_axs[count].set_xticks([1,2,3,4,5],minor=False)
    stacked_axs[count].set_xticklabels(['1','2','3','4','5'],fontsize=12,minor=False)
    stacked_axs[count].set_xlabel('Timepoint',fontsize=12)
    stacked_axs[count].set_title(function,fontsize=14)
    stacked_axs[count].set_ylabel('Expression level (TPM)',fontsize=12)
    count+=1



for blank_plot in [7]:
    stacked_fig.delaxes(stacked_axs.flatten()[blank_plot])

patch_list = []
for CST in CST_color_dict:
    data_key = mpatches.Patch(color=CST_color_dict[CST],label=CST)
    patch_list.append(data_key)    
stacked_fig.legend(handles=patch_list,loc=[0.8,0.2] ,ncol=1,fontsize='large',title='CST')

stacked_fig.savefig("allFunctions_longitudinal.pdf")



    