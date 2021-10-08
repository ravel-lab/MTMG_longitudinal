#!/usr/bin/python

import pandas as pd
import sys

data = pd.read_csv(sys.argv[1],sep="\t")
species_reads = pd.read_csv(sys.argv[2],sep=",",index_col=0)

focal_species = sys.argv[3]

vog_key = pd.read_csv("~/bioinformatics_tools/virgo/0.VOG.VIRGO.tbl.txt",sep="\t",header=None)
vog_key.columns = ['VOG','PC']

species_data = data[data['Taxa']==focal_species]

species_data = species_data.drop(['Length','KEGG','Taxa'],axis=1)
species_data = pd.merge(left=vog_key,right=species_data,on="PC",how='right')
species_data = species_data.drop(['PC'],axis=1)
species_data = species_data.groupby("VOG").sum()

#species_data = species_data.set_index("PC")
keep_samples = list(species_reads[species_reads[focal_species]>=500000].index)

species_data = species_data[keep_samples]

species_data.to_csv("%s_TPM.txt" %(focal_species))
