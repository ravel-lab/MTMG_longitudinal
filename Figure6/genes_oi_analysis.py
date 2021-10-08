#!/usr/bin/python
import pandas as pd 
import sys
import numpy as np

genes_interest = pd.read_csv(sys.argv[1],sep=",")

mt_data = pd.read_csv(sys.argv[2],sep="\t")

mt_data_glcorr = mt_data[mt_data.columns[4:]].div(mt_data['Length'],axis=0)
mt_data_glcorr_rel = mt_data_glcorr.div(mt_data_glcorr.sum(axis=0),axis=1)
mt_data = pd.concat([mt_data[mt_data.columns[0]],mt_data_glcorr_rel],axis=1)

goi_data = pd.merge(left=genes_interest,right=mt_data,on="PC",how="left")

function_abund = goi_data.drop(["KEGG",'PC','Taxa','Length','Hotpep_GH'],axis=1)
function_abund_rel = function_abund.groupby("Type").sum()

function_abund_rel = function_abund_rel.T
function_abund_rel.to_csv("goi_function_relabun.csv",sep=",")

function_options = list(set(goi_data['Type'].values))

for function in function_options:

	function_data = goi_data[goi_data['Type']==function]

	function_data = function_data.drop(['KEGG','Type','PC','Length','Hotpep_GH','Type'],axis=1)

	function_data = function_data.groupby("Taxa").sum()

	function_data_rel =  function_data.div(function_data.sum(axis=0),axis=1)

	function_data_rel.to_csv("%s_taxa_contribs.csv" %(function),sep=",")