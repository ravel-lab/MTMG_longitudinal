#!/usr/bin/python3

import pandas as pd
import numpy as np
from skbio.stats.composition import multiplicative_replacement,clr
import sys

data=pd.read_csv(sys.argv[1],sep=",")

filename = sys.argv[1].split("_TPM")[0]

data = data.set_index("VOG")

data = data.T

data = data.replace(0.0,np.nan)

data_trim = data.loc[:, data.isnull().mean() < .5]

data_trim = data_trim.fillna(0.0)

data_trim.T.to_csv("%s_50cut_TPM.csv" %(filename),sep=",")

data_trim = data_trim.div(1000000)

mult_clr = clr(multiplicative_replacement(data_trim))

clr_out = pd.DataFrame(mult_clr,index=data_trim.index,columns=data_trim.columns)

clr_out.T.to_csv("%s_50cut_CLR.csv" %(filename),sep=",")
