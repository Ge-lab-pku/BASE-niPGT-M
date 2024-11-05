#!/usr/bin/env python
# coding: utf-8

# In[21]:

import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import random
from math import log
from math import exp
pd.options.display.max_columns = 100
pd.set_option('display.max_rows', None)


# In[23]:


# If your phasing data is one file for each chromosome, use this cell.

dfs = []
# Loop through the files and read them into DataFrames
#for i in range(1, 24):
#    file_name = 'Phased_Data_Genome/301F11/Family11_chr' + str(i) + '_genome.phasedHap.txt'
#    df = pd.read_csv(file_name, sep='\t', index_col=0)
#    dfs.append(df)
#merged_df = pd.concat(dfs)
#merged_df.to_csv('Phased_Data_Genome/301F11.txt', sep='\t')


# In[25]:

family = sys.argv[1] 
infile = sys.argv[2]
fa = sys.argv[3]
mo = sys.argv[4]
fchr = sys.argv[5]

# If your phasing data is a whole file, use this cell.
data = pd.read_csv(infile,sep = "\t", index_col=None)


# In[27]:


def row_apply(row,name1,name2):
    return str(row[name1]) + "|" + str(row[name2])

def chrom(row):
    return str(row["CHROM"])

data["fGT_p"] = data.apply(row_apply,axis = 1,
                          name1 = "Phase_F1",name2 = "Phase_F2")

data["mGT_p"] =  data.apply(row_apply,axis = 1,
                          name1 = "Phase_M1",name2 = "Phase_M2")

data = data[["#Chrom","Position","Ref","Alt","fGT_p","mGT_p"]]
data.columns = ["CHROM","POS","REF","ALT",fa,mo] # Change it to the father/mother's name
data["QUAL"] = 50
data["FILTER"] ="PASS"
data["CHROM"] = data.apply(chrom,axis = 1)
# data.replace("chrX","chr23",inplace = True) If you want to change the chrX to chr23 just for simplification in code, use this.


data.to_csv("L"+family + "_chr" +str(fchr)+".hap.txt") # If you want to store the data, use this line.



