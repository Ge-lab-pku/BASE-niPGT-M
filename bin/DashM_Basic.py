#!/usr/bin/env python
# coding: utf-8


'''
This document defines some functions used to preprocess the raw data and phased data.

'''

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import random
from math import log
from math import exp

# Below are two auxiliary functions

def factor(n):
    
    # Compute the factor of n
    
    if n == 0:
        return 1
    n = int(n)
    ans = 1
    for i in range(1,n+1):
        ans *= i
    return ans
    
def C(n,p):
    
    # Compute the Combination Number of n and p
    
    n = int(n)
    p = int(p)
    if n < 0 or p < 0 or p > n:
        return 0
    else:
        return factor(n) / (factor(p) * factor(n-p))


# Below are functions used to extract information in the raw data and formulate a structured data.

def get_GT(sample,name):
    
    '''
    Get Genotype of each SNP from the raw data
    
    sample: a row of the raw data
    
    '''
    
    if "GT" not in sample["FORMAT"].split(":"):
        return None
    else:
        ind = sample["FORMAT"].split(":").index("GT")
        if sample[name].split(":")[ind] not in ["0/0","0/1","1/0","1/1","0|0","0|1","1|0","1|1"]:
            return None
        else:
            GT = sample[name].split(":")[ind]
            if GT in ["0/0","0|0"]:
                return "0/0"
            elif GT in ["0/1","1/0","0|1","1|0"]:
                return "0/1"
            elif GT in ["1|1","1/1"]:
                return "1/1"
            else:
                return None
    
def get_AD(sample,name):
    
    '''
    Get Allele Depth of each SNP from the raw data
    
    sample: a row of the raw data
    
    '''
    
    if "AD" not in sample["FORMAT"].split(":"):
        return None
    else:
        ind = sample["FORMAT"].split(":").index("AD")
        AD = eval("[" + sample[name].split(":")[ind] + "]")
        if len(AD) != 2:
            return None
        else:
            return AD
        
def get_DP(sample,name):
    
    '''
    Get Depth of each SNP from the raw data
    
    sample: a row of the raw data
    
    '''
    
    if "DP" not in sample["FORMAT"].split(":"):
        return None
    else:
        ind = sample["FORMAT"].split(":").index("DP")
        try:
            DP = eval(sample[name].split(":")[ind])
        except:
            return None
        if not isinstance(DP,(int,float,np.float64,np.int64)):
            return None
        else:
            return DP
        
def get_GQ(sample,name):
    
    '''
    Get Genotype Quality of each SNP from the raw data
    
    sample: a row of the raw data
    
    '''
    
    if "GQ" not in sample["FORMAT"].split(":"):
        return None
    else:
        ind = sample["FORMAT"].split(":").index("GQ")
        try:
            GQ = eval(sample[name].split(":")[ind])
        except:
            return None
        return GQ


def data_filter(data,fname,mname,sname,chrom,GQ_thres):
    
    '''
    Filter the raw data transformed from vcf file.
    
    '''
    
    data = data[["CHROM","POS","ID","REF","ALT","QUAL","FORMAT",fname,mname,sname]]
    
    if isinstance(chrom, int):
        data = data[data["CHROM"] == "chr" + str(chrom)]
        
    elif chrom == "X":
        data = data[data["CHROM"].isin(["chrX","chr23"])]
        
    elif chrom == "Y":
        data = data[data["CHROM"].isin(["chrX","chr23"])]
        
    else:
        raise ValueError("The chromosome input in the raw data is invalid ", chrom)
        
    # Extract the genotype of father, mother and sample
    data["fGT"] = data.apply(get_GT,axis = 1,name = fname)
    data["mGT"] = data.apply(get_GT,axis = 1,name = mname)
    data["sGT"] = data.apply(get_GT,axis = 1,name = sname)
    
    # Replacing all | with / because the data here s not phased.
    data.replace("0|0","0/0",inplace = True)
    data.replace("0|1","0/1",inplace = True)
    data.replace("1|0","0/1",inplace = True)
    data.replace("1|1","1/1",inplace = True)
    
    data = data[(data["fGT"].isin(["0/0","0/1","1/0","1/1"]) &
                 data["mGT"].isin(["0/0","0/1","1/0","1/1"]) &
                 data["sGT"].isin(["0/0","0/1","1/0","1/1"]))]
    
    # Extract ADs, DPs, GQs of the father, mother and sample
    data["fDP"] = data.apply(get_DP,axis = 1,name = fname)
    data["mDP"] = data.apply(get_DP,axis = 1,name = mname)
    data["sDP"] = data.apply(get_DP,axis = 1,name = sname)
    
    data["sAD"] = data.apply(get_AD,axis = 1,name = sname)
    
    data["fGQ"] = data.apply(get_GQ,axis = 1,name = fname)
    data["mGQ"] = data.apply(get_GQ,axis = 1,name = mname)
    data["sGQ"] = data.apply(get_GQ,axis = 1,name = sname)
    
    data = data[(data["fDP"] >= 5) & (data["mDP"] >= 5) &
               (data["fGQ"] >= GQ_thres) & (data["mGQ"] >= GQ_thres)]
    
    # We only use SNPs with DP >= 5 and GQ >= certain threshold usually 10).
    
    data = data[["CHROM","POS","ID","QUAL","fGT","mGT","sGT","sAD","sGQ"]]
    
    return data


def data_filter_phased(raw_data,phased_data,fname,mname,sname,chrom,GQ_thres):
    
    '''
    preprocessing the phased data which shows the linkage of SNPs.
    
    '''
    
    data = data_filter(raw_data,fname,mname,sname,chrom,GQ_thres)
    
    if isinstance(chrom, int):
        phased_data = phased_data[phased_data["CHROM"] == "chr" + str(chrom)]
        
    elif chrom == "X":
        phased_data = phased_data[phased_data["CHROM"].isin(["chrX","chr23"])]
        
    elif chrom == "Y":
        phased_data = phased_data[phased_data["CHROM"].isin(["chrX","chr23"])]
        
    else:
        raise ValueError("The chromosome input in the phased data is invalid ", chrom)
        
    phased_data = phased_data[["POS",fname,mname]]
    phased_data.columns = ["POS","fGT_p","mGT_p"] # '_p' means phased
    
    # Merge the raw data and the phased data and drop the NAs
    data = pd.merge(data,phased_data,left_on='POS', right_on='POS')
    data.dropna(axis = 0,subset = np.array(["sAD","fGT","mGT","sGT","fGT_p","mGT_p"]),inplace = True)
    
    def check(row):
        # A local function to check whether the raw data agrees with the phased data
        
        if row["fGT"] == "0/0" and not row["fGT_p"] in ("0|0","0/0"):
            return False
        elif row["fGT"] == "0/1" and not row["fGT_p"] in ("0|1","1|0","0/1","1/0"):
            return False
        elif row["fGT"] == "1/0" and not row["fGT_p"] in ("0|1","1|0","0/1","1/0"):
            return False
        elif row["fGT"] == "1/1" and not row["fGT_p"] in ("1|1","1/1"):
            return False
        
        if row["mGT"] == "0/0" and not row["mGT_p"] in ("0|0","0/0"):
            return False
        elif row["mGT"] == "0/1" and not row["mGT_p"] in ("0|1","1|0","0/1","1/0"):
            return False
        elif row["mGT"] == "1/0" and not row["mGT_p"] in ("0|1","1|0","0/1","1/0"):
            return False
        elif row["mGT"] == "1/1" and not row["mGT_p"] in ("1|1","1/1"):
            return False
        
        return True
    
    data = data.loc[data.apply(check, axis = 1),:]
    
    return data

def extract(data,chrom,fGT,mGT,sGT,sAD,filtered = True):
    
    '''
    Select all SNPs with specific genotype.
    FOr example, select all SNPs with fGT == 0/0 and mGT == 1/1.
    
    '''
    
    if not isinstance(data,pd.DataFrame):
        raise ValueError("Data input must be in the form of pd.DataFrame.")
        
    if isinstance(chrom, int):
        res = data[data["CHROM"] == "chr" + str(chrom)]
    elif chrom == "X":
        res = data[data["CHROM"] == "chrX"]
    elif chrom == "Y":
        res = data[data["CHROM"] == "chrY"]
    else:
        raise ValueError("The chromosome input in the phased data is invalid ")
    
    if fGT == "all":
        fGT = ["0/0","0/1","1/0","1/1"]
    elif isinstance(fGT, str):
        fGT = [fGT]
    
    if mGT == "all":
        mGT = ["0/0","0/1","1/0","1/1"]
    elif isinstance(mGT, str):
        mGT = [mGT]
        
    if sGT == "all":
        sGT = ["0/0","0/1","1/0","1/1"]
    elif isinstance(sGT, str):
        sGT = [sGT]
        
    def extract_single(row,fGT,mGT,sGT,sAD):
        
        '''
        A local function to determine whether an SNP satisfied the condition.
        
        '''
        
        if "sGT" not in row.index:
            return False
        
        if sAD == "all":
            if (row["fGT"] in fGT and row["mGT"] in mGT and
                row["sGT"] in sGT and isinstance(row["sAD"],list) and len(row["sAD"]) == 2):
                return True
            else:
                return False
            
        elif (row["fGT"] in fGT and row["mGT"] in mGT and
            row["sGT"] in sGT and row["sAD"] == sAD):
            return True
        else:
            return False
        
    res = res.loc[res.apply(extract_single,axis = 1,fGT = fGT, mGT = mGT,
                   sGT = sGT, sAD = sAD),:]
    
    return res



