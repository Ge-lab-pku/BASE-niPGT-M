# This ducoment defines functions to label the data after preprocessing.

# version 3.0: no truncation

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats

def label(data, chrom, pos, window_size, ADO_thres, sex = None, output = True):
    
    '''
    
    ADO_thres: The threshoold used to determine which SNPs are classified as 'paternal's DNA only' and 'maternal's DNA only'.
    
    '''
    
    if window_size == "Whole_Genome":
        truncated_data = data
    else:
        truncated_data = data[(data["POS"] <= pos + window_size // 2)&((data["POS"] >= pos - window_size // 2))]
    
    truncated_data = truncated_data[~((truncated_data["fGT_p"] == "0|0")&(truncated_data["mGT_p"] == "0|0")|
                                   (truncated_data["fGT_p"] == "1|1")&(truncated_data["mGT_p"] == "1|1"))]
    # We only use the SNPs with heterogeneous genotype or with different homogeneous genotypes.
    if output:
        print("In the region we require, we detect %d SNPs with either parent's genotype being heterogeneous or their genotypes are homogeneous but different.\n" % truncated_data.shape[0])
    
    total_reads = [truncated_data.iloc[x]["sAD"][0] + truncated_data.iloc[x]["sAD"][1] for x in range(truncated_data.shape[0])]
    
    # Label F or M or FM or F1 or F2
    label = []
    for i in range(truncated_data.shape[0]):
        # Pop A Row
        row = truncated_data.iloc[i]
        sAD = row["sAD"]
        if isinstance(sAD,str):
            sAD = eval(sAD) 
            # For some samples, the format of sAD is string, for others, they are lists.
            
        if row["fGT"] == "0/1" and row["mGT"] == "0/1":
            label.append("Unknown")
            
        elif row["fGT"] == "0/0" and row["mGT"] == "1/1":
            if sAD[0] > ADO_thres and sAD[1] == 0:
                label.append("F")
            elif sAD[1] > ADO_thres and sAD[0] == 0:
                label.append("M")
            elif sAD[1] > ADO_thres // 2 and sAD[0] > ADO_thres // 2:
                label.append("FM")
            else:
                label.append("Unknown")
                
        elif row["fGT"] == "1/1" and row["mGT"] == "0/0":
            if sAD[0] > ADO_thres and sAD[1] == 0:
                label.append("M")
            elif sAD[1] > ADO_thres and sAD[0] == 0:
                label.append("F")
            elif sAD[1] > ADO_thres // 2 and sAD[0] > ADO_thres // 2:
                label.append("FM")
            else:
                label.append("Unknown")
                
        elif row["fGT"] == "0/0" and row["mGT"] == "0/1":
            if sAD[1] > ADO_thres and sAD[0] == 0:
                label.append("M")
            else:
                label.append("Unknown")
                
        elif row["fGT"] == "1/1" and row["mGT"] == "0/1":
            if sAD[0] > ADO_thres and sAD[1] == 0:
                label.append("M")
            else:
                label.append("Unknown")
                
        elif row["fGT"] == "0/1" and row["mGT"] == "0/0":
            if sAD[1] > ADO_thres and sAD[0] == 0:
                if row["fGT_p"] == "0|1":
                    label.append("F2")
                elif row["fGT_p"] == "1|0":
                    label.append("F1")
            elif sAD[1] > ADO_thres // 2 and sAD[0] > ADO_thres // 2:# We  assume there is no paternal contamination.
                if row["fGT_p"] == "0|1":
                    label.append("FM")
                elif row["fGT_p"] == "1|0":
                    label.append("FM")
            else:
                label.append("Unknown")
                
        elif row["fGT"] == "0/1" and row["mGT"] == "1/1":
            if sAD[0] > ADO_thres and sAD[1] == 0:
                if row["fGT_p"] == "0|1":
                    label.append("F1")
                elif row["fGT_p"] == "1|0":
                    label.append("F2")
            elif sAD[1] > ADO_thres // 2 and sAD[0] > ADO_thres // 2: # We assume there is no paternal contamination.
                if row["fGT_p"] == "0|1":
                    label.append("FM")
                elif row["fGT_p"] == "1|0":
                    label.append("FM")
            else:
                label.append("Unknown")
                
    truncated_data["Label"] = label

    return truncated_data