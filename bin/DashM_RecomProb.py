# This document calculates the recombination probablity between SNPs.

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import math

def all_recom_prob(data, chrom, rec_rate_pd, recom_rate, output):
    
    '''
    Calculate the recombination probability within a specific region of the chromosome.
    
    rec_rate_pd: files containing the information for calculation.
    
    Return: A dictionary where the keys are the position pairs of SNPs and the values are the estimated recombination probability.
    '''
    if output:
        print("Calculating Recombination Probability.....")
    recom_dict = {}
    for i in range(1,data.shape[0]):
        large_pos = max(data["POS"].iloc[i],data["POS"].iloc[i-1])
        small_pos = min(data["POS"].iloc[i],data["POS"].iloc[i-1])
        recom_dict[(large_pos,small_pos)] = calculate_recom_prob(chrom, large_pos, small_pos, rec_rate_pd, recom_rate, output)
    if output: 
        print("Finish calculating recombination probability.")
    return recom_dict

def kosambi_theta(d_cM):
    """Calculate recombination probability using Kosambi mapping function"""
    d = d_cM / 100  # Convert cM to Morgans
    return (math.exp(4*d) - 1) / (2*(math.exp(4*d) + 1))

def calculate_recom_prob(chrom, pos1, pos2, rec_rate_pd, recom_rate, output):
    
    '''Using recombination rate mean in the given range
    
    Arg: 
        chrom: number of chrom
        pos1, pos2: position at given chrom
        
    Return estimated recombination probabilty of two sites
    '''
    if pos1 > pos2:
        pos1, pos2 = pos2, pos1
    if (rec_rate_pd['Position(bp)'].iloc[0] > pos1):
        start_phy=rec_rate_pd['Position(bp)'].iloc[0]
    else:
        start_phy=rec_rate_pd.loc[rec_rate_pd['Position(bp)']<=pos1,'Position(bp)'].iloc[-1]
    if(rec_rate_pd['Position(bp)'].iloc[-1]< pos2):
        end_phy=rec_rate_pd['Position(bp)'].iloc[-1]
    else:
        end_phy = rec_rate_pd.loc[rec_rate_pd['Position(bp)']>=pos2,'Position(bp)'].iloc[0]

    start_phy_cm = rec_rate_pd.loc[(rec_rate_pd['Position(bp)']==start_phy),'Map(cM)'].iloc[0]
    end_phy_cm = rec_rate_pd.loc[(rec_rate_pd['Position(bp)']==end_phy),'Map(cM)'].iloc[0]
    d_cM = end_phy_cm - start_phy_cm
    prob = kosambi_theta(d_cM)


        
    #if isinstance(prob,float) and prob > 0 and 0 <= prob <= 0.5: # if we cannot query a recombination probability from the dataset
    if isinstance(prob,float) and 0 <= prob <= 0.5: # if we cannot query a recombination probability from the dataset
        return prob
    else:
        print(pos1,pos2,start_phy,end_phy,start_phy_cm,end_phy_cm,prob)
        #return abs(pos2 - pos1) * recom_rate # this is the average recombination rate in human chromosome.
        raise ValueError("Calculated probability out of valid range [0, 0.5]")

