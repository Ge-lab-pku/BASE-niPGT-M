# This document calculates the recombination probablity between SNPs.

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

def all_recom_prob(data, chrom, loc_pd, rec_rate_pd, recom_rate, output):
    
    '''
    Calculate the recombination probability within a specific region of the chromosome.
    
    loc_pd, rec_rate_pd: files containing the information for calculation.
    
    Return: A dictionary where the keys are the position pairs of SNPs and the values are the estimated recombination probability.
    '''
    if output:
        print("Calculating Recombination Probability.....")
    recom_dict = {}
    for i in range(1,data.shape[0]):
        large_pos = max(data["POS"].iloc[i],data["POS"].iloc[i-1])
        small_pos = min(data["POS"].iloc[i],data["POS"].iloc[i-1])
        recom_dict[(large_pos,small_pos)] = calculate_recom_prob(chrom, large_pos, small_pos, loc_pd, rec_rate_pd, recom_rate, output)
    if output: 
        print("Finish calculating recombination probability.")
    return recom_dict


def calculate_recom_prob(chrom, pos1, pos2, loc_pd, rec_rate_pd, recom_rate, output):
    
    '''Using recombination rate mean in the given range
    
    Arg: 
        chrom: number of chrom
        pos1, pos2: position at given chrom
        
    Return estimated recombination probabilty of two sites
    '''
    
    if pos1 > pos2:
        pos1, pos2 = pos2, pos1
    
    if (chrom != "X" and chrom !="Y"):
        start_loc = loc_pd.loc[loc_pd['phys.loc']<=pos1,'gen.loc'].iloc[-1]
        end_loc = loc_pd.loc[loc_pd['phys.loc']>=pos2,'gen.loc'].iloc[1]
        start_phy=loc_pd.loc[loc_pd['phys.loc']<=pos1,'phys.loc'].iloc[-1]
        end_phy = loc_pd.loc[loc_pd['phys.loc']>=pos2,'phys.loc'].iloc[1]
        distance_in_loc = end_loc - start_loc
        rec_rate_mean = rec_rate_pd.loc[(rec_rate_pd['phys.loc']>=start_phy)&(rec_rate_pd['phys.loc']<=end_phy),'rec.rate'].mean()
        prob = rec_rate_mean * distance_in_loc/100
        
    else:
        start_loc = loc_pd.loc[loc_pd['phys.loc']<=pos1,'gen.loc'].iloc[-1]
        end_loc = loc_pd.loc[loc_pd['phys.loc']>=pos2,'gen.loc'].iloc[1]
        start_phy=loc_pd.loc[loc_pd['phys.loc']<=pos1,'phys.loc'].iloc[-1]
        end_phy = loc_pd.loc[loc_pd['phys.loc']>=pos2,'phys.loc'].iloc[1]
        distance_in_loc = end_loc - start_loc
        rec_rate_mean = rec_rate_pd.loc[(rec_rate_pd['phys.loc']>=start_phy)&(rec_rate_pd['phys.loc']<=end_phy),'rec.rate'].mean()
        prob = rec_rate_mean * distance_in_loc/100
        
    if isinstance(prob,float) and prob > 0 and 0 <= prob <= 1: # if we cannot query a recombination probability from the dataset
        return prob
    else:
        return abs(pos2 - pos1) * recom_rate # this is the average recombination rate in human chromosome.
