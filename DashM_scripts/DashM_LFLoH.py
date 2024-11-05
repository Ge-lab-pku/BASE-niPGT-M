# LFLoH 3.0: This document calbrates the mcc_rate and the LFLoH states

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from DashM_Basic import *
from DashM_RowLikelihood import *
from DashM_Label import *
from DashM_RecomProb import *

def calibrate(data, chrom, true_pos, window_size, sample, parent, err_rate, mcc_rate, ADO_thres, LFLoH_thres, recom_probs, df_save, df_save_path, flag, loc_pd, rec_rate_pd, recom_rate, sex = None):
    
    '''
    Calibrates the MCC Rate and LFLoH(Large Fraction Loss of Hyplotype) in the data.
    
    sample: the name of the sample.
    parent: "father" or "mother".
    recom_probs: the calculated recombination probabilty clculated before.
    LFLoH_thres: used to control the degree of conservation of LFLoH calibration.
    err_rate: estimated error rate.
    mcc_rate: the initially estimated mcc rate.
    df_save: boolean variable, denoting whether we store the result.
    df_save_path: the filename to save the file.
    '''
    
    print("Calbrating the MCC Rates and Hyplotype Loss......\n")

    mcc_flag = flag
    
    if chrom != "X" or (chrom == "X" and sex == "Female"):
        n_update = 0 # Keep record of the number of mcc_rate updates
        mcc_history = [] # Keep record of the hstoric mcc_rate estimates

        while True:
            if n_update >= 10:
                break
                
            new_truncated_data = label(data, chrom, true_pos, window_size, ADO_thres, sex, output = False)
            df = LFLoH(new_truncated_data, LFLoH_thres, mcc_rate, err_rate, recom_probs, chrom, loc_pd, rec_rate_pd, recom_rate)
            new_mcc_rate = recalculate_MCC_Rate(df, mcc_rate)
            print("New MCC Rate: ",new_mcc_rate)
            if abs(new_mcc_rate - mcc_rate) < 0.05:
                break
            else:
                if new_mcc_rate > 0:
                    mcc_rate = new_mcc_rate
                    n_update += 1
                    mcc_history.append(new_mcc_rate)
                    continue
                elif mcc_flag:
                    new_mcc_rate = 0.01
                    break
                else:
                    mcc_rate = 0.01
                    mcc_flag = True
                    continue

        if n_update >= 10:
            new_mcc_rate = np.mean(np.array(mcc_history))
        if new_mcc_rate > 0:
            mcc_rate = new_mcc_rate
        else:
            mcc_rate = 0.01
            
    elif chrom == "X" and sex == "Male":
        if mcc_rate > 0:
            mcc_rate = mcc_rate
        else:
            mcc_rate = 0.01
    print("Final MCC Rate: ",mcc_rate)

    new_truncated_data = label(data,chrom,true_pos,window_size,ADO_thres, output = False)
    df = LFLoH(new_truncated_data, LFLoH_thres, mcc_rate, err_rate, recom_probs, chrom, loc_pd, rec_rate_pd, recom_rate)
    if df_save:
        df.to_csv(df_save_path + sample + "_" + parent + ".txt")
        
    return {"data": df, "mcc_rate_final": mcc_rate}

def LFLoH_probability1(data, mcc_rate, err_rate, recom_dict, states):

    '''
    data contains N (at least 2) SNPs
    Calculate the ratio of the likelihood that SNPs are at each hyplotype state in 'states'.

    mcc_rate, err_rate: the estimated parameter.
    recom_dict: the dictionary of recombination probability
    '''

    assert data.shape[0] >= 2
    prob_dict = {}
    for crt_state in states:
        prob_dict[crt_state] =  calculate_conditional_likelihood(data, mcc_rate, err_rate, recom_dict, crt_state, crt_state)

    maximum = max(list(prob_dict.values()))
    denominator = (np.exp(np.array(list(prob_dict.values())) - maximum)).sum()
    for crt_state in states:
        prob_dict[crt_state] = np.exp(prob_dict[crt_state] - maximum) / denominator
    return prob_dict

def LFLoH_probability2(data, mcc_rate, err_rate, recom_dict, cond_state, states):

    '''
    data contains N (at least 2) SNPs and the hyplotype state of the 1 - (N-1)-th SNPs is cond_state
    Calculate the ratio of the likelihood that the N-th SNP is at different hyplotype state.

    mcc_rate, err_rate: the estimated parameter.
    recom_dict: the dictionary of recombination probability
    cond_state: the hyplotype state of 1-(N-1) -th SNP.
    '''

    assert data.shape[0] >= 2
    prob_dict = {}
    for crt_state in states:
        prob_dict[crt_state] =  calculate_conditional_likelihood(data, mcc_rate, err_rate, recom_dict, cond_state, crt_state)

    maximum = max(list(prob_dict.values()))
    denominator = (np.exp(np.array(list(prob_dict.values())) - maximum)).sum()
    for crt_state in states:
        prob_dict[crt_state] = np.exp(prob_dict[crt_state] - maximum) / denominator
    return prob_dict
    
def LFLoH(data, LFLoH_thres, mcc_rate, err_rate, recom_probs, chrom, loc_pd, rec_rate_pd, recom_rate):

    '''
    Determine the hyplotype state based on the current err_rate and mcc_rate.

    df: the truncated data within the adjacent region of the disease causal mutation site.
    chrom: the number of chromosome
    LFLoH_thres: A dictionary which contains the threshold for determining the LFLoH
    mcc_rate, err_rate: the estimated parameter
    recom_probs: the recombination probability.
    loc_pd, rec_rate_pd: the documents to calculate the recombination probability.
    '''
    
    df = data
    pos_lst = list(df["POS"])
    upper_label = ["Unknown" for x in df["POS"]]
    lower_label = ["Unknown" for x in df["POS"]]
    markers = df[(df["Label"] != "Unknown")] 
    # All SNPs which are not labeled as 'Unknown' in the inital labeling process. We call the markers
    
    for i in range(markers.shape[0]):

        # Upward searching
        marker = markers.iloc[i] # Get the marker
        marker_pos = marker["POS"] # Get the physical position of marker
        crt_index = pos_lst.index(marker_pos) # Get the index of the marker in df.
        crt_label = marker["Label"] # Get the label of the marker
        upper_label[crt_index] = crt_label # The labels of the markers will not change.
        
        if crt_label == "FM":
            continue
        else:
            crt_thres = LFLoH_thres[crt_label]

        upper_index = pos_lst.index(markers.iloc[i-1]["POS"]) if i != 0 else -1
        # The searching region for the current marker is [upper_index+1,lower_index), left closed right open
        
        crt_upper_index = crt_index - 1
        break_reason = 0 
        # to record the reason why the iteration ends. 
        # This is one if the iteration ends due to the likelihood less than the thresholds.
        # This is two if the SNP is beyond the range we consider (we only consider 10 SNPs around for each central SNP).
        
        while True:
            if abs(crt_upper_index - crt_index) >= 10:
                break_reason = 1
                break
            if crt_upper_index == upper_index:
                break_reason = 1
                break
            if crt_label == "F1" or crt_label == "F2":
                states = ["F1","F2","M","FM"]

                if crt_label == "F1": # determine whether the new SNP is "F1"
                    prob_F1 = row_likelihood(df.iloc[crt_upper_index], mcc_rate, "F1", err_rate)[0]
                    prob_F2 = row_likelihood(df.iloc[crt_upper_index], mcc_rate, "F2", err_rate)[2]
                    if prob_F2 > prob_F1:
                        break_reason = 1
                        break
                    else:
                        pass

                if crt_label == "F2": # determine whether the new SNP is "F2"
                    prob_F1 = row_likelihood(df.iloc[crt_upper_index], mcc_rate, "F1", err_rate)[0]
                    prob_F2 = row_likelihood(df.iloc[crt_upper_index], mcc_rate, "F2", err_rate)[2]
                    if prob_F1 > prob_F2:
                        break_reason = 1
                        break
                    else:
                        pass

            elif crt_label == "F" or crt_label == "M":
                states = ["F","M","FM"]
            prob_dict1 = LFLoH_probability1(df.iloc[crt_upper_index:crt_index+1].sort_values(by = "POS", ascending = False), mcc_rate, err_rate, recom_probs, states)
            max_label = max(prob_dict1, key = prob_dict1.get)
            if max_label == crt_label and prob_dict1[max_label] >= crt_thres:
                upper_label[crt_upper_index] = crt_label
                crt_upper_index -= 1 # crt_label and crt_thres do not change, move to the upper SNP.
                continue
            else:
                break_reason = 2
                break
        
        if break_reason == 2 and crt_upper_index > upper_index:   
            if crt_label == "F1" or crt_label == "F2":
                states = ["F1","F2","M","FM"]
                prob_dict2 = LFLoH_probability2(df.iloc[crt_upper_index:crt_index+1].sort_values(by = "POS", ascending = False), mcc_rate, err_rate, recom_probs, crt_label, states)
                max_label = max(prob_dict2, key = prob_dict2.get)

                if max_label == "M" and prob_dict2[max_label] >= LFLoH_thres[max_label]:
                    upper_label[crt_upper_index] = max_label

                    if crt_upper_index - 1 > upper_index:
                        recom_probs.update(all_recom_prob(df.iloc[np.r_[crt_upper_index-1,crt_upper_index+1:crt_index+1]], chrom, loc_pd, rec_rate_pd, recom_rate, output = False))
                        prob_dict2 = LFLoH_probability2(df.iloc[np.r_[crt_upper_index-1,crt_upper_index+1:crt_index+1]].sort_values(by = "POS", ascending = False), mcc_rate, err_rate, recom_probs, crt_label, states)
                        max_label = max(prob_dict2, key = prob_dict2.get)
                        if max_label == "M" and prob_dict2[max_label] >= LFLoH_thres[max_label]:
                            upper_label[crt_upper_index+1] = max_label
                        else:
                            pass
                    else:
                        pass
                else:
                    pass
                
            elif crt_label == "F":
                states = ["F","M","FM"]
                prob_dict2 = LFLoH_probability2(df.iloc[crt_upper_index:crt_index+1].sort_values(by = "POS", ascending = False), mcc_rate, err_rate, recom_probs, crt_label, states)
                max_label = max(prob_dict2, key = prob_dict2.get)

                if max_label == "M" and prob_dict2[max_label] >= LFLoH_thres[max_label]:
                    upper_label[crt_upper_index] = max_label

                    if crt_upper_index - 1 > upper_index:
                        recom_probs.update(all_recom_prob(df.iloc[np.r_[crt_upper_index-1,crt_upper_index+1:crt_index+1]], chrom, loc_pd, rec_rate_pd, recom_rate, output = False))
                        prob_dict2 = LFLoH_probability2(df.iloc[np.r_[crt_upper_index-1,crt_upper_index+1:crt_index+1]].sort_values(by = "POS", ascending = False), mcc_rate, err_rate, recom_probs, crt_label, states)
                        max_label = max(prob_dict2, key = prob_dict2.get)
                        if max_label == "M" and prob_dict2[max_label] >= LFLoH_thres[max_label]:
                            upper_label[crt_upper_index+1] = max_label
                        else:
                            pass
                    else:
                        pass
                else:
                    pass
                
            elif crt_label == "M":
                states = ["F","M","FM"]
                prob_dict2 = LFLoH_probability2(df.iloc[crt_upper_index:crt_index+1].sort_values(by = "POS", ascending = False), mcc_rate, err_rate, recom_probs, crt_label, states)
                max_label = max(prob_dict2, key = prob_dict2.get)

                if max_label == "F" and prob_dict2[max_label] >= LFLoH_thres[max_label]:
                    upper_label[crt_upper_index] = max_label

                    if crt_upper_index - 1 > upper_index:
                        recom_probs.update(all_recom_prob(df.iloc[np.r_[crt_upper_index-1,crt_upper_index+1:crt_index+1]], chrom, loc_pd, rec_rate_pd, recom_rate, output = False))
                        prob_dict2 = LFLoH_probability2(df.iloc[np.r_[crt_upper_index-1,crt_upper_index+1:crt_index+1]].sort_values(by = "POS", ascending = False), mcc_rate, err_rate, recom_probs, crt_label, states)
                        max_label = max(prob_dict2, key = prob_dict2.get)
                        if max_label == "F" and prob_dict2[max_label] >= LFLoH_thres[max_label]:
                            upper_label[crt_upper_index+1] = max_label
                        else:
                            pass
                    else:
                        pass
                else:
                    pass
            else:
                raise ValueError("Label Error:", crt_label)
            

        # Downward searching
        marker = markers.iloc[i] # Get the marker
        marker_pos = marker["POS"] # Get the physical position of marker
        crt_index = pos_lst.index(marker_pos) # Get the index of the marker in df.
        crt_label = marker["Label"] # Get the label of the marker
        lower_label[crt_index] = crt_label # The labels of the markers will not change.
        
        if crt_label == "FM":
            continue
        else:
            crt_thres = LFLoH_thres[crt_label]

        lower_index = pos_lst.index(markers.iloc[i+1]["POS"]) if i != markers.shape[0] - 1 else len(pos_lst)
        # The searching region for the current marker is [upper_index+1,lower_index), left closed right open
        
        crt_lower_index = crt_index + 1
            
        break_reason = 0
        
        while True:
            if abs(crt_lower_index - crt_index) >= 10:
                break_reason = 1
                break
            if crt_lower_index == lower_index:
                break_reason = 1
                break
            if crt_label == "F1" or crt_label == "F2":
                states = ["F1","F2","M","FM"]

                if crt_label == "F1": # determine whether the new SNP is "F1"
                    prob_F1 = row_likelihood(df.iloc[crt_lower_index], mcc_rate, "F1", err_rate)[0]
                    prob_F2 = row_likelihood(df.iloc[crt_lower_index], mcc_rate, "F2", err_rate)[2]
                    if prob_F2 > prob_F1:
                        break_reason = 1
                        break
                    else:
                        pass

                if crt_label == "F2": # determine whether the new SNP is "F2"
                    prob_F1 = row_likelihood(df.iloc[crt_lower_index], mcc_rate, "F1", err_rate)[0]
                    prob_F2 = row_likelihood(df.iloc[crt_lower_index], mcc_rate, "F2", err_rate)[2]
                    if prob_F1 > prob_F2:
                        break_reason = 1
                        break
                    else:
                        pass

            elif crt_label == "F" or crt_label == "M":
                states = ["F","M","FM"]
            prob_dict1 = LFLoH_probability1(df.iloc[crt_index:crt_lower_index+1].sort_values(by = "POS", ascending = True), mcc_rate, err_rate, recom_probs, states)
            max_label = max(prob_dict1, key = prob_dict1.get)
            if max_label == crt_label and prob_dict1[max_label] >= crt_thres:
                lower_label[crt_lower_index] = crt_label
                crt_lower_index += 1 # crt_label and crt_thres do not change, move to the lower SNP.
                continue
            else:
                break_reason = 2
                break
                
        if break_reason == 2 and  crt_lower_index < lower_index:   
                
            if crt_label == "F1" or crt_label == "F2":
                states = ["F1","F2","M","FM"]
                prob_dict2 = LFLoH_probability2(df.iloc[crt_index:crt_lower_index+1].sort_values(by = "POS", ascending = True), mcc_rate, err_rate, recom_probs, crt_label, states)
                max_label = max(prob_dict2, key = prob_dict2.get)

                if max_label == "M" and prob_dict2[max_label] >= LFLoH_thres[max_label]:
                    lower_label[crt_lower_index] = max_label
                    
                    if crt_lower_index + 1 < lower_index:
                        recom_probs.update(all_recom_prob(df.iloc[np.r_[crt_index:crt_lower_index,crt_lower_index+1]], chrom, loc_pd, rec_rate_pd, recom_rate, output = False))
                        prob_dict2 = LFLoH_probability2(df.iloc[np.r_[crt_index:crt_lower_index,crt_lower_index+1]].sort_values(by = "POS", ascending = True), mcc_rate, err_rate, recom_probs, crt_label, states)
                        max_label = max(prob_dict2, key = prob_dict2.get)
                        if max_label == "M" and prob_dict2[max_label] >= LFLoH_thres[max_label]:
                            lower_label[crt_lower_index+1] = max_label
                        else:
                            pass
                    else:
                        pass
                else:
                    pass
                
            elif crt_label == "F":
                states = ["F","M","FM"]
                prob_dict2 = LFLoH_probability2(df.iloc[crt_index:crt_lower_index+1].sort_values(by = "POS", ascending = True), mcc_rate, err_rate, recom_probs, crt_label, states)
                max_label = max(prob_dict2, key = prob_dict2.get)

                if max_label == "M" and prob_dict2[max_label] >= LFLoH_thres[max_label]:
                    lower_label[crt_lower_index] = max_label
                    
                    if crt_lower_index + 1 < lower_index:
                        recom_probs.update(all_recom_prob(df.iloc[np.r_[crt_index:crt_lower_index,crt_lower_index+1]], chrom, loc_pd, rec_rate_pd, recom_rate, output = False))
                        prob_dict2 = LFLoH_probability2(df.iloc[np.r_[crt_index:crt_lower_index,crt_lower_index+1]].sort_values(by = "POS", ascending = True), mcc_rate, err_rate, recom_probs, crt_label, states)
                        max_label = max(prob_dict2, key = prob_dict2.get)
                        if max_label == "M" and prob_dict2[max_label] >= LFLoH_thres[max_label]:
                            lower_label[crt_lower_index+1] = max_label
                        else:
                            pass
                    else:
                        pass
                else:
                    pass
                
            elif crt_label == "M":
                states = ["F","M","FM"]
                prob_dict2 = LFLoH_probability2(df.iloc[crt_index:crt_lower_index+1].sort_values(by = "POS", ascending = True), mcc_rate, err_rate, recom_probs, crt_label, states)
                max_label = max(prob_dict2, key = prob_dict2.get)

                if max_label == "F" and prob_dict2[max_label] >= LFLoH_thres[max_label]:
                    lower_label[crt_lower_index] = max_label
                    
                    if crt_lower_index + 1 < lower_index:
                        recom_probs.update(all_recom_prob(df.iloc[np.r_[crt_index:crt_lower_index,crt_lower_index+1]], chrom, loc_pd, rec_rate_pd, recom_rate, output = False))
                        prob_dict2 = LFLoH_probability2(df.iloc[np.r_[crt_index:crt_lower_index,crt_lower_index+1]].sort_values(by = "POS", ascending = True), mcc_rate, err_rate, recom_probs, crt_label, states)
                        max_label = max(prob_dict2, key = prob_dict2.get)
                        if max_label == "F" and prob_dict2[max_label] >= LFLoH_thres[max_label]:
                            lower_label[crt_lower_index+1] = max_label
                        else:
                            pass
                    else:
                        pass
                else:
                    pass    
            else:
                raise ValueError("Label Error:", crt_label)

    df["Upper_label"] = upper_label
    df["Lower_label"] = lower_label

    def LFLoH_row(row):
        if row["Label"] != "Unknown":
            return row["Label"]
        else:
            label_pair = (row["Upper_label"],row["Lower_label"])
            if label_pair in [("F1","F1"),("F1","F"),("F","F1"),("F1","Unknown"),("Unknown","F1")]:
                return "F1"
            elif label_pair in [("F2","F2"),("F2","F"),("F","F2"),("F2","Unknown"),("Unknown","F2")]:
                return "F2"
            elif label_pair in [("F","F"),("F","Unknown"),("Unknown","F")]:
                return "F"
            elif label_pair in [("M","M"),("M","Unknown"),("Unknown","M")]:
                return "M"
            else:
                return "FM"
                    
    df["Final_label"] = df.apply(LFLoH_row,axis = 1)
            
    return df
    
def recalculate_MCC_Rate(truncated_data, original_mcc_rate):
    '''
    recalculate the mcc rates based on the last iteration and hyplotype states.

    truncated_data: the data used to calculate
    '''

    df1 = truncated_data[(truncated_data["Final_label"] == "FM")]

    df1 = df1[((df1["fGT"] == "0/0")&(df1["mGT"] == "1/1"))|((df1["fGT"] == "1/1")&(df1["mGT"] == "0/0"))]

    df1_reads = [df1.iloc[i]["sAD"][0] + df1.iloc[i]["sAD"][1] for i in range(df1.shape[0])]

    if df1.shape[0] <= 3:
        new_mcc_rate = original_mcc_rate
        return new_mcc_rate
    else:
        df1 = df1[stats.zscore(df1_reads) <= 3] # remove the SNPs with extremely high AD values, like hundreds or thousands of reads.
        print(df1.shape)

    mcc = [0,0]
    mcc_reads = [] # keep record of the (reads from mother, total reads) pairs.
    for i in range(df1.shape[0]):
        row = df1.iloc[i]
        x = row["sAD"]
        if isinstance(x,str):
            x = eval(x)
        if row["fGT"] == "0/0" and row["mGT"] == "1/1":
            mcc[0] += x[0]
            mcc[1] += x[1]
            mcc_reads.append((x[1], x[1] + x[0]))
        elif row["fGT"] == "1/1" and row["mGT"] == "0/0":
            mcc[1] += x[0]
            mcc[0] += x[1]
            mcc_reads.append((x[0], x[1] + x[0]))
        else:
            continue
            
    print("The (reads from mother, total reads) pairs are", mcc_reads)
    
    if sum(mcc) == 0:
        new_mcc_rate = 0
    else:
        new_mcc_rate = mcc[1]/sum(mcc) * 2 - 1
    return new_mcc_rate
