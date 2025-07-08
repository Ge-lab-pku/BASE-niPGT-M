# This document calculates likelihood of each SNP. Version 2.0

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from DashM_Basic import *


def row_likelihood(row, mcc_rate, state, err_rate):
    
    '''
    Calculates the sngle-SNP likelihood.
    
    mcc_rate, err_rate: the parameters estimated before.
    
    state: the hyplotype state. It can be either one of "F", "M", "FM", "F1", "F2"
    '''
    
    sAD = row["sAD"]

    if isinstance(sAD,str):
        sAD = eval(sAD)
        
    reads = sAD[0] + sAD[1]
    AD0 = min(sAD[0],100)
    AD1 = min(sAD[1],100) # truncate to 100 in order to avoid 'math domain error'in the calculation
    fGT = row["fGT_p"]
    mGT = row["mGT_p"]
    
    if state == "F":
        h = np.array([1,1,0,0])
    elif state == "F1":
        h = np.array([1,0,0,0])
    elif state == "F2":
        h = np.array([0,1,0,0])
    elif state == "M":
        h = np.array([0,0,1,1])
    elif state == "FM":
        h = np.array([1,1,1,1])
    else:
        raise ValueError("The hyplotype error", state)
    
    A11 = np.array([(1 - mcc_rate) / 2, 0, 1/2, mcc_rate / 2])
    A12 = np.array([(1 - mcc_rate) / 2, 0, mcc_rate / 2, 1/2])
    A21 = np.array([0, (1 - mcc_rate) / 2, 1/2, mcc_rate / 2])
    A22 = np.array([0, (1 - mcc_rate) / 2, mcc_rate / 2, 1/2])
    
    GT1 = np.array([int(fGT[0]),int(fGT[-1]),int(mGT[0]),int(mGT[-1])])
    GT0 = 1 - GT1
    
    for (i,j) in [(1,1),(1,2),(2,1),(2,2)]:
        i, j = str(i), str(j)

        if np.dot(eval("A" + i + j), h) == 0:
            globals()["rate0_" + i + j] = 0
            globals()["rate1_" + i + j] = 0
        else:
            globals()["rate0_" + i + j] = (np.dot(eval("A" + i + j) * h / np.dot(eval("A" + i + j), h),GT0) * (1 - err_rate) 
                                                    + np.dot(eval("A" + i + j) * h / np.dot(eval("A" + i + j), h),GT1) * (err_rate))

            globals()["rate1_" + i + j] = (np.dot(eval("A" + i + j) * h / np.dot(eval("A" + i + j), h),GT1) * (1 - err_rate) 
                                                    + np.dot(eval("A" + i + j) * h / np.dot(eval("A" + i + j), h),GT0) * (err_rate))

    for (i,j) in [(1,1),(1,2),(2,1),(2,2)]:
        i, j = str(i), str(j)
        globals()["p" + i + j] = (C(reads,AD0) * (eval("rate0_" + i + j)**AD0)
                                              * (eval("rate1_" + i + j)**AD1))

    return p11, p12, p21, p22

#def row_likelihood_slack(row, mcc_rate, f_loss0, m_loss0, err_rate):
#    p_floss = np.array(row_likelihood(row, mcc_rate, "M", err_rate))
#    p_mloss = np.array(row_likelihood(row, mcc_rate, "F", err_rate))
#    p_noloss = np.array(row_likelihood(row, mcc_rate, "FM", err_rate))
#    f_loss, m_loss, no_loss = f_loss0, m_loss0, 1 - f_loss0 - m_loss0
#    return f_loss * p_floss + m_loss * p_mloss + no_loss * p_noloss

def calculate_likelihood_F1(data, mcc_rate, err_rate, recom_dict, separated):
    '''
    calculate the likelihood under the condition that all the hyplotype states are F1
    '''

    log_likelihood = {x:[0,0] for x in data["POS"]} # [F1M1,F1M2]
    row = data.iloc[0]
    start_pos = data.iloc[0]["POS"]
    log_likelihood[start_pos] = [np.log(row_likelihood(row, mcc_rate, "F1", err_rate)[0]),
                             np.log(row_likelihood(row, mcc_rate, "F1", err_rate)[1])] # log-likelihood of [F1M1,F1M2]
    for i in range(1,data.shape[0]):
        row = data.iloc[i]
        current_pos = row["POS"]
        last_pos = data.iloc[i-1]["POS"]

        # calculate the recombnation probability.
        recom_prob = recom_dict[(max(current_pos,last_pos),min(current_pos,last_pos))]
        recom_prob_m = recom_prob

        log_likelihood[current_pos][0] = (np.log(row_likelihood(row,mcc_rate,"F1",err_rate)[0]) + log_likelihood[last_pos][0]
                    + np.log((1 - recom_prob_m)
                    + recom_prob_m * np.exp(log_likelihood[last_pos][1] - log_likelihood[last_pos][0])))
        log_likelihood[current_pos][1] = (np.log(row_likelihood(row,mcc_rate,"F1",err_rate)[1]) + log_likelihood[last_pos][1]
                    + np.log((1 - recom_prob_m)
                    + recom_prob_m * np.exp(log_likelihood[last_pos][0] - log_likelihood[last_pos][1])))
    
    final_pos = data.iloc[data.shape[0]-1]["POS"]    
    maximum = max(log_likelihood[final_pos])

    if separated:
        return log_likelihood[final_pos]
    else:
        ave_log_likelihood = (maximum
                            + np.log(np.exp(log_likelihood[final_pos][0] - maximum)
                                + np.exp(log_likelihood[final_pos][1] - maximum))-np.log(2))

        return ave_log_likelihood

def calculate_likelihood_F2(data, mcc_rate, err_rate, recom_dict, separated):
    '''
    calculate the likelihood under the condition that all the hyplotype states are F2
    '''

    log_likelihood = {x:[0,0] for x in data["POS"]} # [F2M1,F2M2]
    row = data.iloc[0]
    start_pos = data.iloc[0]["POS"]
    log_likelihood[start_pos] = [np.log(row_likelihood(row, mcc_rate, "F2", err_rate)[2]),
                             np.log(row_likelihood(row, mcc_rate, "F2", err_rate)[3])] # log-likelihood of [F2M1,F2M2]
    for i in range(1,data.shape[0]):
        row = data.iloc[i]
        current_pos = row["POS"]
        last_pos = data.iloc[i-1]["POS"]

        # calculate the recombnation probability.
        recom_prob = recom_dict[(max(current_pos,last_pos),min(current_pos,last_pos))]
        recom_prob_m = recom_prob

        log_likelihood[current_pos][0] = (np.log(row_likelihood(row,mcc_rate,"F2",err_rate)[2]) + log_likelihood[last_pos][0]
                    + np.log((1 - recom_prob_m)
                    + recom_prob_m * np.exp(log_likelihood[last_pos][1] - log_likelihood[last_pos][0])))
        log_likelihood[current_pos][1] = (np.log(row_likelihood(row,mcc_rate,"F2",err_rate)[3]) + log_likelihood[last_pos][1]
                    + np.log((1 - recom_prob_m)
                    + recom_prob_m * np.exp(log_likelihood[last_pos][0] - log_likelihood[last_pos][1])))
    
    final_pos = data.iloc[data.shape[0]-1]["POS"]    
    maximum = max(log_likelihood[final_pos])

    if separated:
        return log_likelihood[final_pos]
    else:
        ave_log_likelihood = (maximum
                            + np.log(np.exp(log_likelihood[final_pos][0] - maximum)
                                + np.exp(log_likelihood[final_pos][1] - maximum))-np.log(2))

        return ave_log_likelihood

def calculate_likelihood_others(data, mcc_rate, err_rate, recom_dict, state, separated): 
    
    '''
    calculate the likelihood under the condition that all the hyplotype states are FM or M

    mcc_rate, err_rate: the estimated parameter used to calculate.
    recom_dict: the recombination probability.
    separated: if True, return the log-likelihood of the data conditioanl on the F1M1, F1M2, F2M1, F2M2, where FiMj means that we get the i-th patrnal chain and the j-th maternal chain.
    '''

    log_likelihood = {x:[0,0,0,0] for x in data["POS"]}
    row = data.iloc[0]
    start_pos = data.iloc[0]["POS"]
    log_likelihood[start_pos] = [np.log(row_likelihood(row, mcc_rate, state, err_rate)[0]),
                             np.log(row_likelihood(row, mcc_rate, state, err_rate)[1]),
                             np.log(row_likelihood(row, mcc_rate, state, err_rate)[2]),
                             np.log(row_likelihood(row, mcc_rate, state, err_rate)[3])]

    for i in range(1,data.shape[0]):
        row = data.iloc[i]
        current_pos = row["POS"]
        last_pos = data.iloc[i-1]["POS"]

        # calculate the recombnation probability
        recom_prob = recom_dict[(max(current_pos,last_pos),min(current_pos,last_pos))]
        recom_prob_f = recom_prob
        recom_prob_m = recom_prob
        
        # Iterate
        log_likelihood[current_pos][0] = (np.log(row_likelihood(row,mcc_rate,state,err_rate)[0]) + log_likelihood[last_pos][0]
                    + np.log((1 - recom_prob_f) * (1 - recom_prob_m)
                    + (1 - recom_prob_f) * recom_prob_m * np.exp(log_likelihood[last_pos][1] - log_likelihood[last_pos][0])
                    + (1 - recom_prob_m) * recom_prob_f * np.exp(log_likelihood[last_pos][2] - log_likelihood[last_pos][0])
                    + recom_prob_f * recom_prob_m * np.exp(log_likelihood[last_pos][3] - log_likelihood[last_pos][0])))
        log_likelihood[current_pos][1] = (np.log(row_likelihood(row,mcc_rate,state,err_rate)[1]) + log_likelihood[last_pos][1]
                    + np.log((1 - recom_prob_f) * (1 - recom_prob_m)
                    + (1 - recom_prob_f) * recom_prob_m * np.exp(log_likelihood[last_pos][0] - log_likelihood[last_pos][1])
                    + (1 - recom_prob_m) * recom_prob_f * np.exp(log_likelihood[last_pos][3] - log_likelihood[last_pos][1])
                    + recom_prob_f * recom_prob_m * np.exp(log_likelihood[last_pos][2] - log_likelihood[last_pos][1])))
        log_likelihood[current_pos][2] = (np.log(row_likelihood(row,mcc_rate,state,err_rate)[2]) + log_likelihood[last_pos][2]
                    + np.log((1 - recom_prob_f) * (1 - recom_prob_m)
                    + (1 - recom_prob_m) * recom_prob_f * np.exp(log_likelihood[last_pos][0] - log_likelihood[last_pos][2])
                    + (1 - recom_prob_f) * recom_prob_m * np.exp(log_likelihood[last_pos][3] - log_likelihood[last_pos][2])
                    + recom_prob_f * recom_prob_m * np.exp(log_likelihood[last_pos][1] - log_likelihood[last_pos][2])))
        log_likelihood[current_pos][3] = (np.log(row_likelihood(row,mcc_rate,state,err_rate)[3]) + log_likelihood[last_pos][3]
                    + np.log((1 - recom_prob_f) * (1 - recom_prob_m)
                    + (1 - recom_prob_m) * recom_prob_f * np.exp(log_likelihood[last_pos][1] - log_likelihood[last_pos][3])
                    + (1 - recom_prob_f) * recom_prob_m * np.exp(log_likelihood[last_pos][2] - log_likelihood[last_pos][3])
                    + recom_prob_f * recom_prob_m * np.exp(log_likelihood[last_pos][0] - log_likelihood[last_pos][3]))) 
        
    final_pos = data.iloc[data.shape[0]-1]["POS"]    
    maximum = max(log_likelihood[final_pos])
    
    if separated:
        return log_likelihood[final_pos]
    else:
        ave_log_likelihood = (maximum
                        + np.log(np.exp(log_likelihood[final_pos][0] - maximum)
                            + np.exp(log_likelihood[final_pos][1] - maximum)
                            + np.exp(log_likelihood[final_pos][2] - maximum)
                            + np.exp(log_likelihood[final_pos][3] - maximum))-np.log(4))

        return ave_log_likelihood

def calculate_conditional_likelihood(data, mcc_rate, err_rate, recom_dict, cond_state, crt_state): 
    
    '''
    data: a period of data, which contains at least 2 SNPs.
    Suppose the data contains N SNPs, this calculates the likelihood that the 1 to (N-1) SNPs are at the hyplotype state = cond_state and the last SNP is at hyplotype state = crt_state.

    mcc_rate, err_rate: the estimated parameter used to calculate.
    recom_dict: the recombination probability.

    cond_state: the hyplotype state of the first (N-1) SNPs.
    crt_state: current hyplotype state at the N-th SNP.
    '''

    assert data.shape[0] >= 2
    past_data = data.iloc[:data.shape[0]-1]
    crt_log_likelihood = [0,0,0,0] # crt for "current"
    
    # calculate the recombnation probability.
    current_pos = data.iloc[data.shape[0]-1]["POS"]
    last_pos = data.iloc[data.shape[0]-2]["POS"]
    recom_prob = recom_dict[(max(current_pos,last_pos),min(current_pos,last_pos))]


    if cond_state == "F1": # If the hyplotype state of the first (N-1) SNPs is "F1".
        log_likelihood = calculate_likelihood_F1(past_data, mcc_rate, err_rate, recom_dict, True) # log_likelihood of [F1M1,F1M2]
        recom_prob_m = recom_prob
        recom_prob_f = recom_prob
        row = data.iloc[data.shape[0]-1] # the last SNP (the N-th SNP if the data contains N SNPs)

        if crt_state == "F1":

            crt_log_likelihood[0] = (np.log(row_likelihood(row,mcc_rate,"F1",err_rate)[0]) + log_likelihood[0]
                    + np.log((1 - recom_prob_m) * (1 - recom_prob_f)
                    + recom_prob_m * (1 - recom_prob_f) * np.exp(log_likelihood[1] - log_likelihood[0])))
            crt_log_likelihood[1] = (np.log(row_likelihood(row,mcc_rate,"F1",err_rate)[1]) + log_likelihood[1]
                    + np.log((1 - recom_prob_m) * (1 - recom_prob_f)
                    + recom_prob_m * (1 - recom_prob_f) * np.exp(log_likelihood[0] - log_likelihood[1])))

            maximum = max(crt_log_likelihood[:2])
            return (maximum + np.log(np.exp(crt_log_likelihood[0] - maximum)
                            + np.exp(crt_log_likelihood[1] - maximum))-np.log(2))

        if crt_state == "F2":

            crt_log_likelihood[2] = (np.log(row_likelihood(row,mcc_rate,"F2",err_rate)[2]) + log_likelihood[0]
                    + np.log((1 - recom_prob_m) * (recom_prob_f)
                    + recom_prob_m * (recom_prob_f) * np.exp(log_likelihood[1] - log_likelihood[0])))
            crt_log_likelihood[3] = (np.log(row_likelihood(row,mcc_rate,"F2",err_rate)[3]) + log_likelihood[1]
                    + np.log((1 - recom_prob_m) * (recom_prob_f)
                    + recom_prob_m * (recom_prob_f) * np.exp(log_likelihood[0] - log_likelihood[1])))

            maximum = max(crt_log_likelihood[1:])
            return (maximum + np.log(np.exp(crt_log_likelihood[2] - maximum)
                            + np.exp(crt_log_likelihood[3] - maximum))-np.log(2))

        elif crt_state == "F" or crt_state == "M" or crt_state == "FM":

            crt_log_likelihood[0] = (np.log(row_likelihood(row,mcc_rate,crt_state,err_rate)[0]) + log_likelihood[0]
                    + np.log((1 - recom_prob_m) * (1 - recom_prob_f)
                    + recom_prob_m * (1 - recom_prob_f) * np.exp(log_likelihood[1] - log_likelihood[0])))
            crt_log_likelihood[1] = (np.log(row_likelihood(row,mcc_rate,crt_state,err_rate)[1]) + log_likelihood[1]
                    + np.log((1 - recom_prob_m) * (1 - recom_prob_f)
                    + recom_prob_m * (1 - recom_prob_f) * np.exp(log_likelihood[0] - log_likelihood[1])))
            crt_log_likelihood[2] = (np.log(row_likelihood(row,mcc_rate,crt_state,err_rate)[2]) + log_likelihood[0]
                    + np.log((1 - recom_prob_m) * (recom_prob_f)
                    + recom_prob_m * (recom_prob_f) * np.exp(log_likelihood[1] - log_likelihood[0])))
            crt_log_likelihood[3] = (np.log(row_likelihood(row,mcc_rate,crt_state,err_rate)[3]) + log_likelihood[1]
                    + np.log((1 - recom_prob_m) * (recom_prob_f)
                    + recom_prob_m * (recom_prob_f) * np.exp(log_likelihood[0] - log_likelihood[1])))
            maximum = max(crt_log_likelihood)
            return (maximum + np.log(
                np.exp(crt_log_likelihood[0] - maximum) + np.exp(crt_log_likelihood[1] - maximum)
                + np.exp(crt_log_likelihood[2] - maximum) + np.exp(crt_log_likelihood[3] - maximum)
                )-np.log(4))

    if cond_state == "F2": # If the hyplotype state of the first (N-1) SNPs is "F1".
        log_likelihood = calculate_likelihood_F2(past_data, mcc_rate, err_rate, recom_dict, True) # log_likelihood of [F2M1,F2M2]
        recom_prob_m = recom_prob
        recom_prob_f = recom_prob
        row = data.iloc[data.shape[0]-1] # the last SNP (the N-th SNP if the data contains N SNPs)

        if crt_state == "F1":

            crt_log_likelihood[0] = (np.log(row_likelihood(row,mcc_rate,"F1",err_rate)[0]) + log_likelihood[0]
                    + np.log((1 - recom_prob_m) * (recom_prob_f)
                    + recom_prob_m * (recom_prob_f) * np.exp(log_likelihood[1] - log_likelihood[0])))
            crt_log_likelihood[1] = (np.log(row_likelihood(row,mcc_rate,"F1",err_rate)[1]) + log_likelihood[1]
                    + np.log((1 - recom_prob_m) * (recom_prob_f)
                    + recom_prob_m * (recom_prob_f) * np.exp(log_likelihood[0] - log_likelihood[1])))

            maximum = max(crt_log_likelihood[:2])
            return (maximum + np.log(np.exp(crt_log_likelihood[0] - maximum)
                            + np.exp(crt_log_likelihood[1] - maximum))-np.log(2))

        if crt_state == "F2":

            crt_log_likelihood[2] = (np.log(row_likelihood(row,mcc_rate,"F2",err_rate)[2]) + log_likelihood[0]
                    + np.log((1 - recom_prob_m) * (1 - recom_prob_f)
                    + recom_prob_m * (1 - recom_prob_f) * np.exp(log_likelihood[1] - log_likelihood[0])))
            crt_log_likelihood[3] = (np.log(row_likelihood(row,mcc_rate,"F2",err_rate)[3]) + log_likelihood[1]
                    + np.log((1 - recom_prob_m) * (1 - recom_prob_f)
                    + recom_prob_m * (1 - recom_prob_f) * np.exp(log_likelihood[0] - log_likelihood[1])))

            maximum = max(crt_log_likelihood[1:])
            return (maximum + np.log(np.exp(crt_log_likelihood[2] - maximum)
                            + np.exp(crt_log_likelihood[3] - maximum))-np.log(2))

        elif crt_state == "F" or crt_state == "M" or crt_state == "FM":

            crt_log_likelihood[0] = (np.log(row_likelihood(row,mcc_rate,crt_state,err_rate)[0]) + log_likelihood[0]
                    + np.log((1 - recom_prob_m) * (recom_prob_f)
                    + recom_prob_m * (recom_prob_f) * np.exp(log_likelihood[1] - log_likelihood[0])))
            crt_log_likelihood[1] = (np.log(row_likelihood(row,mcc_rate,crt_state,err_rate)[1]) + log_likelihood[1]
                    + np.log((1 - recom_prob_m) * (recom_prob_f)
                    + recom_prob_m * (recom_prob_f) * np.exp(log_likelihood[0] - log_likelihood[1])))
            crt_log_likelihood[2] = (np.log(row_likelihood(row,mcc_rate,crt_state,err_rate)[2]) + log_likelihood[0]
                    + np.log((1 - recom_prob_m) * (1 - recom_prob_f)
                    + recom_prob_m * (1 - recom_prob_f) * np.exp(log_likelihood[1] - log_likelihood[0])))
            crt_log_likelihood[3] = (np.log(row_likelihood(row,mcc_rate,crt_state,err_rate)[3]) + log_likelihood[1]
                    + np.log((1 - recom_prob_m) * (1 - recom_prob_f)
                    + recom_prob_m * (1 - recom_prob_f) * np.exp(log_likelihood[0] - log_likelihood[1])))
            maximum = max(crt_log_likelihood)
            return (maximum + np.log(
                np.exp(crt_log_likelihood[0] - maximum) + np.exp(crt_log_likelihood[1] - maximum)
                + np.exp(crt_log_likelihood[2] - maximum) + np.exp(crt_log_likelihood[3] - maximum)
                )-np.log(4))

    elif cond_state == "F" or cond_state == "M" or cond_state == "FM":
        log_likelihood = calculate_likelihood_others(past_data, mcc_rate, err_rate, recom_dict, cond_state, True) 
        # log_likelihood of [F1M1, F1M2, F2M1,F2M2]. They should all be non-zero.

        recom_prob_m = recom_prob
        recom_prob_f = recom_prob
        row = data.iloc[data.shape[0]-1] # the last SNP (the N-th SNP if the data contains N SNPs)

        if crt_state == "F1":

            crt_log_likelihood[0] = (np.log(row_likelihood(row,mcc_rate,"F1",err_rate)[0]) + log_likelihood[0]
                    + np.log((1 - recom_prob_m) * (1 - recom_prob_f)
                    + recom_prob_m * (1 - recom_prob_f) * np.exp(log_likelihood[1] - log_likelihood[0])
                    + (1 - recom_prob_m) * (recom_prob_f) * np.exp(log_likelihood[2] - log_likelihood[0])
                    + recom_prob_m * (recom_prob_f) * np.exp(log_likelihood[3] - log_likelihood[0])))
            crt_log_likelihood[1] = (np.log(row_likelihood(row,mcc_rate,"F1",err_rate)[1]) + log_likelihood[1]
                    + np.log((1 - recom_prob_m) * (1 - recom_prob_f)
                    + recom_prob_m * (1 - recom_prob_f) * np.exp(log_likelihood[0] - log_likelihood[1])
                    + (1 - recom_prob_m) * (recom_prob_f) * np.exp(log_likelihood[3] - log_likelihood[1])
                    + recom_prob_m * (recom_prob_f) * np.exp(log_likelihood[2] - log_likelihood[1])))

            maximum = max(crt_log_likelihood[:2])
            return (maximum + np.log(np.exp(crt_log_likelihood[0] - maximum)
                            + np.exp(crt_log_likelihood[1] - maximum))-np.log(2))

        elif crt_state == "F2":

            crt_log_likelihood[2] = (np.log(row_likelihood(row,mcc_rate,"F2",err_rate)[2]) + log_likelihood[2]
                    + np.log((1 - recom_prob_m) * (1 - recom_prob_f)
                    + recom_prob_m * (1 - recom_prob_f) * np.exp(log_likelihood[3] - log_likelihood[2])
                    + (1 - recom_prob_m) * (recom_prob_f) * np.exp(log_likelihood[0] - log_likelihood[2])
                    + recom_prob_m * (recom_prob_f) * np.exp(log_likelihood[1] - log_likelihood[2])))
            crt_log_likelihood[3] = (np.log(row_likelihood(row,mcc_rate,"F2",err_rate)[3]) + log_likelihood[3]
                    + np.log((1 - recom_prob_m) * (1 - recom_prob_f)
                    + recom_prob_m * (1 - recom_prob_f) * np.exp(log_likelihood[2] - log_likelihood[3])
                    + (1 - recom_prob_m) * (recom_prob_f) * np.exp(log_likelihood[1] - log_likelihood[3])
                    + recom_prob_m * (recom_prob_f) * np.exp(log_likelihood[0] - log_likelihood[3])))

            maximum = max(crt_log_likelihood[1:])
            return (maximum + np.log(np.exp(crt_log_likelihood[2] - maximum)
                            + np.exp(crt_log_likelihood[3] - maximum))-np.log(2))

        elif crt_state == "F" or crt_state == "M" or crt_state == "FM":

            crt_log_likelihood[0] = (np.log(row_likelihood(row,mcc_rate,crt_state,err_rate)[0]) + log_likelihood[0]
                    + np.log((1 - recom_prob_m) * (1 - recom_prob_f)
                    + recom_prob_m * (1 - recom_prob_f) * np.exp(log_likelihood[1] - log_likelihood[0])
                    + (1 - recom_prob_m) * (recom_prob_f) * np.exp(log_likelihood[2] - log_likelihood[0])
                    + recom_prob_m * (recom_prob_f) * np.exp(log_likelihood[3] - log_likelihood[0])))
            crt_log_likelihood[1] = (np.log(row_likelihood(row,mcc_rate,crt_state,err_rate)[1]) + log_likelihood[1]
                    + np.log((1 - recom_prob_m) * (1 - recom_prob_f)
                    + recom_prob_m * (1 - recom_prob_f) * np.exp(log_likelihood[0] - log_likelihood[1])
                    + (1 - recom_prob_m) * (recom_prob_f) * np.exp(log_likelihood[3] - log_likelihood[1])
                    + recom_prob_m * (recom_prob_f) * np.exp(log_likelihood[2] - log_likelihood[1])))
            crt_log_likelihood[2] = (np.log(row_likelihood(row,mcc_rate,crt_state,err_rate)[2]) + log_likelihood[2]
                    + np.log((1 - recom_prob_m) * (1 - recom_prob_f)
                    + recom_prob_m * (1 - recom_prob_f) * np.exp(log_likelihood[3] - log_likelihood[2])
                    + (1 - recom_prob_m) * (recom_prob_f) * np.exp(log_likelihood[0] - log_likelihood[2])
                    + recom_prob_m * (recom_prob_f) * np.exp(log_likelihood[1] - log_likelihood[2])))
            crt_log_likelihood[3] = (np.log(row_likelihood(row,mcc_rate,crt_state,err_rate)[3]) + log_likelihood[3]
                    + np.log((1 - recom_prob_m) * (1 - recom_prob_f)
                    + recom_prob_m * (1 - recom_prob_f) * np.exp(log_likelihood[2] - log_likelihood[3])
                    + (1 - recom_prob_m) * (recom_prob_f) * np.exp(log_likelihood[1] - log_likelihood[3])
                    + recom_prob_m * (recom_prob_f) * np.exp(log_likelihood[0] - log_likelihood[3])))
            maximum = max(crt_log_likelihood)
            return (maximum + np.log(
                np.exp(crt_log_likelihood[0] - maximum) + np.exp(crt_log_likelihood[1] - maximum)
                + np.exp(crt_log_likelihood[2] - maximum) + np.exp(crt_log_likelihood[3] - maximum)
                )-np.log(4))