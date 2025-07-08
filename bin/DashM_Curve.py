import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import random
from math import exp
from math import log

from scipy import stats

def stable_average(lst, eps=0.1):
    if len(lst) == 0:
        return None

    for X in [30, 10]:
        for i in range(len(lst) - X + 1,-1,-1):
            sublist = lst[i:i+X]
            if max(sublist) - min(sublist) <= eps:
                return sum(sublist) / len(sublist)

    return None

def calculate_log_likelihood(row,parent):
    if parent == "father":
        maxi1 = max(row["F1M1"], row["F1M2"])
        maxi2 = max(row["F2M1"], row["F2M2"])
        return maxi1 + np.log(np.exp(row["F1M1"] - maxi1) + np.exp(row["F1M2"] - maxi1)) - maxi2 - np.log(np.exp(row["F2M1"] - maxi2) + np.exp(row["F2M2"] - maxi2))
    elif parent == "mother":
        maxi1 = max(row["F1M1"], row["F2M1"])
        maxi2 = max(row["F1M2"], row["F2M2"])
        return maxi1 + np.log(np.exp(row["F1M1"] - maxi1) + np.exp(row["F2M1"] - maxi1)) - maxi2 - np.log(np.exp(row["F1M2"] - maxi2) + np.exp(row["F2M2"] - maxi2))
    else:
        raise NotImplementedError
    
def stable_value_DashM(sample, parent, file_path, eps=0.1):

    # get the data
    data_up = pd.read_csv(file_path + sample + "_" + parent + "_upper.txt", index_col=0)
    data_down = pd.read_csv(file_path + sample + "_" + parent + "_lower.txt", index_col=0)
    data_up.sort_values(by = "POS",ascending = False, inplace = True)
    data_down.sort_values(by = "POS",ascending = True, inplace = True)

    # calculate the stable values
    log_likelihood_up = data_up.apply(calculate_log_likelihood, axis = 1, parent = parent)
    log_likelihood_down = data_down.apply(calculate_log_likelihood, axis = 1, parent = parent)
    stable_up = stable_average(np.round(log_likelihood_up,2).values, eps = eps)
    stable_down = stable_average(np.round(log_likelihood_down,2).values, eps = eps)
    data_up['log_likelihood_up'] = log_likelihood_up ##qg add
    data_down['log_likelihood_down'] = log_likelihood_down ##qg add
    data_up.to_csv(file_path + sample + "_" + parent + "_up_with_log_likelihood.csv" ) ##qg add 
    data_down.to_csv(file_path + sample + "_" + parent + "_down_with_log_likelihood.csv" ) ##qg add
    count_up_greater_than_zero = (data_up['log_likelihood_up'] >= 0).sum()
    count_down_greater_than_zero = (data_down['log_likelihood_down'] >= 0).sum()
    count_up_less_than_zero = (data_up['log_likelihood_up'] <= 0).sum()
    count_down_less_than_zero = (data_down['log_likelihood_down'] <= 0).sum()


    if stable_up is None:
        up_sign = None
    else:
        up_sign = round(stable_up,2)
    if stable_down is None:
        down_sign = None
    else:
        down_sign = round(stable_down,2)

    if data_up.shape[0] < 10:
        up_sign = "Not Enough SNP"
        max_up, min_up = None, None
    else:
        max_up, min_up = round(max(log_likelihood_up),2), round(min(log_likelihood_up),2)
    
    if data_down.shape[0] < 10:
        down_sign = "Not Enough SNP"
        max_down, min_down = None, None
    else: 
        max_down, min_down = round(max(log_likelihood_down),2), round(min(log_likelihood_down),2)

    return (up_sign,down_sign, max_up, min_up, max_down, min_down,count_up_greater_than_zero,count_up_less_than_zero,count_down_greater_than_zero,count_down_less_than_zero)

