import sys
import pandas as pd

def check_values(up_steady_value, up_max_value, up_min_value, down_max_value ,down_min_value):
    check = "Possible"
    check_cond = "None"
    if up_steady_value > 0:
        if up_min_value > -5 and down_min_value > -5:
            check = "Moderate"
        else:
            check_cond = "any(min_up, min_down <= -5)"
    elif up_steady_value < 0:
        if up_max_value < 5 and down_max_value < 5:
            check = "Moderate"
        else:
            check_cond = "any(max_up, max_down >= 5)"
    return check, check_cond


def evaluate_sample(row):
    up_steady, down_steady = row['Up_Steady_Value'], row['Down_Steady_Value']
    max_up, min_up, max_down, min_down = row['Upstream_Max_Value'],row['Upstream_Min_Value'],row['Downstream_Max_Value'],row['Downstream_Min_Value']

    cm_type = row.get('CM_type', 'T-SCM')
    result = "Possible"
    condition = "None"

    if any(row[col] >= threshold for col, threshold in [('Ratio_Low_GQ', 0.4), ('Error_Rate', 0.2),
                                                           ('Original_MCC_Rate', 0.75 if cm_type == 'T-SCM' else 0.85), ('Final_MCC_Rate', 0.75 if cm_type == 'T-SCM' else 0.85)]):
        result = "Rule Out"
        condition= "any(row[col] >= threshold)"
    elif any(row[col] <= threshold for col, threshold in [('SNP_Density', 0.001),('Original_MCC_Rate', -0.5 if cm_type == 'T-SCM' else -0.7), ('Final_MCC_Rate', -0.5 if cm_type == 'T-SCM' else -0.7)]):
        result = "Rule Out"
        condition = "any(row[col] <= threshold)"  ## "Rule Out" determined

    elif pd.isna(up_steady) or pd.isna(down_steady) or up_steady in ['', 'Not Enough SNP'] or down_steady in ['', 'Not Enough SNP']:
        return "Undetermined", " none"

    elif any(abs(steady) < 0.5 for steady in [up_steady, down_steady]):
        return "Undetermined","none and any(steady < 0.5)"

    elif up_steady * down_steady < 0 :
        result = "Undetermined"
        condition = "no none but diff sign" ## "undetermined determined"
    elif any(abs(steady) < 0.5 for steady in [up_steady, down_steady]):
        result = "Undetermined"
        condition = "no none but any(steady < 0.5)" ##none????

    elif any(abs(value) >= 1 for value in [up_steady, down_steady]):
        if abs((up_steady + down_steady)/2) >= 5:
            if any(row[col] == 0 for col in ['Up_positive_number', 'Up_negative_number']) and any(row[col] == 0 for col in ['Down_positive_number', 'Down_negative_number']):
                 result = "High" ##"high" determined
            else:
                 result, condition = check_values(up_steady, max_up, min_up, max_down, min_down) ##"moderate" determined
        elif abs((up_steady + down_steady)/2) >= 3: 
                 result, condition = check_values(up_steady, max_up, min_up, max_down, min_down) ##"moderate" determined
 
    return result , condition



