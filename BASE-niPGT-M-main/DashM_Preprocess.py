#!/usr/bin/env python
# coding: utf-8

'''
    This document defines functions to calculate the error rate, MCC rate and so on.

'''


from DashM_Basic import *


def compute_err_rate(data, chrom, pos, window_size, output = True):
    
    '''
    Computing error rate.
    window_size: the physical distance of the farest SNP we take into consideration from the disease causal mutation site.
    
    '''
    
    if window_size == "Whole_Genome": # Calculate the whole genome error rate.
        f00m00 = extract(data,chrom,"0/0","0/0","all","all",filtered = True)
        f11m11 = extract(data,chrom,"1/1","1/1","all","all",filtered = True)
    else:
        f00m00 = extract(data,chrom,"0/0","0/0","all","all",filtered = True)
        f11m11 = extract(data,chrom,"1/1","1/1","all","all",filtered = True)
        f00m00 = f00m00[(f00m00["POS"] <= pos + window_size // 2)&(f00m00["POS"] >= pos - window_size // 2)]
        f11m11 = f11m11[(f11m11["POS"] <= pos + window_size // 2)&(f11m11["POS"] >= pos - window_size // 2)]  

    n_error = 0 # Keep record of the number of SNPs where at least one erroneous read occurs.
    error_reads = [] # Keep record of the erroneous reads.
    err = [0,0]
    for x in f00m00["sAD"]:
        if isinstance(x,str): # In some raw data, the AD will be a string rather than a list of length 2.
            x = eval(x)
            
        err[0] += x[0]
        err[1] += x[1]
        
        if x[1] >= 1:
            n_error += 1
            error_reads.append([x[1], x[0] + x[1]]) 
            # The first element is the erroneous read and the second is the total read
        
    for x in f11m11["sAD"]:
        if isinstance(x,str):
            x = eval(x)
            
        err[0] += x[1]
        err[1] += x[0]
        
        if x[0] >= 1:
            n_error += 1
            error_reads.append([x[0], x[0] + x[1]]) 
            # The first element is the erroneous read and the second is the total read
            
    if output:
        print("The window size: ", window_size)
        print("Total number of SNP with same homogeneous parental genotype is ", f00m00.shape[0] + f11m11.shape[0])
        print("Total number of SNP with at least one erroneous read is ", n_error)
        #print("The details of the erroneous reads:", error_reads)
        
    if sum(err) == 0:
        err_rate = 0
    else:
        err_rate = err[1] / sum(err)
    print("Error rate estimated is ",err_rate, "\n")
    return err_rate


def calculate_density(data, chrom, pos, window_size, ref_filename, ref_sep, chunksize, output = True):
    
    '''
    Return the average interval between SNPs in the adjacent area of the disease causal mutation site or in the whole genome. Additionally, we return the ratio of the number of SNPs detected versus that of human common SNPs.
    
    ref_filename: the name of the file which contains the reference data of human common SNPs.
    
    ref_sep: the separation signal in the reference data.
    
    chunksize: the chunksize when we read the Common_SNPs.txt. This file is large and hence, we read it by chunks.
    
    '''
    
    result = {}
    
    if window_size == "Whole_Genome":
        new_data = data
        result["n_within"] = new_data.shape[0]
        result["ave_interval"] = (max(new_data["POS"]) - min(new_data["POS"])) / new_data.shape[0]
        
    else:
        new_data = data[(data["POS"] <= pos + window_size // 2)&((data["POS"] >= pos - window_size // 2))]
        
        result["n_within"] = new_data.shape[0]
        if new_data.shape[0] == 0:
            print("In ", window_size, "base pairs round ", pos, " there are no SNPs detected.")
            result["ave_interval"] = np.Inf
        else:
            result["ave_interval"] = window_size / new_data.shape[0]
        
    # Compute the relative ratio of the number of SNPs versus the total number of human common SNPs.
    print("Comparing the SNPs with the human common SNPs.....")
    for chunk in pd.read_csv(ref_filename, sep = ref_sep, chunksize=chunksize):
        if chrom != 23:
            ref_data = chunk[chunk["CHROM"] == chrom]
        else:
            ref_data = chunk[chunk["CHROM"] == "X"]
        if ref_data.shape[0] == 0:
            continue
        else:
            if window_size == "Whole_Genome":
                pass
            else:
                ref_data = ref_data[(ref_data["POS"] <= pos + window_size // 2)&((ref_data["POS"] >= pos - window_size // 2))]
            ref_data = ref_data[ref_data.apply(lambda x:(len(x["ALT"]) == 1) & (len(x["REF"]) == 1), axis = 1)]
            if ref_data.shape[0] == 0:
                continue
                
        ref_pos = list(ref_data["POS"])
        n_common_within = sum([x in ref_pos for x in list(new_data["POS"])])
        result["n_common"] = result.get("n_common", 0) + len(ref_pos)
        result["n_common_within"] = result.get("n_common_within", 0) + n_common_within
        
    result["n_common_ratio"] = result["n_common_within"] / result["n_common"] if result["n_common"] != 0 else 0
    result["n_uncommon_detected"] = result["n_within"] - result["n_common_within"]

    if output:
        print("There are %d SNPs detected within the required region.\nThe average interval between SNPs is %f. \nThere are %d common SNPs in this region in total and %d of them are detected. \nThe ratio of human common SNPs that are detected is %f. \nThere are %d detected SNPs belonging to the uncommon human SNPs, which is %f of the total SNPs detected.\n" % (result["n_within"], result["ave_interval"], result["n_common"], result["n_common_within"], result["n_common_ratio"], result["n_uncommon_detected"], result["n_uncommon_detected"] / result["n_within"]))
    return result

    
def compute_mcc_rate(data, chrom, pos, window_size, sex = None, output = True):
    
    '''
    Computing mcc_rate.
    '''
    
    if (chrom != 23) or (chrom == 23 and sex == "Female"):
        f00m11 = extract(data,chrom,"0/0","1/1","all","all",filtered = True)
        f11m00 = extract(data,chrom,"1/1","0/0","all","all",filtered = True)
        
        if window_size == "Whole_Genome":
            pass
        else:
            f00m11 = f00m11[(f00m11["POS"] <= pos + window_size // 2)&(f00m11["POS"] >= pos - window_size // 2)]
            f11m00 = f11m00[(f11m00["POS"] <= pos + window_size // 2)&(f11m00["POS"] >= pos - window_size // 2)]
            
        mcc = [0,0]
        n_mcc = 0 # Keep record of the number of total number of SNPs with fGT = 00, mGT = 11 or fGT = 11, mGT = 00.
        mcc_reads = [] # keep record of the (reads from mother, total reads) pairs.
        
        for x in f00m11["sAD"]:
            if isinstance(x,str):
                x = eval(x)
            mcc[0] += x[0]
            mcc[1] += x[1]
            n_mcc += 1
            mcc_reads.append((x[1], x[1] + x[0]))
            
        for x in f11m00["sAD"]:
            if isinstance(x,str):
                x = eval(x)
            mcc[0] += x[1]
            mcc[1] += x[0]
            n_mcc += 1
            mcc_reads.append((x[0], x[1] + x[0]))

        if sum(mcc) == 0:
            mcc_rate = 0
        else:
            mcc_rate = mcc[1]/sum(mcc) * 2 - 1
            
        if output:
            print("The disease causal mutation site is on the autosomal chromosomes. \nThe number of SNPs with fGT = 00, mGT = 11 or fGT = 11, mGT = 00 is %d" % n_mcc)
            print("The original estimated mcc_rate is %f .\n" % mcc_rate)
                 
    elif chrom == 23 and sex == "Male" and pos > 5e7:
        print("The disease causal mutatin site is on the chrX and not on the chrY.")
        X_data = extract(data,chrom,"all","0/1","all","all",filtered = True)
        
        if window_size == "Whole_Genome":
            X_data = X_data[X_data["POS"] > 5e7]
            print("We can not compute the mcc rate on chromosome X as a whole. Now we compute t on the part of chrX that has no corresponding part on chrY.")
        else:
            X_data = X_data[((X_data["POS"] <= pos + window_size // 2)&(X_data["POS"] >= pos - window_size // 2))]
        
        mcc = [0,0]
        n_mcc = 0 # Keep record of the number of total number of SNPs with fGT = 00, mGT = 11 or fGT = 11, mGT = 00.
        mcc_reads = [] # keep record of the (reads from mother chain 1, reads from mother chain 2) pairs.
        
        for i in range(X_data.shape[0]):
            x = X_data.iloc[i]["sAD"]
            if isinstance(x,str):
                x = eval(x)
            if X_data.iloc[i]["mGT_p"] == "0|1":
                mcc[0] += x[0]
                mcc[1] += x[1]
                n_mcc += 1
                mcc_reads.append((x[0],x[1]))
            elif X_data.iloc[i]["mGT_p"] == "1|0":
                mcc[0] += x[1]
                mcc[1] += x[0]
                n_mcc += 1
                mcc_reads.append((x[0],x[1]))      
        if sum(mcc) == 0:
            mcc_rate = 0
        else:
            mcc_rate = min(mcc[0],mcc[1]) / sum(mcc)
        if output:
            print("The disease causal mutation site is on the autosomal chromosomes. \nThe number of SNPs with mGT = 01 is %d, \nThe (reads from mother chain 1, reads from mother chain 2) pairs are" % n_mcc, mcc_reads)
            print("The original estimated mcc_rate is %f .\n" % mcc_rate)
            
    else:
        print("This part of function is under construction.")
        
    return {"mcc_rate": mcc_rate, "n_mcc": n_mcc, "mcc_reads": mcc_reads}
 
