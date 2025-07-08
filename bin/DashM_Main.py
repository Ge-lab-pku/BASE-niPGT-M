# Main function of calculating Log-Likelihood Ratio of the data

# Version 2.0

import argparse
import numpy as np
import pandas as pd
import sys
import matplotlib.pyplot as plt
from DashM_Basic import *
from DashM_RowLikelihood import *
from DashM_Preprocess import *
from DashM_Label import *
from DashM_LFLoH import *
from DashM_RecomProb import *
from DashM_Curve import *
from DashM_Judge import *
import os
import warnings
warnings.filterwarnings("ignore")
import openpyxl

results_list = []

def parse_args():
    parser = argparse.ArgumentParser(description="DashM Analysis Pipeline")

    parser.add_argument("--raw_data", required=True, help="Input raw VCF file")
    parser.add_argument("--phased_data", required=True, help="Input phased VCF file")
    parser.add_argument("--chrom", required=True, help="Chromosome name")
    parser.add_argument("--position", type=int, required=True, help="Variant position")
    parser.add_argument("--father", required=True, help="Father sample ID")
    parser.add_argument("--mother", required=True, help="Mother sample ID")
    parser.add_argument("--sample", required=True, help="Proband sample ID")
    parser.add_argument("--cm_type", required=True, help="CM type")
    parser.add_argument("--parent", required=True, help="Parental origin")
    parser.add_argument("--sex", default="unknown", help="Sample sex (male/female/unknown)")

    parser.add_argument("--output_dir", default="./output", help="Output directory")
    
    return parser.parse_args()

def DashM_main(raw_data_name, phased_data_name, chrom, true_pos, fname, mname, sample, parent,
               mode, df_save, df_save_path, figure_save_path, output, 
               ADO_thres, GQ_thres, QUAL_thres, recom_rate, data_sep, cal_density,num_clip,
               LFLoH_thres, rec_rate_pd = None, cm_type = None, 
               ref_filename = "00-common_all.vcf.txt.gz", ref_sep = ",", 
               chunksize = 100000, sex = None, max_length = 250):
    
    print("Start testing...., sample =",sample)
    raw_data = pd.read_csv(raw_data_name,sep = data_sep,index_col = 0)
    raw_data = raw_data[(raw_data["QUAL"] >= QUAL_thres)&(raw_data["FILTER"] == "PASS")]
    phased_data = pd.read_csv(phased_data_name,sep = "\t") #,index_col = 0)
    #phased_data = phased_data[(phased_data["FILTER"] == "PASS")&(phased_data["QUAL"] >= QUAL_thres)]

    pos = true_pos
    if chrom in [f'{i}' for i in range(23)]:
        chrom= int(chrom)
    data = data_filter_phased(raw_data,phased_data,fname,mname,sample,chrom,GQ_thres)
    window_size = 2e7

    # Preprocessing: compute the error rate, density, mcc_rate
    # if mode == 0 or mode == 1 or mode == 2 or mode == 3: this is for removing one or two SNPs closest to the disease-causal mutation site
        
        
    # truncate the sample_GQ.
    print("There are %d SNPs with GQ <= 3" % data[data["sGQ"] <= 3].shape[0])
    Ratio_Low_GQ=data[data["sGQ"] <= 3].shape[0] / data.shape[0]
    print("The ratio of SNPs with low sample GQ is ", Ratio_Low_GQ)

    sample_GQ_thres = max(np.quantile(data["sGQ"],0.2),3)
    # sample_GQ_thres = 3
    print("The sample GQ threshold is %d" % sample_GQ_thres)
    data = data[data["sGQ"] > sample_GQ_thres]

    err_rate = compute_err_rate(data, chrom, pos, window_size, output)
    if err_rate == 0:
        err_rate = 0.001
    
    if cal_density:
        density_info = calculate_density(data, chrom, pos, window_size, ref_filename, ref_sep, chunksize, output)    
        SNP_Density=density_info["n_common_ratio"]
    
    else:
        SNP_Density=1
    mcc_info = compute_mcc_rate(data, chrom, pos, window_size, sex, output)
    # Label the data
    truncated_data = label(data, chrom, true_pos, window_size, ADO_thres, sex, output)


    def is_in_PAR(chrom, pos):
        if(chrom =="X"):
            if( 60000<pos<2699520):
                return 'par1'
            if(154931044<pos<155260560):
                return 'par2'
            return 'non-par'
        elif(chrom =="Y"):
            if(10000<pos<2649520 or 59034050<pos<59363566):
                return 'par'
            return 'non-par'
        else:
            print('check chrom=X or chrom=Y, PAR is the specific regions on chromosome X and Y')
    # Calculate the recombination probability
    if rec_rate_pd is None:
        try:
            if chrom != "X" and chrom !="Y":
                #file_path = f'Recom_Prob/genetic_map_chr{chrom}_combined_b37.20140701.txt'
                file_path = f'Recom_Prob/genetic_map_GRCh37_chr{chrom}.txt'
            elif chrom == 'Y' or is_in_PAR(chrom,true_pos)=='non-par':
                file_path = f'Recom_Prob/genetic_map_GRCh37_chrX.txt'
                #file_path = f'Recom_Prob/genetic_map_chrX_combined_b37.20140701.txt'
            elif is_in_PAR(chrom,true_pos)=='par1':
                file_path = f'Recom_Prob/genetic_map_GRCh37_chrX_par1.txt'
            elif is_in_PAR(chrom,true_pos)=='par2':
                file_path = f'Recom_Prob/genetic_map_GRCh37_chrX_par2.txt'
            else:
                print('ERROR, THE CHROM OR POSITION IS NOT IN REF FILE')
            rec_rate_pd = pd.read_csv(file_path, sep='\s+')
        except FileNotFoundError: 
            print(f"Error: File '{file_path}' does not exist.")

    recom_probs = all_recom_prob(truncated_data, chrom, rec_rate_pd, recom_rate, output)

    # Calibrating the mcc rate
    Original_MCC_Rate = mcc_info["mcc_rate"]

    if Original_MCC_Rate <= 0:
        mcc_rate = 0.01
        flag = True
    else:
        mcc_rate=Original_MCC_Rate
        flag = False   
    calibrate_info = calibrate(data, chrom, true_pos, window_size, sample, parent, err_rate, mcc_rate, ADO_thres, LFLoH_thres, recom_probs, df_save, df_save_path, flag, rec_rate_pd, recom_rate, sex = None)
    df = calibrate_info["data"]
    mcc_rate = calibrate_info["mcc_rate_final"]
    Final_MCC_Rate=mcc_rate
    

    if (chrom != "X" and chrom != "Y"):
        if parent == "father":
            df = df[df["Final_label"] != "M"]
        elif parent == "mother":
            df = df[(df["Final_label"] != "F")&(df["Final_label"] != "F1") & (df["Final_label"] != "F2")]
        elif parent == "both":
            df = df[df["Final_label"] == "FM"]

    elif(is_in_PAR(chrom,true_pos) !='non-par'):
        Original_MCC_Rate, Final_MCC_Rate=0, 0
        if(chrom == "Y"):
            if(sex =="Female"):
                print('The child is healthy')
            elif(sex == 'Male'):
                print('The child is a carrier or affected')
            else:
                print(f'there is no recombination at location {true_pos} on chrom Y, child is sick if it is a boy and its father is sick')
            return 0
        elif(sex == 'Male'):
            if(parent =='father'):
                print('The child is healthy')
                return 0
            df['Final_label']='M'
            print(1111111)
            mcc_rate = 0

        elif(sex == 'Female'):
            if(parent == 'mother'):
                #df['Final_label']='FM'
                pass
            else:
                print('The child is a carrier or affected')
                return 0
            

    pos = true_pos
    data1 = df[(df["POS"] >= pos - window_size//2)&(df["POS"] <= pos)]
    data2 = df[(df["POS"] <= pos + window_size//2)&(df["POS"] >= pos)]
    data1.sort_values(by = "POS", ascending = False, inplace = True)
    data2.sort_values(by = "POS", ascending = True, inplace = True)

    upper_prob = []
    lower_prob = []
    length = max_length

    if parent == "father":
        # Only use the SNPs with paternal genotype is heterozygous (when the disease causal mutation site is on the paternal chromosome)
        data1 = data1[data1["fGT"] == "0/1"]
        data2 = data2[data2["fGT"] == "0/1"]

    elif parent == "mother":
        # Only use the SNPs with maternal genotype is heterozygous (when the disease causal mutation site is on the maternal chromosome)
        data1 = data1[data1["mGT"] == "0/1"]
        data2 = data2[data2["mGT"] == "0/1"]

    else:
        print("Parent error", parent)

    if mode == 0:
        pass
    elif mode == 1:
        data1 = data1.iloc[1:data1.shape[0]]
        data2 = data2.iloc[1:data2.shape[0]]
    elif mode == 2:
        data1 = data1.iloc[2:data1.shape[0]]
        data2 = data2.iloc[2:data2.shape[0]]

    recom_probs.update(all_recom_prob(data1, chrom, rec_rate_pd, recom_rate, output = False))
    recom_probs.update(all_recom_prob(data2, chrom, rec_rate_pd, recom_rate, output = False))
    
    # Begin to iterate calculating the likelihood

    print("Begin to iterately calculating the likelihood.......")
    for j in range(min(data1.shape[0],length)):
        data11 = data1.iloc[:j+1,:]
        data11.sort_values(by = "POS",ascending = True,inplace = True)

        likelihood1 = {x:[0,0,0,0] for x in data11["POS"]}

        row = data11.iloc[0]
        if row["Final_label"] == "FM":
            crt_label = "FM"
        elif row["Final_label"] in ["F","F1","F2"]:
            crt_label = "F"
        elif row["Final_label"] == "M":
            crt_label = "M"
        start_pos = data11.iloc[0]["POS"]
        tmp_likelihood = [log(row_likelihood(row, mcc_rate, crt_label, err_rate)[0]),
                                 log(row_likelihood(row, mcc_rate, crt_label, err_rate)[1]),
                                 log(row_likelihood(row, mcc_rate, crt_label, err_rate)[2]),
                                 log(row_likelihood(row, mcc_rate, crt_label, err_rate)[3])]
        
        # clip
        tmp_max = max(tmp_likelihood)
        likelihood1[start_pos] = [max(tmp_max - num_clip, x) for x in tmp_likelihood]


        for i in range(1,data11.shape[0]):
            row = data11.iloc[i]
            current_pos = row["POS"]
            last_pos = data11.iloc[i-1]["POS"]
            if (current_pos,last_pos) in recom_probs.keys():
                recom_prob = recom_probs[(current_pos,last_pos)]
            else:
                recom_prob = calculate_recom_prob(chrom,current_pos,last_pos)
            no_recom_prob = 1 - recom_prob
            if row["Final_label"] == "FM":
                crt_label = "FM"
            elif row["Final_label"] in ["F","F1","F2"]:
                crt_label = "F"
            elif row["Final_label"] == "M":
                crt_label = "M"
            
            tmp_likelihood = [0,0,0,0]
            tmp_likelihood[0] = (log(row_likelihood(row, mcc_rate, crt_label, err_rate)[0]) + likelihood1[last_pos][0]
                        + log(no_recom_prob**2 
                        + no_recom_prob * recom_prob * np.exp(likelihood1[last_pos][1] - likelihood1[last_pos][0])
                        + no_recom_prob * recom_prob * np.exp(likelihood1[last_pos][2] - likelihood1[last_pos][0])
                        + recom_prob**2 * np.exp(likelihood1[last_pos][3] - likelihood1[last_pos][0])))
            tmp_likelihood[1] = (log(row_likelihood(row, mcc_rate, crt_label, err_rate)[1]) + likelihood1[last_pos][1]
                        + log(no_recom_prob**2 
                        + no_recom_prob * recom_prob * np.exp(likelihood1[last_pos][0] - likelihood1[last_pos][1])
                        + no_recom_prob * recom_prob * np.exp(likelihood1[last_pos][3] - likelihood1[last_pos][1])
                        + recom_prob**2 * np.exp(likelihood1[last_pos][2] - likelihood1[last_pos][1])))
            tmp_likelihood[2] = (log(row_likelihood(row, mcc_rate, crt_label, err_rate)[2]) + likelihood1[last_pos][2]
                        + log(no_recom_prob**2 
                        + no_recom_prob * recom_prob * np.exp(likelihood1[last_pos][0] - likelihood1[last_pos][2])
                        + no_recom_prob * recom_prob * np.exp(likelihood1[last_pos][3] - likelihood1[last_pos][2])
                        + recom_prob**2 * np.exp(likelihood1[last_pos][1] - likelihood1[last_pos][2])))
            tmp_likelihood[3] = (log(row_likelihood(row, mcc_rate, crt_label, err_rate)[3]) + likelihood1[last_pos][3]
                        + log(no_recom_prob**2 
                        + no_recom_prob * recom_prob * np.exp(likelihood1[last_pos][1] - likelihood1[last_pos][3])
                        + no_recom_prob * recom_prob * np.exp(likelihood1[last_pos][2] - likelihood1[last_pos][3])
                        + recom_prob**2 * np.exp(likelihood1[last_pos][0] - likelihood1[last_pos][3])))

            # clip
            tmp_diff = [tmp_likelihood[i] - likelihood1[last_pos][i] for i in range(len(tmp_likelihood))]
            tmp_diff = [max(max(tmp_diff) - num_clip, tmp_diff[i]) for i in range(len(tmp_diff))]
            likelihood1[current_pos] = [likelihood1[last_pos][i] + tmp_diff[i] for i in range(len(tmp_diff))]

        pos1 = data11.iloc[j]["POS"] if data11.shape[0] > 0 else None
        like11 = likelihood1[pos1] if pos1 else [0,0]
        upper_prob.append((j+1,like11))

    for j in range(min(data2.shape[0],length)):
        data22 = data2.iloc[:j+1,:]
        data22.sort_values(by = "POS",ascending = False,inplace = True)

        likelihood2 = {x:[0,0,0,0] for x in data22["POS"]}

        row = data22.iloc[0]
        if row["Final_label"] == "FM":
            crt_label = "FM"
        elif row["Final_label"] in ["F","F1","F2"]:
            crt_label = "F"
        elif row["Final_label"] == "M":
            crt_label = "M"
        start_pos = data22.iloc[0]["POS"]
        tmp_likelihood = [log(row_likelihood(row, mcc_rate, crt_label, err_rate)[0]),
                                 log(row_likelihood(row, mcc_rate, crt_label, err_rate)[1]),
                                 log(row_likelihood(row, mcc_rate, crt_label, err_rate)[2]),
                                 log(row_likelihood(row, mcc_rate, crt_label, err_rate)[3])]
        
        # clip
        tmp_max = max(tmp_likelihood)
        likelihood2[start_pos] = [max(tmp_max - num_clip, x) for x in tmp_likelihood]


        for i in range(1,data22.shape[0]):
            row = data22.iloc[i]
            current_pos = row["POS"]
            last_pos = data22.iloc[i-1]["POS"]
            if (last_pos,current_pos) in recom_probs.keys():
                recom_prob = recom_probs[(last_pos,current_pos)]
            else:
                recom_prob = calculate_recom_prob(chrom,last_pos,current_pos)
            no_recom_prob = 1 - recom_prob
            if row["Final_label"] == "FM":
                crt_label = "FM"
            elif row["Final_label"] in ["F","F1","F2"]:
                crt_label = "F"
            elif row["Final_label"] == "M":
                crt_label = "M"

            tmp_likelihood = [0,0,0,0]
            tmp_likelihood[0] = (log(row_likelihood(row, mcc_rate, crt_label, err_rate)[0]) + likelihood2[last_pos][0]
                        + log(no_recom_prob**2 
                        + no_recom_prob * recom_prob * np.exp(likelihood2[last_pos][1] - likelihood2[last_pos][0])
                        + no_recom_prob * recom_prob * np.exp(likelihood2[last_pos][2] - likelihood2[last_pos][0])
                        + recom_prob**2 * np.exp(likelihood2[last_pos][3] - likelihood2[last_pos][0])))
            tmp_likelihood[1] = (log(row_likelihood(row, mcc_rate, crt_label, err_rate)[1]) + likelihood2[last_pos][1]
                        + log(no_recom_prob**2 
                        + no_recom_prob * recom_prob * np.exp(likelihood2[last_pos][0] - likelihood2[last_pos][1])
                        + no_recom_prob * recom_prob * np.exp(likelihood2[last_pos][3] - likelihood2[last_pos][1])
                        + recom_prob**2 * np.exp(likelihood2[last_pos][2] - likelihood2[last_pos][1])))
            tmp_likelihood[2] = (log(row_likelihood(row, mcc_rate, crt_label, err_rate)[2]) + likelihood2[last_pos][2]
                        + log(no_recom_prob**2 
                        + no_recom_prob * recom_prob * np.exp(likelihood2[last_pos][0] - likelihood2[last_pos][2])
                        + no_recom_prob * recom_prob * np.exp(likelihood2[last_pos][3] - likelihood2[last_pos][2])
                        + recom_prob**2 * np.exp(likelihood2[last_pos][1] - likelihood2[last_pos][2])))
            tmp_likelihood[3] = (log(row_likelihood(row, mcc_rate, crt_label, err_rate)[3]) + likelihood2[last_pos][3]
                        + log(no_recom_prob**2 
                        + no_recom_prob * recom_prob * np.exp(likelihood2[last_pos][1] - likelihood2[last_pos][3])
                        + no_recom_prob * recom_prob * np.exp(likelihood2[last_pos][2] - likelihood2[last_pos][3])
                        + recom_prob**2 * np.exp(likelihood2[last_pos][0] - likelihood2[last_pos][3])))

            # clip
            tmp_diff = [tmp_likelihood[i] - likelihood2[last_pos][i] for i in range(len(tmp_likelihood))]
            tmp_diff = [max(max(tmp_diff) - num_clip, tmp_diff[i]) for i in range(len(tmp_diff))]
            likelihood2[current_pos] = [likelihood2[last_pos][i] + tmp_diff[i] for i in range(len(tmp_diff))]


        pos2 = data22.iloc[j]["POS"] if data22.shape[0] > 0 else None
        like22 = likelihood2[pos2] if pos2 else [0,0]
        lower_prob.append((j+1,like22))

    # Save data and figure
    if df_save:
        data1 = data1.iloc[:min(data1.shape[0],length),:]
        data1["F1M1"] = [x[1][0] for x in upper_prob]
        data1["F1M2"] = [x[1][1] for x in upper_prob]
        data1["F2M1"] = [x[1][2] for x in upper_prob]
        data1["F2M2"] = [x[1][3] for x in upper_prob]

        data1.to_csv(figure_save_path + sample + "_" + parent + "_upper.txt")

        data2 = data2.iloc[:min(data2.shape[0],length),:]
        data2["F1M1"] = [x[1][0] for x in lower_prob]
        data2["F1M2"] = [x[1][1] for x in lower_prob]
        data2["F2M1"] = [x[1][2] for x in lower_prob]
        data2["F2M2"] = [x[1][3] for x in lower_prob]

        data2.to_csv(figure_save_path + sample + "_" + parent + "_lower.txt")
    
    pos = true_pos
    plt.close()
    plt.figure(figsize = (20,10)) ##qg
    if parent == "father":
        #plt.figure(1,figsize = (20,10))
        ax1 = plt.subplot(2,2,1)
        ax2 = plt.subplot(2,2,2)
###qg 3 - 1.5
        if upper_prob:
            plt.sca(ax1)
            plt.plot([(pos - data1.iloc[x[0]-1]["POS"])/1e3 for x in upper_prob],
                 [max(x[1][0],x[1][1]) + 
                  np.log(np.exp(x[1][0] - max(x[1][0],x[1][1]))
                         +np.exp(x[1][1] - max(x[1][0],x[1][1])))
                  - (
                      max(x[1][2],x[1][3]) + 
                  np.log(np.exp(x[1][2] - max(x[1][2],x[1][3]))
                         +np.exp(x[1][3] - max(x[1][2],x[1][3])))
                  )
                  for x in upper_prob],
                linewidth = 3,marker='x',markeredgecolor='black',markersize='2',markeredgewidth=2,
                markerfacecolor = 'none',label = sample)
            plt.legend()
            plt.xlabel("Distance(Kbp) from Upstream Region to \n the Disease Causal Mutation Site.",size = 15,weight = "bold")
            plt.ylabel("Log Likelihood Ratio.",size = 20,weight = "bold")
            #plt.xlim((0,5000))
            plt.xticks(fontsize = 10,weight = "bold")
            plt.yticks(fontsize = 10,weight = "bold")
            ax1.spines['bottom'].set_linewidth(3)
            ax1.spines['left'].set_linewidth(3)

        if lower_prob:
            plt.sca(ax2)
            plt.plot([(data2.iloc[x[0]-1]["POS"] - pos)/1e3 for x in lower_prob],
                     [max(x[1][0],x[1][1]) + 
                  np.log(np.exp(x[1][0] - max(x[1][0],x[1][1]))
                         +np.exp(x[1][1] - max(x[1][0],x[1][1])))
                  - (
                      max(x[1][2],x[1][3]) + 
                  np.log(np.exp(x[1][2] - max(x[1][2],x[1][3]))
                         +np.exp(x[1][3] - max(x[1][2],x[1][3])))
                  )
                  for x in lower_prob],
                    linewidth = 3,marker='x',markeredgecolor='black',markersize="2",markeredgewidth=2,
                    markerfacecolor = 'none',label = sample)
            plt.legend()
            plt.xlabel("Distance(Kbp) from Downstream Region to \n the Disease Causal Mutation Site.",size = 15,weight = "bold")
            plt.ylabel("Log Likelihood Ratio.",size = 20,weight = "bold")
            #plt.xlim((0,5000))
            plt.xticks(fontsize = 10,weight = "bold")
            plt.yticks(fontsize = 10,weight = "bold")
            ax2.spines['bottom'].set_linewidth(3)
            ax2.spines['left'].set_linewidth(3)

        plt.suptitle("Log-likelihood Curve of " + sample + " (Father)" ,fontsize = 20,x=0.5,y=0.92,weight = "bold")
        plt.savefig(figure_save_path + sample + "_"  + str(chrom) + "_" + str(pos) + "_father.svg") #qg
        #plt.show()

    if parent == "mother":
        plt.figure(1,figsize = (20,10))
        ax1 = plt.subplot(2,2,1)
        ax2 = plt.subplot(2,2,2)

        if upper_prob:
            plt.sca(ax1)
            plt.plot([(pos - data1.iloc[x[0]-1]["POS"])/1e3 for x in upper_prob],
                 [max(x[1][0],x[1][2]) + 
                  np.log(np.exp(x[1][0] - max(x[1][0],x[1][2]))
                         +np.exp(x[1][2] - max(x[1][0],x[1][2])))
                  - (
                      max(x[1][1],x[1][3]) + 
                  np.log(np.exp(x[1][1] - max(x[1][1],x[1][3]))
                         +np.exp(x[1][3] - max(x[1][1],x[1][3])))
                  )
                  for x in upper_prob],
                linewidth = 3,marker='x',markeredgecolor='black',markersize='2',markeredgewidth=2,
                markerfacecolor = 'none',label = sample)
            plt.legend()
            plt.xlabel("Distance(Kbp) from Upstream Region to \n the Disease Causal Mutation Site.",size = 15,weight = "bold")
            plt.ylabel("Log Likelihood Ratio.",size = 20,weight = "bold")
            #plt.xlim((0,5000))
            plt.xticks(fontsize = 10,weight = "bold")
            plt.yticks(fontsize = 10,weight = "bold")
            ax1.spines['bottom'].set_linewidth(3)
            ax1.spines['left'].set_linewidth(3)

        if lower_prob:
            plt.sca(ax2)
            plt.plot([(data2.iloc[x[0]-1]["POS"] - pos)/1e3 for x in lower_prob],
                     [max(x[1][0],x[1][2]) + 
                  np.log(np.exp(x[1][0] - max(x[1][0],x[1][2]))
                         +np.exp(x[1][2] - max(x[1][0],x[1][2])))
                  - (
                      max(x[1][1],x[1][3]) + 
                  np.log(np.exp(x[1][1] - max(x[1][1],x[1][3]))
                         +np.exp(x[1][3] - max(x[1][1],x[1][3])))
                  )
                  for x in lower_prob],
                    linewidth = 3,marker='x',markeredgecolor='black',markersize='2',markeredgewidth=2,
                    markerfacecolor = 'none',label = sample)
            plt.legend()
            plt.xlabel("Distance(Kbp) from Downstream Region to \n the Disease Causal Mutation Site.",size = 15,weight = "bold")
            plt.ylabel("Log Likelihood Ratio.",size = 20,weight = "bold")
            #plt.xlim((0,5000))
            plt.xticks(fontsize = 10,weight = "bold")
            plt.yticks(fontsize = 10,weight = "bold")
            ax2.spines['bottom'].set_linewidth(3)
            ax2.spines['left'].set_linewidth(3)

        plt.suptitle("Log-likelihood Curve of " + sample + " (Mother)" ,fontsize = 20,x=0.5,y=0.92,weight = "bold")
        plt.savefig(figure_save_path + sample + "_"  + str(chrom) + "_" + str(pos)  + "_mother.svg") #qg
        #plt.show()

    print(sample + " " + parent + " Finished.")
    #judge
    up_sign,down_sign, max_up, min_up, max_down, min_down,count_up_greater_than_zero,count_up_less_than_zero,count_down_greater_than_zero,count_down_less_than_zero = stable_value_DashM(sample, parent, figure_save_path, eps=0.1)
    print(f'Sample name is {sample}')
    print(f'The steady value of up/down stream are {up_sign}, {down_sign}.')
    print(f'The maximal and minimal value in upstream area are {max_up}, {min_up}.')
    print(f'The maximal and minimal value in downstream area are {max_down}, {min_down}.')
    print(f"Count of values greater than or equal to 0 and less than or equal to 0 of upstream: {count_up_greater_than_zero}, {count_up_less_than_zero}.")
    print(f"Count of values greater than or equal to 0 and less than or equal to 0 of downstream: {count_down_greater_than_zero}, {count_down_less_than_zero}.")
    sample_dict={'name':sample,'Up_Steady_Value':up_sign,'Down_Steady_Value':down_sign,'Upstream_Max_Value':max_up,'Upstream_Min_Value':min_up,'Downstream_Max_Value':max_down, 'Downstream_Min_Value':min_down,'Ratio_Low_GQ':Ratio_Low_GQ,'Error_Rate':err_rate,'Original_MCC_Rate':Original_MCC_Rate,'Final_MCC_Rate':Final_MCC_Rate,'SNP_Density':SNP_Density,'Up_positive_number':count_up_greater_than_zero,'Up_negative_number':count_up_less_than_zero,'Down_positive_number':count_down_greater_than_zero,'Down_negative_number':count_down_less_than_zero,'CM_type':cm_type}
    result,condition=evaluate_sample(sample_dict)
    inheritance_result = ""
    if isinstance(up_sign, (int, float)) and isinstance(down_sign, (int, float)):
        if(up_sign>=0 and down_sign>=0):
            if(parent=="mother"):
                inheritance_result = f"M1 at {true_pos}, chr{chrom}"
                print(f"We suppose the sample inherit M1 at {true_pos}, chr{chrom}")
            elif(parent == "father"):
                inheritance_result = f"F1 at {true_pos}, chr{chrom}"
                print(f"We suppose the sample inherit F1 at {true_pos}, chr{chrom}")
        elif(up_sign<=0 and down_sign<=0):
            if(parent=="mother"):
                inheritance_result = f"M2 at {true_pos}, chr{chrom}"
                print(f"We suppose the sample inherit M2 at {true_pos}, chr{chrom}")
            elif(parent == "father"):
                inheritance_result = f"F2 at {true_pos}, chr{chrom}"
                print(f"We suppose the sample inherit F2 at {true_pos}, chr{chrom}")
        else:
            inheritance_result = "Undetermined"
            print(f"We can't determine the inheritance of the sample at {true_pos}, chr{chrom}")
    else:
        inheritance_result = "Undetermined"
        print(f"We can't determine the inheritance of the sample at {true_pos}, chr{chrom}")    
    print(f'The confidence level is {result}')

    sample_dict = {
        'name': sample,
        'chrom': chrom,
        'true_pos': true_pos,
        'fname': fname,
        'mname': mname,
        'cm_type': cm_type,
        'parent': parent,
        'Up_Steady_Value': up_sign,
        'Down_Steady_Value': down_sign,
        'Upstream_Max_Value': max_up,
        'Upstream_Min_Value': min_up,
        'Downstream_Max_Value': max_down,
        'Downstream_Min_Value': min_down,
        'Ratio_Low_GQ': Ratio_Low_GQ,
        'Error_Rate': err_rate,
        'Original_MCC_Rate': Original_MCC_Rate,
        'Final_MCC_Rate': Final_MCC_Rate,
        'SNP_Density': SNP_Density,
        'Up_positive_number': count_up_greater_than_zero,
        'Up_negative_number': count_up_less_than_zero,
        'Down_positive_number': count_down_greater_than_zero,
        'Down_negative_number': count_down_less_than_zero,
        'confidence_level': result,
        'inheritance_result': inheritance_result
    }

    results_list.append(sample_dict)

if(__name__ == "__main__"):
    args = parse_args()
    print('begin')
    raw_data_name=args.raw_data
    phased_data_name=args.phased_data
    chrom=args.chrom
    true_pos=args.position
    fname=args.father
    mname=args.mother
    sample=args.sample
    parent=args.parent
    cm_type=args.cm_type
    sex = None if args.sex == 'unknown' else args.sex
    mode=0

    output_dir = args.output_dir
    result_dir = output_dir + "/" 
    df_save, df_save_path, figure_save_path= True, result_dir, result_dir
    os.makedirs(df_save_path, exist_ok=True)
    os.makedirs(figure_save_path, exist_ok=True)

    print('savepath=',df_save_path)
    output = True
    GQ_thres = 10
    QUAL_thres = 30
    recom_rate = 1e-8
    data_sep = ","
    cal_density=True
    num_clip, LFLoH_thres=5,{"F":0.9, "F1":0.9, "F2":0.9, "M":0.99, "FM":0.9}
    ADO_thres = 2 
    DashM_main(raw_data_name, phased_data_name, chrom, true_pos, fname, mname, sample, parent,
               mode, df_save, df_save_path, figure_save_path, output, 
               ADO_thres, GQ_thres, QUAL_thres, recom_rate, data_sep, cal_density,num_clip,
               LFLoH_thres, rec_rate_pd = None,cm_type = cm_type, 
               ref_filename = "00-common_all.vcf.txt.gz", ref_sep = ",", 
               chunksize = 100000, sex=sex , max_length = 250)
    if results_list:
        df = pd.DataFrame(results_list)
        output_excel_path = os.path.join(df_save_path, sample + "_"  + str(chrom) + "_" + str(true_pos) + "_" +parent  + '_results.xlsx')
    
        if os.path.exists(output_excel_path):
            existing_df = pd.read_excel(output_excel_path)
            combined_df = pd.concat([existing_df, df], ignore_index=True)
        else:
            combined_df = df
    
        combined_df.to_excel(output_excel_path, index=False)
        print(f"\nResults saved to: {output_excel_path}")

