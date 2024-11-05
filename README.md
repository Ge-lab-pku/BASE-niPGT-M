# **BASE-niPGT-M: BAyesian linkage analySis mEthod for non-invasive Preimplantation Genetic Testing of Monogenic disorders**

# Table of Contents

+ Background
+ Install
+ Overview of the code
+ Prepare datasets
+ Run DashM_Main
+ Outputs

# Background

Preimplantation Genetic Testing(PGT) is of great importance in In Vitro Fertilization. Non-invasive PGTM(ni-PGTM) is a class of PGT techniques that has great potential, but is limited by low accuracy rates. In order to enhance the accuracy of niPGT-M from the data obtained from spent culture medium, we present a computational method BASE-niPGT-M.  This approach addresses the challenges posed by MCC and LFLoH in SCM samples. The Bayesian model takes the sequencing data of SCM and the phased haplotypes of the parents, as input. Samples with extremely low quality are filtered out.

Several parameters were estimated before calculating the likelihood ratio: the sequencing error rate, MCC rate and haplotype status. MCC rate was defined as the ratio of DNA fragments sourcing from the maternal chromosome. The haplotype status was defined as the true parental source of the DNA at each SNP and was divided into ‘Paternal Chromosome Only’, ‘Maternal Chromosome Only’, or ‘Both Parental Chromosome’. Parental chromosome only, for instance, suggested that the DNA detected at this SNP solely came from the father. MCC rate and haplotype status were repeatedly calibrated until convergence.

To calculate the likelihood, the recombination probability between adjacent SNPs is determined and the likelihood at a single SNP is calculated. The recombination probability can be obtained from a dataset, while the single-SNP likelihood was attained from a binomial model, whose parameters were determined by the estimated sequencing error rate, MCC rate and haplotype status.

Unlike previous approaches for high-quality sequencing data, which often predetermined the number of SNPs to be considered, our model incrementally incorporates SNPs, starting from the disease-causing mutation site. This process allows us to calculate the log-likelihood ratio of inheriting disease-carrying chromosomes versus disease-free chromosomes for each SNP subset, ultimately resulting in a log-likelihood ratio curve. Typically, the curve started at a point where the ordinate was close to zero and eventually stabilized at a plateau as sufficient SNPs were added. This curve, with its distinct characteristics, enables us to classify the results into four categories: Highly Confident, Moderately Confident, Possibly confident, and Undetermined. In the first three categories, we can pinpoint the inherited chromosome.  This information can then be used to execute the disease carrying analysis for the monogenic disease.

# Install

The project mainly use python codes. The Require python packages are as follows:

## Python Packages

+ `scipy`
+ `matplotlib`
+ `numpy`
+ `pandas`
+ `seaborn`

Besides, there are three softwares that may be helpful:

## External Software

+ `bcftools`
+ `tabix`
+ `R`

# Overview of the code

The main code is in DashM_Main.py,  which takes family sequencing information, parents phasing information, and sample information as input, and outputs an image of the log-likelihood curve and a judgment of whether the sample is pathogenic or not. The image will be saved in **figure** file and the judgement will be print in terminal.

# Prepare datasets

The two datasets we need are:

+ `joint_calling_ann.txt`: the joint analysis of the entire family
+ `genome.phasedHap.txt`: the phased data of the parents of the family

## Preprocess for joint_calling_ann

Usually the joint analysis file is of vcf form and zipped:

+ `joint_calling_ann_vcf.gz`

The following steps can do a quick decompression.

```bash
bcftools view -r {diseased_chr} {joint_calling_ann.vcf.gz} -O z -o {diseased_chr}.vcf.gz
tabix {diseased_chr}.vcf.gz
Rscript prepare_vcf.R {diseased_chr}.vcf.gz
```

> diseased_chr: the chromosome where the disease is located, `chr1`,`chrX` for instance

## Preprocess for genome.phasedHap

As for `genome.phasedHap.txt`, a few preprocess is needed for follow-up analysis.

```bash
python Phasing_data_preprocess.py {family_id} {genome.phasedHap.txt} {father_name} {mother_name} {diseased_chr}
```

> family_id: the id of the family ready for the ni-PGTM. `1`,`2`,`3` for instance. If there is only one family, name it whatever you like
> father_name/mother_name: father's id and mother's id in the family

# Run DashM_Main

Before running **DashM_Main.py** , the following parameters is needed to be specified:

```python
raw_data_name           # Output file generated by prepare_vcf.R
phased_data_name        # Output file generated by Phasing_data_preprocess.py
chrom                   # The chromosome where the disease-causing variant is located. This can be an autosome (e.g., 1, 2, etc.) or a sex chromosome ('X', 'Y')
true_pos                # The exact position of the disease-causing variant on the chromosome
fname, mname            # The names of the father and mother in the dataset
sample                  # The sample name corresponding to the individual being analyzed
cm_type                 # The sample type: T-SCM (cultural media) or D-SCM (cultural media of discarded embryos)
parent                  # Indicates the parent who carries the disease-causing variant. Set this to either father or mother
sex                     # The sex of the sample. Set this to "Male" or "Female" for sex chromosomes, and None for autosomes
```

Running **DashM_Main.py** also need the common SNP file from https://ftp.ncbi.nih.gov/snp/organisms/human_9606_b151_GRCh37p13/VCF/00-common_all.vcf.gz

Then run DashM_Main.py

```bash
cd DashM_scripts
python DashM_Main.py {raw_data_name} {phase_data_name} {chrom} {true_pos} {fname} {mname} {sample} {cm_type} {parent} {sex}
```

# Outputs

Our result will be shown in the terminal. For Further research, you can use the loglikelihood ratio curve saved in 'figure'
The confidence levels are categorized based on the following standards:

## Highly Confident

+ The average steady value, which is the average of upstream and downstream steady values, is ≥ 5,  and there are no differing signs.
+ At least one steady value ≥ 1.

## Moderately Confident

+ The maximum value with differing signs, defined as the maximum number of differing signs in the curve (considering absolute values), is < 5, while the average steady value is ≥ 3.
+ At least one steady value ≥ 1.

## Possibly Confident

+ All other samples that are not categorized as Highly Confident, Moderately Confident, or Undetermined

## Undetermined

1. Steady values have differing signs.
2. Either side has a steady value < 0.5, indicating low confidence evidence.
3. **Rule Out** conditions:
    + Error rate ≥ 0.2.
    + Low GQ  rate ≥ 0.4.
    + Original MCC rate or final MCC rate:
        + ≥ 0.75 or ≤ -0.5 for cultural media.
        + ≥ 0.85 or ≤ -0.7 for cultural media of discarded embryos.
    + SNP density ≤ 0.001.

# Disclosure
L.H., H.G., R.Z. and X.S.X. are coinventors on patent application PCT/CN2023/125605 that includes BASE-niPGT-M.
