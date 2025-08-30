# **BASE-niPGT-M: BAyesian linkage analySis mEthod for non-invasive Preimplantation Genetic Testing of Monogenic disorders**

# Table of Contents

+ Background
+ Install
+ Pipeline Overview

# Background

Preimplantation Genetic Testing(PGT) is of great importance in In Vitro Fertilization. Non-invasive PGTM(ni-PGTM) is a class of PGT techniques that has great potential, but is limited by low accuracy rates. In order to enhance the accuracy of niPGT-M from the data obtained from spent culture medium, we present a computational method BASE-niPGT-M.  This approach addresses the challenges posed by MCC and LFLoH in SCM samples. The Bayesian model takes the sequencing data of SCM and the phased haplotypes of the parents, as input. Samples with extremely low quality are filtered out.

Several parameters were estimated before calculating the likelihood ratio: the sequencing error rate, MCC rate and haplotype status. MCC rate was defined as the ratio of DNA fragments sourcing from the maternal chromosome. The haplotype status was defined as the true parental source of the DNA at each SNP and was divided into ‘Paternal Chromosome Only’, ‘Maternal Chromosome Only’, or ‘Both Parental Chromosome’. Parental chromosome only, for instance, suggested that the DNA detected at this SNP solely came from the father. MCC rate and haplotype status were repeatedly calibrated until convergence.

To calculate the likelihood, the recombination probability between adjacent SNPs is determined and the likelihood at a single SNP is calculated. The recombination probability can be obtained from a dataset, while the single-SNP likelihood was attained from a binomial model, whose parameters were determined by the estimated sequencing error rate, MCC rate and haplotype status.

Unlike previous approaches for high-quality sequencing data, which often predetermined the number of SNPs to be considered, our model incrementally incorporates SNPs, starting from the disease-causing mutation site. This process allows us to calculate the log-likelihood ratio of inheriting disease-carrying chromosomes versus disease-free chromosomes for each SNP subset, ultimately resulting in a log-likelihood ratio curve. Typically, the curve started at a point where the ordinate was close to zero and eventually stabilized at a plateau as sufficient SNPs were added. This curve, with its distinct characteristics, enables us to classify the results into four categories: Highly Confident, Moderately Confident, Possibly confident, and Undetermined. In the first three categories, we can pinpoint the inherited chromosome.  This information can then be used to execute the disease carrying analysis for the monogenic disease.

# Install

The following external tools are required:

## External software and package versions tested successfully
+ `nextflow` (v24.04.4.5917)
+ `bcftools` (v1.9)
+ `tabix` (v1.15.1)
+ `python` (v3.11.2) 
+ `R` (v4.2.0) with the `vcfR` package 

This project is primarily written in Python. The required Python packages are:

## Python Packages

+ `scipy`
+ `matplotlib`
+ `numpy`
+ `pandas`
+ `seaborn`
+ `openpyxl`

#### Installation Time

The full installation process (including downloading dependencies) typically takes **5 to 10 minutes** on a system.

# Pipeline Overview

This nextflow pipeline  takes family joint variant  information, Parental phasing information, and sample information as input.

It outputs:
- A log-likelihood curve visualization (**SVG file**)
- A pathogenicity assessment **report** (**XLSX file**)

Core functionality is implemented in `bin/DashM_Main.py`.


## Dataset Preparation
sample_sheet.tsv : 

|  familyID | father  |  mother | jointvcf  | phasedvcf  | chrom  | position  | sample  | CMtype  | parent  | sex |
| ------------ | ------------ | ------------ | ------------ | ------------ | ------------ | ------------ | ------------ | ------------ | ------------ |------------ |
|FAM001|FATHER  | MOTHER  | /data/joint_calling_ann.vcf.gz  | /data/pahsed.vcf.gz  | 1  | 436546  | SAMPLE|  T-SCM |  mother | unknown  |
|FAM001| FATHER  |   MOTHER |  /data/joint_calling_ann_X.vcf.gz |  /data/pahsed_X.vcf.gz |   X  | 32716070  | SAMPLE|  T-SCM |  mother | Female|


**Column Descriptions**:
- `familyID`: Target family identifier
- `father`/`mother`: Parental sample names
- `jointvcf`: Joint variant calling results (VCF)
- `phasedvcf`: Parental phased haplotypes (VCF)
- `chrom`: Chromosome containing variant (autosomes 1-22, X, Y)
- `position`: Variant genomic position
- `sample`: Analyzed sample name
- `CMtype`: Sample type (`T-SCM`  or `D-SCM` )
- `parent`: Carrier parent (`father` or `mother`)
- `sex`: Sample sex (`Male`, `Female`, or `unknown` for autosomes)


**Required Resource**:  
The pipeline requires `data/00-common_all.vcf.txt.gz`, derived from:  
https://ftp.ncbi.nih.gov/snp/organisms/human_9606_b151_GRCh37p13/VCF/00-common_all.vcf.gz  
*(Processed to extract SNP information)*


## Execution Command
```bash
nextflow run main.nf \
  --sample_sheet sample_sheet.tsv \
  -resume \
  --outdir results
```

#### Runtime Performance

*   **Per-sample runtime:** ~6 minutes
*   **Tested system:** Red Hat Enterprise Linux 7 (RHEL 7 / CentOS 7)

#### The confidence levels are categorized based on the following standards:

#### *Highly Confident*

+ The average steady value, which is the average of upstream and downstream steady values, is ≥ 5,  and there are no differing signs.
+ At least one steady value ≥ 1.

#### *Moderately Confident*

+ The maximum value with differing signs, defined as the maximum number of differing signs in the curve (considering absolute values), is < 5, while the average steady value is ≥ 3.
+ At least one steady value ≥ 1.

#### *Possibly Confident*

+ All other samples that are not categorized as Highly Confident, Moderately Confident, or Undetermined

#### *Undetermined*

1. Steady values have differing signs.
2. Either side has a steady value < 0.5, indicating low confidence evidence.
3. ***Rule Out*** conditions:
    + Error rate ≥ 0.2.
    + Low GQ  rate ≥ 0.4.
    + Original MCC rate or final MCC rate:
        + ≥ 0.75 or ≤ -0.5 for cultural media.
        + ≥ 0.85 or ≤ -0.7 for cultural media (T-SCM) of discarded embryos(D-SCM).
    + SNP density ≤ 0.001.

# Disclosure
L.H., H.G., R.Z. and X.S.X. are coinventors on patent application PCT/CN2023/125605 that includes BASE-niPGT-M.
