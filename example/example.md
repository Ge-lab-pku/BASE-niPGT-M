
# This is an example of our model

## Dataset

Let's demonstrate with simulated data, which you can find in the **example** folder.We suppose disease lies on chr1, pos=398036 and chromosome F2 carries the disease.Few SNPs arround that position is selected to demonstrate the data.

+ L1_chr1.vcf.txt: the joint analysis data

```plain
,CHROM,POS,ID,REF,ALT,QUAL,FILTER,FORMAT,SAMPLE,MOTHER,FATHER
34,chr1,378320,35,T,G,1000,PASS,GT:AD:DP:GQ:PL,"1/1:0,7:7:21:251,21,0","1/1:0,26:26:78:782,78,34","1/0:0,12:12:36:373,36,0"  
35,chr1,395536,36,G,A,1000,PASS,GT:AD:DP:GQ:PL,"1/1:0,11:11:33:398,33,0","1/1:0,26:26:78:782,78,35","1/0:0,12:12:36:373,36,0"  
36,chr1,396064,37,A,G,1000,PASS,GT:AD:DP:GQ:PGT:PID:PL:PS,"1/1:1,5:6:19:1|0:80396058_CA_C:19,0,20:80396058","1/1:0,26:26:78:782,78,36","1/0:0,12:12:36:373,36,0"  
37,chr1,398036,38,G,A,1000,PASS,GT:AD:DP:GQ:PGT:PID:PL:PS,"1/1:0,1:1:30:.:.:0,30,343:.","1/1:0,26:26:78:782,78,37","1/0:0,12:12:36:373,36,0"  
38,chr1,401569,39,G,A,1000,PASS,GT:AD:DP:GQ:PL,"1/1:0,11:11:33:365,33,0","1/1:0,26:26:78:782,78,38","1/0:0,12:12:36:373,36,0"  
39,chr1,404748,40,G,C,1000,PASS,GT:AD:DP:GQ:PL,"1/1:0,7:7:21:249,21,0","1/1:0,26:26:78:782,78,39","1/0:0,12:12:36:373,36,0"  
40,chr1,406271,41,A,G,1000,PASS,GT:AD:DP:GQ:PL,"1/1:0,5:5:15:155,15,0","1/1:0,26:26:78:782,78,40","1/0:0,12:12:36:373,36,0"  
41,chr1,413342,42,A,C,1000,PASS,GT:AD:DP:GQ:PL,"0/1:10,8:18:27:0,27,405","1/1:0,26:26:78:782,78,41","1/0:0,12:12:36:373,36,0"
```

+ L1_chr1.hap.txt: the phased data

```plain
,CHROM,POS,REF,ALT,FATHER,MOTHER,QUAL,FILTER
27,chr1,368740,G,C,1|0,1|1,50,PASS
28,chr1,371446,T,C,1|0,1|1,50,PASS
29,chr1,378320,T,G,1|0,1|1,50,PASS
30,chr1,395536,G,A,1|0,1|1,50,PASS
31,chr1,398036,G,A,1|0,1|1,50,PASS
32,chr1,401569,G,A,1|0,1|1,50,PASS
33,chr1,404748,G,C,1|0,1|1,50,PASS
34,chr1,406271,A,G,1|0,1|1,50,PASS
35,chr1,413342,A,C,1|0,1|1,50,PASS
36,chr1,417103,G,A,1|0,1|1,50,PASS
37,chr1,419910,C,G,1|0,0|0,50,PASS
```

From the data we suppose F1 is inherented at pos 398036. But SNP at 374776 and 413342 seems to stand against that suppose.  

Let's use our model to solve that problem, make sure that all the package is installed before running our model.  

```bash
cd DashM_scripts
python DashM_Main.py xxx\DashM\example\L1_chr1.vcf.txt xxx\DashM\example\L1_chr1.hap.txt 1 398036 FATHER MOTHER SAMPLE T-SCM father
```

After running DashM_Main.py, the loglikelihood ratio can be found in  **SAMPLE_father.svg**  
