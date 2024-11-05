##env 4.2
library(vcfR)

args <- commandArgs(trailingOnly = TRUE)

vcf <- read.vcfR(args[1])

data1 <- vcf@fix
data2 <- vcf@gt

data <- cbind(data1, data2)

write.table(data, paste0(strsplit(args[1], '\\.')[[1]][1],".vcf.txt"), sep=",")

