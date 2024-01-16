#!/bin/bash

library(vcfR)
library(gaston)
library(hierfstat)
library(SNPRelate)

args <- commandArgs(trailingOnly = TRUE)
replicateN <- args[[1]]

vcf_file <- read.VCF(paste0("/work/FAC/FBM/DEE/jgoudet/default/pkergoat/2pop_1sel_500snps_100overlap/data/Simu_testFST_", replicateN, ".vcf"))
print(vcf_file)
vcf_mat <- as.matrix(vcf_file)
#beta_mat <- beta.dosage(vcf_mat)

#extarct snp id for the selected mutation
position = vcf_file@snps$id[vcf_file@snps$pos == 5e6]

#check count of alternate allele per pop
genosPOP1selmut = vcf_mat[1:1000,position]
genosPOP2selmut = vcf_mat[1001:2000,position]

#estimate freq sel mut per pop
frqPOP1 = sum(genosPOP1selmut)/2000
frqPOP2 = sum(genosPOP2selmut)/2000
freq_check <- c(frqPOP1, frqPOP2)

#FST analysis
showfile.gds(closeall=TRUE)
write.bed.matrix(vcf_file, "./BED_object")
snpgdsBED2GDS("BED_object", out.gdsfn = "./GDS_object")
genomic_data = snpgdsOpen("GDS_object")

Fst_calc <- snpgdsFst(gdsobj = genomic_data, population = as.factor(c(rep(1, 1000), rep(2,1000))), method="W&C84", remove.monosnp = FALSE, autosome.only=FALSE)
Fst_calc_maf_01 <- snpgdsFst(gdsobj = genomic_data, population = as.factor(c(rep(1, 1000), rep(2,1000))), method="W&C84", remove.monosnp = FALSE, autosome.only=FALSE, maf=0.01)

save(Fst_calc_maf_01, Fst_calc, freq_check, vcf_file, file = paste0("/work/FAC/FBM/DEE/jgoudet/default/pkergoat/2pop_1sel_500snps_100overlap/data/FST_calc_500", replicateN, ".Rdata"))

