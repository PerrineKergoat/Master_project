#!/bin/bash

library(vcfR)
library(hierfstat)
library(gaston)
library(stats)
library(StatMatch)
library(dgof)
library(outliers)
library(nortest)

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

vcf_mat_pop1 <- vcf_mat[1:1000,]
vcf_mat_pop2 <- vcf_mat[1001:2000,]

### Get SNPs and parameters to cut the genomes

data_snps <- vcf_file@snps
nb_snps <- nrow(data_snps)
l_window <- 500
overlap <- 100
nb_windows <- (nb_snps-l_window) %/% (overlap) + 1
start_bound <- seq(from = 1, to = nb_snps-l_window, by = overlap)

print("NB of windows")
print(nb_windows)
print("NB of SNPs in window")
print(nb_snps)
print("Window parameters caculated")

### Create the windows and calculate beta kinship for each and for the rest of the windows for the whole population

list_vcf <- list()
list_beta <- list()
list_beta_prime <- list()

list_vcf_pop1 <- list()
list_beta_pop1 <- list()
list_beta_prime_pop1 <- list()

list_vcf_pop2 <- list()
list_beta_pop2 <- list()
list_beta_prime_pop2 <- list()

columns <- c("window", "mean Beta", "max Beta", "min Beta", "variance Beta", "pi", "Tajima's D", "Fst")
pop_summary <- data.frame(matrix(nrow = nb_windows, ncol = length(columns)))
pop1_summary <- data.frame(matrix(nrow = nb_windows, ncol = length(columns)))
pop2_summary <- data.frame(matrix(nrow = nb_windows, ncol = length(columns)))
colnames(pop_summary) = columns
colnames(pop1_summary) = columns
colnames(pop2_summary) = columns

i <- 0

#par(mfrow = c(2,2))
for (start in start_bound){

  i <- i + 1
  end <- start + l_window
  
  print(i)

  if (end > nb_snps) {end = nb_snps}
  if (start == nb_snps) {break}
  if (end-start < l_window) {break}

  else {
    print(paste0("Starting window ", i))
    
    reduced_vcf <- vcf_mat[,start:end]
    #list_vcf[[i]] <- reduced_vcf
    reduced_beta <- beta.dosage(reduced_vcf)
    #list_beta[[i]] <- reduced_beta
    
    reduced_vcf_pop1 <- vcf_mat_pop1[,start:end]
    #list_vcf_pop1[[i]] <- reduced_vcf_pop1
    reduced_beta_pop1 <- beta.dosage(reduced_vcf_pop1)
    #list_beta_pop1[[i]] <- reduced_beta_pop1
    
    reduced_vcf_pop2 <- vcf_mat_pop2[,start:end]
    #list_vcf_pop2[[i]] <- reduced_vcf_pop2
    reduced_beta_pop2 <- beta.dosage(reduced_vcf_pop2)
    #list_beta_pop2[[i]] <- reduced_beta_pop2

    if (start == 1){
      rest_vcf <- vcf_mat[,end:nb_snps]
      rest_vcf_pop1 <- vcf_mat_pop1[,end:nb_snps]
      rest_vcf_pop2 <- vcf_mat_pop2[,end:nb_snps]
    }
    if (end == nb_snps){
      rest_vcf <- vcf_mat[,1:start]
      rest_vcf_pop1 <- vcf_mat_pop1[,1:start]
      rest_vcf_pop2 <- vcf_mat_pop2[,1:start]
    }
    else {
      rest_vcf <- vcf_mat[,c(1:start,end:nb_snps)]
      rest_vcf_pop1 <- vcf_mat_pop1[,c(1:start,end:nb_snps)]
      rest_vcf_pop2 <- vcf_mat_pop2[,c(1:start,end:nb_snps)]
    }
    
    #rest_beta_prime <- beta.dosage(rest_vcf)
    ##list_beta_prime[[i]] <- rest_beta_prime
    #rest_beta_prime_pop1 <- beta.dosage(rest_vcf_pop1)
    ##list_beta_prime_pop1[[i]] <- rest_beta_prime_pop1
    #rest_beta_prime_pop2 <- beta.dosage(rest_vcf_pop2)
    ##list_beta_prime_pop2[[i]] <- rest_beta_prime_pop2
    
    filtered_beta <- reduced_beta
    filtered_beta <- filtered_beta[upper.tri(filtered_beta, diag = FALSE)]
    
    pop_summary[i, "window"] <- i
    pop_summary[i, "mean Beta"] <- mean(filtered_beta)
    pop_summary[i, "max Beta"] <- max(filtered_beta)
    pop_summary[i, "min Beta"] <- min(filtered_beta)
    pop_summary[i, "variance Beta"] <- var(filtered_beta)
    pop_summary[i, "pi"] <- pi.dosage(reduced_vcf)
    pop_summary[i, "Tajima's D"] <- TajimaD.dosage(reduced_vcf)
    
    filtered_beta_pop1 <- reduced_beta_pop1
    filtered_beta_pop1 <- filtered_beta_pop1[upper.tri(filtered_beta_pop1, diag = FALSE)]
    
    pop1_summary[i, "window"] <- i
    pop1_summary[i, "mean Beta"] <- mean(filtered_beta_pop1)
    pop1_summary[i, "max Beta"] <- max(filtered_beta_pop1)
    pop1_summary[i, "min Beta"] <- min(filtered_beta_pop1)
    pop1_summary[i, "variance Beta"] <- var(filtered_beta_pop1)
    pop1_summary[i, "pi"] <- pi.dosage(reduced_vcf_pop1)
    pop1_summary[i, "Tajima's D"] <- TajimaD.dosage(reduced_vcf_pop1)
    
    filtered_beta_pop2 <- reduced_beta_pop2
    filtered_beta_pop2 <- filtered_beta_pop2[upper.tri(filtered_beta_pop2, diag = FALSE)]
    
    pop2_summary[i, "window"] <- i
    pop2_summary[i, "mean Beta"] <- mean(filtered_beta_pop2)
    pop2_summary[i, "max Beta"] <- max(filtered_beta_pop2)
    pop2_summary[i, "min Beta"] <- min(filtered_beta_pop2)
    pop2_summary[i, "variance Beta"] <- var(filtered_beta_pop2)
    pop2_summary[i, "pi"] <- pi.dosage(reduced_vcf_pop2)
    pop2_summary[i, "Tajima's D"] <- TajimaD.dosage(reduced_vcf_pop2)
  }
}

print("Metrics caculated")

### Calculate eigenvalues for the windows anf the rest beta kinship matrix 

#List_betaDIFF = vector(mode = "list", length = nb_windows)
#eigenval_betaDIFF = vector(mode = "list", length = nb_windows)


### Graph plotting 

position <- (match(5e6, vcf_file@snps$pos))

start_bound <- as.list(start_bound)
list_start <- list()
ind <- 0

for (lim in start_bound) {
  if (position > lim) {
    if (position < lim + l_window) {
      ind <- ind + 1
      print(match(lim, start_bound))
      list_start[[ind]] <- match(lim, start_bound)
    }
  }
}


save.image(paste0("/work/FAC/FBM/DEE/jgoudet/default/pkergoat/2pop_1sel_500snps_100overlap/data/TestFST_500_", replicateN, ".RData"))
save(pop_summary, pop1_summary, pop2_summary, nb_windows, start_bound, list_start, freq_check, file = paste0("~/WORK/2pop_1sel_500snps_100overlap/data/Graph_testFST_500_", replicateN, ".RData"))

print("Script Comparison.R done")
