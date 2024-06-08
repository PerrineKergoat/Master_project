#!/bin/bash

library(vcfR)
library(gaston)
library(hierfstat)
library(SNPRelate)

args <- commandArgs(trailingOnly = TRUE)
job_id <- args[[1]]
rep_nb <- args[[2]]
nb_pop <- args[[3]]
gen_nb <- args[[4]]
path_data <- args[[5]]

vcf_file <- read.VCF(paste0(path_data, "Recap_", gen_nb,"_", job_id, "_", rep_nb, ".vcf"))

# Convert VCF into matrix 
vcf_mat <- as.matrix(vcf_file)

# Extract SNP id for the beneficial mutation (selection)
position = vcf_file@snps$id[vcf_file@snps$pos == 5e6]

# Estimate frequency of beneficial mutation per subpopulations
i <- 1
freq_check <- c()
for (pop1 in seq(from = 1, to = nb_pop, by = 1)){
  first_ind_subpop <- i
  last_ind_subpop <- 1000 * pop1
  i <- i + 1000

  mut_sel_pop <- vcf_mat[first_ind_subpop:last_ind_subpop, position]
  freq_pop <- sum(mut_sel_pop)/2000
  freq_check[pop1] <- freq_pop
}

# Get SNPs and parameters to cut the genomes
data_snps <- vcf_file@snps
nb_snps <- nrow(data_snps)

list_start_position <- seq(from = 1, to = 1e7-56000, by = 14000)
list_end_position <- seq(from = 70000, to = 1e7 + 14000, by =14000)
data_lim <- data.frame(list_start_position, list_end_position)
data_lim$list_SNPs <- rep(c(NA), length(list_start_position))

# Creation of the boundaries to cut the genome by bp position
for (id_SNP in vcf_file@snps$id) {
  pos_SNP <- vcf_file@snps$pos[vcf_file@snps$id == id_SNP]

  if(as.integer((pos_SNP %% 14000) * 14000) != as.integer(pos_SNP)) {
    mult14_lim <- pos_SNP %/% 14000

    pos_min <- (mult14_lim + 1) * 14000
    pos_max <- mult14_lim * 14000 +1

    if (pos_min < 70000) {
      pos_min <- 70000
    }

    if (pos_min > 10010000) {
      pos_min <- 10010000
    }

    if (pos_max > 9940001) {
      pos_max <- 9940001
    }

    ind_min <- match(pos_min,list_end_position)
    ind_max <- match(pos_max, list_start_position)

    wind_snp <-  seq(from = ind_min, to = ind_max, by = 1)

    for (window in wind_snp) {
      prov_list <- data_lim[window, 3]
      prov_list <- as.list(prov_list)
      prov_list <- append(prov_list, id_SNP)
      data_lim[window, 3] <- toString(prov_list)
    }

  } else {
    mult14_lim <- pos_SNP %/% 14000

    pos_min <- (mult14_lim + 1) * 14000
    pos_max <- mult14_lim * 14000 +1

    if (pos_min < 70000) {
      pos_min <- 70000
    }

    if (pos_max > 10010000) {
      pos_max <- 10010000
    }

    ind_min <- match(pos_min,list_end_position) -1
    ind_max <- match(pos_max, list_start_position) -1

    wind_snp <-  seq(from = ind_min, to = ind_max, by = 1)

    for (window in wind_snp) {
      prov_list <- data_lim[window, 3]
      prov_list <- as.list(prov_list)
      prov_list <- append(prov_list, id_SNP)
      data_lim[window, 3] <- toString(prov_list)
    }
  }
}

# Translation of the windows boundaries in SNPs number / id
data_lim$nb_SNPs <- 0
data_lim$start_SNP <- 0
data_lim$end_SNP <- 0

for (nb_row in c(1:length(list_start_position))) {

  rem_NA <- data_lim$list_SNPs[nb_row]
  rem_NA <- strsplit(rem_NA, ", ")[[1]]
  rem_NA <- rem_NA[-1]
  nb_SNPs <- length(rem_NA)

  start_snp <- as.numeric(head(rem_NA, n=1))
  end_snp <- as.numeric(tail(rem_NA, n=1))
  nb_snp <- end_snp - start_snp

  rem_NA <- toString(rem_NA)
  data_lim$list_SNPs[nb_row] <- rem_NA
  data_lim$nb_SNPs[nb_row] <- nb_snp
  data_lim$start_SNP[nb_row] <- start_snp
  data_lim$end_SNP[nb_row] <- end_snp
}

start_bound <- data_lim$start_SNP
end_bound <- data_lim$end_SNP
nb_windows <- length(start_bound)

print(paste0("NB of windows: ", nb_windows))
print(paste0("NB of SNPs: ", nb_snps))
print("Windows parameters caculated")

# Selection of individuals to compare the same number of individuals between the different simulations
j <- 1
subset_vcf_mat <- c()
for (pop2 in seq(from = 1, to = nb_pop, by = 1)){
  start_sel <- j
  end_sel <- (1000/as.numeric(nb_pop)) +j -1
  j <- j + 1000

  subset_vcf_mat <- append(subset_vcf_mat, seq(from = start_sel, to = end_sel, by = 1))
}

sub_vcf_mat <- vcf_mat[subset_vcf_mat,]

# Beta kinship and pi calculation for each window (whole population scale)
columns <- c("window", "beta_variance", "pi")
pop_summary <- data.frame(matrix(nrow = nb_windows, ncol = length(columns)))
colnames(pop_summary) = columns

i <- 0
for (start in start_bound){
  i <- i + 1
  end <- end_bound[i]
  start <- as.integer(start)
  end <- as.integer(end)
  l_window <- end - start

  if (end > nb_snps) {end = nb_snps}
  if (start == nb_snps) {break}
  if (end-start < l_window) {break}

  else {
    print(paste0("Calculating beta and pi in window: ", i))

    reduced_vcf <- sub_vcf_mat[,start:end]  # reduced_vcf only contains the SNPs on the window
    reduced_beta <- beta.dosage(reduced_vcf)

    filtered_beta <- reduced_beta
    filtered_beta <- filtered_beta[upper.tri(filtered_beta, diag = FALSE)]  # beta matrix is symmetrical
                                                                            # so the variance of one half
                                                                            # equals the variance of the
                                                                            # whole matrix

    pop_summary[i, "window"] <- i
    pop_summary[i, "beta_variance"] <- var(filtered_beta)
    pop_summary[i, "pi"] <- pi.dosage(reduced_vcf)
  }
}
print("Beta and pi metrics calculated")

# FST analysis
showfile.gds(closeall=TRUE)
write.bed.matrix(vcf_file, paste0("./BED/BED_object_", job_id, "_", rep_nb))
snpgdsBED2GDS(paste0("./BED/BED_object_", job_id, "_", rep_nb), out.gdsfn = paste0("./GDS/GDS_object_", job_id, "_", rep_nb))
genomic_data = snpgdsOpen(paste0("./GDS/GDS_object_", job_id, "_", rep_nb))

pop_vector <- c() # for the pop-specific Fst we are calculating we need to explain the function what subpop. ind. belong to
for (pop3 in seq(from = 1, to = nb_pop, by = 1)){
  pop_vector <- append(pop_vector, rep(pop3, 1000))
}

Fst_calc <- snpgdsFst(gdsobj = genomic_data,
                      population = as.factor(pop_vector),
                      method="W&C84", remove.monosnp = FALSE, autosome.only=FALSE)
Fst_calc_maf <- snpgdsFst(gdsobj = genomic_data,
                             population = as.factor(pop_vector),
                             method="W&C84", remove.monosnp = FALSE, autosome.only=FALSE,
                             maf=0.005)

# Fst filtered by window
k <- 0
snps_window <- list()
snps_fst <- as.vector(na.omit(Fst_calc$FstSNP))

for (bound in start_bound) {
  k <- k + 1
  start_mean <- bound
  end_mean <- end_bound[k]
  snps_window[k] <- mean(snps_fst[start_mean : end_mean])
}

print("Fst calculated")

# Save everything we need for the Replicates. R script and for the graphs
save(pop_summary, nb_windows, snps_window, data_lim, position, freq_check, Fst_calc, Fst_calc_maf,
     file = paste0(path_data, "Graph_env_", gen_nb,"_", job_id, "_", rep_nb, ".RData"))

print(paste0("Script Comparison.R done for ", gen_nb," generations"))
