#!/bin/bash

args = commandArgs()

job_id <- args[6]
nb_rep <- args[7]
nb_rep <- as.integer(nb_rep)
percent <- args[8]
percent <- as.integer(percent)
path_data <- args[[9]]
path_rec_out <- args[[10]]

################################### FUNCTIONS ###################################

# Lists all the files obtained for a generation number for a simulation job id
List_files <- function(gen_nb, job_id) {  # gen_nb is a string
  list_files <- list.files(path = path_data,
                               pattern = paste0("Graph_env_", gen_nb, "_", job_id))
  print(paste0("List files ", gen_nb, " generations obtained"))

  return(list_files)  # returns a list with all the file names
}



# Creates the df where interesting windows and all the parameters are saved
Summary_df_creation <- function(metric, gen_list, nb_rep) {  # metric is a string, gen_list is a list of strings
  col_names <- c(paste0("Threshold_", metric,"_value"),
                 paste0("p_value_", metric),
                 "Windows_number",
                 paste0("Windows_index_", metric),
                 paste0("Top_", metric,"_values"))

  for (gen in gen_list) {
    name <- paste0("summary_replicates_", metric,"_", gen)
    df <- data.frame(matrix(nrow = nb_rep, ncol = 5))
    colnames(df) <- col_names
    assign(name, df, envir = parent.frame())
  }
}

# Detects Fst "positive" signals (ie. windows potentially under selection)
Fst_norm_trans_percentile <- function(snps_window, percent) {
  data_fst <- as.numeric(snps_window)

  # Normal transformation of fst data using log transformation
  trans_data_fst <- log10(data_fst)
  trans_data_fst <- as.list(trans_data_fst)

  # Establish the parameters of the normal distribution of Fst data
  mean_trans_fst <- mean(na.omit(as.numeric(trans_data_fst)))
  sd_trans_fst <- sd(na.omit(as.numeric(trans_data_fst)))

  # Get the percentile-threshold
  trans_lim_fst_95 <- qnorm((100-percent)/100, mean = mean_trans_fst, sd = sd_trans_fst)
  pval_fst_95 <- pnorm(trans_lim_fst_95, mean = mean_trans_fst, sd = sd_trans_fst, lower.tail = FALSE)
  lim_fst_95 <- 10^trans_lim_fst_95

  # Find the values over the threshold
  index_top5_fst <- which(data_fst >= lim_fst_95)
  top5_fst <- data_fst[index_top5_fst]

  # Sum up everything and return the list
  summary_fst <- list(lim_fst_95, pval_fst_95, nb_windows, toString(index_top5_fst), toString(top5_fst))

  return(summary_fst)
}

# Detects Pi "positive" signals (ie. windows potentially under selection)
Pi_norm_trans_percentile <- function(pop_summary, percent) {
  data_pi <- pop_summary$pi

  # Establish the parameters of the normal distribution of pi data (no transformation needed)
  mean_pi <- mean(data_pi)
  sd_pi <- sd(data_pi)

  # Get the percentile-threshold
  lim_pi_5 <- qnorm((percent)/100, mean = mean_pi, sd = sd_pi)
  pval_pi_5 <- pnorm(lim_pi_5, mean = mean_pi, sd = sd_pi, lower.tail = FALSE)

  # Find the values below the threshold
  index_top5_pi <- which(data_pi <= lim_pi_5)
  top5_pi <- data_pi[index_top5_pi]

  # Sum up everything and return the list
  summary_pi <- list(lim_pi_5, pval_pi_5, nb_windows, toString(index_top5_pi), toString(top5_pi))

  return(summary_pi)
}

# Detects Beta "positive" signals (ie. windows potentially under selection)
Beta_norm_trans_percentile <- function(pop_summary, percent) {
  data_beta <- pop_summary$beta_variance

  # Normal transformation of beta data using log transformation
  max_beta <- max(data_beta)
  min_beta <- min(data_beta)
  trans_data_beta <- (data_beta - max_beta)/(max_beta - min_beta)
  trans_data_beta <- log10(trans_data_beta + 1)
  trans_data_beta[is.infinite(trans_data_beta)] <- NA

  # Establish the parameters of the normal distribution of beta data
  mean_trans_beta <- mean(na.omit(trans_data_beta))
  sd_trans_beta <- sd(na.omit(trans_data_beta))

  # Get the percentile-threshold
  trans_lim_beta_95 <- qnorm((100-percent)/100, mean = mean_trans_beta, sd = sd_trans_beta)
  pval_beta_95 <- pnorm(trans_lim_beta_95, mean = mean_trans_beta, sd = sd_trans_beta, lower.tail = FALSE)

  # Find the values over the threshold
  index_top5_beta <- which(trans_data_beta >= trans_lim_beta_95)
  top5_beta <- data_beta[index_top5_beta]
  lim_beta_95 <- "NC"

  # Sum up everything and return the list
  summary_beta <- list(lim_beta_95, pval_beta_95, nb_windows, toString(index_top5_beta), toString(top5_beta))

  return(summary_beta)
}

# Extract the recombination rate of the portion of the genome containing the beneficial mutation
Rec_rate_analysis <- function(path_rec_out, file_name) {
  print(paste0(path_rec_out, file_name))
  
  out_data <- read.csv(paste0(path_rec_out, file_name))
  list_rate <- as.vector(out_data$position)
  candidate_pos <- 0
  
  for (rate in list_rate) {
    if (rate <= 5e6) {
      candidate_pos <- rate
    }
  }
  
  index_rate <- match(candidate_pos, list_rate) + 1
  print(out_data$rate[index_rate])
  return(out_data$rate[index_rate])
}

#################################################################################
job_id <- 44679626
nb_rep <- 3
percent <- 1
path_data <- "/Users/perrinekergoat/Desktop/CheckBeta/data/"
path_rec_out <- "/Users/perrinekergoat/Desktop/CheckBeta/rec_rate_out/"

### List files for a simulation
list_files_600 <- List_files("600", job_id)
list_files_800 <- List_files("800", job_id)
list_files_1000 <- List_files("1000", job_id)
print("Files listed")

### Create df where info will be saved for fst, pi and beta ###
Summary_df_creation("fst", c("600", "800", "1000"), nb_rep)
Summary_df_creation("pi", c("600", "800", "1000"), nb_rep)
Summary_df_creation("beta", c("600", "800", "1000"), nb_rep)
print("DF created")

#List files for the recombination rate analysis
list_files_recomb_rate_600 <- List_files_recomb_rate("600", job_id)
list_files_recomb_rate_800 <- List_files_recomb_rate("800", job_id)
list_files_recomb_rate_1000 <- List_files_recomb_rate("1000", job_id)

### Create dataframe where info will be saved for the recomb. rate ###
summary_recomb_rates_600 <- data.frame()
summary_recomb_rates_800 <- data.frame()
summary_recomb_rates_1000 <- data.frame()


### Analyse replicates of 600 generations ###
i<- 0
for (file_name in list_files_600){
  load(paste0(path_data, file_name))
  i <- i+1
  print(paste0("File 600 gen (fst, pi, beta): ", i))

  summary_replicates_fst_600[i, ] <- Fst_norm_trans_percentile(snps_window, percent)
  summary_replicates_pi_600[i, ] <- Pi_norm_trans_percentile(pop_summary, percent)
  summary_replicates_beta_600[i, ] <- Beta_norm_trans_percentile(pop_summary, percent)
}
print("Fst, pi and beta calculated for 600 generations")

### Analysis of the changing recombination rate 
i<- 0
summary_recomb_rate_600 <- data.frame()
for (file_name in list_files_recomb_rate_600){
  i <- i+1
  print(paste0("File 600 gen (recombination rate): ", i))
  
  summary_recomb_rate_600[i, 1] <- Rec_rate_analysis(path_rec_out, file_name)
}
print("Recombination rate analysed for 600 generations")

save(summary_replicates_fst_600, summary_replicates_pi_600, summary_replicates_beta_600, summary_recomb_rate_600,
     data_lim, position, freq_check, Fst_calc, Fst_calc_maf,
     file = paste0(path_data, "Graph_replicates_600_", job_id, ".RData"))

print("Analysis of 600 generations saved")



# Lists all the files obtained for a generation number for a simulation job id with changing recombination rate 
List_files_recomb_rate <- function(gen_nb, job_id) {  # gen_nb is a string
  list_files <- list.files(path = path_rec_out,
                           pattern = paste0("out_", job_id))
  print(paste0("List files varying recomb. rate ", gen_nb, " generations obtained"))
  
  return(list_files)  # returns a list with all the file names
}



### Analyse replicates of 800 generations #
i<- 0
for (file_name in list_files_800){
  load(paste0(path_data, file_name))
  i <- i+1
  print(paste0("File 800 gen (fst, pi, beta): ", i))

  summary_replicates_fst_800[i, ] <- Fst_norm_trans_percentile(snps_window, percent)
  summary_replicates_pi_800[i, ] <- Pi_norm_trans_percentile(pop_summary, percent)
  summary_replicates_beta_800[i, ] <- Beta_norm_trans_percentile(pop_summary, percent)
}
print("Fst, pi and beta calculated for 800 generations")

### Analysis of the changing recombination rate for 800 generations
i<- 0
summary_recomb_rate_800 <- data.frame()
for (file_name in list_files_recomb_rate_800){
  i <- i+1
  summary_recomb_rate_800[i, 1] <- Rec_rate_analysis(path_rec_out, file_name)
}
print("Recombination rate analysed for 800 generations")

save(summary_replicates_fst_800, summary_replicates_pi_800, summary_replicates_beta_800, summary_recomb_rate_800,
     data_lim, position, freq_check, Fst_calc, Fst_calc_maf,
     file = paste0(path_data, "Graph_replicates_800_", job_id, ".RData"))

print("Analysis of 800 generations saved")

### Analyse replicates of 1000 generations
i<- 0
for (file_name in list_files_1000){
  load(paste0(path_data, file_name))
  i <- i+1
  print(paste0("File 1000 gen (fst, pi, beta): ", i))

  summary_replicates_fst_1000[i, ] <- Fst_norm_trans_percentile(snps_window, percent)
  summary_replicates_pi_1000[i, ] <- Pi_norm_trans_percentile(pop_summary, percent)
  summary_replicates_beta_1000[i, ] <- Beta_norm_trans_percentile(pop_summary, percent)
}
print("Fst, pi and beta calculated for 1000 generations")

### Analysis of the changing recombination rate for 1000 generations
i<- 0
summary_recomb_rate_1000 <- data.frame()
for (file_name in list_files_recomb_rate_1000){
  i <- i+1
  summary_recomb_rate_1000[i, 1] <- Rec_rate_analysis(path_rec_out, file_name)
}
print("Recombination rate analysed for 1000 generations")

save(summary_replicates_fst_1000, summary_replicates_pi_1000, summary_replicates_beta_1000, summary_recomb_rate_1000,
     data_lim, position, freq_check, Fst_calc, Fst_calc_maf,
     file = paste0(path_data, "Graph_replicates_1000_", job_id, ".RData"))

print("Analysis of 1000 generations done")

print("Script Replicates.R done")
