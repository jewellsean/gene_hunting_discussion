# Reproduce Figure 2 of Jewell and Witten (2018) discussion of Sesia et al. (2018).
# Figure 2: Effect of knockoff model misspecification on FDR (top row) and Power (bottom row) for
# hidden Markov (left column) and Markov chain (right column) covariates
#
# This file creates all of the experiment files to produce Figure 2
# Command line arguments:
# signal_amplitude -- scalar value of beta
# covariate_type --
#    'dmc' -- discrete markov chain
#    'hmm' -- hidden markov chain
# knockoff_type (same options as covariate_type + gaussian)
# variable_importance_measure (coefdiff or tree)
# config_file -- location of config file
# returns:
#    no return value, writes fdps and powers to csv file in save directory specified in config file
required_packages <- c("yaml")
new_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]
if(length(new_packages)) install.packages(new_packages, repos = "http://cran.us.r-project.org")

library(yaml)
source("utils.R")
## Load in settings via CLI and config file
args <- commandArgs(trailingOnly = TRUE)
signal_amplitude <- as.numeric(args[1])
covariate_type <- args[2]
knockoff_type <- args[3]
variable_importance_measure <- args[4]
config_file <- args[5]

configs <- yaml.load_file(config_file)

n <- configs$n_obs
p <- configs$p
number_simulations <- configs$number_simulations
active_support_size <- configs$active_support_size
write_directory <- configs$write_directory
fdr <- configs$fdr
knockoff_offset <- 1 # always want to have conservative offset

# change state space size based on type of covariate model
if (covariate_type == 'hmm') {
  K <- 9
  covariate_state_offset <- -4
} else {
  K <- 5
  covariate_state_offset <- -2
}

out <- simulate_full_independence(number_simulations, n, p, K, active_support_size, signal_amplitude,
                                  covariate_state_offset, fdr, knockoff_offset, covariate_type, knockoff_type, variable_importance_measure)


df <- data.frame(covariate_type = covariate_type, signal_amplitude = signal_amplitude,
                 n = n, p = p, K = K, number_simulations = number_simulations, active_support_size = active_support_size, 
                 fdr = fdr, knockoff_offset = knockoff_offset, fdps = out$fdps, pows = out$pows,
                 knockoff_type = knockoff_type, variable_importance_measure)

dir.create(write_directory, showWarnings = FALSE)
# possibility for name clashes, but should be fairly low chance.
file_name <- paste0(write_directory, sample.int(1e10, size = 1, replace = F), "-", knockoff_type, covariate_type, ".csv")
write.csv(df, file = file_name, row.names = FALSE)
