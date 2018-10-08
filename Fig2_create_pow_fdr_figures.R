# Reproduce Figure 2 of Jewell and Witten (2018) discussion of Sesia et al. (2018).
# Figure 2: Effect of knockoff model misspecification on FDR (top row) and Power (bottom row) for
# hidden Markov (left column) and Markov chain (right column) covariates
#
# Before running this plotting script
# 1. Configure config.yml
# 2. Run run_experiments.sh
# these two steps will generate required simualtions
required_packages <- c("gridExtra", "tidyverse", "yaml")
new_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]
if(length(new_packages)) install.packages(new_packages)

args <- commandArgs(trailingOnly = TRUE)
library(gridExtra)
library(tidyverse)
library(yaml)


config_file <- args[1]
configs <- yaml.load_file(config_file)
experiment_directory <- configs$load_experiments

df <- NULL
experiment_files <- list.files(experiment_directory, pattern = ".csv")

if (is.null(experiment_files)) {
  stop(paste0('This script analyzes results from a large scale simualation study. Results, generated from main.R,
must be saved in experiment_directory: ', experiment_directory))
}

## collect resutls from all experiments
for (experiment_file in experiment_files) {
  file_i <- paste0(experiment_directory, "/", experiment_file)
  df <- rbind(df, read.csv(file_i))
}  

# only need subset of stored information
df <- df %>% mutate(signal_amplitude = as.factor(signal_amplitude)) %>%
  select(signal_amplitude, fdps, pows, covariate_type, knockoff_type) 

# rename and reorder factors for plotting
df$knockoff_type <- factor(df$knockoff_type, levels = c("hmm", "dmc", "gaussian"),
                           labels = c("Hidden Markov Model", "Markov Chain", "Gaussian Approximation"))

df$covariate_type <- factor(df$covariate_type, levels = c("hmm", "dmc"),
                           labels = c("Covariate Generation: Hidden Markov Model", "Covariate Generation: Markov Chain"))

p_fdr <- df %>% 
  ggplot(aes(x = signal_amplitude, y = fdps)) + 
  geom_boxplot(aes(fill = knockoff_type)) +
  geom_hline(yintercept = 0.1, col = "red") +
  ylim(0, 0.5) +
  ylab('FDR') + 
  xlab('') +
  facet_wrap(~covariate_type) +
  theme_bw() +
  scale_fill_grey(start = 0.5, end = .9) +
  theme(legend.position="none")
  
p_pow <- df %>% 
  ggplot(aes(x = signal_amplitude, y = pows)) + 
  geom_boxplot(aes(fill = knockoff_type)) +
  ylim(0, 1) +
  ylab('Power') + 
  xlab('Signal Amplitude') +
  facet_wrap(~covariate_type) +
  theme_bw() +
  theme(legend.position="bottom") + 
  scale_fill_grey(start = 0.5, end = .9, name = "Knockoff Generated Under: ")


dir.create(configs$figure_directory, showWarnings = FALSE)
pdf(paste0(configs$figure_directory, "pow_fdr_simulated.pdf"), width = 8.5, height = 5)
grid.arrange(p_fdr, p_pow, nrow = 2, heights = c(3, 4))
dev.off()