# Reproduce Figure 1 of Jewell and Witten (2018) discussion of Sesia et al. (2018).
# Figure 1: Validity of constructed knockoffs for Markov chain covariates: pairwise correlations of covariates
# versus pairwise correlations of knockoffs to check exchangeability condition (3.1) of Candès et al (2018).
# Motivated from the tutorial GWAS tutorial from the author's webpage
# https://web.stanford.edu/group/candes/knockoffs/tutorials/gwas_tutorial.html#knockoff_diagnostics

required_packages <- c("SNPknock", "knockoff", "tidyverse", "yaml", "gridExtra", "grid")
new_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]
if(length(new_packages)) install.packages(new_packages)


args <- commandArgs(trailingOnly = TRUE)

library(SNPknock)
library(knockoff)
library(tidyverse)
library(yaml)
library(gridExtra)
library(grid)
source("utils.R")

## simualtion parameters
config_file <- args[1]
configs <- yaml.load_file(config_file)
n <- configs$n_obs
p <- configs$p
K <- 5

set.seed(1)

# simulate discrete markov chain and associated discrete markov chain and gaussian knockoffs
markov_states <- simulate_discrete_markov_states(n, p, K)
knockoffs_dmc <- SNPknock.knockoffDMC(markov_states$covariates, markov_states$prob_initial, markov_states$transition_matrices, display_progress = FALSE)
knockoffs_gaussian <- create.second_order(markov_states$covariates)

# calculate pairwise correlations for covariates and knockoffs
df_markov <- compute_covariate_knockoff_corr(markov_states$covariates, knockoffs_dmc)
df_markov$knockoff_type <- "Markov"

df_gaussian <- compute_covariate_knockoff_corr(markov_states$covariates, knockoffs_gaussian)
df_gaussian$knockoff_type <- "Gaussian"

df <- rbind(df_markov, df_gaussian)
df$knockoff_type <- factor(df$knockoff_type, levels = c("Markov", "Gaussian"), 
                           labels = c("Knockoffs generated under Markov chain",
                                      "Knockoffs generated under Gaussian approximation"))

# compare the pairwise correlations under the correct and misspecified models
p1 <- df %>% filter(knockoff_type == "Knockoffs generated under Markov chain") %>% 
  ggplot(aes(x = covariate_corr, y = knockoff_corr)) + 
  geom_point(alpha = 0.25) + 
  facet_wrap(~knockoff_type) + 
  xlim(c(-0.1, 0.6)) +
  ylim(c(-0.1, 0.6)) +
  xlab('Pairwise covariate correlation') + 
  ylab('Pairwise knockoff correlation') + 
  geom_abline(slope = 1, intercept = 0, col = "red") +
  theme_bw()

p2 <- df %>% filter(knockoff_type == "Knockoffs generated under Gaussian approximation") %>% 
  ggplot(aes(x = covariate_corr, y = knockoff_corr)) + 
  geom_point(alpha = 0.25) + 
  xlim(c(-0.1, 0.6)) +
  ylim(c(-0.1, 0.6)) +
  facet_wrap(~knockoff_type) + 
  xlab('Pairwise covariate correlation') + 
  ylab('Pairwise knockoff correlation') + 
  geom_abline(slope = 1, intercept = 0, col = "red") +
  theme_bw()

## https://stackoverflow.com/questions/17576381/corner-labels-in-ggplot2
p1 <- arrangeGrob(p1, top = textGrob("(a)", x = unit(0, "npc")
                                               , y   = unit(0.5, "npc"), just=c("left","top"),
                                               gp=gpar(col="black", fontsize=12)))


p2 <- arrangeGrob(p2, top = textGrob("(b)", x = unit(0, "npc")
                                     , y   = unit(0.5, "npc"), just=c("left","top"),
                                     gp=gpar(col="black", fontsize=12)))


dir.create(configs$figure_directory, showWarnings = FALSE)

pdf(paste0(configs$figure_directory, "diagonistics.pdf"),height = 4, width = 8.5)
grid.arrange(p1, p2, ncol = 2)
dev.off()



