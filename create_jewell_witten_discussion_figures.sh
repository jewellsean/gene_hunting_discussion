#!/usr/bin/env bash

echo "Generating Figure 1 of Jewell and Witten (2018) discussion"

RScript Fig1_create_diag_corr_figures.R configs/config_desktop.yml


echo "Running simulation experiments"
echo "NB: Ensure that configuration file is updated."
echo "For desktop computation considerations the default simulation study is smaller than that used in Figure 2!"
source run_experiments.sh


echo "Generating Figure 2 of Jewell and Witten (2018) discussion"
RScript Fig2_create_pow_fdr_figures.R configs/config_desktop.yml
