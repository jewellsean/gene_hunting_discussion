## Generate gammas for MC model 
library(yaml)
seed <- 123
config_file <- "configs/config_desktop.yml"
configs <- yaml.load_file(config_file)
p <- configs$p
set.seed(seed)

gammas <- runif(p - 1, 0, 0.5)

readr::write_csv(data.frame(gam = gammas), path = "gamma_hyperparameters.csv")
