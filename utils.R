required_packages <- c("SNPknock", "knockoff")
new_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]
if(length(new_packages)) install.packages(new_packages, repos = "http://cran.us.r-project.org")

library(SNPknock)
library(knockoff)

expit <- function(x) { exp(x) / (1 + exp(x))}

simulate_discrete_markov_states <- function(n, p, K) {
    prob_initial <- rep(1 / K, K) # uniform inital distribution
    gammas <- readr::read_csv("gamma_hyperparameters.csv")$gam
    transition_matrices <- array(data = NA, dim = c((p - 1), K, K)) # tensor to store each transition matrix
    for (transition_i in 1 : (p - 1)) {
        transition_matrices[transition_i, ,] <- (1 / (K - 1)) * (1 -
            (1 / K) -
            gammas[transition_i] * (1 - (1 / K)))
        diag(transition_matrices[transition_i, ,]) <- (1 / K) + gammas[transition_i] * (1 - (1 / K))
    }
    X <- SNPknock.models.sampleDMC(prob_initial, transition_matrices, n)
    return(list(covariates = X, prob_initial = prob_initial, transition_matrices = transition_matrices))
}

simulate_hidden_markov_states <- function(n, p, K) {
    stopifnot(K == 9)
    prob_initial <- c(1, rep(0, K - 1)) # start at inital state 1
    gam <- 0.35 # to match ln 223-225 p9 supplement of gene hunting

    ## base transition matrix construction
    base_transition_matrix <- matrix(0, nrow = K, ncol = K)
    for (l in 1 : K) {
        for (k in 1 : K) {
            if ((k - 1) == l %% K) {
                base_transition_matrix[l, k] <- 0.1
            }
        }
    }
    diag(base_transition_matrix) = 0.9

    stopifnot(rowSums(base_transition_matrix) == rep(1, K))

    ## base emissions matrix construction
    base_hidden_emission_matrix <- matrix(0, nrow = K, ncol = K)
    xs <- - 4 : 4
    zs <- 0 : (K - 1)
    for (z in 1 : K) {
        for (x in 1 : K) {
            if (((xs[x] + 4) == zs[z]) || ((xs[x] + 4) == zs[z] + 1)) {
                base_hidden_emission_matrix[z, x] <- gam / 2
            } else if (((xs[x] + 4) == 0) && (zs[z] == (K - 1))) {
                base_hidden_emission_matrix[z, x] <- gam / 2
            } else {
                base_hidden_emission_matrix[z, x] <- (1 - gam) / (K - 2)
            }
        }
    }

    stopifnot(rowSums(base_hidden_emission_matrix) == rep(1, K))

    transition_matrices <- array(data = NA, dim = c((p - 1), K, K)) # tensor to store each transition matrix
    hidden_emission_matrices <- array(data = NA, dim = c(p, K, K)) # tensor to store each transition matrix
    for (transition_i in 1 : (p - 1)) {
        transition_matrices[transition_i, ,] <- base_transition_matrix
        hidden_emission_matrices[transition_i, ,] <- base_hidden_emission_matrix
    }
    hidden_emission_matrices[p, ,] <- base_hidden_emission_matrix

    X = SNPknock.models.sampleHMM(prob_initial, transition_matrices, hidden_emission_matrices, n)

    return(list(covariates = X, prob_initial = prob_initial, transition_matrices = transition_matrices,
    hidden_emission_matrices = hidden_emission_matrices))
}

simulate_emissions <- function(n, p, signal_amplitude, active_support_size, covariates, covariate_state_offset) {
    active_support <- sample.int(p, size = active_support_size)
    model_coefs <- numeric(p)
    model_coefs[active_support] <- signal_amplitude / sqrt(n)
    # transform states to by offset
    prob_success <- expit((covariates + covariate_state_offset) %*% model_coefs)
    responses <- rbinom(n, 1, prob = prob_success)
    return(list(responses = responses, active_support = active_support))
}

compute_variable_importance <- function(responses, covariates, covariate_state_offset, covariate_knockoffs, variable_variable_importance_measure) {
    if (variable_variable_importance_measure == 'coefdiff') {
        variable_importance = stat.glmnet_coefdiff(covariates + covariate_state_offset,
        covariate_knockoffs + covariate_state_offset,
        responses,
        family = "binomial")
    } else if (variable_importance_measure == 'tree') {
        variable_importance = stat.random_forest(covariates + covariate_state_offset,
        covariate_knockoffs + covariate_state_offset,
        as.factor(responses))
    } else {
        stop(paste0("variable importance measure ", variable_importance_measure, " not implemented"))
    }
    return(variable_importance)
}

compute_knockoff_filter_dmc <- function(responses, covariates, prob_initial, transition_matrices,
covariate_state_offset = - 2, fdr = 0.1, offset = 0, variable_importance_measure){

    covariate_knockoffs <- SNPknock.knockoffDMC(covariates, prob_initial, transition_matrices, display_progress = FALSE)
    variable_importance <- compute_variable_importance(responses, covariates, covariate_state_offset, covariate_knockoffs, variable_importance_measure)
    adaptive_threshold <- knockoff.threshold(variable_importance, fdr, offset)
    return(list(variable_importance = variable_importance, adaptive_threshold = adaptive_threshold, covariate_knockoffs = covariate_knockoffs))
}

compute_knockoff_filter_gaussian <- function(responses, covariates, covariate_state_offset = - 2, fdr = 0.1, offset = 0, variable_importance_measure) {
    covariate_knockoffs <- create.second_order(covariates)
    variable_importance <- compute_variable_importance(responses, covariates, covariate_state_offset, covariate_knockoffs, variable_importance_measure)
    adaptive_threshold <- knockoff.threshold(variable_importance, fdr, offset)
    return(list(variable_importance = variable_importance, adaptive_threshold = adaptive_threshold))
}

compute_knockoff_filter_hmm <- function(responses, covariates, prob_initial, transition_matrices, hidden_emission_matrices,
covariate_state_offset = - 4, fdr = 0.1, offset = 0, variable_importance_measure) {

    covariate_knockoffs <- SNPknock.knockoffHMM(covariates, prob_initial, transition_matrices,
    hidden_emission_matrices, display_progress = F)
    variable_importance <- compute_variable_importance(responses, covariates, covariate_state_offset, covariate_knockoffs, variable_importance_measure)
    adaptive_threshold <- knockoff.threshold(variable_importance, fdr, offset)
    return(list(variable_importance = variable_importance, adaptive_threshold = adaptive_threshold, covariate_knockoffs = covariate_knockoffs))
}

analyze_knockoffs <- function(variable_importance, adaptive_threshold, active_support) {
    discoveries <- which(variable_importance >= adaptive_threshold)
    n_discoveries <- length(discoveries)
    n_true_discoveries <- length(intersect(discoveries, active_support))
    n_false_discoveries <- n_discoveries - n_true_discoveries
    FDP <- n_false_discoveries / max(c(1, n_discoveries))
    Power <- n_true_discoveries / length(active_support)
    return(list(fdp = FDP, pow = Power))
}

simulate_full_independence <- function(number_simulations, n, p, K, active_support_size, signal_amplitude,
covariate_state_offset, fdr = 0.1, offset = 1, covariate_type,
knockoff_type, variable_importance_measure) {

    fdps <- numeric(number_simulations)
    pows <- numeric(number_simulations)

    for (experiment_i in 1 : number_simulations) {

        if (covariate_type == 'hmm') {
            markov_states <- simulate_hidden_markov_states(n, p, K)
        } else {
            markov_states <- simulate_discrete_markov_states(n, p, K)
        }

        emiss <- simulate_emissions(n, p, signal_amplitude, active_support_size, markov_states$covariates, covariate_state_offset)

        if (knockoff_type == 'dmc') {
            kfilt <- compute_knockoff_filter_dmc(emiss$responses, markov_states$covariates, markov_states$prob_initial,
            markov_states$transition_matrices, covariate_state_offset, fdr, offset, variable_importance_measure)
        } else if (knockoff_type == 'gaussian') {
            kfilt <- compute_knockoff_filter_gaussian(emiss$responses, markov_states$covariates,
            covariate_state_offset, fdr, offset, variable_importance_measure)
        } else if (knockoff_type == 'hmm') {
            kfilt <- compute_knockoff_filter_hmm(emiss$responses, markov_states$covariates, markov_states$prob_initial,
            markov_states$transition_matrices, markov_states$hidden_emission_matrices,
            covariate_state_offset, fdr, offset, variable_importance_measure)
        } else {
            stop(paste0("Knockoff type must be either dmc, hmm, or gaussian. Current value is: ", knockoff_type))
        }
        out <- analyze_knockoffs(kfilt$variable_importance, kfilt$adaptive_threshold, emiss$active_support)
        fdps[experiment_i] <- out$fdp
        pows[experiment_i] <- out$pow
    }

    return(list(fdps = fdps, pows = pows))
}

## knockoff diag. to check eqn (3.1) of Candès et al. 2018
## From tutorial https://web.stanford.edu/group/candes/knockoffs/tutorials/gwas_tutorial.html#knockoff_diagnostics
compute_covariate_knockoff_corr <- function(covariate, knockoff) {
    covariate_corr = sapply(2:p, function(j) cor(covariate[, j - 1],covariate[, j]))
    knockoff_corr = sapply(2:p, function(j) cor(knockoff[, j - 1],knockoff[, j]))
    return(data.frame(covariate_corr = covariate_corr, knockoff_corr = knockoff_corr, feature = 1:length(knockoff_corr)))
}