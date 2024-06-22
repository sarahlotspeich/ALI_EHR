# Load library --------------------------------------------
## Need to use the dev version of sleev: 
## devtools::install_github("dragontaoran/sleev")
library(sleev, lib.loc = "/home/lotspes/packages/") ## for SMLE logistic regression
library(magrittr) ## for les pipes %>%
library(auditDesignR, lib.loc = "/home/lotspes/packages/") ## for audit sampling functions

# Set parameters that won't be varied in the loop
## These values will be set as the defaults in the sim_data() function for convenience
lambda_age = 45.66 / 10 ## mean of Poisson for Z
beta0 = -1.93 ## intercept in model of Y|X,Z
beta1 = 1.88 ## coefficient on X in model of Y|X,Z
beta2 = 0.10 ## coefficient on Z in model of Y|X,Z
pS = c(0.2500000, 0.9870130, 0.4549098, 0.1450000, 0.0580000, 
       0.2490119, 0.3138501, 0.3316391, 0.3111111, 0.0000000) ## probability of stressor = YES
pM = c(0.996, 0.153, 0.002, 0.000, 0.000, 
       0.494, 0.213, 0.213, 0.955, 0.983) ## probability of stressor = NA
N = 1000 ## total sample size (Phase I)
pV = 0.1 ## proportion of patients to be validated (Phase II)
n = ceiling(pV * N) ## number of patients to be validated (Phase II)
nsieve = 16 ## number of B-spline sieves

# Random seed to be used for each simulation setting
args = commandArgs(TRUE)
## When running on the cluster, give each array a unique seed by adding the array ID to 11422
sim_seed = 11422 #+ as.integer(args)

reps = 1000 ## Number of replications

# Function to simulate data
sim_data = function(audit_recovery, tpr, fpr) {
  ## Simulate continuous error-free covariate: age at first encounter 
  ### from Poisson(lambda_age) 
  Z = rpois(n = N, 
            lambda = lambda_age)
  
  ## Simulate error-free (validated) version of error-prone covariate: ALI 
  ### Begin with stress indicators (50 per person) from Bernoulli (pS)
  S = rbinom(n = N * 10,
             size = 1, 
             prob = rep(x = pS, times = N))
  S_mat = matrix(data = S, 
                 nrow = N, 
                 ncol = 10, 
                 byrow = TRUE)
  X = rowMeans(S_mat)
  
  ## Simulate error-prone (EHR) version of the covariate: ALI* 
  ### Begin with stress indicators (10 per person) from Bernoulli (1 / (1 + exp(-(gamma0 + gamma1 S))))
  gamma0 = - log((1 - fpr) / fpr) ### define intercept such that P(S* = 0|S = 0) = FPR
  gamma1 = - log((1 - tpr) / tpr) - gamma0 ### define slope such that P(S* = 1|S = 1) = TPR
  Sstar = rbinom(n = N * 10, 
                 size = 1, 
                 prob = 1 / (1 + exp(- (gamma0 + gamma1 * S))))
  
  ### Simulate missingness in EHR version of the covariate
  Sstar_miss = rbinom(n = N * 10, 
                      size = 1, 
                      prob = rep(x = pM, times = N))
  Sstar[which(Sstar_miss == 1)] = NA
  Sstar_mat = matrix(data = Sstar, 
                     nrow = N, 
                     ncol = 10, 
                     byrow = TRUE)
  Xstar = rowMeans(Sstar_mat, 
                    na.rm = TRUE)
  
  ## Simulate outcome: healthcare utilization
  Y = rbinom(n = N, 
             size = 1, 
             prob = (1 + exp(-(beta0 + beta1 * X - beta2 * Z))) ^ (- 1))
  
  ## Simulate imperfect audit recovery
  recovered = sample(x = which(Sstar_miss == 1), 
                     size = audit_recovery * length(which(Sstar_miss == 1)), 
                     replace = FALSE)
  not_recovered = setdiff(x = which(Sstar_miss == 1), 
                          y = recovered)
  S[not_recovered] = NA ## components not recovered in audit
  S_mat = matrix(data = S, 
                 nrow = N, 
                 ncol = 10, 
                 byrow = TRUE)
  Xval = rowMeans(S_mat, 
                    na.rm = TRUE) ## ALI based on recovered components
  
  ## Create dataset
  dat = data.frame(id = 1:N, X, Xstar, Xval, Y, Z)
  
  # Return dataset
  return(dat)
}

# Function to simulate data and then fit all models
sim_data_fit = function(id, audit_recovery, tpr = 0.95, fpr = 0.05, q = 0.5) {
  results = data.frame(sim = id, 
                       optMLE2_beta0 = NA, optMLE2_beta1 = NA, optMLE2_beta2 = NA, 
                       optMLE2_conv_msg = NA, optMLE_wave1_est = "MLE", 
                       optMLE_data_resampled = FALSE, optMLE_sing_inf = FALSE)
  
  # Simulate data 
  temp = sim_data(tpr = tpr, 
                  fpr = fpr, 
                  audit_recovery = audit_recovery) 
  
  # Setup B-splines
  B = splines::bs(x = temp$Xstar, 
                  df = nsieve, 
                  Boundary.knots = range(temp$Xstar), 
                  intercept = TRUE)
  colnames(B) = paste0("bs", seq(1, nsieve))
  temp = cbind(temp, B)
  
  # Draw a BCC* of n / 2 in the first wave 
  ## Stratify X* at the median
  temp$Xstar_strat = as.numeric(temp$Xstar <= median(temp$Xstar))
  ## Create Wave 1 validation indicator
  temp$V1 = sample_bcc(dat = temp,
                       phI = nrow(temp), 
                       phII = (n / 2), 
                       sample_on = c("Y", "Xstar_strat"), 
                       wave1_Validated = NULL)
  
  # Create stratified versions of Phase I variables  -------------------------
  temp$Z_strat = as.numeric(temp$Z <= lambda_age)
  temp$X_strat = as.numeric(temp$X <= quantile(x = temp$Xstar, probs = q) )
  temp$Xstar_strat = as.numeric(temp$Xstar <= quantile(x = temp$X, probs = q))
  
  # Create matrix of complete data  ------------------------------------------
  complete_data = expand.grid(Y = c(0, 1), 
                              X_strat = c(0, 1), 
                              Xstar_strat = c(0, 1),
                              Z_strat = c(0, 1),
                              V = c(0, 1))
  
  # Cross-tabulate Phase I strata --------------------------------------------
  N00 = with(temp, table(Y, Xstar_strat))[1,1]
  N01 = with(temp, table(Y, Xstar_strat))[1,2]
  N10 = with(temp, table(Y, Xstar_strat))[2,1]
  N11 = with(temp, table(Y, Xstar_strat))[2,2]
  stratN = list(N00 = N00, N01 = N01, N10 = N10, N11 = N11)
  
  # Sample wave 1 audits and estimate parameters and score  ------------------
  wave1_strat = data.frame(n00 = with(temp[temp$V1 == 1, ], table(Y, Xstar_strat))[1, 1],
                           n01 = with(temp[temp$V1 == 1, ], table(Y, Xstar_strat))[1, 2],
                           n10 = with(temp[temp$V1 == 1, ], table(Y, Xstar_strat))[2, 1],
                           n11 = with(temp[temp$V1 == 1, ], table(Y, Xstar_strat))[2, 2])
  mle_wave1 = twophase_mle(dat = temp, 
                           Y_val = "Y", 
                           Y_unval = NULL, 
                           X_val = "X_strat", 
                           X_unval = "Xstar_strat", 
                           addl_covar = "Z_strat",
                           Validated = "V1", 
                           nondiff_X_unval = TRUE)
  
  if(mle_wave1$conv) {
    beta_hat = mle_wave1$mod_Y_val$Est[3]
    eta_hat = with(mle_wave1, 
                   c(mod_Y_val$Est[1:2], 
                     mod_Y_unval$Est, 
                     mod_X_unval$Est, 
                     mod_X_val$Est))
  } else {
    ## If MLEs don't converge, use complete case estimates -------------------
    mod_Y_val = glm(formula = Y ~ Z_strat + X_strat, 
                    data = temp, 
                    family = "binomial", 
                    subset = V1 == 1)
    mod_X_unval = glm(formula = Xstar_strat ~ X_strat + Z_strat, 
                      data = temp, 
                      family = "binomial", 
                      subset = V1 == 1)
    mod_X_val = glm(formula =  X_strat ~ Z_strat, 
                    data = temp, 
                    family = "binomial", 
                    subset = V1 == 1)
    beta_hat = mod_Y_val$coefficients[3]
    eta_hat = c(mod_Y_val$coefficients[1:2], 
                mod_X_unval$coefficients, 
                mod_X_val$coefficients)
    results[1, "optMLE_wave1_est"] = "CC"
  }
  s_hat = score(comp_dat = complete_data, 
                Y_val = "Y", 
                Y_unval = NULL, 
                X_val = "X_strat", 
                X_unval = "Xstar_strat", 
                addl_covar = "Z_strat",
                Validated = "V", 
                beta = beta_hat, 
                eta = eta_hat, 
                nondiff_X_unval = TRUE)
  
  grid_search = optMLE_grid(phI = nrow(temp), 
                            phII = (n / 2), 
                            phI_strat = stratN, 
                            phIIa_strat = wave1_strat, 
                            min_n = 0, 
                            sample_on = c("Y", "Xstar_strat"), 
                            indiv_score = s_hat)
  
  ## Create Wave 2 validation indicator 
  if (grid_search$findOptimal) {
    opt_des2 = grid_search$min_var_design
    opt_des2[, c("n00", "n01", "n10", "n11")] = opt_des2[, c("n00", "n01", "n10", "n11")] - 
      with(wave1_strat, c(n00, n01, n10, n11))
    temp$V2 = sample_optMLE(dat = temp, 
                            sample_on = c("Y", "Xstar_strat"), 
                            des = opt_des2, 
                            wave1_Validated = "V1")
    
    ## Create any wave validation indicator
    temp$V = pmax(temp$V1, temp$V2)
    
    ## Check for empty sieves in validated data
    sieve_sums = temp %>%
      dplyr::filter(V == 1) %>% 
      dplyr::select(dplyr::starts_with("bs")) %>%
      colSums()
    
    ## If any sieves are empty, re-sample waves 1 and 2
    while(any(sieve_sums == 0)) {
      ## Save logical indicator that validation study was re-sampled
      results[1, "optMLE_data_resampled"] = TRUE
      
      # Simulate data 
      temp = sim_data(tpr = tpr, 
                      fpr = fpr, 
                      audit_recovery = audit_recovery) 
      
      # Setup B-splines
      B = splines::bs(x = temp$Xstar, 
                      df = nsieve, 
                      Boundary.knots = range(temp$Xstar), 
                      intercept = TRUE)
      colnames(B) = paste0("bs", seq(1, nsieve))
      temp = cbind(temp, B)
      
      # Draw a BCC* of n / 2 in the first wave 
      ## Stratify X* at the median
      temp$Xstar_strat = as.numeric(temp$Xstar <= median(temp$Xstar))
      ## Create Wave 1 validation indicator
      temp$V1 = sample_bcc(dat = temp,
                           phI = nrow(temp), 
                           phII = (n / 2), 
                           sample_on = c("Y", "Xstar_strat"), 
                           wave1_Validated = NULL)
      
      # Create stratified versions of Phase I variables  -------------------------
      temp$Z_strat = as.numeric(temp$Z <= lambda_age)
      temp$X_strat = as.numeric(temp$X <= quantile(x = temp$Xstar, probs = q) )
      temp$Xstar_strat = as.numeric(temp$Xstar <= quantile(x = temp$X, probs = q))
      
      # Create matrix of complete data  ------------------------------------------
      complete_data = expand.grid(Y = c(0, 1), 
                                  X_strat = c(0, 1), 
                                  Xstar_strat = c(0, 1),
                                  Z_strat = c(0, 1),
                                  V = c(0, 1))
      
      # Cross-tabulate Phase I strata --------------------------------------------
      N00 = with(temp, table(Y, Xstar_strat))[1,1]
      N01 = with(temp, table(Y, Xstar_strat))[1,2]
      N10 = with(temp, table(Y, Xstar_strat))[2,1]
      N11 = with(temp, table(Y, Xstar_strat))[2,2]
      stratN = list(N00 = N00, N01 = N01, N10 = N10, N11 = N11)
      
      # Sample wave 1 audits and estimate parameters and score  ------------------
      wave1_strat = data.frame(n00 = with(temp[temp$V1 == 1, ], table(Y, Xstar_strat))[1, 1],
                               n01 = with(temp[temp$V1 == 1, ], table(Y, Xstar_strat))[1, 2],
                               n10 = with(temp[temp$V1 == 1, ], table(Y, Xstar_strat))[2, 1],
                               n11 = with(temp[temp$V1 == 1, ], table(Y, Xstar_strat))[2, 2])
      mle_wave1 = twophase_mle(dat = temp, 
                               Y_val = "Y", 
                               Y_unval = NULL, 
                               X_val = "X_strat", 
                               X_unval = "Xstar_strat", 
                               addl_covar = "Z_strat",
                               Validated = "V1", 
                               nondiff_X_unval = TRUE)
      
      if(mle_wave1$conv) {
        beta_hat = mle_wave1$mod_Y_val$Est[3]
        eta_hat = with(mle_wave1, 
                       c(mod_Y_val$Est[1:2], 
                         mod_Y_unval$Est, 
                         mod_X_unval$Est, 
                         mod_X_val$Est))
      } else {
        ## If MLEs don't converge, use complete case estimates -------------------
        mod_Y_val = glm(formula = Y ~ Z_strat + X_strat, 
                        data = temp, 
                        family = "binomial", 
                        subset = V1 == 1)
        mod_X_unval = glm(formula = Xstar_strat ~ X_strat + Z_strat, 
                          data = temp, 
                          family = "binomial", 
                          subset = V1 == 1)
        mod_X_val = glm(formula =  X_strat ~ Z_strat, 
                        data = temp, 
                        family = "binomial", 
                        subset = V1 == 1)
        beta_hat = mod_Y_val$coefficients[3]
        eta_hat = c(mod_Y_val$coefficients[1:2], 
                    mod_X_unval$coefficients, 
                    mod_X_val$coefficients)
        results[1, "optMLE_wave1_est"] = "CC"
      }
      s_hat = score(comp_dat = complete_data, 
                    Y_val = "Y", 
                    Y_unval = NULL, 
                    X_val = "X_strat", 
                    X_unval = "Xstar_strat", 
                    addl_covar = "Z_strat",
                    Validated = "V", 
                    beta = beta_hat, 
                    eta = eta_hat, 
                    nondiff_X_unval = TRUE)
      
      grid_search = optMLE_grid(phI = nrow(temp), 
                                phII = (n / 2), 
                                phI_strat = stratN, 
                                phIIa_strat = wave1_strat, 
                                min_n = 0, 
                                sample_on = c("Y", "Xstar_strat"), 
                                indiv_score = s_hat)
      
      ## Create Wave 2 validation indicator 
      if (grid_search$findOptimal) {
        opt_des2 = grid_search$min_var_design
        opt_des2[, c("n00", "n01", "n10", "n11")] = opt_des2[, c("n00", "n01", "n10", "n11")] - 
          with(wave1_strat, c(n00, n01, n10, n11))
        temp$V2 = sample_optMLE(dat = temp, 
                                sample_on = c("Y", "Xstar_strat"), 
                                des = opt_des2, 
                                wave1_Validated = "V1")
        
        ## Create any wave validation indicator
        temp$V = pmax(temp$V1, temp$V2)
        
        ## Check for empty sieves in validated data
        sieve_sums = temp %>%
          dplyr::filter(V == 1) %>% 
          dplyr::select(dplyr::starts_with("bs")) %>%
          colSums()
      } else {
        results[1, "optMLE_sing_inf"] = TRUE
      }
    }
  }
  
  if (grid_search$findOptimal) {
    ## Create Xmiss to be NA if V = 0
    temp = temp %>% 
      dplyr::mutate(Xmiss = dplyr::if_else(condition = V == 0, 
                                           true = NA, 
                                           false = Xval))
    
    ## Fit SMLE
    fit = logistic2ph(
      Y_unval = NULL,
      Y = "Y",
      X_unval = "Xstar",
      X = "Xmiss",
      Z = "Z",
      Bspline = colnames(B),
      data = temp,
      hn_scale = 1,
      noSE = TRUE,
      TOL = 1e-04,
      MAX_ITER = 1000
    )
    results[1, paste0("optMLE2_beta", 0:2)] = fit$coefficients$Estimate
    results[1, "optMLE2_conv_msg"] = fit$converge_msg
  } else {
    results[1, "optMLE_sing_inf"] = TRUE
  }
  
  # Return results
  return(results)
}

# Loop over different audit recovery rates (high, moderate, low)
for (r in c(0.9, 0.5, 0.25)) {
  set.seed(sim_seed) # Be reproducible
  sims = do.call(what = rbind,
                 args = sapply(X = 1:reps,
                               FUN = sim_data_fit, 
                               simplify = FALSE,
                               audit_recovery = r)) %>% 
    dplyr::mutate(audit_recovery = r)
  sims %>% 
    write.csv(file = paste0("vary_wave2_design/wave2_optMLE2_recover", 100 * r, ".csv"), 
              row.names = FALSE)
}