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
audit_recovery = 1 ## proportion of missing data recovered through the audit

# Random seed to be used for each simulation setting
args = commandArgs(TRUE)
## When running on the cluster, give each array a unique seed by adding the array ID to 11422
sim_seed = 11422 #+ as.integer(args)

reps = 1000 ## Number of replications

# Function to simulate data
sim_data = function(audit_recovery = 1, tpr, fpr) {
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
sim_data_fit = function(id, tpr, fpr) {
  results = data.frame(sim = id, 
                       bcc_beta0 = NA, bcc_beta1 = NA, bcc_beta2 = NA, 
                       bcc_conv_msg = NA, bcc_data_resampled = FALSE)
  
  # Simulate data 
  temp = sim_data(tpr = tpr, 
                  fpr = fpr) 
  
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

  ## Create Wave 2 validation indicator 
  temp$V2 = sample_bcc(dat = temp,
                       phI = nrow(temp), 
                       phII = (n / 2), 
                       sample_on = c("Y", "Xstar_strat"), 
                       wave1_Validated = temp$V1 == 1)
  
  ## Create any wave validation indicator
  temp$V = pmax(temp$V1, temp$V2)
  
  ## Check for empty sieves in validated data
  sieve_sums = temp %>%
    dplyr::filter(V == 1) %>% 
    dplyr::select(dplyr::starts_with("bs")) %>%
    colSums()
  
  ## If any sieves are empty, generate new data re-sample waves 1 and 2
  while(any(sieve_sums == 0)) {
    ## Save logical indicator that validation study was re-sampled
    results[1, "bcc_data_resampled"] = TRUE
    
    # Simulate data 
    temp = sim_data(tpr = tpr, 
                    fpr = fpr) 
    
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
    
    ## Create Wave 2 validation indicator 
    temp$V2 = sample_bcc(dat = temp,
                         phI = nrow(temp), 
                         phII = (n / 2), 
                         sample_on = c("Y", "Xstar_strat"), 
                         wave1_Validated = temp$V1 == 1)
    
    ## Create any wave validation indicator
    temp$V = pmax(temp$V1, temp$V2)
    
    ## Check for empty sieves in validated data
    sieve_sums = temp %>%
      dplyr::filter(V == 1) %>% 
      dplyr::select(dplyr::starts_with("bs")) %>%
      colSums()
  }
  
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
  results[1, paste0("bcc_beta", 0:2)] = fit$coefficients$Estimate
  results[1, "bcc_conv_msg"] = fit$converge_msg
  
  # Return results
  return(results)
}

# Error setting 1 (extra light): 99% TPR/1% FPR
set.seed(sim_seed) # Be reproducible
tpr99_fpr01 = do.call(what = rbind,
                      args = sapply(X = 1:reps,
                                    FUN = sim_data_fit, 
                                    simplify = FALSE,
                                    tpr = 0.99,
                                    fpr = 0.01)) %>% 
  dplyr::mutate(tpr = 0.99, 
                fpr = 0.01)
tpr99_fpr01 %>% 
  write.csv(file = "vary_wave2_design/wave2_bcc_tpr99_fpr01.csv", 
            row.names = FALSE)

# Error setting 2 (light): 95% TPR/5% FPR
set.seed(sim_seed) # Be reproducible
tpr95_fpr05 = do.call(what = rbind,
                      args = sapply(X = 1:reps,
                                    FUN = sim_data_fit, 
                                    simplify = FALSE,
                                    tpr = 0.95,
                                    fpr = 0.05)) %>% 
  dplyr::mutate(tpr = 0.95, 
                fpr = 0.05)
tpr95_fpr05 %>% 
  write.csv(file = "vary_wave2_design/wave2_bcc_tpr95_fpr05_seed.csv", 
            row.names = FALSE)

# Error setting 3 (moderate): 80% TPR/20% FPR
set.seed(sim_seed) # Be reproducible
tpr80_fpr20 = do.call(what = rbind,
                      args = sapply(X = 1:reps,
                                    FUN = sim_data_fit, 
                                    simplify = FALSE,
                                    tpr = 0.80,
                                    fpr = 0.20)) %>% 
  dplyr::mutate(tpr = 0.8, 
                fpr = 0.2)
tpr80_fpr20 %>% 
  write.csv(file = "vary_wave2_design/wave2_bcc_tpr80_fpr20_seed.csv", 
            row.names = FALSE)

# Error setting 4 (heavy): 50% TPR/50% FPR
set.seed(sim_seed) # Be reproducible
tpr50_fpr50 = do.call(what = rbind,
                      args = sapply(X = 1:reps,
                                    FUN = sim_data_fit, 
                                    simplify = FALSE,
                                    tpr = 0.50,
                                    fpr = 0.50)) %>% 
  dplyr::mutate(tpr = 0.5, 
                fpr = 0.5)
tpr50_fpr50 %>% 
  write.csv(file = "vary_wave2_design/wave2_bcc_tpr50_fpr50_seed.csv", 
            row.names = FALSE)
