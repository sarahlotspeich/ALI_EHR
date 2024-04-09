# Load library --------------------------------------------
# library(sleev) ## for SMLE logistic regression
source("~/Documents/logreg2phRonly/R/logreg2ph.R")
source("~/Documents/logreg2phRonly/R/hessian_row.R")
source("~/Documents/logreg2phRonly/R/observed_data_loglik.R")
source("~/Documents/logreg2phRonly/R/pl_theta.R")
source("~/Documents/logreg2phRonly/R/profile_out.R")

# Load library --------------------------------------------
library(auditDesignR) ## for audit sampling functions

# Set parameters that won't be varied in the loop
## These values will be set as the defaults in the sim_data() function for convenience
lambda_age = 45.66 ## mean of Poisson for Z
beta0 = -1.57 ## intercept in model of Y|X,Z
beta1 = 0.95 ## coefficient on X in model of Y|X,Z
beta2 = 0.01 ## coefficient on Z in model of Y|X,Z
pS = c(0.2500000, 0.9870130, 0.4549098, 0.1450000, 0.0580000, 
       0.2490119, 0.3138501, 0.3316391, 0.3111111, 0.0000000) ## probability of stressor = YES
pM = c(0.996, 0.153, 0.002, 0.000, 0.000, 
       0.494, 0.213, 0.213, 0.955, 0.983) ## probability of stressor = NA
N = 1000 ## total sample size (Phase I)
pV = 0.1 ## proportion of patients to be validated (Phase II)
nsieve = 16 ## number of B-spline sieves
# tpr = 0.95 ## true positive rate of EHR components
# fpr = 0.05 ## false positive rate of EHR components
audit_recovery = 1 ## proportion of missing data recovered through the audit

# Function to simulate data
sim_data = function(tpr, fpr) {
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
  Xstar1 = rowMeans(Sstar_mat, 
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
  Xstar2 = rowMeans(S_mat, 
                    na.rm = TRUE) ## ALI based on recovered components

  ## Create dataset
  dat = data.frame(id = 1:N, X, Xstar1, Xstar2, Y, Z)
  
  # Return dataset
  return(dat)
}

# Function to get validation indicators
get_V = function(design, data) {
  if (design == "SRS") {
    sample_srs(phI = N, 
               phII = ceiling(pV * N))
  } else if (design == "CC") {
    sample_cc(dat = data, 
              phI = N, 
              phII = ceiling(pV * N), 
              sample_on = "Y")
  } else if (design == "BCC") {
    data$Xstar_strat = as.numeric(data$Xstar1 <= 0.5)
    sample_bcc(dat = data,
               phI = N, 
               phII = ceiling(pV * N), 
               sample_on = c("Y", "Xstar_strat"))
  } else if (design == "RESID") {
    sample_resid(formula = Y ~ Xstar1 + Z, 
                 family = "binomial", 
                 dat = data, 
                 phI = N, 
                 phII = ceiling(pV * N))
  } else if (design == "SFS") {
    sample_sfs(formula = Y ~ Xstar1 + Z, 
               family = "binomial", 
               dat = data, 
               phI = N, 
               phII = ceiling(pV * N),
               X = "Xstar1")
  }
}

# Function to simulate data and then fit all models
sim_data_fit = function(id, tpr, fpr) {
  results = data.frame(sim = id, 
                       srs_beta0 = NA, srs_beta1 = NA, srs_beta2 = NA, srs_conv_msg = NA, srs_resampled_V = FALSE, 
                       cc_beta0 = NA, cc_beta1 = NA, cc_beta2 = NA, cc_conv_msg = NA, cc_resampled_V = FALSE, 
                       bcc_beta0 = NA, bcc_beta1 = NA, bcc_beta2 = NA, bcc_conv_msg = NA, bcc_resampled_V = FALSE, 
                       resid_beta0 = NA, resid_beta1 = NA, resid_beta2 = NA, resid_conv_msg = NA, resid_resampled_V = FALSE)
  
  # Simulate data 
  temp = sim_data(tpr = tpr, 
                  fpr = fpr) 
  
  # Setup B-splines
  B = splines::bs(x = temp$Xstar1, 
                  df = nsieve, 
                  Boundary.knots = range(temp$Xstar1), 
                  intercept = TRUE)
  colnames(B) = paste0("bs", seq(1, nsieve))
  temp = cbind(temp, B)
  
  # Loop over the different designs
  for (des in c("srs", "cc", "bcc", "resid", "sfs")) {
    #print(des)
    
    ## Create validation indicator 
    temp$V = get_V(design = toupper(des), 
                   data = temp)
    
    ## Check for empty sieves in B-splines of validated data
    # while(0 %in% colSums(temp[which(temp$V == 1), paste0("bs", 1:nsieve)])) {
    #   ## Re-create validation indicator 
    #   temp$V = get_V(design = toupper(des), 
    #                  data = temp)
    #   
    #   ## Validation indicators re-sampled
    #   results[1, paste0(des, "_resampled_V")] = TRUE
    # }
    
    ## Create Xmiss to be NA if V = 0
    temp = temp |> 
      dplyr::mutate(Xmiss = dplyr::if_else(condition = V == 0, 
                                           true = NA, 
                                           false = Xstar2))
    
    ## Fit SMLE
    fit = sleev::log
    
    suppressMessages(fit <- logreg2ph(
      Y_unval = NULL,
      Y_val = "Y",
      X_unval = "Xstar1",
      X_val = "Xmiss",
      C = "Z",
      Validated = "V",
      Bspline = colnames(B),
      data = temp,
      h_N_scale = 1,
      noSE = FALSE,
      TOL = 1e-04,
      MAX_ITER = 1000
    ))
    results[1, paste0(des, "_beta", 0:2)] = fit$model_coeff$coeff
    results[1, paste0(des, "_conv_msg")] = fit$converged_msg
  }
  
  # Return results
  return(results)
}

# Error setting 1 (extra light): 99% TPR/1% FPR
set.seed(918) # Be reproducible
tpr99_fpr01 = do.call(what = rbind,
                      args = pbapply::pbsapply(X = 1:1000,
                                               FUN = sim_data_fit, 
                                               simplify = FALSE,
                                               tpr = 0.99,
                                               fpr = 0.01)) |> 
  dplyr::mutate(tpr = 0.99, 
                fpr = 0.01)
tpr99_fpr01 |> 
  write.csv("ALI_EHR/sim-data/vary_designs_tpr99_fpr01.csv", 
            row.names = FALSE)

# Error setting 2 (light): 95% TPR/5% FPR
set.seed(918) # Be reproducible
tpr95_fpr05 = do.call(what = rbind,
                      args = pbapply::pbsapply(X = 1:1000,
                                               FUN = sim_data_fit, 
                                               simplify = FALSE,
                                               tpr = 0.95,
                                               fpr = 0.05)) |> 
  dplyr::mutate(tpr = 0.95, 
                fpr = 0.05)
tpr95_fpr05 |> 
  write.csv("ALI_EHR/sim-data/vary_designs_tpr95_fpr05.csv", 
            row.names = FALSE)

# Error setting 3 (moderate): 80% TPR/20% FPR
set.seed(918) # Be reproducible
tpr80_fpr20 = do.call(what = rbind,
                      args = pbapply::pbsapply(X = 1:1000,
                                               FUN = sim_data_fit, 
                                               simplify = FALSE,
                                               tpr = 0.80,
                                               fpr = 0.20)) |> 
  dplyr::mutate(tpr = 0.8, 
                fpr = 0.2)
tpr80_fpr20 |> 
  write.csv("ALI_EHR/sim-data/vary_designs_tpr80_fpr20.csv", 
            row.names = FALSE)

# Error setting 4 (heavy): 50% TPR/50% FPR
set.seed(918) # Be reproducible
tpr50_fpr50 = do.call(what = rbind,
                      args = pbapply::pbsapply(X = 1:1000,
                                               FUN = sim_data_fit, 
                                               simplify = FALSE,
                                               tpr = 0.50,
                                               fpr = 0.50)) |> 
  dplyr::mutate(tpr = 0.5, 
                fpr = 0.5)
tpr50_fpr50 |> 
  write.csv("ALI_EHR/sim-data/vary_designs_tpr50_fpr50.csv", 
            row.names = FALSE)
