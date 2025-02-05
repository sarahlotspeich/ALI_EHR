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
n = ceiling(pV * N) ## number of patients to be validated (Phase II)
nsieve = 16 ## number of B-spline sieves
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
get_V = function(design, data, Xstar_cut = 0.5, X_cut = 0.5, Z_cut = 46) {
  if (design == "SRS") {
    sample_srs(phI = N, 
               phII = n)
  } else if (design == "CC") {
    sample_cc(dat = data, 
              phI = N, 
              phII = n, 
              sample_on = "Y")
  } else if (design == "BCC") {
    data$Xstar_strat = as.numeric(data$Xstar1 <= Xstar_cut)
    sample_bcc(dat = data,
               phI = N, 
               phII = n, 
               sample_on = c("Y", "Xstar_strat"))
  } else if (design == "RESID") {
    sample_resid(formula = Y ~ Xstar1 + Z, 
                 family = "binomial", 
                 dat = data, 
                 phI = N, 
                 phII = n)
  } else if (design == "SFS") {
    sample_sfs(formula = Y ~ Z, 
               family = "binomial", 
               dat = data, 
               phI = N, 
               phII = n,
               X = "Xstar1")
  } else if (design == "OPTMLE") {
    # Create stratified versions of Phase I variables  -------------------------
    data$Z_strat = as.numeric(data$Z <= Z_cut)
    data$X_strat = as.numeric(data$X <= X_cut)
    data$Xstar_strat = as.numeric(data$Xstar1 <= Xstar_cut)
    
    # Obtain full cohort estimates of parameters -------------------------------
    mle_wave1 = twophase_mle(dat = cbind(V = 1, data), 
                             Y_val = "Y", 
                             Y_unval = NULL, 
                             X_val = "X_strat", 
                             X_unval = "Xstar_strat", 
                             addl_covar = "Z_strat",
                             Validated = "V", 
                             nondiff_X_unval = TRUE)
    if(mle_wave1$conv) {
      beta = mle_wave1$mod_Y_val$Est[3]
      eta = with(mle_wave1, 
                 c(mod_Y_val$Est[1:2], 
                   mod_Y_unval$Est, 
                   mod_X_unval$Est, 
                   mod_X_val$Est))
      
      # Create matrix of complete data  ------------------------------------------
      complete_data = expand.grid(Y = c(0, 1), 
                                  X_strat = c(0, 1), 
                                  Xstar_strat = c(0, 1),
                                  Z_strat = c(0, 1),
                                  V = c(0, 1))
      
      # Calculate score based on full cohort estimates  --------------------------
      s = score(comp_dat = complete_data, 
                Y_val = "Y", 
                X_val = "X_strat", 
                X_unval = "Xstar_strat", 
                addl_covar = "Z_strat",
                Validated = "V", 
                beta = beta, 
                eta = eta, 
                nondiff_X_unval = TRUE)
      
      # Cross-tabulate Phase I strata --------------------------------------------
      N00 = with(data, table(Y, Xstar_strat))[1,1]
      N01 = with(data, table(Y, Xstar_strat))[1,2]
      N10 = with(data, table(Y, Xstar_strat))[2,1]
      N11 = with(data, table(Y, Xstar_strat))[2,2]
      stratN = list(N00 = N00, N01 = N01, N10 = N10, N11 = N11)
      
      # Grid search to find optimal Phase II strata  -----------------------------
      grid_search = optMLE_grid(phI = N, 
                                phII = n, 
                                phI_strat = stratN, 
                                min_n = 10, 
                                sample_on = c("Y", "Xstar_strat"), 
                                indiv_score = s)
      if (grid_search$findOptimal) {
        opt_des = grid_search$min_var_design
        sample_optMLE(dat = data, 
                      sample_on = c("Y", "Xstar_strat"), 
                      des = opt_des)
      } 
    } else {
      rep(NA, N)
    }
  } else if (design == "OPTMLE2") {
    # Create stratified versions of Phase I variables  -------------------------
    data$Z_strat = as.numeric(data$Z <= Z_cut)
    data$X_strat = as.numeric(data$X <= X_cut)
    data$Xstar_strat = as.numeric(data$Xstar1 <= Xstar_cut)
    
    # Create matrix of complete data  ------------------------------------------
    complete_data = expand.grid(Y = c(0, 1), 
                                X_strat = c(0, 1), 
                                Xstar_strat = c(0, 1),
                                Z_strat = c(0, 1),
                                V = c(0, 1))
    
    # Cross-tabulate Phase I strata --------------------------------------------
    N00 = with(data, table(Y, Xstar_strat))[1,1]
    N01 = with(data, table(Y, Xstar_strat))[1,2]
    N10 = with(data, table(Y, Xstar_strat))[2,1]
    N11 = with(data, table(Y, Xstar_strat))[2,2]
    stratN = list(N00 = N00, N01 = N01, N10 = N10, N11 = N11)
    
    # Sample wave 1 audits and estimate parameters and score  ------------------
    V_wave1 = sample_bcc(dat = data, 
                         phI = N, 
                         phII = (n / 2), 
                         sample_on = c("Y", "Xstar_strat"))
    wave1_strat = data.frame(n00 = with(data[V_wave1 == 1, ], table(Y, Xstar_strat))[1, 1],
                             n01 = with(data[V_wave1 == 1, ], table(Y, Xstar_strat))[1, 2],
                             n10 = with(data[V_wave1 == 1, ], table(Y, Xstar_strat))[2, 1],
                             n11 = with(data[V_wave1 == 1, ], table(Y, Xstar_strat))[2, 2])
    mle_wave1 = twophase_mle(dat = cbind(V = V_wave1, data), 
                             Y_val = "Y", 
                             Y_unval = NULL, 
                             X_val = "X_strat", 
                             X_unval = "Xstar_strat", 
                             addl_covar = "Z_strat",
                             Validated = "V", 
                             nondiff_X_unval = TRUE)
    
    if(mle_wave1$conv) {
      beta_hat = mle_wave1$mod_Y_val$Est[3]
      eta_hat = with(mle_wave1, 
                     c(mod_Y_val$Est[1:2], 
                       mod_Y_unval$Est, 
                       mod_X_unval$Est, 
                       mod_X_val$Est))
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
      
      grid_search = optMLE_grid(phI = N, 
                                phII = (n / 2), 
                                phI_strat = stratN, 
                                phIIa_strat = wave1_strat, 
                                min_n = 0, 
                                sample_on = c("Y", "Xstar_strat"), 
                                indiv_score = s_hat)
      
      if (grid_search$findOptimal) {
        opt_des2 = grid_search$min_var_design
        opt_des2[, c("n00", "n01", "n10", "n11")] = opt_des2[, c("n00", "n01", "n10", "n11")] - with(wave1_strat, c(n00, n01, n10, n11))
        pmax(V_wave1, 
             sample_optMLE(dat = cbind(V = V_wave1, data), 
                           sample_on = c("Y", "Xstar_strat"), 
                           des = opt_des2, 
                           wave1_Validated = "V"))
      }
    } else {
      rep(NA, N)
    }
  }
}

# Function to simulate data and then fit all models
sim_data_fit = function(id, tpr, fpr) {
  results = data.frame(sim = id, 
                       srs_beta0 = NA, srs_beta1 = NA, srs_beta2 = NA, srs_conv_msg = NA, srs_resampled_V = FALSE, 
                       cc_beta0 = NA, cc_beta1 = NA, cc_beta2 = NA, cc_conv_msg = NA, cc_resampled_V = FALSE, 
                       bcc_beta0 = NA, bcc_beta1 = NA, bcc_beta2 = NA, bcc_conv_msg = NA, bcc_resampled_V = FALSE, 
                       resid_beta0 = NA, resid_beta1 = NA, resid_beta2 = NA, resid_conv_msg = NA, resid_resampled_V = FALSE,
                       sfs_beta0 = NA, sfs_beta1 = NA, sfs_beta2 = NA, sfs_conv_msg = NA, sfs_resampled_V = FALSE,
                       optmle_beta0 = NA, optmle_beta1 = NA, optmle_beta2 = NA, optmle_conv_msg = NA, optmle_resampled_V = FALSE,
                       optmle2_beta0 = NA, optmle2_beta1 = NA, optmle2_beta2 = NA, optmle2_conv_msg = NA, optmle2_resampled_V = FALSE)
  
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
  for (des in c("srs", "cc", "bcc", "resid", "sfs", "optmle", "optmle2")) {
    #print(des)
    
    ## Create validation indicator 
    temp$V = get_V(design = toupper(des), 
                   data = temp, 
                   Xstar_cut = median(temp$Xstar1), 
                   X_cut = median(temp$X))
    
    ## Check for empty sieves in B-splines of validated data
    # while(0 %in% colSums(temp[which(temp$V == 1), paste0("bs", 1:nsieve)])) {
    #   ## Re-create validation indicator 
    #   temp$V = get_V(design = toupper(des), 
    #                  data = temp)
    #   
    #   ## Validation indicators re-sampled
    #   results[1, paste0(des, "_resampled_V")] = TRUE
    # }
    
    if (!any(is.na(temp$V))) {
      ## Create Xmiss to be NA if V = 0
      temp = temp |> 
        dplyr::mutate(Xmiss = dplyr::if_else(condition = V == 0, 
                                             true = NA, 
                                             false = Xstar2))
      
      ## Fit SMLE
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
    } else {
      results[1, paste0(des, "_conv_msg")] = "MLEs Didn't Converge"
    }
  }
  
  # Return results
  return(results)
}

# Error setting 1 (extra light): 99% TPR/1% FPR
set.seed(918) # Be reproducible
tpr99_fpr01 = do.call(what = rbind,
                      args = pbapply::pbsapply(X = 1:100,
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
                      args = pbapply::pbsapply(X = 1:100,
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
                      args = pbapply::pbsapply(X = 1:100,
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
                      args = pbapply::pbsapply(X = 1:100,
                                               FUN = sim_data_fit, 
                                               simplify = FALSE,
                                               tpr = 0.50,
                                               fpr = 0.50)) |> 
  dplyr::mutate(tpr = 0.5, 
                fpr = 0.5)
tpr50_fpr50 |> 
  write.csv("ALI_EHR/sim-data/vary_designs_tpr50_fpr50.csv", 
            row.names = FALSE)
