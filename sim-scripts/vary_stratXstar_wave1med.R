# Load libraries ---------------------------------------------------------------
library(sleev) ## for SMLE logistic regression
library(auditDesignR) ## for audit sampling functions
library(magrittr) ## for pipes

# Random seed to be used for each simulation setting ---------------------------
args = commandArgs(TRUE) # 0 
## When running on the cluster, give each array a unique seed by adding the array ID to 11422
sim_seed = 11422 + as.integer(args)

# Set parameters that won't be varied in the loop
## These values will be set as the defaults in the sim_data() function for convenience
lambda_age = 45.66 / 5 ## mean of Poisson for Z
beta0 = -1.57 ## intercept in model of Y|X,Z
beta1 = 0.95 ## coefficient on X in model of Y|X,Z
beta2 = 0.1 ## coefficient on Z in model of Y|X,Z
pS = c(0.2500000, 0.9870130, 0.4549098, 0.1450000, 0.0580000, 
       0.2490119, 0.3138501, 0.3316391, 0.3111111, 0.0000000) ## probability of stressor = YES
pM = c(0.996, 0.153, 0.002, 0.000, 0.000, 
       0.494, 0.213, 0.213, 0.955, 0.983) ## probability of stressor = NA
N = 1000 ## total sample size (Phase I)
pV = 0.1 ## proportion of patients to be validated (Phase II)
n = ceiling(pV * N) ## number of patients to be validated (Phase II)
nsieve = 14 ## number of B-spline sieves
audit_recovery = 1 ## proportion of missing data recovered through the audit

# Function to simulate data
sim_data = function(tpr, fpr) {
  ## Simulate continuous error-free covariate: age at first encounter (in 5-year increments)
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
             prob = (1 + exp(-(beta0 + beta1 * X + beta2 * Z))) ^ (- 1))
  
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
get_V = function(data, Xstar_cut = 0.5, X_cut = 0.5, Z_cut = lambda_age / 2) {
  # Create stratified versions of Phase I variables  -------------------------
  ## For wave 1 (medians)
  data$Z_strat = as.numeric(data$Z <= Z_cut)
  data$X_strat1 = as.numeric(data$X <= median(data$X))
  data$Xstar_strat1 = as.numeric(data$Xstar1 <= median(data$Xstar1))
  ## For wave 2 (varied)
  data$Xstar_strat2 = as.numeric(data$Xstar1 <= quantile(x = data$Xstar1, probs = Xstar_cut))
  data$X_strat2 = as.numeric(data$X <= quantile(x = data$X, probs = X_cut))
  #  ------------------------- Create stratified versions of Phase I variables
  ############################################################################
  # Sample wave 1 audits using BCC* on Y and X* (median)  --------------------
  data$V = sample_bcc(dat = data, 
                      phI = N, 
                      phII = (n / 2), 
                      sample_on = c("Y", "Xstar_strat1"))
  #  -------------------- Sample wave 1 audits using BCC* on Y and X* (median)
  ############################################################################
  # Cross-tabulate strata sizes ----------------------------------------------
  ## Phase I strata based on X* (varied) -------------------------------------
  N00 = with(data, sum(Y == 0 & Xstar_strat2 == 0))
  N01 = with(data, sum(Y == 0 & Xstar_strat2 == 1))
  N10 = with(data, sum(Y == 1 & Xstar_strat2 == 0))
  N11 = with(data, sum(Y == 1 & Xstar_strat2 == 1))
  stratN = list(N00 = N00, N01 = N01, N10 = N10, N11 = N11)
  ## Cross-tabulate Phase IIa strata based on X* (varied) --------------------
  wave1_strat = data.frame(n00 = with(data, sum(V == 1 & Y == 0 & Xstar_strat2 == 0)),
                           n01 = with(data, sum(V == 1 & Y == 0 & Xstar_strat2 == 1)),
                           n10 = with(data, sum(V == 1 & Y == 1 & Xstar_strat2 == 0)),
                           n11 = with(data, sum(V == 1 & Y == 1 & Xstar_strat2 == 1)))
  # ---------------------------------------------- Cross-tabulate strata sizes
  ############################################################################
  # Estimate parameters using Phase IIa audits + the rest of Phase I ---------
  mle_wave1 = twophase_mle(dat = data, 
                           Y_val = "Y", 
                           Y_unval = NULL, 
                           X_val = "X_strat2", 
                           X_unval = "Xstar_strat2", 
                           addl_covar = "Z_strat",
                           Validated = "V", 
                           nondiff_X_unval = TRUE)
  if(mle_wave1$conv) {
    ## If MLEs converge, use them! -------------------------------------------
    beta_hat = mle_wave1$mod_Y_val$Est[3]
    eta_hat = with(mle_wave1, 
                   c(mod_Y_val$Est[1:2], 
                     mod_Y_unval$Est, 
                     mod_X_unval$Est, 
                     mod_X_val$Est))
    wave1est = "MLE"
    ## ------------------------------------------- If MLEs converge, use them!
  } else { 
    ## If MLEs don't converge, use complete case estimates -------------------
    mod_Y_val = glm(formula = Y ~ Z_strat + X_strat2, 
                    data = data, 
                    family = "binomial", 
                    subset = V == 1)
    mod_X_unval = glm(formula = Xstar_strat2 ~ X_strat2 + Z_strat, 
                      data = data, 
                      family = "binomial", 
                      subset = V == 1)
    mod_X_val = glm(formula =  X_strat2 ~ Z_strat, 
                    data = data, 
                    family = "binomial", 
                    subset = V == 1)
    beta_hat = mod_Y_val$coefficients[3]
    eta_hat = c(mod_Y_val$coefficients[1:2], 
                mod_X_unval$coefficients, 
                mod_X_val$coefficients)
    wave1est = "CC"
    ## ------------------- If MLEs don't converge, use complete case estimates
  }
  # --------- Estimate parameters using Phase IIa audits + the rest of Phase I
  ############################################################################
  # Using parameter estimates, find optimal wave 2 strata sizes --------------
  ## Create matrix of complete data  -----------------------------------------
  complete_data = expand.grid(Y = c(0, 1), 
                              X_strat2 = c(0, 1), 
                              Xstar_strat2 = c(0, 1),
                              Z_strat = c(0, 1),
                              V = c(0, 1))
  ##  ----------------------------------------- Create matrix of complete data
  ############################################################################
  ## Calculate score vectors for each row in complete_data -------------------
  s_hat = score(comp_dat = complete_data, 
                Y_val = "Y", 
                Y_unval = NULL, 
                X_val = "X_strat2", 
                X_unval = "Xstar_strat2", 
                addl_covar = "Z_strat",
                Validated = "V", 
                beta = beta_hat, 
                eta = eta_hat, 
                nondiff_X_unval = TRUE)
  ## ------------------- Calculate score vectors for each row in complete_data
  ############################################################################
  ## Run grid search to find strata sizes to minimize expected variance ------
  grid_search = optMLE_grid(phI = N, 
                            phII = (n / 2), 
                            phI_strat = stratN, 
                            phIIa_strat = wave1_strat, 
                            min_n = 0, 
                            sample_on = c("Y", "Xstar_strat2"), 
                            indiv_score = s_hat)
  ## ------ Run grid search to find strata sizes to minimize expected variance
  ############################################################################
  ## If grid search found a minimum, sample according to the optimal design --
  if (grid_search$findOptimal) {
    opt_des2 = grid_search$min_var_design ### the optimal design (to be sampled in wave 1 + wave 2)
    opt_des2[, c("n00", "n01", "n10", "n11")] = opt_des2[, c("n00", "n01", "n10", "n11")] - 
      with(wave1_strat, c(n00, n01, n10, n11)) ### the optimal design (just to be sampled in wave 2)
    V_des = list(
      V = pmax(data$V, 
               sample_optMLE(dat = data, 
                             sample_on = c("Y", "Xstar_strat2"), 
                             des = opt_des2, 
                             wave1_Validated = "V")), ### validation indicator (wave 1 or 2)
      opt_des = grid_search$min_var_design,
      opt_des_msg = grid_search$message,
      wave1params = wave1est
    )
  } else {
    V_des = list(
      V = rep(0, nrow(data)), ### validation indicator (wave 1 or 2)
      opt_des = grid_search$min_var_design,
      opt_des_msg = grid_search$message,
      wave1params = wave1est
    )
  }
  return(V_des)
}

# Function to simulate data and then fit all models
sim_data_fit = function(id, tpr, fpr) {
  results = data.frame(sim = id, 
                       q = c(0.1, 0.25, 0.5, 0.75, 0.9),
                       beta0 = NA, beta1 = NA, beta2 = NA, 
                       conv_msg = NA, resampled_V = FALSE, 
                       des_n00 = NA, des_n01 = NA, des_n10 = NA, des_n11 = NA, 
                       wave1est = NA)
  
  # Simulate data 
  data = sim_data(tpr = tpr, 
                  fpr = fpr) 
  
  # Setup B-splines
  B = splines::bs(x = data$Xstar1, 
                  df = nsieve, 
                  Boundary.knots = range(data$Xstar1), 
                  intercept = TRUE, 
                  degree = 3)
  colnames(B) = paste0("bs", seq(1, nsieve))
  data = cbind(data, B)
  
  # Loop over the different designs
  for (r in 1:nrow(results)) {
    V_des = get_V(data = data, 
                  Xstar_cut = quantile(x = data$Xstar1, probs = results$q[r]), 
                  X_cut = quantile(x = data$X, probs = results$q[r]))
    ## If optimal design was located, sample and estimate
    if (sum(V_des$V) != 0) {
      ### Append the validation indicators to the data (wave 1 or 2) -----------
      data$V = V_des$V 
      
      ### Create Xmiss to be NA if V = 0 ---------------------------------------
      data = data %>%
        dplyr::mutate(Xmiss = dplyr::if_else(condition = V == 0, 
                                             true = NA, 
                                             false = Xstar2))
      
      ### Fit SMLE -------------------------------------------------------------
      suppressMessages(fit <- logistic2ph(
        Y_unval = NULL,
        Y = "Y",
        X_unval = "Xstar1",
        X = "Xmiss",
        Z = "Z",
        Bspline = colnames(B),
        data = data,
        hn_scale = 1,
        noSE = FALSE,
        TOL = 1e-04,
        MAX_ITER = 1000
      ))
      
      # while (fit$converge_msg == "B-spline error") {
      #   # V_des = get_V(data = data, 
      #   #               Xstar_cut = quantile(x = data$Xstar1, probs = results$q[r]), 
      #   #               X_cut = quantile(x = data$X, probs = results$q[r]))
      #   data$V = V_des$V ## Append the validation indicators to the data
      #   ## Create Xmiss to be NA if V = 0
      #   data = data %>%
      #     dplyr::mutate(Xmiss = dplyr::if_else(condition = V == 0, 
      #                                          true = NA, 
      #                                          false = Xstar2))
      #   
      #   ## Fit SMLE
      #   suppressMessages(fit <- logistic2ph(
      #     Y_unval = NULL,
      #     Y = "Y",
      #     X_unval = "Xstar1",
      #     X = "Xmiss",
      #     Z = "Z",
      #     Bspline = colnames(B),
      #     data = data,
      #     hn_scale = 1,
      #     noSE = FALSE,
      #     TOL = 1e-04,
      #     MAX_ITER = 1000
      #   ))
      #   results[r, "resampled_V"] = TRUE
      # }
      results[r, "wave1est"] = V_des$wave1params
      results[r, paste0("beta", 0:2)] = fit$coefficients$Estimate
      results[r, "conv_msg"] = fit$converge_msg
      results[r, paste0("des_", c("n00", "n01", "n10", "n11"))] = V_des$opt_des[, c("n00", "n01", "n10", "n11")]
    } else {
      results[r, "conv_msg"] = V_des$opt_des_msg
    }
  }
  
  # Return results
  return(results)
}

# Error setting 1 (extra light): 99% TPR/1% FPR
set.seed(sim_seed) # Be reproducible
for (x in 1:100) {
  sim_data_fit(id = x, tpr = 0.99, fpr = 0.01)
  print(x)
}
tpr99_fpr01 = do.call(what = rbind,
                      args = sapply(X = 1:100,
                                    FUN = sim_data_fit, 
                                    simplify = FALSE,
                                    tpr = 0.99,
                                    fpr = 0.01)) %>% 
  dplyr::mutate(tpr = 0.99, 
                fpr = 0.01)
tpr99_fpr01 %>%
  write.csv(paste0("vary_stratXstar_wave1med/vary_stratXstar_tpr99_fpr01_seed", sim_seed, ".csv"), 
            row.names = FALSE)

# Error setting 2 (light): 95% TPR/5% FPR
set.seed(sim_seed) # Be reproducible
tpr95_fpr05 = do.call(what = rbind,
                      args = sapply(X = 1:100,
                                    FUN = sim_data_fit, 
                                    simplify = FALSE,
                                    tpr = 0.95,
                                    fpr = 0.05)) %>% 
  dplyr::mutate(tpr = 0.95, 
                fpr = 0.05)
tpr95_fpr05 %>%
  write.csv(paste0("vary_stratXstar_wave1med/vary_stratXstar_tpr95_fpr05_seed", sim_seed, ".csv"), 
            row.names = FALSE)

# Error setting 3 (moderate): 80% TPR/20% FPR
set.seed(sim_seed) # Be reproducible
tpr80_fpr20 = do.call(what = rbind,
                      args = sapply(X = 1:100,
                                    FUN = sim_data_fit, 
                                    simplify = FALSE,
                                    tpr = 0.80,
                                    fpr = 0.20)) %>%
  dplyr::mutate(tpr = 0.8, 
                fpr = 0.2)
tpr80_fpr20 %>%
  write.csv(paste0("vary_stratXstar_wave1med/vary_stratXstar_tpr80_fpr20_seed", sim_seed, ".csv"),
            row.names = FALSE)

# Error setting 4 (heavy): 50% TPR/50% FPR
set.seed(sim_seed) # Be reproducible
tpr50_fpr50 = do.call(what = rbind,
                      args = sapply(X = 1:100,
                                    FUN = sim_data_fit, 
                                    simplify = FALSE,
                                    tpr = 0.50,
                                    fpr = 0.50)) %>% 
  dplyr::mutate(tpr = 0.5, 
                fpr = 0.5)
tpr50_fpr50 %>%
  write.csv(paste0("vary_stratXstar_wave1med/vary_stratXstar_tpr50_fpr50_seed", sim_seed, ".csv"),
            row.names = FALSE)

