# Load library --------------------------------------------
# library(sleev) ## for SMLE logistic regression
source("~/Documents/logreg2phRonly/R/logreg2ph.R")
source("~/Documents/logreg2phRonly/R/hessian_row.R")
source("~/Documents/logreg2phRonly/R/observed_data_loglik.R")
source("~/Documents/logreg2phRonly/R/pl_theta.R")
source("~/Documents/logreg2phRonly/R/profile_out.R")

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
pV = 0.1 ## proportion of patients to be validated
nsieve = 16 ## number of B-spline sieves
tpr = 0.95 ## true positive rate of EHR components
fpr = 0.05 ## false positive rate of EHR components

# Function to simulate data
sim_data = function(N = 1000, audit_recovery = 1) {
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
  X = rowMeans(S_mat, 
               na.rm = TRUE) ## ALI based on recovered components
  
  ## Select subset for validation
  audit_rows = sample(x = 1:N, 
                      size = ceiling(pV * N), 
                      replace = FALSE)
  
  ## Create dataset
  dat = data.frame(id = 1:N, X, Xstar, Y, Z)
  
  ## Create validation indicator 
  dat$V = as.numeric(dat$id %in% audit_rows)
  
  # Return dataset
  return(dat)
}

# Function to simulate data and then fit all models
sim_data_fit = function(id, N = 1000, audit_recovery = 1) {
  results = data.frame(sim = id, 
                       gs_beta0 = NA, gs_beta1 = NA, gs_beta2 = NA,
                       naive_beta0 = NA, naive_beta1 = NA, naive_beta2 = NA,
                       cc_beta0 = NA, cc_beta1 = NA, cc_beta2 = NA,
                       smle_beta0 = NA, smle_beta1 = NA, smle_beta2 = NA)
  
  # Simulate data 
  temp = sim_data(N = N, 
                  audit_recovery = audit_recovery) |> 
    dplyr::mutate(Xmiss = dplyr::if_else(condition = V == 0, 
                                         true = NA, 
                                         false = X))
  
  # 1. Gold standard
  fit = glm(formula = Y ~ X + Z, 
            family = "binomial", 
            data = temp)
  results[1, c("gs_beta0", "gs_beta1", "gs_beta2")] = coefficients(fit)
  
  # 2. Naive
  fit = glm(formula = Y ~ Xstar + Z, 
            family = "binomial", 
            data = temp)
  results[1, c("naive_beta0", "naive_beta1", "naive_beta2")] = coefficients(fit)
  
  # 3. Gold standard
  fit = glm(formula = Y ~ X + Z, 
            family = "binomial", 
            data = temp, 
            subset = V == 1)
  results[1, c("cc_beta0", "cc_beta1", "cc_beta2")] = coefficients(fit)
  
  # 4. SMLE
  ## Setup B-splines
  B = splines::bs(x = temp$Xstar, 
                  df = nsieve, 
                  Boundary.knots = range(temp$Xstar), 
                  intercept = TRUE)
  colnames(B) = paste0("bs", seq(1, nsieve))
  temp = cbind(temp, B)
  
  ## Fit SMLE
  suppressMessages(fit <- logreg2ph(
    Y_unval = NULL,
    Y_val = "Y",
    X_unval = "Xstar",
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
  
  results[1, c("smle_beta0", "smle_beta1", "smle_beta2")] = fit$model_coeff$coeff
  
  # Return results
  return(results)
}

# Setting 1 (full recovery): 100% of missing components found in EHR
set.seed(918) # Be reproducible
recover100 = do.call(what = rbind,
                     args = pbapply::pbsapply(X = 1:1000,
                                              FUN = sim_data_fit, 
                                              simplify = FALSE,
                                              N = 1000,
                                              audit_recovery = 1)) |> 
  dplyr::mutate(prop_recovered = 1)
recover100 |> 
  write.csv("ALI_EHR/sim-data/recover100.csv", 
            row.names = FALSE)

# Setting 2 (high recovery): 90% of missing components found in EHR
set.seed(918) # Be reproducible
recover90 = do.call(what = rbind,
                    args = pbapply::pbsapply(X = 1:1000,
                                             FUN = sim_data_fit, 
                                             simplify = FALSE,
                                             N = 1000,
                                             audit_recovery = 0.9)) |> 
  dplyr::mutate(prop_recovered = 0.9)
recover90 |> 
  write.csv("ALI_EHR/sim-data/recover90.csv", 
            row.names = FALSE)

# Setting 3 (moderate recovery): 50% of missing components found in EHR
set.seed(918) # Be reproducible
recover50 = do.call(what = rbind,
                    args = pbapply::pbsapply(X = 1:1000,
                                             FUN = sim_data_fit, 
                                             simplify = FALSE,
                                             N = 1000,
                                             audit_recovery = 0.5)) |> 
  dplyr::mutate(prop_recovered = 0.5)
recover50 |> 
  write.csv("ALI_EHR/sim-data/recover50.csv", 
            row.names = FALSE)

# Setting 4 (low recovery): 25% of missing components found in EHR
set.seed(918) # Be reproducible
recover25 = do.call(what = rbind,
                    args = pbapply::pbsapply(X = 1:1000,
                                             FUN = sim_data_fit, 
                                             simplify = FALSE,
                                             N = 1000,
                                             audit_recovery = 0.25)) |> 
  dplyr::mutate(prop_recovered = 0.25)
recover25 |> 
  write.csv("ALI_EHR/sim-data/recover25.csv", 
            row.names = FALSE)