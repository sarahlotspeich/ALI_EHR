# //////////////////////////////////////////////////////////////////////
# Replicate Figure S1 in Supplementary Materials  //////////////////////
# Caption begins "Estimated log prevalence ratios for food access $X_P$ 
# on health. The five possible ways to include the analysis model... ///
# //////////////////////////////////////////////////////////////////////

# Load packages
library(dplyr) # To wrangle data
library(tidyr) # To transform data
library(ggplot2) # To create plots
library(latex2exp) # To create LaTex labels for plots

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

# Function to simulate data
sim_data = function(N = 1000, tpr = 0.95, fpr = 0.05, audit_recovery = 1) {
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
  
  ## Select subset for validation
  audit_rows = sample(x = 1:N, 
                      size = ceiling(pV * N), 
                      replace = FALSE)
  
  ## Create dataset
  dat = data.frame(id = 1:N, X, Xstar1, Xstar2, Y, Z)
  
  ## Create validation indicator 
  dat$V = as.numeric(dat$id %in% audit_rows)
  
  # Return dataset
  return(dat)
}

save_dat = data.frame()
for(rec in c(0.25, 0.5, 0.9, 1)) {
  save_dat = sim_data(N = 10000, 
                      audit_recovery = rec) |> 
    dplyr::mutate(perc_recovered = rec) |> 
    dplyr::bind_rows(save_dat)
}

save_dat |>
  dplyr::mutate(row_id = 1:dplyr::n(),
                perc_recovered = factor(x = perc_recovered, 
                                        levels = c(0.25, 0.5, 0.9, 1),
                                        labels = paste0(c(25, 50, 90, 100), "% of Missing Data Recovered"))) |> 
  ggplot(aes(x = X, y = Xstar2)) + 
  geom_point(alpha = 0.5) + 
  geom_smooth(method = "lm", se = FALSE) + 
  theme_minimal() + 
  facet_wrap(~ perc_recovered) + 
  xlab("True Allostatic Load Index") + 
  ylab("Validated Allostatic Load Index (After Chart Review)") + 
  theme(strip.background = element_rect(fill = "black"), 
        strip.text = element_text(colour = "white"))
ggsave(filename = "~/Documents/ALI_EHR/figures/sim_settings_audit_recovery.png", 
       device = "png", width = 7, height = 5, units = "in")

save_dat = data.frame()
error_sett = data.frame(tpr = c(0.99, 0.95, 0.8, 0.5), 
                        fpr = c(0.01, 0.05, 0.2, 0.5))
for(s in 1:nrow(error_sett)) {
  save_dat = sim_data(N = 10000, 
                      tpr = error_sett$tpr[s], 
                      fpr = error_sett$fpr[s]) |> 
    dplyr::mutate(tpr_fpr = paste0("TPR = ", 
                                   round(100 * error_sett$tpr[s]), 
                                   "%, FPR = ", 
                                   round(100 * error_sett$fpr[s]), "%")) |> 
    dplyr::bind_rows(save_dat)
}

save_dat |>
  dplyr::mutate(row_id = 1:dplyr::n(),
                tpr_fpr = factor(x = tpr_fpr, 
                                 levels = c("TPR = 99%, FPR = 1%",
                                            "TPR = 95%, FPR = 5%",
                                            "TPR = 80%, FPR = 20%",
                                            "TPR = 50%, FPR = 50%"))) |> 
  ggplot(aes(x = X, y = Xstar1)) + 
  geom_point(alpha = 0.5) + 
  geom_smooth(method = "lm", se = FALSE) + 
  theme_minimal() + 
  facet_wrap(~tpr_fpr) + 
  xlab("True Allostatic Load Index") + 
  ylab("Validated Allostatic Load Index (After Chart Review)") + 
  theme(strip.background = element_rect(fill = "black"), 
        strip.text = element_text(colour = "white"))
ggsave(filename = "~/Documents/ALI_EHR/figures/sim_settings_vary_tpr_fpr.png", 
       device = "png", width = 7, height = 5, units = "in")
