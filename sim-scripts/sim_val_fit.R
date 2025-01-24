# Load packages ----------------------------------------------------------------
### RUN ONCE: devtools::install_github("sarahlotspeich/auditDesignR", ref = "main")
library(auditDesignR)

# Create matrix of complete data for the optimal design ------------------------
complete_data = expand.grid(Y = c(0, 1), ## binary outcome
                            Xval_strat = c(0, 1), ## stratified ALI (validated)
                            Xstar_strat = c(0, 1), ## stratified ALI (unvalidated)
                            V1 = c(0, 1)) ## validation indicators

# Function to sample Wave 2 of validation --------------------------------------
wave2_val = function(data, val_design) {
  if (val_design == "SRS") {
    data$V2 = sample_srs(phI = nrow(data), 
                         phII = 48, 
                         wave1_Validated = data$V1 == 1)  
  } else if (val_design == "CC") {
    data$V2 = sample_cc(dat = data,
                        phI = nrow(data), 
                        phII = 48, 
                        sample_on = "Y", 
                        wave1_Validated = data$V1 == 1)
  } else if (val_design == "BCC*") {
    data$V2 = sample_bcc(dat = data,
                         phI = nrow(data), 
                         phII = 48, 
                         sample_on = c("Y", "Xstar_strat"), 
                         wave1_Validated = data$V1 == 1)
  } else if (val_design == "RS") {
    data$V2 = sample_resid(formula = Y ~ Xstar + Z, 
                           family = "binomial", 
                           dat = data, 
                           phI = nrow(data), 
                           phII = 48, 
                           wave1_Validated = data$V1 == 1)
  } else if (val_design == "optMLE") {
    ## Create stratified versions of additional Phase I variables  -------------
    data$Z_strat = as.numeric(data$Z <= median(data$Z))
    data$Xval_strat = as.numeric(data$Xval <= median(data$Xval))

    # Create matrix of complete data  ------------------------------------------
    complete_data = expand.grid(Y = c(0, 1), 
                                Xval_strat = c(0, 1), 
                                Xstar_strat = c(0, 1),
                                Z_strat = c(0, 1),
                                V = c(0, 1))
    
    # Cross-tabulate Phase I strata --------------------------------------------
    N00 = with(data, sum(Y == 0 & Xstar_strat == 0))
    N01 = with(data, sum(Y == 0 & Xstar_strat == 1))
    N10 = with(data, sum(Y == 1 & Xstar_strat == 0))
    N11 = with(data, sum(Y == 1 & Xstar_strat == 1))
    stratN = list(N00 = N00, 
                  N01 = N01, 
                  N10 = N10, 
                  N11 = N11)
    
    # Sample wave 1 audits and estimate parameters and score  ------------------
    wave1_strat = data.frame(n00 = with(data[data$V1 == 1, ], sum(Y == 0 & Xstar_strat == 0)),
                             n01 = with(data[data$V1 == 1, ], sum(Y == 0 & Xstar_strat == 1)),
                             n10 = with(data[data$V1 == 1, ], sum(Y == 1 & Xstar_strat == 0)),
                             n11 = with(data[data$V1 == 1, ], sum(Y == 1 & Xstar_strat == 1)))
    mle_wave1 = twophase_mle(dat = data, 
                             Y_val = "Y", 
                             Y_unval = NULL, 
                             X_val = "Xval_strat", 
                             X_unval = "Xstar_strat", 
                             addl_covar = "Z_strat",
                             Validated = "V1", 
                             nondiff_X_unval = TRUE)
    
    if (mle_wave1$conv) {
      beta_hat = mle_wave1$mod_Y_val$Est[3]
      eta_hat = with(mle_wave1, 
                     c(mod_Y_val$Est[1:2], 
                       mod_Y_unval$Est, 
                       mod_X_unval$Est, 
                       mod_X_val$Est))
    } else {
      ## If MLEs don't converge, use complete case estimates -------------------
      mod_Y_val = glm(formula = Y ~ Z_strat + Xval_strat, 
                      data = data, 
                      family = "binomial", 
                      subset = V1 == 1)
      mod_X_unval = glm(formula = Xstar_strat ~ Xval_strat + Z_strat, 
                        data = data, 
                        family = "binomial", 
                        subset = V1 == 1)
      mod_X_val = glm(formula =  Xval_strat ~ Z_strat, 
                      data = data, 
                      family = "binomial", 
                      subset = V1 == 1)
      beta_hat = mod_Y_val$coefficients[3]
      eta_hat = c(mod_Y_val$coefficients[1:2], 
                  mod_X_unval$coefficients, 
                  mod_X_val$coefficients)
      print("Used Wave I complete-case estimates instead.")
    }
    
    s_hat = score(comp_dat = complete_data, 
                  Y_val = "Y", 
                  Y_unval = NULL, 
                  X_val = "Xval_strat", 
                  X_unval = "Xstar_strat", 
                  addl_covar = "Z_strat",
                  Validated = "V", 
                  beta = beta_hat, 
                  eta = eta_hat, 
                  nondiff_X_unval = TRUE)
    
    grid_search = optMLE_grid(phI = nrow(data), 
                              phII = 48, 
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
      data$V2 = sample_optMLE(dat = data, 
                              sample_on = c("Y", "Xstar_strat"), 
                              des = opt_des2, 
                              wave1_Validated = "V1")
    } else {
      data$V2 = NA
    }
  } 
  
  return(data)
}

# Function to simulate data and then fit all models ----------------------------
sim_val_fit = function(id, tpr = 0.95, fpr = 0.05, audit_recovery = 1, val_design) {
  results = data.frame(sim = id, val_design, 
                       gs_beta0 = NA, gs_beta1 = NA, gs_beta2 = NA,
                       naive_beta0 = NA, naive_beta1 = NA, naive_beta2 = NA,
                       cc_beta0 = NA, cc_beta1 = NA, cc_beta2 = NA,
                       smle_beta0 = NA, smle_beta1 = NA, smle_beta2 = NA, 
                       smle_conv_msg = NA, empty_sieve = FALSE, nsieve = NA)
  
  # Simulate data 
  temp = sim_ali_data(tpr = tpr, 
                      fpr = fpr, 
                      audit_recovery = audit_recovery)
  
  # Setup B-splines
  B = bs(x = temp$Xstar, ## Error-prone ALI (from EHR)
         df = nsieve, 
         Boundary.knots = range(temp$Xstar), 
         intercept = TRUE)
  colnames(B) = paste0("bs", seq(1, nsieve))
  temp = cbind(temp, B)

  # Wave 1 validation: Draw a BCC* of 52 from (Y, X*) strata
  ## Stratify X* at the median
  temp$Xstar_strat = as.numeric(temp$Xstar <= median(temp$Xstar))
  
  ## Create Wave 1 validation indicator
  temp$V1 = sample_bcc(dat = temp,
                       phI = nrow(temp), 
                       phII = 52, 
                       sample_on = c("Y", "Xstar_strat"), 
                       wave1_Validated = NULL)
  
  # Wave 2 validation: Draw remaining 48 from specified design
  ## Create Wave 2 validation indicator 
  temp = temp |> 
    wave2_val(val_design = val_design)

  # Check that the Wave2 sampling was successful
  if(!any(is.na(temp$V2))) {
    # Final validation: Create indicator for Wave 1 *or* Wave 2
    temp$V = pmax(temp$V1, temp$V2)
    
    # Check for empty sieves in B-splines of validated data
    results[1, "empty_sieve"] = 0 %in% colSums(temp[which(temp$V == 1), colnames(B)])
    ## If present, decrease number of B-splines one at a time 
    temp_nsieve = nsieve ### Initialize at current value
    while(0 %in% colSums(temp[which(temp$V == 1), colnames(B)]) & temp_nsieve >= 5) {
      ### Decrease number of B-splines by 1
      temp_nsieve = temp_nsieve - 1
      
      ### Delete existing B-splines
      temp[, grep(pattern = "bs", x = colnames(B), value = TRUE)] = NULL
      
      ### Setup new B-splines
      B = bs(x = temp$Xstar, ## Error-prone ALI (from EHR)
             df = temp_nsieve, 
             Boundary.knots = range(temp$Xstar), 
             intercept = TRUE)
      colnames(B) = paste0("bs", seq(1, temp_nsieve))
      temp = cbind(temp, B)
    }
    
    if (!0 %in% colSums(temp[which(temp$V == 1), colnames(B)])) {
      # Save the number of B-splines used 
      results[1, "nsieve"] = temp_nsieve
      
      # Create Xmiss to be NA if V = 0 (for unvalidated patients)
      temp = temp |> 
        mutate(Xmiss = if_else(condition = V == 0, 
                               true = NA, 
                               false = Xval))
      
      # 1. Gold standard analysis
      fit = glm(formula = Y ~ X + Z, 
                family = "binomial", 
                data = temp)
      results[1, c("gs_beta0", "gs_beta1", "gs_beta2")] = coefficients(fit)
      
      # 2. Naive analysis
      fit = glm(formula = Y ~ Xstar + Z, 
                family = "binomial", 
                data = temp)
      results[1, c("naive_beta0", "naive_beta1", "naive_beta2")] = coefficients(fit)
      
      # 3. Complete case analysis
      fit = glm(formula = Y ~ Xmiss + Z, 
                family = "binomial", 
                data = temp, 
                subset = V == 1)
      results[1, c("cc_beta0", "cc_beta1", "cc_beta2")] = coefficients(fit)
      
      # 4. SMLE analysis
      fit <- logistic2ph(
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
      results[1, c("smle_beta0", "smle_beta1", "smle_beta2")] = fit$coefficients$Estimate
      results[1, "smle_conv_msg"] = fit$converge
    }
  } 
  
  # Return results
  return(results)
}