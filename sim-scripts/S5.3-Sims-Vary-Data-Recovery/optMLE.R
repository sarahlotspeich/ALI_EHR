################################################################################
# Simulations for Supplemental Section S.5.2 -----------------------------------
## Held fixed: Sample size, validation proportion, error rates -----------------
## Varied: Data recovery rate and validation study design for chart reviews ----
################################################################################

# Load packages ----------------------------------------------------------------
### RUN ONCE: devtools::install_github("dragontaoran/sleev", ref = "main")
library(splines) ## for B-splines 
library(sleev) ## for SMLE logistic regression
library(dplyr) ## for data wrangling
library(pbapply) ## for apply functions with progress bar

# Load data generating function sim_ali_data() from GitHub ---------------------
library(devtools) # To source an R script from GitHub
source_url("https://raw.githubusercontent.com/sarahlotspeich/ALI_EHR/refs/heads/main/sim-scripts/sim_ali_data.R")

# Load data generating function sim_ali_data() from GitHub ---------------------
source_url("https://raw.githubusercontent.com/sarahlotspeich/ALI_EHR/refs/heads/main/sim-scripts/sim_val_fit.R")

# Set parameters that won't be varied in the loop ------------------------------
val_design = "optMLE" ## validation study design
N = 1000 ## total sample size
pV = 0.1 ## proportion of patients to be validated
nsieve = 16 ## number of B-spline sieves
nsim = 1000 ## number of replications per setting

# Consider multiple settings with 0 - 100% missing components recovered --------
all_recover = data.frame() ## Initialize empty data frame to hold results 
for (audit_recovery in c(1, 0.9, 0.5, 0.25, 0)) {
  set.seed(918) ## Be reproducible
  all_recover = do.call(what = rbind,
                 args = pbsapply(X = 1:nsim,
                                 FUN = sim_val_fit, 
                                 simplify = FALSE, 
                                 audit_recovery = audit_recovery, 
                                 val_design = val_design)) |> 
    mutate(prop_recovered = audit_recovery) |> 
    bind_rows(all_recover)
  all_recover |> 
    write.csv(paste0(val_design, ".csv"), 
              row.names = FALSE)
}