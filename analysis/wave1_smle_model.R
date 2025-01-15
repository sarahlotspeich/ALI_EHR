# Load packages
library(dplyr) ## for data wrangling
### RUN ONCE: devtools::install_github("dragontaoran/sleev", ref = "main")
library(sleev) ## for the SMLE

# Load data
## ALI components after Pilot + Wave I validation 
val_data = read.csv("~/Documents/Allostatic_load_audits/all_ali_dat.csv") |> 
  filter(DATA == "Pilot + Wave I Validation", 
         VALIDATED)
nrow(val_data) ## CHECK: n = 52 validated 
## ALI before validation and age / hospitalizations
unval_data = read.csv("~/Documents/Allostatic_load_audits/all_ali_dat.csv") |> 
  filter(DATA == "EHR (Before Validation)") |> 
  rename(ALI_STAR = ALI)
nrow(unval_data) ## CHECK: N = 1000 unvalidated 
## Merge validated + unvalidated 
data = val_data |> 
  select(PAT_MRN_ID, ALI) |> 
  right_join(
    unval_data |> 
      select(PAT_MRN_ID, ANY_ENCOUNTERS, AGE_AT_ENCOUNTER, ALI_STAR)
  )

# Recenter age at 18 and rescale to be in 10-year increments and
## ALI/ALI* to be in 0.1-point increments. 
data = data |> 
  mutate(AGE_AT_ENCOUNTER_10 = (AGE_AT_ENCOUNTER - 18) / 10,
         ALI_01 = ALI / 0.1, 
         ALI_STAR_01 = ALI_STAR / 0.1)

# Estimate parameters using Phase IIa audits + the rest of Phase I -------------
## Setup B-splines
B = splines::bs(x = data$ALI_STAR, 
                df = 8, 
                Boundary.knots = range(data$ALI_STAR), 
                intercept = TRUE, 
                degree = 3)
colnames(B) = paste0("bs", seq(1, 8))
data = data |> 
  bind_cols(B)

### Fit SMLE model Y ~ X + Z----------------------------------------------------
suppressMessages(fit <- logistic2ph(
  Y_unval = NULL,
  Y = "ANY_ENCOUNTERS",
  X_unval = "ALI_STAR",
  X = "ALI",
  Z = "AGE_AT_ENCOUNTER_10",
  Bspline = colnames(B),
  data = data,
  hn_scale = 1,
  noSE = FALSE,
  TOL = 1e-04,
  MAX_ITER = 1000
))

# Transform log odds ratio for ALI from 1-point to 0.1-point scale
beta1 = fit$coefficients$Estimate[2] ## original 
beta1_t = beta1 * 0.1 ## transformed
se_beta1 = fit$coefficients$SE[2] ## original
se_beta1_t = se_beta1 * 0.1 ## transformed

# Create vectors of log odds ratios and standard errors
logORs = c(fit$coefficients$Estimate[1], beta1_t, fit$coefficients$Estimate[3])
se_logORS = c(fit$coefficients$SE[1], se_beta1_t, fit$coefficients$SE[3])

# Estimates and 95% confidence intervals for the odds ratios
## Odds ratios
exp(logORs)
## Lower bounds 
exp(logORs - 1.96 * se_logORS)
## Upper bounds 
exp(logORs + 1.96 * se_logORS)
