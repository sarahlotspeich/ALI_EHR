# Load packages
library(dplyr) ## for data wrangling
### RUN ONCE: devtools::install_github("sarahlotspeich/logiSieve", ref = "main")
library(logiSieve) ## for the SMLEs

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

# Recenter age at 18 and rescale to be in 10-year increments
data = data |> 
  mutate(AGE_AT_ENCOUNTER_10 = (AGE_AT_ENCOUNTER - 18) / 10)

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

## Check for empty B-splines among validated rows 
colSums(B[which(!is.na(data$ALI)), ]) ### None! 

## Fit SMLE model Y ~ X + Z
fit = logiSieve(
  analysis_formula = ANY_ENCOUNTERS ~ ALI + AGE_AT_ENCOUNTER_10, 
  error_formula = paste("ALI ~", paste(colnames(B), collapse = "+")), 
  data = data, 
  pert_scale = 1 / 4) #### needed a smaller perturbation since \beta_3 is close to 0 

# Transform log odds ratio for ALI from 1-point to 0.1-point scale
beta1 = fit$model_coeff$coeff[2] ## original 
beta1_t = beta1 * 0.1 ## transformed
se_beta1 = fit$model_coeff$se[2] ## original
se_beta1_t = se_beta1 * 0.1 ## transformed

# Create vectors of log odds ratios and standard errors
logORs = c(fit$model_coeff$coeff[1], beta1_t, fit$model_coeff$coeff[3])
se_logORS = c(fit$model_coeff$se[1], se_beta1_t, fit$model_coeff$se[3])

# Estimates and 95% confidence intervals for the odds ratios
exp(logORs) ## Odds ratios: 0.173, 1.207, 1.104
exp(logORs - 1.96 * se_logORS) ## Lower bounds: 0.103, 1.104 1.007 
exp(logORs + 1.96 * se_logORS) ## Upper bounds: 0.292, 1.319, 1.210