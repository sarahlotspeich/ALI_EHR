# Load data
## ALI components before validation 
unval_data = read.csv("~/Documents/Allostatic_load_audits/all_ali_dat.csv") |> 
  filter(DATA == "EHR (Before Validation)")

# Recenter age at 18 and rescale to be in 10-year increments and
## ALI/ALI* to be in 0.1-point increments. 
unval_data = unval_data |> 
  mutate(AGE_AT_ENCOUNTER_10 = AGE_AT_ENCOUNTER / 10,
         ALI_01 = ALI / 0.1)

# Fit naive model Y ~ X* + Z
naive_mod = glm(formula = ANY_ENCOUNTERS ~ ALI + AGE_AT_ENCOUNTER_10, 
                family = "binomial", 
                data = unval_data)
coefficients(summary(naive_mod))

# Transform log odds ratio for ALI from 1-point to 0.1-point scale
beta1 = naive_mod$coefficients[2] ## original 
beta1_t = beta1 * 0.1 ## transformed
se_beta1 = sqrt(diag(vcov(naive_mod)))[2] ## original
se_beta1_t = se_beta1 * 0.1 ## transformed

# Create vectors of log odds ratios and standard errors
logORs = c(naive_mod$coefficients[1], beta1_t, naive_mod$coefficients[3])
se_logORS = c(sqrt(diag(vcov(naive_mod)))[1], se_beta1_t, sqrt(diag(vcov(naive_mod)))[3])

# Estimates and 95% confidence intervals for the odds ratios
## Odds ratios
exp(logORs)
## Lower bounds 
exp(logORs - 1.96 * se_logORS)
## Upper bounds 
exp(logORs + 1.96 * se_logORS)