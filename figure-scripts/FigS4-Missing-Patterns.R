# Load packages 
library(dplyr) ## for data wrangling
library(naniar) ## for missing data plots

# Load data
## Error-Prone Data from the EHR
summ_data = read.csv("~/Documents/ALI_EHR/analysis/patient-data/summary_data.csv") |> 
  mutate(ANY_ENCOUNTERS = as.numeric(NUM_ENCOUNTERS > 0)) ### Create indicator of >= 1 healthcare encounter (ED visit or hospital admission)

## Subset to only the variables needed for the analysis 
summ_data = summ_data |> 
  select(PAT_MRN_ID, ANY_ENCOUNTERS, AGE_AT_ENCOUNTER, SEX, 
         CREAT_C, ALB, BMI, BP_SYSTOLIC, BP_DIASTOLIC, A1C, CHOL, TRIG, CRP, HCST)

## Plot missingness
summ_data |> 
  magrittr::set_colnames(c("PAT_MRN_ID", "ANY_ENCOUNTERS", "AGE", "SEX",
                           "CREATININE CLEARANCE",  "SERUM ALBUMIN",
                           "BODY MASS INDEX", "SYSTOLIC BLOOD PRESSURE",
                           "DIASTOLIC BLOOD PRESSURE", "HEMOGLOBIN A1C", 
                           "TOTAL CHOLESTEROL", "TRIGLYCERIDES",
                           "C-REACTIVE PROTEIN", "HOMOCYSTEINE")) |> 
  gg_miss_upset()

## Save it 
ggsave(filename = "~/Documents/ALI_EHR/figures/FigS4_Missing_Patterns.png", 
       device = "png", width = 8, height = 5, units = "in")
