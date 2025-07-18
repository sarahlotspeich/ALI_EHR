# Load packages 
library(dplyr) ## for data wrangling
library(tidyr) ## for other data wrangling
library(kableExtra) ## for table making

# Load data (wave I + pilot) ---------------------------------------------------
## Read in long audit data (one row per patient per audited variable)
long_audit_dat1 = dat = read.csv("~/Documents/Allostatic_load_audits/Pilot_First4/pilot_audits_for_analysis.csv") |> 
  select(-contains("Rabeya"), -contains("Aidan"), -contains("Amelia")) |> 
  rename(Reviewed_Value = Reviewed_Value_Sheetal, 
         Reviewed_Value_Num = Reviewed_Value_Sheetal_Num, 
         Notes = Notes_Sheetal, 
         Finding = Finding_Sheetal) |> 
  bind_rows(
    read.csv("~/Documents/Allostatic_load_audits/Wave1_Next44_StartingMarch7/wave1_audits_for_analysis.csv")
  ) |> 
  filter(!(Variable_Name %in% c("HEIGHT", "WEIGHT")))

# Load data ( wave II) -----------------------------------------
## Wave II audits 
long_audit_dat2 = read.csv("~/Documents/Allostatic_load_audits/Wave2_Last48_StartingJune/wave2_audits_for_analysis.csv") |> 
  filter(!(Variable_Name %in% c("HEIGHT", "WEIGHT")))

# Combine them -------------------------------------------------
long_audit_dat = long_audit_dat1 |> 
  bind_rows(long_audit_dat2) |> 
  mutate(Variable_Name = factor(x = Variable_Name, 
                                levels = c("BP_SYSTOLIC", 
                                           "BP_DIASTOLIC", 
                                           "BMI", 
                                           "ALB", 
                                           "CHOL", 
                                           "TRIG", 
                                           "A1C", 
                                           "CREAT_C", 
                                           "CRP", 
                                           "HCST"), 
                                labels = c("Systolic Blood Pressure", 
                                           "Diastolic Blood Pressure", 
                                           "Body Mass Index", 
                                           "Serum Albumin", 
                                           "Cholesterol", 
                                           "Triglycerides", 
                                           "Hemoglobin A1C", 
                                           "Creatinine Clearance", 
                                           "C-Reactive Protein", 
                                           "Homocysteine"
                                )), 
         Finding = factor(x = Finding, 
                          levels = c("Extracted Value Correct", 
                                     "Extracted Value Incorrect", 
                                     "Extracted Value Not Found", 
                                     "No Auxiliary Information Found",
                                     "Auxiliary Information Found" 
                          )))

# Calculate percent of each finding ----------------------------
summ_audit_dat = long_audit_dat |> 
  group_by(Category, Variable_Name, Finding) |> 
  summarize(Num = n()) |> 
  group_by(Category, Variable_Name) |> 
  mutate(Total_Num = sum(Num), 
         Perc = Num / Total_Num, 
         Num_Perc = paste0(Num, " (", 
                           ifelse(test = round(100 * Perc) >= 1, 
                                  yes = round(100 * Perc), 
                                  no = "$<1$"),
                           "\\%)")) |> 
  select(Category, Variable_Name, Total_Num, Finding, Num_Perc) |> 
  spread(key = Finding, value = Num_Perc, fill = "0 (0\\%)")
summ_audit_dat |> 
  ungroup() |> 
  select(-Category) |> 
  kable(format = "latex", digits = 0, booktabs = T, escape = FALSE)
### Collapsed by Category 
long_audit_dat |> 
  group_by(Category, Finding) |> 
  summarize(Num = n()) |> 
  group_by(Category) |> 
  mutate(Total_Num = sum(Num), 
         Perc = Num / Total_Num, 
         Num_Perc = paste0(Num, " (", 
                           ifelse(test = round(100 * Perc) >= 1, 
                                  yes = round(100 * Perc), 
                                  no = "$<1$"),
                           "\\%)")) |> 
  select(Category, Total_Num, Finding, Num_Perc) |> 
  spread(key = Finding, value = Num_Perc, fill = "0 (0\\%)") |> 
  kable(format = "latex", digits = 0, booktabs = T, escape = FALSE) |> 
  row_spec(row = 1:2, bold = TRUE)