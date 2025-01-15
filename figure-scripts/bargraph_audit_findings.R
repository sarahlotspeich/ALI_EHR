# Load packages 
library(dplyr) ## for data wrangling
library(ggplot2) ## for pretty plots

# Define colors 
cols = c("#ff99ff", "#8bdddb", "#787ff6", "#ffbd59", "#7dd5f6")

# Load data (wave I + pilot) ---------------------------------------------------
## Read in long audit data (one row per patient per audited variable)
long_audit_dat = dat = read.csv("~/Documents/Allostatic_load_audits/Pilot_First4/pilot_audits_for_analysis.csv") |> 
  select(-contains("Rabeya"), -contains("Aidan"), -contains("Amelia")) |> 
  rename(Reviewed_Value = Reviewed_Value_Sheetal, 
         Reviewed_Value_Num = Reviewed_Value_Sheetal_Num, 
         Notes = Notes_Sheetal, 
         Finding = Finding_Sheetal) |> 
  bind_rows(
    read.csv("~/Documents/Allostatic_load_audits/Wave1_Next44_StartingMarch7/wave1_audits_for_analysis.csv")
  ) |> 
  filter(!(Variable_Name %in% c("HEIGHT", "WEIGHT")))
nrow(long_audit_dat) ## 7605 audited data points 
table(long_audit_dat$Category) ## 1049 labs, 6556 vitals
long_audit_dat |> 
  group_by(PAT_MRN_ID) |> 
  summarize(NUM_AUDITED = dplyr::n()) |> 
  pull(NUM_AUDITED) |> 
  summary() ## min = 13, median = 65.5, max = 763 data points per patient

## ALI components after wave 1 validation
order_levels = long_audit_dat |> 
  group_by(Variable_Name) |> 
  summarize(num = n()) |> 
  arrange(desc(num)) |> 
  pull(Variable_Name)
long_audit_dat |> 
  mutate(
    COMP = factor(x = Variable_Name, 
                  levels = order_levels, 
                  labels = c("SYSTOLIC BLOOD PRESSURE", 
                             "DIASTOLIC BLOOD PRESSURE", 
                             "BODY MASS INDEX", 
                             "SERUM ALBUMIN", 
                             "CHOLEST-\nEROL", 
                             "TRIGLY-\nCERIDES", 
                             "HEMOGLOBIN A1C", 
                             "CREATININE CLEARANCE", 
                             "C-REACTIVE PROTEIN", 
                             "HOMO-\nCYSTEINE"
                             )), 
    Finding = factor(x = Finding, 
                     levels = c("Extracted Value Correct", 
                                "Extracted Value Incorrect", 
                                "Extracted Value Not Found", 
                                "No Auxiliary Information Found",
                                "Auxiliary Information Found" 
                                ))
  ) |> 
  ggplot(aes(x = COMP, fill = Finding)) + 
  geom_bar(position="fill") + 
  scale_fill_manual(values = cols, 
                    name = "Auditor Finding:", 
                    labels = function(x) stringr::str_wrap(x, width = 5)) + 
  theme_minimal() + 
  theme(legend.title = element_text(face = "bold"),
        strip.text = element_text(face = "bold"),
        axis.title = element_text(face = "bold"),
        title = element_text(face = "bold"),
        legend.position = "right") + 
  scale_x_discrete(labels = function(x) stringr::str_wrap(x, width = 8)) +
  scale_y_continuous(labels = scales::percent) + 
  labs(x = "Component of the Allostatic Load Index",
       y = "Proportion of Validated Measurements", 
       title = "A) Findings from Pilot + Wave I Validation") + 
  coord_flip()
long_audit_dat |>
  dplyr::group_by(Variable_Name, Finding) |> 
  dplyr::summarize(Num = dplyr::n()) |> 
  dplyr::group_by(Variable_Name) |> 
  dplyr::mutate(Prop = Num / sum(Num)) |> 
  write.csv("~/Documents/Allostatic_load_audits/wave1_findings.csv", 
            row.names = FALSE)

## Save it
ggsave(filename = "~/Documents/ALI_EHR/figures/bargraph_wave1_audit_findings.png", 
       device = "png", width = 8, height = 4.5, units = "in")

# Load data ( wave II) -----------------------------------------
## Wave II audits 
long_audit_dat2 = read.csv("~/Documents/Allostatic_load_audits/Wave2_Last48_StartingJune/wave2_audits_for_analysis.csv") |> 
  filter(!(Variable_Name %in% c("HEIGHT", "WEIGHT")))
nrow(long_audit_dat2) ## 3867 audited data points 
table(long_audit_dat2$Category) ## 718 labs, 3149 vitals
long_audit_dat2 |> 
  group_by(PAT_MRN_ID) |> 
  summarize(NUM_AUDITED = dplyr::n()) |> 
  pull(NUM_AUDITED) |> 
  summary() ## min = 16, median = 65.5, max = 453 data points per patient

## ALI components after wave 1 validation
long_audit_dat2 |> 
  mutate(
    COMP = factor(x = Variable_Name, 
                  levels = order_levels, 
                  labels = c("SYSTOLIC BLOOD PRESSURE", 
                             "DIASTOLIC BLOOD PRESSURE", 
                             "BODY MASS INDEX", 
                             "SERUM ALBUMIN", 
                             "CHOLEST-\nEROL", 
                             "TRIGLY-\nCERIDES", 
                             "HEMOGLOBIN A1C", 
                             "CREATININE CLEARANCE", 
                             "C-REACTIVE PROTEIN", 
                             "HOMO-\nCYSTEINE"
                  )), 
    Finding = factor(x = Finding, 
                     levels = c("Extracted Value Correct", 
                                "Extracted Value Incorrect", 
                                "Extracted Value Not Found", 
                                "No Auxiliary Information Found",
                                "Auxiliary Information Found" 
                     ))
  ) |> 
  ggplot(aes(x = COMP, fill = Finding)) + 
  geom_bar(position="fill") + 
  scale_fill_manual(values = cols, 
                    name = "Auditor Finding:", 
                    labels = function(x) stringr::str_wrap(x, width = 5)) + 
  theme_minimal() + 
  theme(legend.title = element_text(face = "bold"),
        strip.text = element_text(face = "bold"),
        axis.title = element_text(face = "bold"),
        title = element_text(face = "bold"),
        legend.position = "right") + 
  scale_x_discrete(labels = function(x) stringr::str_wrap(x, width = 8)) +
  scale_y_continuous(labels = scales::percent) + 
  labs(x = "Component of the Allostatic Load Index",
       y = "Proportion of Validated Measurements", 
       title = "B) Wave II Validation") + 
       #title = "B) Pilot + Wave I + Wave II Validation") + 
  coord_flip()

## Save it
ggsave(filename = "~/Documents/ALI_EHR/figures/bargraph_wave2only_audit_findings.png", 
       device = "png", width = 8, height = 4.5, units = "in")

## Read in long audit data (one row per patient per audited variable)
all_audit_dat = long_audit_dat |> 
  bind_rows(long_audit_dat2)

## ALI components after wave 1 validation
all_audit_dat |> 
  mutate(
    COMP = factor(x = Variable_Name, 
                  levels = order_levels, 
                  labels = c("SYSTOLIC BLOOD PRESSURE", 
                             "DIASTOLIC BLOOD PRESSURE", 
                             "BODY MASS INDEX", 
                             "SERUM ALBUMIN", 
                             "CHOLEST-\nEROL", 
                             "TRIGLY-\nCERIDES", 
                             "HEMOGLOBIN A1C", 
                             "CREATININE CLEARANCE", 
                             "C-REACTIVE PROTEIN", 
                             "HOMO-\nCYSTEINE"
                  )), 
    Finding = factor(x = Finding, 
                     levels = c("Extracted Value Correct", 
                                "Extracted Value Incorrect", 
                                "Extracted Value Not Found", 
                                "No Auxiliary Information Found",
                                "Auxiliary Information Found" 
                     ))
  ) |> 
  ggplot(aes(x = COMP, fill = Finding)) + 
  geom_bar(position="fill") + 
  scale_fill_manual(values = cols, 
                    name = "Auditor Finding:", 
                    labels = function(x) stringr::str_wrap(x, width = 5)) + 
  theme_minimal() + 
  theme(legend.title = element_text(face = "bold"),
        strip.text = element_text(face = "bold"),
        axis.title = element_text(face = "bold"),
        title = element_text(face = "bold"),
        legend.position = "right") + 
  scale_x_discrete(labels = function(x) stringr::str_wrap(x, width = 8)) +
  scale_y_continuous(labels = scales::percent) + 
  labs(x = "Component of the Allostatic Load Index",
       y = "Proportion of Validated Measurements") + 
  coord_flip()

## Save it
ggsave(filename = "~/Documents/ALI_EHR/figures/bargraph_bothwaves_audit_findings.png", 
       device = "png", width = 8, height = 4.5, units = "in")
