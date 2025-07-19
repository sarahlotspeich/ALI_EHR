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
  filter(!(Variable_Name %in% c("HEIGHT", "WEIGHT"))) |> 
  mutate(Audit = "Pilot + Wave I Validation") |> 
  select(PAT_MRN_ID, Audit, Category, Variable_Name, Finding)
nrow(long_audit_dat) ## 7605 audited data points 
table(long_audit_dat$Category) ## 1049 labs, 6556 vitals
long_audit_dat |> 
  group_by(PAT_MRN_ID) |> 
  summarize(NUM_AUDITED = dplyr::n()) |> 
  pull(NUM_AUDITED) |> 
  summary() ## min = 13, median = 65.5, max = 763 data points per patient

# Load data (wave II) -----------------------------------------
## Wave II audits 
long_audit_dat2 = read.csv("~/Documents/Allostatic_load_audits/Wave2_Last48_StartingJune/wave2_audits_for_analysis.csv") |> 
  filter(!(Variable_Name %in% c("HEIGHT", "WEIGHT"))) |> 
  mutate(Audit = "Wave II Validation") |> 
  select(PAT_MRN_ID, Audit, Category, Variable_Name, Finding)
nrow(long_audit_dat2) ## 3867 audited data points 
table(long_audit_dat2$Category) ## 718 labs, 3149 vitals
long_audit_dat2 |> 
  group_by(PAT_MRN_ID) |> 
  summarize(NUM_AUDITED = dplyr::n()) |> 
  pull(NUM_AUDITED) |> 
  summary() ## min = 16, median = 65.5, max = 453 data points per patient

## All Waves (Pilot + Wave I + Wave II)
long_audit_dat_all = long_audit_dat |> 
  bind_rows(long_audit_dat2) |> 
  mutate(Audit = "All Waves of Validation") 

## ALI components after wave 1 validation
order_levels = long_audit_dat |> 
  group_by(Variable_Name) |> 
  summarize(num = n()) |> 
  arrange(desc(num)) |> 
  pull(Variable_Name)

## Make bar graph of audit findings for numeric measurements
long_audit_dat |> 
  bind_rows(long_audit_dat2) |> 
  bind_rows(long_audit_dat_all) |> 
  mutate(
    COMP = factor(x = Variable_Name, 
                  levels = order_levels, 
                  labels = c("Systolic Blood Pressure", 
                             "Diastolic Blood Pressure", 
                             "Body Mass Index", 
                             "Serum Albumin", 
                             "Cholest-\nerol", 
                             "Trigly-\ncerides", 
                             "Hemoglobin A1C", 
                             "Creatinine Clearance", 
                             "C-Reactive Protein", 
                             "Homo-\ncysteine"
                             )), 
    Finding = factor(x = Finding, 
                     levels = rev(c("Extracted Value Correct", 
                                "Extracted Value Incorrect", 
                                "Extracted Value Not Found", 
                                "No Auxiliary Information Found",
                                "Auxiliary Information Found" 
                                ))), 
    Audit = factor(x = Audit, 
                   levels = c("Pilot + Wave I Validation", 
                              "Wave II Validation", 
                              "All Waves of Validation"))) |> 
  ggplot(aes(x = COMP, fill = Finding)) + 
  geom_bar(position="fill", color = "black") + 
  scale_fill_manual(values = rev(cols), 
                    name = "Auditor\nFinding:", 
                    labels = function(x) stringr::str_wrap(x, width = 5)) + 
  guides(fill = guide_legend(reverse=TRUE)) + 
  theme_minimal(base_size = 14) + 
  theme(legend.title = element_text(face = "bold"),
        strip.text = element_text(face = "bold", color = "white"),
        strip.background = element_rect(fill = "black"),
        axis.title = element_text(face = "bold"),
        title = element_text(face = "bold"),
        legend.position = "right", 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.key.spacing.y = unit(1/2, "cm")) + 
  scale_x_discrete(labels = function(x) stringr::str_wrap(x, width = 8)) +
  scale_y_continuous(labels = scales::percent) + 
  labs(x = "Allostatic Load Index Component",
       y = "Proportion of Validated Measurements") + 
  coord_flip() + 
  facet_wrap(~Audit)

## Save it
ggsave(filename = "~/Documents/ALI_EHR/figures/Fig5_Findings_Bargraph.png", 
       device = "png", width = 10, height = 6, units = "in")
