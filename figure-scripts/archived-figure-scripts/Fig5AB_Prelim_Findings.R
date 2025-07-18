# Load packages
library(dplyr) ## To wrangle data
library(tidyr) ## To transform data
library(ggplot2) ## To create plots
library(ggtext) ## To italicize part of titles/labels
library(stringr) ## To wrap plot titles/text
library(ggpubr) ## To combine plots

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

## Make bar graph of audit findings for numeric measurements
plot_a = long_audit_dat |> 
  mutate(
    COMP = factor(x = Variable_Name, 
                  levels = order_levels, 
                  labels = c("Systolic Blood Pressure", 
                             "Diastolic Blood Pressure", 
                             "Body Mass Index", 
                             "Serum Albumin", 
                             "Cholest-\nerol", 
                             "Trigly-\nceries", 
                             "Hemoglobin A1C", 
                             "Creatinine Clearance", 
                             "C-Reactive Protein", 
                             "Homo-\ncysteine"
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
  labs(x = "Numeric Measurement in the Allostatic Load Index",
       y = "Proportion of Validated Measurements", 
       title = "A) Pilot + Wave I Validation") + 
  coord_flip()

## Save audit findings 
long_audit_dat |>
  dplyr::group_by(Variable_Name, Finding) |> 
  dplyr::summarize(Num = dplyr::n()) |> 
  dplyr::group_by(Variable_Name) |> 
  dplyr::mutate(Prop = Num / sum(Num)) |> 
  write.csv("~/Documents/Allostatic_load_audits/wave1_findings.csv", 
            row.names = FALSE)

# Load data
## ALI components before validation
unval_data = read.csv("~/Documents/Allostatic_load_audits/all_ali_dat.csv") |> 
  filter(DATA == "EHR (Before Validation)") |> 
  select(-DATA, -VALIDATED, -starts_with("ALI"), -ANY_ENCOUNTERS, -AGE_AT_ENCOUNTER, -SEX) |> 
  gather(key = "COMPONENT", value = "UNVAL", -1)

## ALI components after Pilot + Wave I validation 
val_data = read.csv("~/Documents/Allostatic_load_audits/all_ali_dat.csv") |> 
  filter(DATA == "Pilot + Wave I Validation", 
         VALIDATED) |> 
  select(-DATA, -VALIDATED, -starts_with("ALI"), -ANY_ENCOUNTERS, -AGE_AT_ENCOUNTER, -SEX) |> 
  gather(key = "COMPONENT", value = "VAL", -1)

## Merge validated + unvalidated 
data = val_data |> 
  left_join(unval_data) |> 
  mutate(VAL = factor(x = VAL, 
                      levels = c(1, 0, NA), 
                      labels = c("Yes", "No", "Missing"), exclude = NULL), 
         UNVAL = factor(x = UNVAL, 
                        levels = c(1, 0, NA), 
                        labels = c("Yes", "No", "Missing"), exclude = NULL))

## Calculate true positive rate (TPR) = P(UNVAL = 1 | VAL = 1)
with(data, sum(UNVAL == "Yes" & VAL == "Yes") / sum(VAL == "Yes" & UNVAL != "Missing")) #### 100%

## Calculate true positive rate (FPR) = P(UNVAL = 1 | VAL = 0)
with(data, sum(UNVAL == "Yes" & VAL == "No") / sum(VAL == "No" & UNVAL != "Missing")) #### <1%

## Calculate missing data recovery rate 
### i.e., percent of missing data from EHR that were nonmissing after validation
with(data, sum(UNVAL == "Missing" & VAL != "Missing") / sum(UNVAL == "Missing")) #### 27% 

## Plot boxplot of coefficient estimates
all_combo = expand.grid(VAL = c("Yes", "No", "Missing"), 
                        UNVAL = c("Yes", "No", "Missing"))

## Make heatmap of unval vs. val
plot_b = data |> 
  group_by(VAL, UNVAL) |> 
  summarize(num = n()) |>
  full_join(all_combo) |> 
  mutate(num = if_else(condition = is.na(num),
                       true = 0, 
                       false = num)) |> 
  ggplot(aes(x = UNVAL, y = VAL, fill = num)) +
  geom_tile(color = "black",
            lwd = 0.5,
            linetype = 1) + 
  geom_text(aes(label = num), 
            color = "black", 
            size = 5) + 
  scale_fill_gradientn(colors = cols, 
                       name = "Data\nPoints:", 
                       guide = guide_colorbar(frame.colour = "black", 
                                              ticks.colour = "black", 
                                              barwidth = 1, 
                                              barheight = 10)) + 
  theme_minimal(base_size = 14) +
  theme(legend.position = "right",
        title = element_text(face = "bold"),
        axis.title = element_text(face = "bold"),
        axis.text = element_markdown(),
        legend.title = element_text(face = "bold")) + 
  labs(x = "Unvalidated Component from EHR",
       y = "Validated Component from Chart Review", 
       title = "A) Pilot + Wave I Validation")

## Combine them 
ggarrange(plot_a, plot_b)

## Save it 
ggsave(filename = "~/Documents/ALI_EHR/figures/Fig6a_Prelim_Heatmap.png",
       device = "png", width = 8, height = 4, units = "in")
