################################################################################
## Setup ///////////////////////////////////////////////////////////////////////
################################################################################
library(dplyr) ## for data wrangling
library(tidyr) ## for data wrangling
library(irr) ## for Fleiss' kappa 

################################################################################
## Read in data ////////////////////////////////////////////////////////////////
################################################################################
pilot_audited = read.csv("~/Documents/Allostatic_load_audits/Pilot_First4/pilot_audits_for_analysis.csv")

################################################################################
## Make a table of frequencies of each finding per data point //////////////////
################################################################################
unique(pilot_audited$Finding_Aidan) ## there are 5 audit findings
# [1] "Extracted Value Correct"        "Extracted Value Not Found"     
# [3] "Extracted Value Incorrect"      "No Auxiliary Information Found"
# [5] "Auxiliary Information Found"   

## Define variables for labs and vitals 
lab_vars = c("CRP", "HCST", "CREAT_C", "CHOL", "TRIG", "A1C", "ALB")
vitals_vars = c("BMI", "BP_SYSTOLIC", "BP_DIASTOLIC")

## Subset to them (exclude height and weight)
pilot_audited = pilot_audited |> 
  dplyr::filter(Variable_Name %in% c(lab_vars, vitals_vars))

pilot_audited |> 
  dplyr::filter(Finding_Rabeya != "Not Audited") |> 
  dplyr::mutate(
    Same_Finding = dplyr::case_when(
      Category == "Labs" & 
        Finding_Sheetal == Finding_Rabeya & 
        Finding_Rabeya == Finding_Aidan ~ TRUE, 
      Category == "Vitals" & 
        Finding_Sheetal == Finding_Rabeya & 
        Finding_Rabeya == Finding_Aidan & 
        Finding_Aidan == Finding_Amelia ~ TRUE, 
      .default = FALSE
    ),
    Same_Finding = factor(x = Same_Finding, 
                          levels = c(FALSE, TRUE), 
                          labels = c("Not in Agreement", "In Agreement")), 
    Variable_Name = factor(x = Variable_Name, 
                           levels = c(lab_vars, vitals_vars), 
                           labels = c("C-Reactive Protein", "Homocysteine", 
                                      "Creatinine Clearance", "Total Cholesterol", 
                                      "Trigycerides", "Hemoglobin A1C", 
                                      "Serum Albumin", "Body Mass Index", 
                                      "Systolic Blood Pressure", "Diastolic Blood Pressure"))) |> 
  ggplot(aes(x = Variable_Name, fill = Same_Finding)) + 
  geom_bar(color = "black") + 
  scale_fill_manual(values = c("#ff99ff", "#8bdddb"), 
                    name = "Auditors' Findings") + 
  #facet_wrap(~Category, scales = "free") + 
  labs(x = "Variable", 
       y = "Data Points") + 
  #coord_flip() + 
  theme_minimal() + 
  theme(legend.position = "top") + 
  scale_x_discrete(labels = function(x) stringr::str_wrap(x, width = 10))

## Save the plot 
ggsave("~/Documents/ALI_EHR/figures/pilot_agreement_box.png", 
       device = "png", width = 8, height = 5)

## Four auditors reviewed the lab variables
vitals_freq = pilot_audited |> 
  filter(Variable_Name %in% vitals_vars, 
         Finding_Rabeya != "Not Audited") |> 
  select(Finding_Sheetal, Finding_Rabeya, Finding_Aidan, Finding_Amelia)
dim(vitals_freq) ## N = 457 vitals variables were validated by n = 4 auditors 

### Test agreement in vitals variables 
kappam.fleiss(vitals_freq)               # Fleiss' Kappa (0.81)
kappam.fleiss(vitals_freq, exact=TRUE)   # Exact Kappa (0.81)

## Three auditors reviewed the lab variables
lab_freq = pilot_audited |> 
  filter(Variable_Name %in% lab_vars) |> 
  select(Finding_Sheetal, Finding_Rabeya, Finding_Aidan)
dim(lab_freq) ## N = 148 lab variables were validated by n = 3 auditors 

### Test agreement in lab variables 
kappam.fleiss(lab_freq)               # Fleiss' Kappa (0.90)
kappam.fleiss(lab_freq, exact=TRUE)   # Exact Kappa (0.90)

