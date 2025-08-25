# Load packages 
library(tidyr) ## to unpivot data
library(dplyr) ## for data wrangling
library(ggplot2) ## for pretty plots

# Define colors 
cols = c("#ff99ff", "#8bdddb", "#787ff6", "#ffbd59", "#7dd5f6")

# Load data
## ALI components before and after validation (all waves)
all_data = read.csv("~/Documents/ALI_EHR/analysis/patient-data/all_ali_dat.csv") |> 
  mutate(DATA = case_when(DATA == "Wave II Validation" ~ "All Waves of Validation", 
                          .default = DATA))

# Create new dataframe with number of missing values per ALI component
num_miss = all_data |> 
  group_by(DATA) |> 
  select(-PAT_MRN_ID, -ANY_ENCOUNTERS, -AGE_AT_ENCOUNTER, -SEX, 
         -starts_with("ALI"), -VALIDATED) |> 
  summarize_all(function(x) sum(is.na(x))) |> 
  gather(key = "Variable", value = "Num_Missing", -1)

# Create bar plot of missing values per component (colored by wave)
order_levels = num_miss |> 
  filter(DATA == "EHR (Before Validation)") |> 
  arrange(desc(Num_Missing)) |> 
  pull(Variable)
num_miss |> 
  mutate(Variable = factor(x = Variable, 
                           levels = order_levels, 
                           labels = c("Creatinine Clearance", "Homo-\ncysteine",
                                      "C-Reactive Protein", "Hemoglobin A1C", 
                                      "Cholest-\nerol", "Trigly-\ncerides", 
                                      "Serum Albumin", "Body Mass Index", 
                                      "Systolic Blood Pressure", 
                                      "Diastolic Blood Pressure")), 
         DATA = factor(x = DATA, 
                       levels = c("EHR (Before Validation)", 
                                  "Pilot + Wave I Validation", 
                                  "Wave II Validation", 
                                  "All Waves of Validation"), 
                       labels = c("EHR (Before Validation)", 
                                  "Pilot and Wave I Validation", 
                                  "Wave II Validation", 
                                  "All Waves of Validation"))) |> 
  ggplot(aes(x = Variable, 
             y = Num_Missing, 
             fill = DATA)) + 
  geom_bar(stat = "identity", 
           position = "dodge", 
           color = "black") + 
  geom_text(aes(label=Num_Missing), 
            vjust = -1, 
            size = 3, 
            position = position_dodge(width = 1)) + 
  theme_minimal(base_size = 12) + 
  labs(x = "Allostatic Load Index Component", 
       y = "Number of Patients", 
       title = "A) Missing Values by Component") +
  theme(title = element_text(face = "bold"), 
        legend.position = "inside", 
        legend.position.inside = c(1, 0.9),
        legend.text = element_text(size = 12), 
        legend.title = element_text(size = 12, face = "bold"), 
        legend.justification = "right", 
        legend.background = element_rect(fill = "white")) + 
  scale_fill_manual(values = cols, name = "Data:") + 
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) + 
  scale_x_discrete(labels = function(x) stringr::str_wrap(x, width = 8))

## Save it 
ggsave(filename = "~/Documents/ALI_EHR/figures/FigS5A_missing_by_component.png", 
       device = "png", width = 10, height = 6, units = "in")
