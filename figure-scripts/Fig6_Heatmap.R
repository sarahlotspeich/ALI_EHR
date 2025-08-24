# //////////////////////////////////////////////////////////////////////
# Replicate Figure XX //////////////////////////////////////////////////
# Caption begins "XXX ..." /////////////////////////////////////////////
# //////////////////////////////////////////////////////////////////////

# Load packages
library(dplyr) ## To wrangle data
library(tidyr) ## To transform data
library(ggplot2) ## To create plots
library(ggtext) ## To italicize part of titles/labels
library(stringr) ## To wrap plot titles/text

# Define colors
paper_colors = c("#ff99ff", "#787ff6", "#8bdddb", "#7dd5f6", "#ffbd59")

# Load data
## ALI components before validation
unval_data = read.csv("~/Documents/ALI_EHR/analysis/patient-data/all_ali_dat.csv") |> 
  filter(DATA == "EHR (Before Validation)") |> 
  select(-DATA, -VALIDATED, -starts_with("ALI"), 
         -ANY_ENCOUNTERS, -AGE_AT_ENCOUNTER, -SEX) |> 
  gather(key = "COMPONENT", value = "UNVAL", -1)

## ALI components after Pilot and Wave I validation 
val_data = read.csv("~/Documents/ALI_EHR/analysis/patient-data/all_ali_dat.csv") |> 
  filter(VALIDATED) |> 
  select(-VALIDATED, -starts_with("ALI"), 
         -ANY_ENCOUNTERS, -AGE_AT_ENCOUNTER, -SEX) |> 
  gather(key = "COMPONENT", value = "VAL", -c(1, 12))
comb_val_data = val_data |> 
  mutate(DATA = "All Waves of Validation") |> 
  bind_rows(val_data)

## Merge validated + unvalidated 
data = comb_val_data |> 
  left_join(unval_data) |> 
  mutate(VAL = factor(x = VAL, 
                      levels = c(1, 0, NA), 
                      labels = c("Yes", "No", "Missing"), exclude = NULL), 
         UNVAL = factor(x = UNVAL, 
                        levels = c(1, 0, NA), 
                        labels = c("Yes", "No", "Missing"), exclude = NULL), 
         DATA = factor(x = DATA, 
                       levels = c("Pilot + Wave I Validation", 
                                  "Wave II Validation", 
                                  "All Waves of Validation"), 
                       labels = c("Pilot and Wave I Validation", 
                                  "Wave II Validation", 
                                  "All Waves of Validation")))

## Calculate error rates and data recovery
### True positive rate (TPR) = P(UNVAL = 1 | VAL = 1)
### False positive rate (FPR) = P(UNVAL = 1 | VAL = 0)
data |> 
  group_by(DATA) |> 
  summarize(TPR = sum(UNVAL == "Yes" & VAL == "Yes") / sum(VAL == "Yes" & UNVAL != "Missing"), 
            FPR = sum(UNVAL == "Yes" & VAL == "Yes") / sum(VAL == "Yes" & UNVAL != "Missing"), 
            Recovery = sum(UNVAL == "Missing" & VAL != "Missing") / sum(UNVAL == "Missing"))

## Plot boxplot of coefficient estimates
all_combo = expand.grid(DATA = c("Pilot and Wave I Validation", 
                                 "Wave II Validation", 
                                 "All Waves of Validation"), 
                        VAL = c("Yes", "No", "Missing"), 
                        UNVAL = c("Yes", "No", "Missing"))

data |> 
  group_by(DATA, VAL, UNVAL) |> 
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
  scale_fill_gradientn(colors = paper_colors, 
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
        legend.title = element_text(face = "bold"), 
        strip.background = element_rect(fill = "black"), 
        strip.text = element_text(face = "bold", color = "white")) + 
  labs(x = "Unvalidated Allostatic Load Index\nComponent (from the EHR)",
       y = "Validated Allostatic Load Index\nComponent (from Chart Review)") + 
  facet_wrap(~DATA)
ggsave(filename = "~/Documents/ALI_EHR/figures/Fig6_Heatmap.png",
       device = "png", width = 9, height = 4, units = "in")
