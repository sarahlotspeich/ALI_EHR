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
unval_data = read.csv("~/Documents/Allostatic_load_audits/all_ali_dat.csv") |> 
  filter(DATA == "EHR (Before Validation)") |> 
  select(-DATA, -VALIDATED, -starts_with("ALI"), -ANY_ENCOUNTERS, -AGE_AT_ENCOUNTER, -SEX) |> 
  gather(key = "COMPONENT", value = "UNVAL", -1)

## ALI components after Wave II validation 
val_data = read.csv("~/Documents/Allostatic_load_audits/all_ali_dat.csv") |> 
  filter(DATA == "Wave II Validation", 
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
with(data, sum(UNVAL == "Yes" & VAL == "Yes") / sum(VAL == "Yes" & UNVAL != "Missing")) #### 99%

## Calculate true positive rate (FPR) = P(UNVAL = 1 | VAL = 0)
with(data, sum(UNVAL == "Yes" & VAL == "No") / sum(VAL == "No" & UNVAL != "Missing")) #### 3%

## Calculate missing data recovery rate 
### i.e., percent of missing data from EHR that were nonmissing after validation
with(data, sum(UNVAL == "Missing" & VAL != "Missing") / sum(UNVAL == "Missing")) #### 17% 

## Plot boxplot of coefficient estimates
all_combo = expand.grid(VAL = c("Yes", "No", "Missing"), 
                        UNVAL = c("Yes", "No", "Missing"))

data |> 
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
        legend.title = element_text(face = "bold")) + 
  labs(x = "Unvalidated Component from EHR",
       y = "Validated Component from Chart Review", 
       title = "B) Wave II Validation")
ggsave(filename = "~/Documents/ALI_EHR/figures/Fig6B_WaveII_Heatmap.png",
       device = "png", width = 8, height = 4, units = "in")
