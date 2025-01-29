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
library(latex2exp) ## To include equations

# Define colors
paper_colors = c("#ff99ff", "#787ff6", "#8bdddb", "#7dd5f6", "#ffbd59")

# Define true parameter values from the sims
beta0 = -1.75 ## intercept in model of Y|X,Z
beta1 = 1.88 ## coefficient on X in model of Y|X,Z
beta2 = 0.10 ## coefficient on Z in model of Y|X,Z

# Read in simulation results
url_stem = "https://raw.githubusercontent.com/sarahlotspeich/ALI_EHR/main/sim-data/S5.3-Sims-Vary-Data-Recovery/"
designs = c("SRS", "CC", "BCC", "optMLE", "RS")
sim_res = do.call(what = dplyr::bind_rows, 
                  args = lapply(X = paste0(url_stem, designs, ".csv"), 
                                FUN = read.csv))

## Describe non-convergence due to max iterations reached (<= 3.4% of reps)
sim_res |> 
  group_by(prop_recovered, val_design, smle_conv_msg) |> 
  summarize(num = dplyr::n()) |> 
  ungroup() |> 
  filter(!smle_conv_msg) |> 
  arrange(desc(num))

## Describe resampling due to empty B-spline sieve
### This was most common among score function and optimal sample
sim_res |> 
  group_by(prop_recovered, val_design, empty_sieve) |> 
  summarize(num = dplyr::n(), 
            median_nsieve = median(nsieve)) |> 
  ungroup() |> 
  filter(empty_sieve) |> 
  arrange(desc(num))

## Create factor variables 
sim_res = sim_res |> 
  mutate(prop_recovered = factor(x = prop_recovered,
                                 levels = c(1, 0.9, 0.5, 0.25, 0),
                                 labels = c("*Full:*<br>100%",
                                            "*High:*<br>90%",
                                            "*Moderate:*<br>50%",
                                            "*Low:*<br>25%", 
                                            "*None:*<br>0%")), 
         val_design = factor(x = val_design, 
                             levels = c("SRS", "CC", "BCC*", "optMLE", "RS"), 
                             labels = c("SRS", "CC", "BCC", "OPT", "RS")
                             # labels = c("Simple Random Sampling", "Case-\nControl Sampling", 
                             #            "Balanced Case-\nControl Sampling", "Optimal Sampling", 
                             #            "Residual Sampling")
                             )
         ) 

# Plot line graph of efficiency
srs_eff = sim_res |>
  filter(val_design == "SRS") |> 
  group_by(prop_recovered) |>
  summarize(srs_eff = 1 / var(smle_beta1, na.rm = TRUE))
sim_res |>
  group_by(prop_recovered, val_design) |>
  summarize(eff = 1 / var(smle_beta1, na.rm = TRUE))  |> 
  dplyr::left_join(srs_eff) |> 
  dplyr::mutate(re = eff / srs_eff) |> 
  ggplot(aes(x = prop_recovered, y = re, 
             color = val_design, group = val_design)) +
  geom_point(size = 2) +
  geom_line(linewidth = 1.2, alpha = 1) +
  scale_color_manual(values = paper_colors,
                     name = "Validation\nStudy\nDesign:", 
                     labels = function(x) str_wrap(x, width = 4)) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "right",
        title = element_text(face = "bold"),
        axis.title = element_text(face = "bold"),
        axis.text = element_markdown(),
        legend.title = element_text(face = "bold")) + 
  labs(x = "Recovery Rate for Missing Allostatic Load Index (ALI) Components in the EHR",
       #y = "Relative Efficiency (to Simple Random Sampling)",
       y = "Relative Efficiency to SRS")
ggsave(filename = "~/Documents/ALI_EHR/figures/FigS8_Recovery_RE_Linegraph.png",
       device = "png", width = 8, height = 5, units = "in")
