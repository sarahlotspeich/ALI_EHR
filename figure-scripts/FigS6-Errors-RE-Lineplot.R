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

# Define true parameter values from the sims
beta0 = -1.75 ## intercept in model of Y|X,Z
beta1 = 1.88 ## coefficient on X in model of Y|X,Z
beta2 = 0.10 ## coefficient on Z in model of Y|X,Z

# Read in simulation results
url_stem = "https://raw.githubusercontent.com/sarahlotspeich/ALI_EHR/main/sim-data/S5.2-Sims-Vary-Error-Rates/"
designs = c("SRS", "CC", "BCC", "optMLE", "RS")
sim_res = do.call(what = dplyr::bind_rows, 
                  args = lapply(X = paste0(url_stem, designs, ".csv"), 
                                FUN = read.csv))

## Create factor variables 
sim_res = sim_res |> 
  mutate(error_sett = paste0("TPR = ", round(100 * tpr), "%, FPR = ", round(100 * fpr), "%"), ## Create factor for error setting (TPR, FPR)
         error_sett = factor(x = error_sett,
                             levels = c("TPR = 99%, FPR = 1%",
                                        "TPR = 95%, FPR = 5%",
                                        "TPR = 80%, FPR = 20%",
                                        "TPR = 50%, FPR = 50%"),
                             labels = c("*Extra Low:*<br>TPR = 99%,<br>FPR = 1%",
                                        "*Low:*<br>TPR = 95%,<br>FPR = 5%",
                                        "*Moderate:*<br>TPR = 80%,<br>FPR = 20%",
                                        "*Heavy:*<br>TPR = 50%,<br>FPR = 50%")), 
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
  #filter(val_design == "Simple Random Sampling") |> 
  filter(val_design == "SRS") |> 
  group_by(error_sett) |>
  summarize(srs_eff = 1 / var(smle_beta1, na.rm = TRUE))
sim_res |>
  group_by(error_sett, val_design) |>
  summarize(eff = 1 / var(smle_beta1, na.rm = TRUE))  |> 
  dplyr::left_join(srs_eff) |> 
  dplyr::mutate(re = eff / srs_eff) |> 
  ggplot(aes(x = error_sett, y = re, 
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
  labs(x = "Error Rates in Allostatic Load Index (ALI) Components from EHR",
       #y = "Relative Efficiency (to Simple Random Sampling)",
       y = "Relative Efficiency to SRS")
ggsave(filename = "~/Documents/ALI_EHR/figures/FigS6_Errors_RE_Linegraph.png",
       device = "png", width = 8, height = 5, units = "in")
