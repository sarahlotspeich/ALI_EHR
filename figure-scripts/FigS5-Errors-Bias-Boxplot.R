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

## Describe non-convergence due to max iterations reached (<= 1.2% of reps)
sim_res |> 
  group_by(tpr, fpr, val_design, smle_conv_msg) |> 
  summarize(num = dplyr::n()) |> 
  ungroup() |> 
  filter(!smle_conv_msg) |> 
  arrange(desc(num))

## Describe resampling due to empty B-spline sieve
### This was most common among optimal sample
sim_res |> 
  group_by(tpr, fpr, val_design, empty_sieve) |> 
  summarize(num = dplyr::n(), 
            median_nsieve = median(nsieve)) |> 
  ungroup() |> 
  filter(empty_sieve) |> 
  arrange(desc(num))

## Calculate average bias for the SMLE analysis 
sim_res |> 
  group_by(tpr, fpr) |> 
  summarize(bias = mean(smle_beta1 - beta1, na.rm = TRUE) / beta1)

## Calculate average bias for the naive analysis 
sim_res |> 
  group_by(tpr, fpr) |> 
  summarize(bias = mean(naive_beta1 - beta1, na.rm = TRUE) / beta1)

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

## Plot boxplot of coefficient estimates
sim_res |> 
  ggplot(aes(x = error_sett, y = smle_beta1, fill = val_design)) +
  geom_hline(yintercept = beta1,
             linetype = 2,
             color = "black") +
  annotate(geom = "text", 
           x = 0.3, 
           y = beta1 + 0.05, 
           label = "Truth =\n1.88", size = 3) + 
  geom_boxplot() +
  scale_fill_manual(values = paper_colors,
                    name = "Validation\nStudy\nDesign:", 
                    labels = function(x) str_wrap(x, width = 4)) +
  scale_x_discrete(expand = expand_scale(add = c(1,1))) + 
  theme_minimal(base_size = 14) +
  theme(legend.position = "right",
        title = element_text(face = "bold"),
        axis.title = element_text(face = "bold"),
        axis.text = element_markdown(),
        legend.title = element_text(face = "bold")) + 
  labs(x = "Error Rates in Allostatic Load Index (ALI) Components from EHR",
       y = "Estimated Log Odds Ratio on ALI")
ggsave(filename = "~/Documents/ALI_EHR/figures/FigS5_Errors_logOR_Boxplot.png",
       device = "png", width = 8, height = 5, units = "in")