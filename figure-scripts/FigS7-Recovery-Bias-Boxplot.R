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
### This was most common among optimal sample
sim_res |> 
  group_by(prop_recovered, val_design, empty_sieve) |> 
  summarize(num = dplyr::n(), 
            median_nsieve = median(nsieve)) |> 
  ungroup() |> 
  filter(empty_sieve) |> 
  arrange(desc(num))

## There were also a few replications where the optMLE didn't converge 
### So no models were fit (all were reps with r = 0)
sim_res |> 
  filter(is.na(smle_conv_msg)) |> 
  group_by(prop_recovered, val_design) |> 
  summarize(num = dplyr::n()) |> 
  ungroup() |> 
  filter(num > 0)

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

## Plot boxplot of coefficient estimates
sim_res |> 
  ggplot(aes(x = prop_recovered, y = smle_beta1, fill = val_design)) +
  geom_hline(yintercept = beta1,
             linetype = 2,
             color = "black") +
  annotate(geom = "text", 
           x = 0.3, 
           y = beta1 + 0.05, 
           label = "Truth =\n1.88", size = 3) + 
  geom_hline(yintercept = mean(sim_res$naive_beta1),
             linetype = "dotted",
             color = "black") +
  annotate(geom = "text", 
           x = 5.7, 
           y = mean(sim_res$naive_beta1) + 0.05, 
           label = "Naive =\n0.83", size = 3) + 
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
  labs(x = "Recovery Rate for Missing Allostatic Load Index (ALI) Components in the EHR",
       y = "Estimated Log Odds Ratio on ALI")
ggsave(filename = "~/Documents/ALI_EHR/figures/FigS7_Recovery_logOR_Boxplot.png",
       device = "png", width = 8, height = 5, units = "in")