# //////////////////////////////////////////////////////////////////////
# Replicate Figure XX //////////////////////////////////////////////////
# Caption begins "XXX ..." /////////////////////////////////////////////
# //////////////////////////////////////////////////////////////////////

# Load packages
library(dplyr) ## To wrangle data
library(tidyr) ## To transform data
library(ggplot2) ## To create plots

# Define colors
slide_colors = c("#C3CFFA", "#C9B0B0", "#DD6D53", "#60CAD6", "#F0D290", "#fac3cf") 

# Plot results using wave-1 SMLE parameter estimates to simulate data ////////
## Define true parameter values from the sims
beta0 = -1.93 ## intercept in model of Y|X,Z
beta1 = 1.88 ## coefficient on X in model of Y|X,Z
beta2 = 0.10 ## coefficient on Z in model of Y|X,Z

## Read in other designs
url_stem = "https://raw.githubusercontent.com/sarahlotspeich/ALI_EHR/main/sim-data/smle_other_designs/"
sim_urls = c(paste0(url_stem, "tpr50_fpr50_seed", 11422:11431, ".csv"),
             paste0(url_stem, "tpr80_fpr20_seed", 11422:11431, ".csv"),
             paste0(url_stem, "tpr95_fpr05_seed", 11422:11431, ".csv"),
             paste0(url_stem, "tpr99_fpr01_seed", 11422:11431, ".csv")
)
sim_res = do.call(what = bind_rows,
                  args = lapply(X = sim_urls,
                                FUN = read.csv)) |>
  dplyr::mutate(error_sett = paste0("TPR = ", round(100 * tpr), "%, FPR = ", round(100 * fpr), "%"),
                error_sett = factor(x = error_sett,
                                    levels = c("TPR = 99%, FPR = 1%",
                                               "TPR = 95%, FPR = 5%",
                                               "TPR = 80%, FPR = 20%",
                                               "TPR = 50%, FPR = 50%")))

# Read in optimal design
url_stem = "https://raw.githubusercontent.com/sarahlotspeich/ALI_EHR/main/sim-data/vary_stratXstar_wave1med/smle/"
sim_urls = c(paste0(url_stem, "vary_stratXstar_tpr50_fpr50_seed", 11422:11431, ".csv"),
             paste0(url_stem, "vary_stratXstar_tpr80_fpr20_seed", 11422:11431, ".csv"),
             paste0(url_stem, "vary_stratXstar_tpr95_fpr05_seed", 11422:11431, ".csv"),
             paste0(url_stem, "vary_stratXstar_tpr99_fpr01_seed", 11422:11431, ".csv")
)
sim_res = do.call(what = bind_rows,
                  args = lapply(X = sim_urls,
                                FUN = read.csv)) |>
  dplyr::filter(q == 0.9) |> 
  dplyr::mutate(design = "OPTMLE2",
                error_sett = paste0("TPR = ", round(100 * tpr), "%, FPR = ", round(100 * fpr), "%"),
                error_sett = factor(x = error_sett,
                                    levels = c("TPR = 99%, FPR = 1%",
                                               "TPR = 95%, FPR = 5%",
                                               "TPR = 80%, FPR = 20%",
                                               "TPR = 50%, FPR = 50%"))) |> 
  dplyr::bind_rows(sim_res)

## Exclude SFS (empty sieves) and optMLE
sim_res = sim_res |> 
  filter(design != "OPTMLE", 
         design != "SFS")

## Define levels for design
sim_res = sim_res |> 
  mutate(design = factor(x = design, 
                         levels = c("SRS", "CC", "BCC", "RESID", "OPTMLE2"), 
                         labels = c("SRS", "CC", "BCC*", "Residual", "optMLE-2")))

## Describe B-spline errors and non-convergence
sim_res |> 
  group_by(design, error_sett, conv_msg) |> 
  summarize(num = dplyr::n()) |> 
  ungroup() |> 
  filter(conv_msg != "Converged") |> 
  arrange(desc(num))

## Plot boxplot of coefficient estimates
sim_res |>
  ggplot(aes(x = error_sett, y = beta1, fill = design)) +
  geom_boxplot() +
  scale_fill_manual(values = slide_colors,
                    name = "Stratified X* At:") +
  geom_hline(yintercept = beta1,
             linetype = 2,
             color = "black") +
  xlab("Error Rates in Allostatic Load Index Components") +
  ylab("Estimated Coefficient") +
  theme_minimal(base_size = 18) +
  theme(legend.position = "top",
        axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold"))
ggsave(filename = "~/Documents/ALI_EHR/figures/smle_all_tpr_fpr_boxplot.png",
       device = "png", width = 12, height = 7, units = "in")

# Plot average stratum sizes per design
avg_des = sim_res |> 
  group_by(design, error_sett) |> 
  select(starts_with("des_")) |> 
  summarize_all(mean, na.rm = TRUE) |> 
  gather("stratum", "size", -c(1:2)) |> 
  mutate(stratum = factor(stratum, 
                          levels = c("des_n00", "des_n01", 
                                     "des_n10", "des_n11"), 
                          labels = c("(Y=0,X*=0)", "(Y=0,X*=1)", "(Y=1,X*=0)", "(Y=1,X*=1)")))
avg_des |> 
  ggplot(aes(x = design, y = stratum, fill = size)) + 
  geom_tile(color = "black", lwd = 1.2) + 
  facet_wrap(~error_sett) + 
  geom_text(aes(label = round(size)), size = 12) + 
  theme_minimal(base_size = 20) + 
  scale_fill_gradientn(
    colours = colorRampPalette(c("#C3CFFA", "#60CAD6", "#F0D290", "#DD6D53"))(100),
    guide = guide_colourbar(direction = "horizontal",
                            barwidth = 15,
                            barheight = 2),
    name = "Number of Patients Validated from Stratum:") + 
  theme(legend.position = "top", 
        strip.background = element_rect(fill = "#C3CFFA", color = "black", linewidth = 1.2),
        strip.text = element_text(color = "black")) + 
  ylab("Stratum")
ggsave(filename = "~/Documents/ALI_EHR/figures/smle_all_designs_heatmap.png",
       device = "png", width = 12, height = 7, units = "in")

# Plot line graph of efficiency
sim_res = sim_res |>
  dplyr::group_by(error_sett, design) |>
  dplyr::summarize(eff = 1 / var(beta1, na.rm = TRUE))
sim_res |>
  dplyr::filter(design != "SFS") |> 
  ggplot(aes(x = error_sett, y = eff, color = design, group = design)) +
  geom_point(size = 2) +
  geom_line(linewidth = 1.2) +
  scale_color_manual(values = slide_colors,
                     name = "Design:") +
  xlab("Error Rates in Allostatic Load Index Components") +
  ylab("Efficiency") +
  theme_minimal(base_size = 18) +
  theme(legend.position = "top",
        axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold"))
ggsave(filename = "~/Documents/ALI_EHR/figures/smle_all_tpr_fpr_line_eff.png",
       device = "png", width = 12, height = 7, units = "in")

# Plot line graph of efficiency
srs_eff = sim_res |> 
  dplyr::filter(design == "SRS") |>
  dplyr::rename(srs_eff = eff) |> 
  dplyr::select(-design)
sim_res = sim_res |> 
  dplyr::left_join(srs_eff) |> 
  dplyr::mutate(re = eff / srs_eff)
sim_res |>
  dplyr::filter(design != "SFS") |> 
  ggplot(aes(x = error_sett, y = re, color = design, group = design)) +
  geom_point(size = 2) +
  geom_line(linewidth = 1.2) +
  scale_color_manual(values = slide_colors,
                     name = "Design:") +
  xlab("Error Rates in Allostatic Load Index Components") +
  ylab("Efficiency") +
  theme_minimal(base_size = 18) +
  theme(legend.position = "top",
        axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold"))
ggsave(filename = "~/Documents/ALI_EHR/figures/smle_all_tpr_fpr_line_re.png",
       device = "png", width = 12, height = 7, units = "in")
