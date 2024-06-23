# //////////////////////////////////////////////////////////////////////
# Replicate Figure XX //////////////////////////////////////////////////
# Caption begins "XXX ..." /////////////////////////////////////////////
# //////////////////////////////////////////////////////////////////////

# Load packages
library(dplyr) ## To wrangle data
library(tidyr) ## To transform data
library(ggplot2) ## To create plots

# Define colors
slide_colors = c("#38544f", "#78c0b5", "#dd6682", "#662708", "#ffcdd9")

# Define true parameter values from the sims
beta0 = -1.93 ## intercept in model of Y|X,Z
beta1 = 1.88 ## coefficient on X in model of Y|X,Z
beta2 = 0.10 ## coefficient on Z in model of Y|X,Z

# Read in simulation results
url_stem = "https://raw.githubusercontent.com/sarahlotspeich/ALI_EHR/main/sim-data/vary_wave2_design/vary_tpr_fpr/"
error_setts = c("tpr50_fpr50", "tpr80_fpr20", "tpr95_fpr05", "tpr99_fpr01")

sim_res = data.frame()
for(des in c("wave2_srs_", "wave2_cc_", "wave2_bcc_", "wave2_score_", "wave2_residual_", "wave2_optMLE2_")) {
  des_urls = paste0(url_stem, des, error_setts, ".csv") ## URLs to results on GitHub
  des_res = do.call(what = bind_rows,
                    args = lapply(X = des_urls,
                                  FUN = read.csv)) |> ## Read them in and stack them 
    dplyr::mutate(error_sett = paste0("TPR = ", round(100 * tpr), "%, FPR = ", round(100 * fpr), "%"), ## Create factor for error setting (TPR, FPR)
                  error_sett = factor(x = error_sett,
                                      levels = c("TPR = 99%, FPR = 1%",
                                                 "TPR = 95%, FPR = 5%",
                                                 "TPR = 80%, FPR = 20%",
                                                 "TPR = 50%, FPR = 50%"),
                                      labels = c("TPR = 99%,\nFPR = 1%",
                                                 "TPR = 95%,\nFPR = 5%",
                                                 "TPR = 80%,\nFPR = 20%",
                                                 "TPR = 50%,\nFPR = 50%")))
  if (des == "wave2_srs_") {
    sim_res = des_res
  } else {
    sim_res = sim_res |> 
      left_join(des_res)
  }
}

## Describe non-convergence due to max iterations reached (<= 2.1% of reps)
sim_res |> 
  select(error_sett, ends_with("conv_msg")) |> 
  gather(key = "design", value = "conv_msg", -1) |> 
  group_by(design, error_sett, conv_msg) |> 
  summarize(num = dplyr::n()) |> 
  ungroup() |> 
  filter(conv_msg != "Converged") |> 
  arrange(desc(num))

## Describe resampling due to empty B-spline sieve
### This was most common among score function and optimal sample
sim_res |> 
  select(error_sett, ends_with("resampled")) |> 
  gather(key = "design", value = "resampled", -1) |> 
  group_by(design, error_sett, resampled) |> 
  summarize(num = dplyr::n()) |> 
  ungroup() |> 
  filter(resampled) |> 
  arrange(desc(num))

# Create long data with just beta 1
long_beta1 = sim_res |> 
  select(error_sett, ends_with("beta1")) |> 
  gather(key = "design", value = "est", -1) |> 
  mutate(design = factor(x = design, 
                         levels = c("srs_beta1", "cc_beta1", "bcc_beta1", "score_beta1", "resid_beta1", "optMLE2_beta1"), 
                         labels = c("SRS", "CC", "BCC*", "Weighted Score", "Residual", "optMLE-2")))

## Plot boxplot of coefficient estimates
long_beta1 |>
  dplyr::filter(design != "Weighted Score") |> 
  ggplot(aes(x = error_sett, y = est, fill = design)) +
  geom_boxplot() +
  scale_fill_manual(values = slide_colors,
                    name = "Validation Study Design:") +
  geom_hline(yintercept = beta1,
             linetype = 2,
             color = "black") +
  xlab("Error Rates in Allostatic Load Index Components") +
  ylab("Estimated Log Odds Ratio") +
  theme_minimal(base_size = 14) +
  theme(legend.position = "right",
        axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold"))
ggsave(filename = "~/Documents/ALI_EHR/figures/vary_wave2_design_boxplot.png",
       device = "png", width = 8, height = 5, units = "in")

# Plot line graph of efficiency
long_beta1 |>
  dplyr::filter(design != "Weighted Score") |> 
  group_by(error_sett, design) |>
  summarize(eff = 1 / var(est, na.rm = TRUE)) |> 
  ggplot(aes(x = error_sett, y = eff, color = design, group = design)) +
  geom_point(size = 2) +
  geom_line(linewidth = 1.2) +
  scale_color_manual(values = slide_colors,
                     name = "Design:") +
  xlab("Error Rates in Allostatic Load Index Components") +
  ylab("Efficiency") +
  theme_minimal(base_size = 18) +
  theme(legend.position = "right",
        axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold"))
ggsave(filename = "~/Documents/ALI_EHR/figures/vary_wave2_design_line_efficiency.png",
       device = "png", width = 8, height = 5, units = "in")

# Plot line graph of efficiency
srs_eff = long_beta1 |>
  dplyr::filter(design != "Weighted Score") |> 
  group_by(error_sett, design) |>
  summarize(eff = 1 / var(est, na.rm = TRUE)) |> 
  dplyr::filter(design == "SRS") |>
  dplyr::rename(srs_eff = eff) |> 
  dplyr::select(-design)
long_beta1 |>
  dplyr::filter(design != "Weighted Score") |> 
  group_by(error_sett, design) |>
  summarize(eff = 1 / var(est, na.rm = TRUE))  |> 
  dplyr::left_join(srs_eff) |> 
  dplyr::mutate(re = eff / srs_eff) |> 
  ggplot(aes(x = error_sett, y = re, color = design, group = design)) +
  geom_point(size = 2) +
  geom_line(linewidth = 1.2) +
  scale_color_manual(values = slide_colors,
                     name = "Design:") +
  xlab("Error Rates in Allostatic Load Index Components") +
  ylab("Relative Efficiency (to SRS)") +
  theme_minimal(base_size = 18) +
  theme(legend.position = "right",
        axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold"))
ggsave(filename = "~/Documents/ALI_EHR/figures/vary_wave2_design_line_relative_efficiency.png",
       device = "png", width = 8, height = 5, units = "in")
