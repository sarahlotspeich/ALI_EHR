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
url_stem = "https://raw.githubusercontent.com/sarahlotspeich/ALI_EHR/main/sim-data/vary_wave2_design/vary_audit_recovery/"
recov_setts = paste0("recover", c(25, 50, 90, 100))

sim_res = data.frame()
for(des in c("wave2_srs_", "wave2_cc_", "wave2_bcc_", "wave2_score_", "wave2_residual_", "wave2_optMLE2_")) {
  des_urls = paste0(url_stem, des, recov_setts, ".csv") ## URLs to results on GitHub
  des_res = do.call(what = bind_rows,
                    args = lapply(X = des_urls,
                                  FUN = read.csv)) ## Read them in and stack them 
  if (des == "wave2_srs_") {
    sim_res = des_res
  } else {
    sim_res = sim_res |> 
      left_join(des_res)
  }
}

## The 100% audit recovery sims need to have audit_recovery defined
sim_res = sim_res |> 
  mutate(audit_recovery = replace_na(data = audit_recovery, 
                                            replace = 1))

## Describe non-convergence due to max iterations reached (<= 3% of reps)
sim_res |> 
  select(audit_recovery, ends_with("conv_msg")) |> 
  gather(key = "design", value = "conv_msg", -1) |> 
  group_by(design, audit_recovery, conv_msg) |> 
  summarize(num = n()) |> 
  ungroup() |> 
  filter(conv_msg != "Converged") |> 
  arrange(desc(num))

## Describe resampling due to empty B-spline sieve
### This was most common among score function and optMLE-2 sample
sim_res |> 
  select(audit_recovery, ends_with("resampled")) |> 
  gather(key = "design", value = "resampled", -1) |> 
  group_by(design, audit_recovery, resampled) |> 
  summarize(num = n()) |> 
  ungroup() |> 
  filter(resampled) |> 
  arrange(desc(num))

# Create long data with just beta 1
long_beta1 = sim_res |> 
  select(audit_recovery, ends_with("beta1")) |> 
  gather(key = "design", value = "est", -1) |> 
  mutate(design = factor(x = design, 
                         levels = c("srs_beta1", "cc_beta1", "bcc_beta1", "score_beta1", "resid_beta1", "optMLE2_beta1"), 
                         labels = c("SRS", "CC", "BCC*", "Weighted Score", "Residual", "optMLE-2")))

## Plot boxplot of coefficient estimates
long_beta1 |>
  filter(design != "Weighted Score") |> 
  mutate(audit_recovery = factor(x = audit_recovery, 
                                        levels = c(1, 0.9, 0.5, 0.25), 
                                        labels = c("100%", "90%", "50%", "25%"))) |> 
  ggplot(aes(x = audit_recovery, y = est, fill = design)) +
  geom_boxplot() +
  scale_fill_manual(values = slide_colors,
                    name = "Validation Study Design:") +
  geom_hline(yintercept = beta1,
             linetype = 2,
             color = "black") +
  xlab("Audit Recovery Rate for Missing Allostatic Load Index Components") +
  ylab("Estimated Log Odds Ratio") +
  theme_minimal(base_size = 14) +
  theme(legend.position = "right",
        axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold"))
ggsave(filename = "~/Documents/ALI_EHR/figures/vary_wave2_design_recover_boxplot.png",
       device = "png", width = 8, height = 5, units = "in")

# Plot line graph of efficiency
long_beta1 |>
  filter(design != "Weighted Score") |> 
  mutate(audit_recovery = factor(x = audit_recovery, 
                                        levels = c(1, 0.9, 0.5, 0.25), 
                                        labels = c("100%", "90%", "50%", "25%"))) |> 
  group_by(audit_recovery, design) |>
  summarize(eff = 1 / var(est, na.rm = TRUE)) |> 
  ggplot(aes(x = audit_recovery, y = eff, color = design, group = design)) +
  geom_point(size = 2) +
  geom_line(linewidth = 1.2) +
  scale_color_manual(values = slide_colors,
                     name = "Design:") +
  xlab("Missing Allostatic Load Index Components Recovered") +
  ylab("Efficiency") +
  theme_minimal(base_size = 18) +
  theme(legend.position = "right",
        axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold"))
ggsave(filename = "~/Documents/ALI_EHR/figures/vary_wave2_design_recover_line_efficiency.png",
       device = "png", width = 8, height = 5, units = "in")

# Plot line graph of efficiency
srs_eff = long_beta1 |>
  filter(design != "Weighted Score") |> 
  group_by(audit_recovery, design) |>
  summarize(eff = 1 / var(est, na.rm = TRUE)) |> 
  filter(design == "SRS") |>
  rename(srs_eff = eff) |> 
  select(-design)
long_beta1 |>
  filter(design != "Weighted Score") |> 
  group_by(audit_recovery, design) |>
  summarize(eff = 1 / var(est, na.rm = TRUE))  |> 
  left_join(srs_eff) |> 
  mutate(re = eff / srs_eff, 
         audit_recovery = factor(x = audit_recovery, 
                                 levels = c(1, 0.9, 0.5, 0.25), 
                                 labels = c("100%", "90%", "50%", "25%"))) |> 
  
  ggplot(aes(x = audit_recovery, y = re, color = design, group = design)) +
  geom_point(size = 2) +
  geom_line(linewidth = 1.2) +
  scale_color_manual(values = slide_colors,
                     name = "Design:") +
  xlab("Missing Allostatic Load Index Components Recovered") +
  ylab("Relative Efficiency (to SRS)") +
  theme_minimal(base_size = 18) +
  theme(legend.position = "right",
        axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold"))
ggsave(filename = "~/Documents/ALI_EHR/figures/vary_wave2_design_recover_line_relative_efficiency.png",
       device = "png", width = 8, height = 5, units = "in")
