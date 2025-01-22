# //////////////////////////////////////////////////////////////////////
# Replicate Figure S1 in Supplementary Materials  //////////////////////
# Caption begins "Estimated log prevalence ratios for food access $X_P$
# on health. The five possible ways to include the analysis model... ///
# //////////////////////////////////////////////////////////////////////

# Load packages
library(dplyr) # To wrangle data
library(tidyr) # To transform data
library(ggplot2) # To create plots
library(ggthemes) # To get colorblind palletes
library(latex2exp) # To create LaTex labels for plots

# Read in data
sim_dir = "~/Documents/ALI_EHR/sim-data/"
sim_files = paste0(sim_dir, list.files(path = sim_dir, pattern = "tpr"))
sim_res = do.call(what = dplyr::bind_rows,
                  args = lapply(X = sim_files,
                                FUN = read.csv))

# True parameter values from the sims
beta0 = -1.57 ## intercept in model of Y|X,Z
beta1 = 0.95 ## coefficient on X in model of Y|X,Z
beta2 = 0.01 ## coefficient on Z in model of Y|X,Z

# Reshape the data
sim_res = sim_res |>
  dplyr::mutate(error_sett = paste0("TPR = ", round(100 * tpr), "%, FPR = ", round(100 * fpr), "%"),
                error_sett = factor(x = error_sett,
                                    levels = c("TPR = 99%, FPR = 1%",
                                               "TPR = 95%, FPR = 5%",
                                               "TPR = 80%, FPR = 20%",
                                               "TPR = 50%, FPR = 50%"))) |>
  dplyr::select(sim, error_sett, dplyr::ends_with("beta1")) |>
  tidyr::gather("method", "beta", -c(1:2)) |>
  dplyr::mutate(method = factor(x = method,
                                levels = c("naive_beta1", "gs_beta1", "smle_beta1", "cc_beta1"),
                                labels = c("Naive", "Gold Standard", "SMLE", "Complete Case")))

# Plot results
slide_colors = c("#C3CFFA", "#C9B0B0", "#DD6D53", "#60CAD6", "#F0D290") ## colorblind_pal()(5)[-1]
sim_res |>
  ggplot(aes(x = error_sett, y = beta, fill = method)) +
  geom_boxplot() +
  scale_fill_manual(values = slide_colors,
                    name = "Method:") +
  geom_hline(yintercept = beta1,
             linetype = 2,
             color = "black") +
  xlab("Error Rates in Allostatic Load Index Components") +
  ylab("Estimated Coefficient") +
  theme_minimal(base_size = 18) +
  theme(legend.position = "top",
        axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold"))
ggsave(filename = "~/Dropbox (Wake Forest University)/5 - CONFERENCES/2 - Slides/2024/ALI-EHR-ENAR-Mar2024/vary_fpr_tpr.png",
       device = "png", width = 12, height = 7, units = "in")

# Read in data
sim_files = paste0(sim_dir, list.files(path = sim_dir, pattern = "recover"))
sim_res = do.call(what = dplyr::bind_rows,
                  args = lapply(X = sim_files,
                                FUN = read.csv)) |>
  dplyr::select(sim, prop_recovered, dplyr::ends_with("beta1")) |>
  tidyr::gather("method", "beta", -c(1:2)) |>
  dplyr::mutate(method = factor(x = method,
                                levels = c("naive_beta1", "gs_beta1", "smle_beta1", "cc_beta1"),
                                labels = c("Naive", "Gold Standard", "SMLE", "Complete Case")),
                prop_recovered = factor(x = prop_recovered,
                                        levels = c(1, 0.9, 0.5, 0.25),
                                        labels = c("100%", "90%", "50%", "25%")))

# Plot results
sim_res |>
  ggplot(aes(x = prop_recovered, y = beta, fill = method)) +
  geom_boxplot() +
  scale_fill_manual(values = slide_colors,
                    name = "Method:") +
  geom_hline(yintercept = beta1,
             linetype = 2,
             color = "black") +
  xlab("Percent of Missing Allostatic Load Index Components Recovered from Audit") +
  ylab("Estimated Coefficient") +
  theme_minimal(base_size = 18) +
  theme(legend.position = "top",
        axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold")) +
  ylim(c(-10, 10))
ggsave(filename = "~/Dropbox (Wake Forest University)/5 - CONFERENCES/2 - Slides/2024/ALI-EHR-ENAR-Mar2024/vary_audit_recovery.png",
       device = "png", width = 12, height = 7, units = "in")

# Read in data
sim_files = paste0(sim_dir, list.files(path = sim_dir, pattern = "vary_designs"))
sim_res = do.call(what = dplyr::bind_rows,
                  args = lapply(X = sim_files,
                                FUN = read.csv)) |>
  dplyr::mutate(error_sett = paste0("TPR = ", round(100 * tpr), "%, FPR = ", round(100 * fpr), "%"),
                error_sett = factor(x = error_sett,
                                    levels = c("TPR = 99%, FPR = 1%",
                                               "TPR = 95%, FPR = 5%",
                                               "TPR = 80%, FPR = 20%",
                                               "TPR = 50%, FPR = 50%"))) |>
  dplyr::select(sim, error_sett, dplyr::ends_with("beta1")) |>
  tidyr::gather("design", "beta", -c(1:2)) |>
  dplyr::mutate(design = factor(x = design,
                                levels = c("srs_beta1", "cc_beta1", "bcc_beta1", "resid_beta1", "sfs_beta1"),
                                labels = c("SRS", "CC", "BCC*", "Residual", "Score Function"))) |>
  dplyr::filter(!is.na(beta)) |>
  dplyr::group_by(error_sett, design) |>
  dplyr::mutate(group_id = 1:dplyr::n()) |>
  dplyr::filter(group_id <= 500)

# Plot results
sim_res |>
  ggplot(aes(x = error_sett, y = beta, fill = design)) +
  geom_boxplot() +
  scale_fill_manual(values = slide_colors,
                    name = "Method:") +
  geom_hline(yintercept = beta1,
             linetype = 2,
             color = "black") +
  xlab("Percent of Missing Allostatic Load Index Components Recovered from Audit") +
  ylab("Estimated Coefficient") +
  theme_minimal(base_size = 18) +
  theme(legend.position = "top",
        axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold")) +
  ylim(c(-10, 10))
ggsave(filename = "~/Dropbox (Wake Forest University)/5 - CONFERENCES/2 - Slides/2024/ALI-EHR-ENAR-Mar2024/vary_audit_designs_tpr_fpr.png",
       device = "png", width = 12, height = 7, units = "in")

# Read in data
sim_res = sim_res |>
  dplyr::group_by(error_sett, design) |>
  dplyr::summarize(eff = 1 / var(beta, na.rm = TRUE))

srs = sim_res |>
  dplyr::filter(design == "SRS") |>
  dplyr::select(error_sett, eff) |>
  dplyr::rename(eff_srs = eff)

sim_res = sim_res |>
  dplyr::filter(design != "SRS") |>
  dplyr::left_join(srs) |>
  dplyr::mutate(re = eff / eff_srs)

# Plot results
sim_res |>
  ggplot(aes(x = error_sett, y = re, color = design, group = design)) +
  geom_point(size = 2) +
  geom_line(linewidth = 1.2) +
  scale_color_manual(values = slide_colors,
                    name = "Design:") +
  geom_hline(yintercept = 1,
             linetype = 2,
             color = "black") +
  xlab("Error Rates in Allostatic Load Index Components") +
  ylab("Relative Efficiency to SRS") +
  theme_minimal(base_size = 18) +
  theme(legend.position = "top",
        axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold"))
ggsave(filename = "~/Dropbox (Wake Forest University)/5 - CONFERENCES/2 - Slides/2024/ALI-EHR-ENAR-Mar2024/vary_audit_designs_tpr_fpr_line.png",
       device = "png", width = 12, height = 7, units = "in")

# Read in data
sim_dir = "~/Documents/ALI_EHR/sim-data/vary_stratXstar_optMLE2/"
sim_files = paste0(sim_dir, list.files(path = sim_dir, pattern = "tpr"))
sim_res = do.call(what = dplyr::bind_rows,
                  args = lapply(X = sim_files,
                                FUN = read.csv)) |>
  dplyr::select(sim, tpr, fpr, dplyr::ends_with("beta1")) |>
  tidyr::gather("strat", "beta1", -c(1:3)) |>
  dplyr::mutate(error_sett = paste0("TPR = ", round(100 * tpr), "%, FPR = ", round(100 * fpr), "%"),
                error_sett = factor(x = error_sett,
                                    levels = c("TPR = 99%, FPR = 1%",
                                               "TPR = 95%, FPR = 5%",
                                               "TPR = 80%, FPR = 20%",
                                               "TPR = 50%, FPR = 50%")),
                strat = factor(x = strat,
                               levels = c("strat10_beta1",
                                          "strat25_beta1",
                                          "strat50_beta1",
                                          "strat75_beta1",
                                          "strat90_beta1" ),
                               labels = c("10th Percentile",
                                          "25th Percentile",
                                          "50th Percentile",
                                          "75th Percentile",
                                          "90th Percentile")))

sim_res |>
  ggplot(aes(x = error_sett, y = beta1, fill = strat)) +
  geom_boxplot() +
  scale_fill_manual(values = slide_colors,
                     name = "Design:") +
  geom_hline(yintercept = beta1,
             linetype = 2,
             color = "black") +
  xlab("Error Rates in Allostatic Load Index Components") +
  ylab("Relative Efficiency to SRS") +
  theme_minimal(base_size = 18) +
  theme(legend.position = "top",
        axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold"))
ggsave(filename = "~/Dropbox (Wake Forest University)/5 - CONFERENCES/2 - Slides/2024/ALI-EHR-ENAR-Mar2024/vary_stratXstar_optMLE2_tpr_fpr_boxplot.png",
       device = "png", width = 12, height = 7, units = "in")

sim_res = sim_res |>
  dplyr::group_by(error_sett, strat) |>
  dplyr::summarize(eff = 1 / var(beta1, na.rm = TRUE))

# Plot results
sim_res |>
  ggplot(aes(x = error_sett, y = eff, color = strat, group = strat)) +
  geom_point(size = 2) +
  geom_line(linewidth = 1.2) +
  scale_color_manual(values = slide_colors,
                     name = "Design:") +
  geom_hline(yintercept = 1,
             linetype = 2,
             color = "black") +
  xlab("Error Rates in Allostatic Load Index Components") +
  ylab("Relative Efficiency to SRS") +
  theme_minimal(base_size = 18) +
  theme(legend.position = "top",
        axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold"))
ggsave(filename = "~/Dropbox (Wake Forest University)/5 - CONFERENCES/2 - Slides/2024/ALI-EHR-ENAR-Mar2024/vary_stratXstar_optMLE2_tpr_fpr_line.png",
       device = "png", width = 12, height = 7, units = "in")

# Read in data
sim_dir = "~/Documents/ALI_EHR/sim-data/vary_stratXstar_optMLE/"
sim_files = paste0(sim_dir, list.files(path = sim_dir, pattern = "tpr"))
sim_res = do.call(what = dplyr::bind_rows,
                  args = lapply(X = sim_files,
                                FUN = read.csv)) |>
  dplyr::select(sim, tpr, fpr, dplyr::ends_with("beta1")) |>
  tidyr::gather("strat", "beta1", -c(1:3)) |>
  dplyr::mutate(error_sett = paste0("TPR = ", round(100 * tpr), "%, FPR = ", round(100 * fpr), "%"),
                error_sett = factor(x = error_sett,
                                    levels = c("TPR = 99%, FPR = 1%",
                                               "TPR = 95%, FPR = 5%",
                                               "TPR = 80%, FPR = 20%",
                                               "TPR = 50%, FPR = 50%")),
                strat = factor(x = strat,
                               levels = c("strat10_beta1",
                                          "strat25_beta1",
                                          "strat50_beta1",
                                          "strat75_beta1",
                                          "strat90_beta1" ),
                               labels = c("10th Percentile",
                                          "25th Percentile",
                                          "50th Percentile",
                                          "75th Percentile",
                                          "90th Percentile")))

sim_res |>
  ggplot(aes(x = error_sett, y = beta1, fill = strat)) +
  geom_boxplot() +
  scale_fill_manual(values = slide_colors,
                    name = "Design:") +
  geom_hline(yintercept = beta1,
             linetype = 2,
             color = "black") +
  xlab("Error Rates in Allostatic Load Index Components") +
  ylab("Relative Efficiency to SRS") +
  theme_minimal(base_size = 18) +
  theme(legend.position = "top",
        axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold"))
ggsave(filename = "~/Dropbox (Wake Forest University)/5 - CONFERENCES/2 - Slides/2024/ALI-EHR-ENAR-Mar2024/vary_stratXstar_optMLE_tpr_fpr_boxplot.png",
       device = "png", width = 12, height = 7, units = "in")

sim_res = sim_res |>
  dplyr::group_by(error_sett, strat) |>
  dplyr::summarize(eff = 1 / var(beta1, na.rm = TRUE))

# Plot results
sim_res |>
  ggplot(aes(x = error_sett, y = eff, color = strat, group = strat)) +
  geom_point(size = 2) +
  geom_line(linewidth = 1.2) +
  scale_color_manual(values = slide_colors,
                     name = "Design:") +
  geom_hline(yintercept = 1,
             linetype = 2,
             color = "black") +
  xlab("Error Rates in Allostatic Load Index Components") +
  ylab("Relative Efficiency to SRS") +
  theme_minimal(base_size = 18) +
  theme(legend.position = "top",
        axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold"))
ggsave(filename = "~/Dropbox (Wake Forest University)/5 - CONFERENCES/2 - Slides/2024/ALI-EHR-ENAR-Mar2024/vary_stratXstar_optMLE_tpr_fpr_line.png",
       device = "png", width = 12, height = 7, units = "in")
