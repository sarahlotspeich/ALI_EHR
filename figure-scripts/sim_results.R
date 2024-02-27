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
sim_res |> 
  ggplot(aes(x = error_sett, y = beta, fill = method)) + 
  geom_boxplot() + 
  scale_fill_manual(values = colorblind_pal()(5)[-1], 
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

# Read in data 
sim_dir = "~/Documents/ALI_EHR/sim-data/"
sim_files = paste0(sim_dir, list.files(path = sim_dir, pattern = "recover"))
sim_res = do.call(what = dplyr::bind_rows, 
                  args = lapply(X = sim_files, 
                                FUN = read.csv))

# Plot results 
sim_res |> 
  dplyr::mutate(prop_recovered = factor(x = prop_recovered, 
                                        levels = c(1, 0.9, 0.5, 0.25), 
                                        labels = c("100% Recovered", "90% Recovered", "50% Recovered", "25% Recovered"))) |> 
  dplyr::select(sim, prop_recovered, dplyr::ends_with("beta1")) |> 
  tidyr::gather("method", "beta", -c(1:2)) |> 
  ggplot(aes(x = prop_recovered, y = beta, fill = method)) + 
  geom_boxplot() + 
  scale_fill_manual(values = colorblind_pal()(5)[-1], 
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
