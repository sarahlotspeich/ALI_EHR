############################################################################################
## SETUP ///////////////////////////////////////////////////////////////////////////////////
############################################################################################
## Load libraries
library(ggplot2) ## to make maps
library(dplyr) ## for data wrangling
library(stringr) ## to wrap axis labels

## Read in forest plot data (all models)
plot_dat = read.csv("https://raw.githubusercontent.com/sarahlotspeich/ALI_EHR/refs/heads/main/analysis/forest_plot_data.csv")

## Define color scheme
cols = c("#ff99ff", "#8bdddb", "#787ff6", "#ffbd59", "#7dd5f6")

# Make the forest plot 
forest = plot_dat |> 
  mutate(Coefficient = factor(x = Coefficient, 
                              levels = c("Intercept", 
                                         "Age", 
                                         "ALI"),
                              labels = c("Intercept", 
                                         "Age",
                                         "Association Between Whole Person Health\nand Healthcare Utilization (Adjusting for Age)")),
         Analysis = factor(x = Analysis, 
                           levels = c("Unvalidated Data", 
                                      "Preliminary Validation Data", 
                                      "Final Validation Data"), 
                           labels = c("Original EHR (0 Patients Validated)", 
                                      "Wave I (52 Patients Validated)", 
                                      "Wave II (100 Patients Validated)"))) |> 
  filter(Coefficient == "Association Between Whole Person Health\nand Healthcare Utilization (Adjusting for Age)") |> 
  ggplot(aes(x = Analysis, y = Estimate)) + 
  geom_point(position = position_dodge(width = 0.5), 
             size = 3, color = cols[3]) + 
  geom_errorbar(aes(ymin = LB, ymax = UB), 
                linewidth = 1.2, color = cols[3]) + 
  geom_hline(yintercept = 1, linetype = 2, color = slide_cols[2]) + 
  scale_color_manual(values = slide_cols[c(1, 2, 5)], 
                     name = "Stage of Validation:") + 
  scale_x_discrete(labels = function(x) str_wrap(x, width = 15)) + 
  theme_minimal(base_size = 14) + 
  facet_wrap(~Coefficient) + 
  theme(axis.title = element_text(face = "bold"), 
        strip.background = element_rect(fill = "black"),
        strip.text = element_text(face = "bold", color = "white")) + 
  xlab("Stage of Validation") + 
  ylab("Odds Ratio (95% Confidence Interval)")
forest

ggsave(filename = "~/Documents/ALI_EHR/figures/forest_plot_ALI_only.png", 
       plot = forest, device = "png", width = 6, height = 6, units = "in")
