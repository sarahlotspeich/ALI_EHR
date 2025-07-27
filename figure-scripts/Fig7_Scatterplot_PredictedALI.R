# Load packages ----------------------------------------------------------------
library(dplyr) ## for data wrangling
### RUN ONCE: devtools::install_github("sarahlotspeich/logiSieve", ref = "main")
library(logiSieve) ## for the SMLEs

# Load data --------------------------------------------------------------------
## ALI components after all waves of validation 
val_data = read.csv("~/Documents/Allostatic_load_audits/all_ali_dat.csv") |> 
  filter(DATA == "Pilot + Wave I Validation", 
         VALIDATED) |> 
  bind_rows(
    read.csv("~/Documents/Allostatic_load_audits/all_ali_dat.csv") |> 
      filter(DATA == "Wave II Validation", 
             VALIDATED)
  )
nrow(val_data) ## CHECK: n = 100 validated 
## ALI before validation and age / hospitalizations
unval_data = read.csv("~/Documents/Allostatic_load_audits/all_ali_dat.csv") |> 
  filter(DATA == "EHR (Before Validation)") |> 
  rename(ALI_STAR = ALI)
nrow(unval_data) ## CHECK: N = 1000 unvalidated 
## Merge validated + unvalidated 
data = val_data |> 
  select(PAT_MRN_ID, ALI) |> 
  right_join(
    unval_data |> 
      select(PAT_MRN_ID, ANY_ENCOUNTERS, AGE_AT_ENCOUNTER, ALI_STAR)
  )

## Recenter age at 18 and rescale to be in 10-year increments
data = data |> 
  mutate(AGE_AT_ENCOUNTER_10 = (AGE_AT_ENCOUNTER - 18) / 10)

# Estimate parameters using all audits + the rest of Phase I -------------------
## Setup B-splines
B = splines::bs(x = data$ALI_STAR, 
                df = 16, 
                Boundary.knots = range(data$ALI_STAR), 
                intercept = TRUE, 
                degree = 3)
colnames(B) = paste0("bs", seq(1, 16))
data = data |> 
  bind_cols(B)

## Fit SMLE model Y ~ X + Z
fit = logiSieve(
  analysis_formula = ANY_ENCOUNTERS ~ ALI + AGE_AT_ENCOUNTER_10, 
  error_formula = paste("ALI ~", paste(colnames(B), collapse = "+")), 
  data = data, 
  output = "all") #### request "all" output to get B-spline coefficients 

## Extract B-spline coefficients 
data$ALI_HAT = fit$predicted

## Make plot 
paper_colors = c("#ff99ff", "#787ff6", "#8bdddb", "#7dd5f6", "#ffbd59")
library(ggplot2)
data |> 
  dplyr::mutate(Validated = factor(x = !is.na(ALI), 
                                   levels = c(TRUE, FALSE), 
                                   labels = c("Chart Review", "Predicted"))) |> 
  ggplot(aes(x = ALI_STAR, y = ALI_HAT)) + 
  geom_point(aes(color = Validated, shape = Validated), size = 2) + 
  geom_smooth(color = paper_colors[3], linewidth = 1.2, se = FALSE) + 
  geom_abline(intercept = 0, slope = 1, linetype = 2, 
              linewidth = 1.2) + 
  scale_color_manual(values = paper_colors, 
                     name = "Source of\nValidated:", 
                     labels = function(x) stringr::str_wrap(x, width = 4)) + 
  scale_shape_manual(values = c(17, 15), 
                     name = "Source of\nValidated:", 
                     labels = function(x) stringr::str_wrap(x, width = 4)) + 
  theme_minimal(base_size = 14) +
  theme(legend.position = "right",
        title = element_text(face = "bold"),
        axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold")) + 
  labs(x = "Unvalidated Allostatic Load Index\n(from the EHR)",
       y = "Validated Allostatic Load Index\n(from Chart Reviews / Predictions)")
ggsave(filename = "~/Documents/ALI_EHR/figures/Fig7_Scatterplot_PredictedALI.png",
       device = "png", width = 8, height = 5, units = "in")