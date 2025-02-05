slide_cols = c("#ff99ff", "#8bdddb", "#787ff6", "#ffbd59", "#7dd5f6")

# Read in data 
orig_data = read.csv("~/Documents/Allostatic_load_audits/ali_dat.csv") |> 
  dplyr::rename(ALI_STAR = ALI) |> 
  dplyr::select(PAT_MRN_ID, ANY_ENCOUNTERS, AGE_AT_ENCOUNTER, ALI_STAR)
wave1_data = read.csv("~/Documents/Allostatic_load_audits/ali_dat_audit1.csv") |> 
  dplyr::select(PAT_MRN_ID, ANY_ENCOUNTERS, AGE_AT_ENCOUNTER, ALI)
summ_data = orig_data |> 
  dplyr::left_join(wave1_data) |> 
  dplyr::mutate(V = as.numeric(!is.na(ALI)))

# Rescale age at first encounter to be in 10-year increments and
# ALI/ALI* to be in 0.1-point increments. 
summ_data = summ_data |> 
  mutate(AGE_AT_ENCOUNTER10 = AGE_AT_ENCOUNTER / 10,
         ALI_01 = ALI / 0.1, 
         ALI_STAR_01 = ALI_STAR / 0.1)

# Naive model 
naive_mod = glm(formula = ANY_ENCOUNTERS ~ ALI_STAR_01 + AGE_AT_ENCOUNTER10, 
                family = "binomial", 
                data = summ_data)
knitr::kable(coefficients(summary(naive_mod)), digits = 3)
knitr::kable(round(exp(confint(object = naive_mod, level = 0.95)), 3))

plot_dat = data.frame(Method = "Naive", 
                      Coefficient = c("Intercept", "ALI/ALI*", "Age"),
                      OR = exp(naive_mod$coefficients), 
                      LB = exp(naive_mod$coefficients - 1.96 * sqrt(diag(vcov(naive_mod)))), 
                      UB = exp(naive_mod$coefficients + 1.96 * sqrt(diag(vcov(naive_mod)))))

# Estimate parameters using Phase IIa audits + the rest of Phase I -------------
## Setup B-splines as a 100 x 8 matrix 
B = splines::bs(x = summ_data$ALI_STAR, 
                df = 8, 
                Boundary.knots = range(summ_data$ALI_STAR), 
                intercept = TRUE, 
                degree = 3)
## Rename columns bs1, ..., bs8
colnames(B) = paste0("bs", seq(1, 8))
## Append the B-splines to the original dataset
summ_data = summ_data |> 
  dplyr::bind_cols(B)

### Fit SMLE -------------------------------------------------------------
library(sleev)
suppressMessages(fit <- logistic2ph(
  Y_unval = NULL,
  Y = "ANY_ENCOUNTERS",
  X_unval = "ALI_STAR",
  X = "ALI",
  Z = "AGE_AT_ENCOUNTER10",
  Bspline = colnames(B),
  data = summ_data,
  hn_scale = 1,
  noSE = FALSE,
  TOL = 1e-04,
  MAX_ITER = 1000
))

## Rescale the coefficient on ALI to be in 0.1-point increments
fit$coefficients[2, 1:2] = 0.1 * fit$coefficients[2, 1:2]

fit$coefficients |> 
  kableExtra::kable(format = "html", digits = 7) |> 
  kableExtra::kable_styling()

data.frame(Method = "SMLE", 
           Coefficient = c("Intercept", "ALI/ALI*", "Age"),
           OR = exp(fit$coefficients$Estimate), 
           LB = exp(fit$coefficients$Estimate - 
                      1.96 * fit$coefficients$SE), 
           UB = exp(fit$coefficients$Estimate + 
                      1.96 * fit$coefficients$SE)) |> 
  kableExtra::kable(format = "html", digits = 7) |> 
  kableExtra::kable_styling()

plot_dat = data.frame(Method = "SMLE", 
                      Coefficient = c("Intercept", "ALI/ALI*", "Age"),
                      OR = exp(fit$coefficients$Estimate), 
                      LB = exp(fit$coefficients$Estimate - 1.96 * fit$coefficients$SE), 
                      UB = exp(fit$coefficients$Estimate + 1.96 * fit$coefficients$SE)) |> 
  dplyr::bind_rows(plot_dat)

# Make the forest plot 
forest = plot_dat |> 
  mutate(Coefficient = factor(x = Coefficient, 
                              levels = c("Intercept", 
                                         "Age",
                                         "ALI/ALI*")),
         Method = factor(x = Method, 
                         levels = c("Naive", "SMLE", "Complete Case"), 
                         labels = c("Naive", "SMLE", "Complete\nCase"))) |> 
  ggplot(aes(x = Method, y = OR, col = Method)) + 
  geom_point(position = position_dodge(width = 0.5), 
             size = 2) + 
  geom_errorbar(aes(ymin = LB, ymax = UB), 
                linewidth = 1) + 
  geom_hline(yintercept = 1, linetype = 2, color = slide_cols[4]) + 
  scale_color_manual(values = slide_cols[c(1, 2, 5)], 
                     name = "Method:") + 
  theme_minimal(base_size = 14) + 
  facet_wrap(~Coefficient) + 
  theme(legend.position = "left", 
        legend.title = element_text(face = "bold"),
        strip.background = element_rect(fill = slide_cols[3]),
        strip.text = element_text(color = "white"),
        axis.text.x = element_blank(), 
        axis.title.x = element_blank()) + 
  ylab("Odds Ratio (95% Confidence Interval)")
forest
ggsave(filename = "~/Documents/ALI_EHR/figures/forest_plot.png", 
       plot = forest, device = "png", width = 8, height = 5, units = "in")