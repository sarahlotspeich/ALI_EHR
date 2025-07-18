# Load packages 
library(dplyr) ## for data wrangling
library(ggplot2) ## for pretty plots

# Define colors 
cols = c("#ff99ff", "#8bdddb", "#787ff6", "#ffbd59", "#7dd5f6")

# Load data
## ALI components before and after validation (all waves)
all_data = read.csv("~/Documents/Allostatic_load_audits/all_ali_dat.csv")
### Create vector of IDS for validated patients
val_ids = unique(all_data$PAT_MRN_ID[all_data$VALIDATED])
### Subset to validated patients 
all_data = all_data |> 
  filter(PAT_MRN_ID %in% val_ids)

# Create scatterplot 
all_data |> 
  group_by(PAT_MRN_ID) |> 
  summarize(UNVAL_ALI = min(if_else(condition = !VALIDATED, 
                                    true = ALI, 
                                    false = 9999)), 
            VAL_ALI = min(if_else(condition = VALIDATED, 
                                  true = ALI, 
                                  false = 9999))) |>
  ggplot(aes(x = UNVAL_ALI, 
             y = VAL_ALI)) + 
  geom_abline(slope = 1, 
              intercept = 0, 
              linetype = 2, 
              size = 1, 
              color = "gray") + 
  geom_point(color = cols[4], size = 3) + 
  theme_minimal(base_size = 12) + 
  labs(x = "Unvalidated ALI (from EHR)", 
       y = "Validated ALI (from Chart Review)") +
  theme(title = element_text(face = "bold")) + 
  scale_x_continuous(expand = expansion(mult = c(0.05, 0.05))) + 
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.05))) + 
  coord_equal()

## Save it 
ggsave(filename = "~/Documents/ALI_EHR/figures/scatterplot_unval_vs_val_ali.png", 
       device = "png", width = 6, height = 6, units = "in")