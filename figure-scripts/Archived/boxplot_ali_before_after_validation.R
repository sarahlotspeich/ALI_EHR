# Load packages 
library(tidyr) ## to unpivot data
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

## Create columns for unvalidated/validated ALI
plot_data = all_data |> 
  group_by(PAT_MRN_ID) |> 
  summarize(UNVAL_ALI = min(if_else(condition = !VALIDATED, 
                                    true = ALI, 
                                    false = 9999)), 
            VAL_ALI = min(if_else(condition = VALIDATED, 
                                  true = ALI, 
                                  false = 9999))) |>
  gather("DATA", "ALI", -1) |>
  mutate(DATA = factor(x = DATA, 
                       levels = c("UNVAL_ALI", "VAL_ALI"), 
                       labels = c("Unvalidated (from EHR)", "Validated (from Chart Review)"))) 

# Save medians unvalidated/validated
meds = plot_data |> 
  group_by(DATA) |> 
  summarize(MEDIAN = median(ALI))

# Create histogram of allostatic load index
plot_data |> 
  ggplot(aes(x = DATA, 
             y = ALI, 
             fill = DATA)) + 
  geom_boxplot(position = "dodge", 
               color = "black") + 
  geom_text(data = meds, 
            aes(x = DATA, y = MEDIAN, 
                label = formatC(round(MEDIAN, 2), 
                                format = 'f', flag='0', digits = 2)), 
            size = 5, hjust = 0.5, vjust = 5.5) + 
  theme_minimal(base_size = 12) + 
  labs(x = "Data", 
       y = "Allostatic Load Index (ALI)") +
  theme(title = element_text(face = "bold"), 
        legend.position = "inside", 
        legend.position.inside = c(1, 0.9),
        legend.text = element_text(size = 12), 
        legend.title = element_text(size = 12, face = "bold"), 
        legend.justification = "right", 
        legend.background = element_rect(fill = "white")) + 
  scale_fill_manual(values = cols, guide = "none") + 
  scale_x_discrete(labels = function(x) stringr::str_wrap(x, width = 12)) + 
  coord_flip()

## Save it 
ggsave(filename = "~/Documents/ALI_EHR/figures/boxplot_ali_before_after_validation.png", 
       device = "png", width = 8, height = 4, units = "in")