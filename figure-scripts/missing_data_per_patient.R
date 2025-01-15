# Load packages 
library(tidyr) ## to unpivot data
library(dplyr) ## for data wrangling
library(ggplot2) ## for pretty plots

# Define colors 
cols = c("#ff99ff", "#8bdddb", "#787ff6", "#ffbd59", "#7dd5f6")

# Load data
## ALI components before and after validation (all waves)
all_data = read.csv("~/Documents/Allostatic_load_audits/all_ali_dat.csv")

## Add a column with the number of missing values per patient 
all_data$NUM_MISSING = apply(X = all_data[, -c(1:7, 18)], 
                             MARGIN = 1, 
                             FUN = function(x) sum(is.na(x)))

# Create bar plot of missing values per patient (colored by wave)
all_data |> 
  ggplot(aes(x = NUM_MISSING, 
             fill = DATA)) + 
  geom_bar(stat = "count", 
           position = "dodge", 
           color = "black") + 
  geom_text(stat='count', 
            aes(label=..count..), 
            vjust = -1, 
            size = 3, 
            position = position_dodge(width = 1)) + 
  theme_minimal(base_size = 12) + 
  labs(x = "Missing ALI Components", 
       y = "Number of Patients", 
       title = "B) Missingness in ALI Components by Patient") +
  theme(title = element_text(face = "bold"), 
        legend.position = "inside", 
        legend.position.inside = c(1, 0.9),
        legend.text = element_text(size = 12), 
        legend.title = element_text(size = 12, face = "bold"), 
        legend.justification = "right", 
        legend.background = element_rect(fill = "white")) + 
  scale_fill_manual(values = cols, name = "Data:") + 
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) + 
  scale_x_continuous(breaks = 0:10)

## Save it 
ggsave(filename = "~/Documents/ALI_EHR/figures/missing_data_per_patient.png", 
       device = "png", width = 10, height = 6, units = "in")
