# Load packages 
library(dplyr) ## for data wrangling
library(ggplot2) ## for pretty plots

# Define colors 
cols = c("#ff99ff", "#8bdddb", "#787ff6", "#ffbd59", "#7dd5f6")

# Load data
## ALI components before 
unval_data = read.csv("~/Documents/Allostatic_load_audits/all_ali_dat.csv") |> 
  filter(DATA == "EHR (Before Validation)")

# Inspect median and interquartile range
unval_data |> 
  pull(ALI) |> 
  summary()
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.0000  0.1667  0.3333  0.3407  0.5000  1.0000 

# Number of unique values 
unval_data |> 
  pull(ALI) |> 
  unique() |> 
  length()

# Create histogram of allostatic load index
unval_data |> 
  ggplot(aes(x = ALI)) + 
  geom_histogram(fill = cols[3], color = "black", bins = 15) + 
  geom_vline(xintercept = 0.333, linetype = 2, color = "black") + 
  theme_minimal(base_size = 12) + 
  labs(x = "Unvalidated Allostatic Load Index (from the EHR)", 
       y = "Number of Patients") +
  theme(title = element_text(face = "bold")) 

## Save it 
ggsave(filename = "~/Documents/ALI_EHR/figures/histogram_ali_before_validation.png", 
       device = "png", width = 8, height = 5, units = "in")