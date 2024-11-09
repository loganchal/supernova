library(here)
library(tidyverse)
library(ggplot2)

#import data
count_data <- readxl::read_excel(here('raw_data', 'zooplankton_count_data.xlsx'))

# Filter data for the two sites of interest, excluding the first subsamples
compare_data <- count_data %>%
  filter(site %in% c("trevor_middle_margin", "trevor_middle_center")) %>%
  filter(subsample_id != "trevor_middle_center_1" & subsample_id != "trevor_middle_margin_1")



#Amphipod Testing

#normality test
shapiro_amphipod_center <- shapiro.test(compare_data$amphipod_count[compare_data$site == "trevor_middle_center"])
shapiro_amphipod_margin <- shapiro.test(compare_data$amphipod_count[compare_data$site == "trevor_middle_margin"])

#homogeneity test 
bartlett_amphipod <- bartlett.test(amphipod_count ~ site, data = compare_data)

#results
shapiro_amphipod_margin
shapiro_amphipod_center
bartlett_amphipod

# Run t-tests for amphipod and copepod counts between the two sites
amphipod_t_test <- t.test(compare_data$amphipod_count ~ compare_data$site)
copepod_t_test <- t.test(compare_data$copepod_count ~ compare_data$site)

# Display t-test results
amphipod_t_test
copepod_t_test








# Combine data for plotting by species
compare_data_long <- compare_data %>%
  select(site, amphipod_count, copepod_count) %>%
  pivot_longer(cols = ends_with("count"), names_to = "species", values_to = "count")

# Plot counts for both amphipods and copepods between the two sites
ggplot(compare_data_long, aes(x = site, y = count, fill = species)) +
  geom_boxplot(alpha = 0.7) +  # Boxplot for distribution
  geom_jitter(size = 2, color = "black", alpha = 0.5, position = position_dodge(0.8)) +  # Individual points without width
  labs(title = "Comparison of Amphipod and Copepod Counts Between Sites (Excluding First Subsamples)",
       x = "Site",
       y = "Count",
       fill = "Species") +
  theme_minimal()