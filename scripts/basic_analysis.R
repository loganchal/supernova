library(here)
library(tidyverse)
library(ggplot2)

#import data
count_data <- readxl::read_excel(here('raw_data', 'zooplankton_count_data.xlsx'))

# Filter data for the two sites of interest, excluding the first subsamples
compare_data <- count_data %>%
  filter(site %in% c("trevor_middle_margin", "trevor_middle_center")) %>%
  filter(subsample_id != "trevor_middle_center_1" & subsample_id != "trevor_middle_margin_1")