rm(list = ls())

library(tidyverse)
library(readxl)
library(robustbase)

# import data 
sim_data <- read_csv("1_Data/RawData/EPIC2/Riccardo_EPIC.csv") 

# filter to winter wheat
df_sim_data <- sim_data |> 
  filter(CROP == "WWHT", 
         IRR == "rf",
         WTH == "hist"             
        ) |> 
  # convert dry matter yield to economic yield (12% moisture content)
  mutate(YLD = YLD_DM / 0.88) |>
  mutate(FTN = round(FTN)) 

# summary statistics 
df_sim_data |> 
  group_by(FTN) |> 
  summarise(
    mean_YLD = mean(YLD),
    median_YLD = median(YLD),
    sd_YLD = sd(YLD),
    min_YLD = min(YLD),
    max_YLD = max(YLD),
    n = n()
  ) |> 
  mutate(across(where(is.numeric), ~ round(., 2)))

# difference between 2 grid locations 
df_sim_data |> 
  group_by(SimUID) |> 
  summarise(
    mean_FTN = mean(FTN),
    sd_FTN = sd(FTN),
    min_FTN = min(FTN),
    max_FTN = max(FTN),
    mean_YLD = mean(YLD),
    median_YLD = median(YLD),
    sd_YLD = sd(YLD),
    min_YLD = min(YLD),
    max_YLD = max(YLD),
    n = n()
  ) |> 
  mutate(across(where(is.numeric), ~ round(., 2)))

# create separate datasets for each grid location
grid1_data <- df_sim_data |> 
  filter(SimUID == 52338)

grid2_data <- df_sim_data |> 
  filter(SimUID == 52339)

# write_csv(grid1_data, "1_Data/Grid1_YieldData.csv")
# write_csv(grid2_data, "1_Data/Grid2_YieldData.csv")

### plot ###########################################################################

ggplot(grid1_data, aes(x = FTN, y = YLD)) +
  geom_point(alpha = 0.5)


# plot yield function with data points
ggplot(grid1_data, aes(x = FTN, y = YLD)) +
  geom_point(alpha = 0.5) +
  stat_smooth(method = "lm", 
              formula = y ~ I(sqrt(x)) + x, 
              color = "blue", se = FALSE) +
  stat_smooth(method = "lm", 
              formula = y ~ x, 
              color = "red", se = FALSE) +
  labs(x = "N Fertilization (kgN/ha)",
       y = "Yield (t/ha)") +
  theme_minimal()

# plot yield level by fertilization treatment
box_plot <- ggplot(grid1_data, aes(x = as.factor(FTN), y = YLD)) +
  geom_boxplot() +
  labs(x = "Fertilization Treatment (kgN/ha)",
       y = "Yield (t/ha)") +
  theme_minimal() +
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14))

box_plot


# plot yield level by year
box_plot1 <- ggplot(grid1_data, aes(x = as.factor(YR), y = YLD)) +
  geom_boxplot() +
  labs(x = NULL,
       y = "Yield (t/ha)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.y = element_text(size = 14))

box_plot1
ggsave("3_Outputs/Fig_Yield_by_Year.png", box_plot1, width = 8, height = 6)


