rm(list = ls())

library(tidyverse)
library(readxl)

##############################################################################
# EPIC Crop Simulation Data
##############################################################################

sim_data <- read_csv("1_Data/RawData/EPIC2/Riccardo_EPIC.csv") 

df <- sim_data |> 
  filter(CROP == "WWHT", 
         SimUID == 52338) |> # select only one grid location  
  select(-SimUID)
  
ggplot(df, aes(x = FTN, y = YLD_DM)) +
  geom_point(alpha = 0.5)




reg_yield_function <- lmrob(YLD_DM ~ I(sqrt(FTN)) + FTN, data = df, method = "MM")
summary(reg_yield_function)

# append residuals to dataframe
df <- df %>% 
  mutate(residuals = reg_yield_function$residuals,
         abs_residuals = abs(residuals))

reg_variation_function <- lmrob(abs_residuals ~ I(sqrt(FTN)), 
                                  data = df, method = "MM")
summary(reg_variation_function)
