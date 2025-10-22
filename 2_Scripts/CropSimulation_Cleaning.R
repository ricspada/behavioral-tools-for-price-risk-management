rm(list = ls())

library(tidyverse)
library(readxl)
library(robustbase)

##############################################################################
# EPIC Crop Simulation Data
##############################################################################

# SimUID
# -----------------------------
# SimUID	simulation grid ID


# CROP (crop abbreviations)
# -----------------------------
# CORN	corn
# CSIL	maize silage
# FPEA	field peas
# OATS	oats
# POTA	potatoes
# RAPE	winter rape
# SBAR	barley
# SGBT	sugar beet
# SOYB	soybeans
# SUNF	sunflower
# WRYE	winter rye
# WWHT	winter wheat

# ROT (cropping system)
# -----------------------------
# CRS1	crop rotation #1
# CRS2	crop rotation #2
# CRS3	crop rotation #3
# CRS4	crop rotation #4
# CRS5	crop rotation #5
# CRS6	crop rotation #6
# CRS7	crop rotation #7
# MONO	nonocrop

# FER (N-fertilization scenario)
# -----------------------------
# BAU	business as usual
# N10	max 10 kgN/ha/yr
# N50	max 50 kgN/ha/yr
# N100	max 100 kgN/ha/yr
# N250	max 250 kgN/ha/yr


# IRR (irrigation scenario)
# -----------------------------
# ir	irrigated
# rf	rainfed


# WTH (meteorological forcing)
# -----------------------------
# cmip6_ipsl_ssp126	climate change (CMIP6, IPSL, ssp 126)
# cmip6_ipsl_ssp585	climate change (CMIP6, IPSL, ssp 585)
# cmip6_mpi_ssp126	climate change (CMIP6, MPI, ssp 126)
# cmip6_mpi_ssp585	climate change (CMIP6, MPI, ssp 585)
# hist			historical data


# TIL (crop residue management scenario)
# -----------------------------
# contill_bau	conventional tillage, business-as-usual crop residue harvest
# contill_r00	conventional tillage, crop residue harvest 100%
# contill_r30	conventional tillage, crop residue harvest 70%
# contill_r60	conventional tillage, crop residue harvest 40%
# contill_r90	conventional tillage, crop residue harvest 10%
# mintill_cons	minimum tillage, low crop residue harvest + mulching


# YLD_DM* (crop yield)
# -----------------------------
# YLD_DM	dry-matter crop yield in tDM/ha

# *the following water contents can be used to convert from dry matter to economic yield:
# CORN	15%
# CSIL	70%
# FPEA	12%
# OATS	10%
# POTA	80%
# RAPE	8%
# SBAR	12%
# SGBT	84%
# SOYB	13%
# SUNF	6%
# WRYE	12%
# WWHT	12%

# FTN (N fertilization)
# -----------------------------
# FTN	N application rate in kgN/ha


sim_data <- read_csv("1_Data/RawData/EPIC2/Riccardo_EPIC.csv") 

df <- sim_data |> 
  filter(CROP == "WWHT", 
#         SimUID == 52338, # select only one grid location
         IRR == "rf" # rainfed only
        ) 

ggplot(df, aes(x = FTN, y = YLD_DM)) +
  geom_point(alpha = 0.5)




reg_yield_function <- lmrob(YLD_DM ~ sqrt(FTN) + FTN, data = df, method = "MM")
summary(reg_yield_function)

# append residuals to dataframe
df <- df %>% 
  mutate(residuals = reg_yield_function$residuals,
         abs_residuals = abs(residuals))

reg_variation_function <- lmrob(abs_residuals ~ I(sqrt(FTN)), 
                                  data = df, method = "MM")
summary(reg_variation_function)

# plot yield function with data points
ggplot(df, aes(x = FTN, y = YLD_DM)) +
  geom_point(alpha = 0.5) +
  stat_smooth(method = "lm", 
              formula = y ~ sqrt(x) + x, 
              color = "blue", se = FALSE) +
  labs(title = "Yield Function",
       x = "N Fertilization (kgN/ha)",
       y = "Yield (tDM/ha)") +
  theme_minimal()
