rm(list = ls()) # clean environment

# Load required libraries
library(tidyverse)
library(readxl)
library(stargazer)
library(zoo)

# Set seed for reproducibility
set.seed(1234)

################################################################################
# 1. NITROGEN-YIELD RELATIONSHIP
################################################################################

# 1.1 Get data from crop simulation model 
sim_data <- read_csv("1_Data/RawData/EPIC2/Riccardo_EPIC.csv") 

df <- sim_data |> 
  filter(CROP == "WWHT", 
         SimUID == 52338) # select only one grid location  
  
# 1.2 Just & Pope production function estimation

reg_yield_function <- lmrob(YLD_DM ~ I(sqrt(FTN)) + FTN, data = df, method = "MM")
summary(reg_yield_function)

# append absolute value of residuals to dataframe
df <- df |> 
  mutate(abs_residuals = abs(reg_yield_function$residuals))

reg_variation_function <- lmrob(abs_residuals ~ I(sqrt(FTN)), 
                                  data = df, method = "MM")
summary(reg_variation_function)

# summary table of regression results 
stargazer(reg_yield_function, reg_variation_function,
          type = "text",
          title = "Yield and Variation Function Estimates",
          dep.var.labels = c("Yield Function", "Variation Function"),
          covariate.labels = c("sqrt(N)", "N", "Intercept"),
          omit.stat = c("f", "ser"),
          no.space = TRUE)

################################################################################
# 2. NITROGEN PRICES - WHEAT PRICES RELATIONSHIP
################################################################################

### Abbreviations and explanations ##############################################

### Prices refer to the Cologne Lowland

# Fertilizers (prices are in €/100 kg)

# CAN = Calcium Ammonium Nitrate (CAN), 27% N
# UAN = Urea Ammonium Nitrate (UAN), 28% N
# ASN = Ammonium Sulphate Nitrate (ASN),26% N + 13% S
# AS  = Ammonium Sulphate (AS),21% N + 24% S
# urea_P = Urea 46% N, granulated, protected (Urea_P), 46% N
# urea   = Urea 46% N, granulated, (Urea), 46% N
# DAP = Diammonium phosphate, (DAP) 18% N + 46% P2O5
# TSP = Triple superphosphate (TSP) 45 % P2O5
# Kornkali = Kornkali + Mg, 38 % K2O + 6 % MgO

# Wheat (prices are in €/100 kg)

# Bread_Wheat = Bread wheat (higher quality)
# Feed_Wheat  = Feed wheat  (lower quality)


### 2.1 Import Data ###############################################################

### Fertilizer prices 
f_data <- read_csv("1_Data/FertilizerPrices.csv")

ggplot(f_data, aes(x = Date, color = Product, fill = Product)) +
  geom_ribbon(aes(ymin = Min_Price, ymax = Max_Price), alpha = 0.2, color = NA) +
  geom_line(aes(y = Avg_Price), linewidth = 0.8) +
  facet_wrap(~ Product, scales = "fixed", ncol = 3) +
  labs(y = "Price (€/100 kg)",
       x = "Date") +
  theme_minimal() +
  theme(legend.position = "none")

### Wheat prices 
w_data <- read_csv("1_Data/Wheatprices.csv")

ggplot(w_data, aes(x = Date, color = Product, fill = Product)) +
  geom_ribbon(aes(ymin = Min_Price, ymax = Max_Price), alpha = 0.2, color = NA) +
  geom_line(aes(y = Avg_Price), linewidth = 0.8) +
  scale_color_brewer(palette = "Dark2") +
  scale_fill_brewer(palette = "Set2") +
  labs(y = "Price (€/100 kg)",
       x = "Date") +
  theme_minimal() +
  theme(legend.position = "bottom")

### 2.2 QVAR and copula
f_data <- f_data[f_data[[2]] == "CAN", ] # Keep only CAN data
w_data <- w_data[w_data[[2]] == "Bread_Wheat", ]

# keep only common date datapoints
dat <- inner_join(w_data, f_data, by = "Date")

# Create variables
nn <- nrow(dat)-1 #number of observations
nL <- 4 # lags
n <- nn-nL # trim for lags

w_p <- dat[ ,5]
w_p <- na.approx(dat[ ,5]) # Linear approximation of missing values (obs. 671 672 673 675 747)
w_p1 <- dat[(nL-1):(nn-1), 5] #wheat lags
w_p2 <- dat[(nL-2):(nn-2), 5]
w_p3 <- dat[(nL-3):(nn-3), 5]
w_p4 <- dat[(nL-4):(nn-4), 5]

f_p <- dat[ ,9]
f_p1 <- dat[(nL-1):(nn-1), 9] #fertilizer lags
f_p2 <- dat[(nL-2):(nn-2), 9]
f_p3 <- dat[(nL-3):(nn-3), 9]
f_p4 <- dat[(nL-4):(nn-4), 9]

yr <- dat[nL:nn, 1]

# Combine data
vardata <- cbind(w_p, f_p)
View(vardata)

################################################################################
# 3. ADJUSTED CONTRACT PRICING
################################################################################

# Burn rate 



################################################################################
# 4. MONTE CARLO SIMULATION OF UTILITY FUNCTIONS
################################################################################

### 4.1 Broad framing ###########################################################

# 4.1.1 Expected Utility Theory (EUT) 

# 4.1.2 Cumulative Prospect Theory (CPT) 

# 4.1.3 Ambiguity Aversion (Alpha model) 


### 4.2 Narrow Framing ####################################################

# 4.2.1 Expected Utility Theory (EUT) 

# 4.2.2 Cumulative Prospect Theory (CPT) 





################################################################################
# END OF SCRIPT
################################################################################