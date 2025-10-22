rm(list = ls()) # clean environment

# Load required libraries
library(quantreg)
library(readxl)
library(stargazer)
library(robustbase)
library(tidyverse)
library(vars)

# Set seed for reproducibility
set.seed(1234)

################################################################################
# 1. NITROGEN-YIELD RELATIONSHIP
################################################################################

# 1.1 Get data from crop simulation model 
  
# 1.2 Just & Pope production function estimation


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
w_data <- w_data[w_data[[2]] == "Bread_Wheat", ] # alt. 'Feed_Wheat'

dat <- inner_join(w_data, f_data, by = "Date") # keep only 'common date'-datapoints

dat <- dat %>%
  rename(w_p = 5, f_p = 9) %>%
  mutate(
    w_p1 = lag(w_p, 1), # wheat lags
    w_p2 = lag(w_p, 2),
    w_p3 = lag(w_p, 3),
    w_p4 = lag(w_p, 4),
    f_p1 = lag(f_p, 1), # fert lags
    f_p2 = lag(f_p, 2),
    f_p3 = lag(f_p, 3),
    f_p4 = lag(f_p, 4),
    w_ps1 = w_p1^2,     # squares
    f_ps1 = f_p1^2
  ) %>% filter(!is.na(w_p4)) # remove NAs

# VAR select
vardat <- dat %>% dplyr::select(w_p, f_p)
colnames(vardat) <- c("Wheat price", "Fertilizer price")

a <- VARselect(vardat, lag.max = 10, type = "const")
a$selection #SC (BIC): 2 lags
# summary(nvar <- VAR(vardat, p=2))

# add exo variables
VARselect(vardat, lag.max =4, exogen = cbind(dat$w_ps1))
VARselect(vardat, lag.max =4, exogen = cbind(dat$f_ps1))
VARselect(vardat, lag.max =4, exogen = cbind(dat$w_ps1, dat$w_ps1)) #SC: 2

#marginal spec
mw = w_p ~ w_p1 + f_p1 + w_p2 + f_p2 #wheat
mf = f_p ~ w_p1 + f_p1 + w_p2 + f_p2 #fert

#set quantiles
taus <- c(.1, .3, .5, .7, .9)

psu_r2s <- data.frame()
rho <- function(u, tau) sum(u * (tau - (u < 0))) #rho calc from JP code

for (tau in taus) {
  fit <- rq(mw, tau = tau, data = dat) #fit the q reg for current tau
  
  y <- dat$w_p  #actual values
  yhat <- predict(fit) #predicted by model
  
  # compute rho and rho0
  rho <- rho(y - yhat, tau)
  rho0  <- rho(y - quantile(y, probs = tau), tau)
  
  
  r2 <- 1 - rho / rho0 #pseudo-R2
  
  psu_r2s <- rbind(psu_r2s, data.frame(tau = tau, pseudoR2 = r2)) #store
}

psu_r2s

# QVARs
rq(mw, tau = taus, data = dat)
rq(mf, tau = taus, data = dat)
################################################################################
# 3. ADJUSTED CONTRACT PRICING
################################################################################

# Burn rate 



################################################################################
# 4. MONTE CARLO SIMULATION OF UTILITY OF PROFITS FUNCTIONS
################################################################################

### 4.1 Broad framing ###########################################################

# 4.1.1 Expected Utility Theory (EUT) 

# 4.1.2 Cumulative Prospect Theory (CPT) 

# 4.1.3 Ambiguity Aversion (Alpha model) 


### 4.2 Narrow Framing ####################################################

# 4.2.1 Expected Utility Theory (EUT) 

# 4.2.2 Cumulative Prospect Theory (CPT) 







