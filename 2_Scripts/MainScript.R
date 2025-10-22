rm(list = ls()) # clean environment

# Load required libraries
library(quantreg)
library(readxl)
library(stargazer)
library(tidyverse)
library(vars)

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
w_data <- w_data[w_data[[2]] == "Bread_Wheat", ] # alt. 'Feed_Wheat'

dat <- inner_join(w_data, f_data, by = "Date") # keep only 'common date'-datapoints

dat <- dat %>%
  rename(w_p = Avg_Price.x, f_p = Avg_Price.y) %>%
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
  ) %>% filter(!is.na(w_p4)) # remove NAs, only 2 lags are used -> over filtered?

# VAR select
vardat <- dat %>% dplyr::select(w_p, f_p)

a <- VARselect(vardat, lag.max = 10, type = "const")
a$selection #SC (BIC): 2 lags
# summary(nvar <- VAR(vardat, p=2))

# add exo variables
VARselect(vardat, lag.max =4, exogen = cbind(dat$w_ps1))
VARselect(vardat, lag.max =4, exogen = cbind(dat$f_ps1))
VARselect(vardat, lag.max =4, exogen = cbind(dat$w_ps1, dat$f_ps1)) #SC: 2

#marginal spec
mw = w_p ~ w_p1 + f_p1 + w_p2 + f_p2 #wheat
mf = f_p ~ w_p1 + f_p1 + w_p2 + f_p2 #fert

#psuedo R2s
taus_sparse <- c(.1, .3, .5, .7, .9)
psu_r2s <- data.frame()
rho <- function(u, tau) sum(u * (tau - (u < 0))) #rho calc from JP code

for (tau in taus_sparse) {
  fit <- rq(mw, tau = tau, data = dat) #fit the q reg for current tau
  
  y <- dat$w_p  #actual values
  yhat <- predict(fit) #predicted by model
  
  # compute rho and rho0
  rho_m <- rho(y - yhat, tau)
  rho0  <- rho(y - quantile(y, probs = tau), tau)
  
  
  r2 <- 1 - rho_m / rho0 #pseudo-R2
  
  psu_r2s <- rbind(psu_r2s, data.frame(tau = tau, pseudoR2 = r2)) #store
}
psu_r2s

# QVARs
taus_dense <- seq(0.01, 0.99, by = 0.01)
fit_w <- rq(mw, tau = taus_dense, data = dat)
fit_f <- rq(mf, tau = taus_dense, data = dat)
coeff_w <- fit_w$coeff
coeff_f <- fit_f$coeff

# Independent variables
w_p <- dat$w_p
w_p1 <- dat$w_p1
w_p2 <- dat$w_p2

f_p <- dat$f_p
f_p1 <- dat$f_p1
f_p2 <- dat$f_p2

X <- cbind(1, w_p1, f_p1, w_p2, f_p2) #col of 1s for b0

# Inverse distributions
y_w <- X %*% coeff_w
y_w <- taus_dense(apply(y_w, 1, cummax))

y_f <- X %*% coeff_f
y_f <- taus_dense(apply(y_f, 1, cummax))

# Copula
n = nrow(vardat)

F_w <- array(0, dim=c(n,1))
F_f <- array(0, dim=c(n,1))

for (ii in 1:n) {
  F_w[ii] <- taus_dense[min(which(min(abs(w_p[ii] - y_w[ii,])) == abs(w_p[ii] - y_w[ii,])))]
  F_f[ii] <- taus_dense[min(which(min(abs(f_p[ii] - y_f[ii,])) == abs(f_p[ii] - y_f[ii,])))]
}

C <- F_w ~ F_f
summary(lm(C))
rC <- rq(C, tau = taus_sparse)
summary(rC, se = "boot", brmethod = "xy")
rC_full <- rq(C, tau = taus_dense)
coeff_rC <- rC_full$coeff

plot(taus_dense, coeff_rC[2,], type="l", col="steelblue", lwd=2,
     xlab="tau of fertilizer price",
     ylab="slope")
abline(h=0, lty=2)

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