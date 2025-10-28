rm(list = ls()) # clean environment

# Load required libraries
library(quantreg)
library(readxl)
library(stargazer)
library(robustbase)
library(tidyverse)
library(vars)
library(tseries)
library(patchwork)

# Set seed for reproducibility
set.seed(1234)

################################################################################
# 1. NITROGEN-YIELD RELATIONSHIP
################################################################################

### 1.1 Get data from crop simulation model 
sim_yield <- read_csv("1_Data/Grid1_YieldData.csv") 

# summary statistics 
summary <- sim_yield |> 
  group_by(FTN) |> 
  summarise(
    mean_YLD = mean(YLD),
    sd_YLD = sd(YLD),
    min_YLD = min(YLD),
    max_YLD = max(YLD),
    n = n()
  ) |> 
  mutate(across(where(is.numeric), ~ round(., 2)))

summary
# write_csv(summary, "3_Outputs/Table_YieldSummary.csv")

### 1.2 Just & Pope production function estimation
yield_function <- lmrob(YLD_DM ~ sqrt(FTN) + FTN, data = sim_yield, method = "MM")

# append residuals
sim_yield <- sim_yield |> 
  mutate(residuals = yield_function$residuals,
         abs_residuals = abs(residuals))

variation_function <- lmrob(abs_residuals ~ sqrt(FTN), data = sim_yield, method = "MM")

# table of results
stargazer(yield_function, variation_function, 
          type = "text",
          dep.var.labels = c("Yield function", "Variation function"),
          covariate.labels = c("sqrt(N)", "N", "Intercept")
         )


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


### 2.2 QVAR and copula #############################################################

f_data <- f_data[f_data[[2]] == "CAN", ] # Keep only CAN data
w_data <- w_data[w_data[[2]] == "Bread_Wheat", ] # alt. 'Feed_Wheat'

# Check time ranges
range(f_data$Date)
range(w_data$Date)

# Check for matching times
common_times <- intersect(f_data$Date, w_data$Date)
length(common_times)

# Check time intervals (see if regularly spaced)
diff(f_data$Date) 
diff(w_data$Date)  

# Merge datasets on common dates
dat <- inner_join(w_data, f_data, by = "Date") # keep only 'common date'-datapoints
dat <- dat |> 
  rename(w_p = Avg_Price.x, f_p = Avg_Price.y) |> 
  dplyr::select(Date, w_p, f_p) 

range(dat$Date)
diff(dat$Date) 

# check for time intervals that are not 14 days
sum(diff(dat$Date) != 14)

# we dont have regular 14-day intervals, data is recorder every 1st and 3rd monday of the month (sometimes 3 week gaps are introduced)
# do we need to ensure regular time intervals?

# plot prices
price_plot <- ggplot(dat, aes(x = Date)) +  
  geom_line(aes(y = f_p,  linetype = "Calcium Ammonium Nitrate price"), linewidth = 0.5) +
  geom_line(aes(y = w_p,  linetype = "Wheat price"), linewidth = 0.5) +
  scale_linetype_manual(values = c("Calcium Ammonium Nitrate price" = "dashed", "Wheat price" = "solid")) +
  scale_y_continuous(limits = c(0, NA), expand = c(0, 0)) +
  labs(y = "Price (€/100 kg)",
       x = "Date",
       linetype = NULL) +
  theme_classic() +
  theme(legend.position = c(0.1, 0.9),
        legend.justification = c(0, 1),
        legend.background = element_rect(fill = "white", color = "black"))

price_plot
#ggsave("3_Outputs/Fig_Prices.pdf", plot = price_plot, width = 6, height = 4)


# check for stationarity (p-value < 0.01 -> stationary)
adf.test(dat$w_p, k=2) 
adf.test(dat$f_p, k=2) 

pp.test(dat$w_p)
pp.test(dat$f_p)

# time series are non-stationary: take first differences
dif_dat <- dat |> 
  mutate(
    w_p = c(NA, diff(w_p)),
    f_p = c(NA, diff(f_p))
  ) |> 
  drop_na()

# check for stationarity again
adf.test(dif_dat$w_p, k=2)
adf.test(dif_dat$f_p, k=2)

pp.test(dif_dat$w_p)
pp.test(dif_dat$f_p)


# plot differenced prices
plot_dif_w <- ggplot(dif_dat, aes(x = Date)) +
  geom_line(aes(y = w_p), linewidth = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(y = "Wheat", x = NULL) +
  theme_classic()

plot_dif_f <- ggplot(dif_dat, aes(x = Date)) +
  geom_line(aes(y = f_p), linewidth = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(y = "Fertilizer", x = NULL) +
  theme_classic()

# Stack vertically
combined_dif_plots <- plot_dif_f / plot_dif_w
combined_dif_plots
# ggsave("3_Outputs/Fig_DifferencedPrices.pdf", plot = combined_dif_plots, width = 8, height = 6)


# Create lags and squares
dat <- dif_dat %>%
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
  ) |> 
  drop_na()
  
# VAR select
vardat <- dat %>% dplyr::select(w_p, f_p)

# VAR Specification
# select var order based on SBIC
a <- VARselect(vardat, lag.max = 10, type = "const")
a$selection #SC (BIC): 1 lag
summary(nvar <- VAR(vardat, p=1))
summary(a)
a

# add exogenous variables
VARselect(vardat, lag.max =4, exogen = cbind(dat$w_ps1))
VARselect(vardat, lag.max =4, exogen = cbind(dat$f_ps1))
VARselect(vardat, lag.max =4, exogen = cbind(dat$w_ps1, dat$f_ps1)) #SC: 2

# Specify and Estimate QVAR

# define individual marginal specifications
mw <- w_p ~ w_p1 + f_p1  #wheat
mf <- f_p ~ w_p1 + f_p1  #fert


# Calculate pseudo R^2s
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

f_p <- dat$f_p
f_p1 <- dat$f_p1

X <- cbind(1, w_p1, f_p1) #col of 1s for b0

# Inverse distributions
y_w <- X %*% coeff_w
y_w <- t(apply(y_w, 1, cummax))

y_f <- X %*% coeff_f
y_f <- t(apply(y_f, 1, cummax))

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

### Price path simulation ###########################################################

Nt <- 180   # number of periods to simulate (months)
Ns <- 1000 # number of draws
Ys <- array(0, dim=c(2, Nt, Ns))

# Set initial conditions 
Ys[1, ,] <- dat$f_p[n] # get last price of fertilizer 
Ys[2, ,] <- dat$w_p[n] # get last price of wheat

# QVAR Simulation
for (it in (3:Nt)) {
  for (is in (1:Ns)) {
    
    f_pi1 <- Ys[1, it-1, is]
    w_pi1 <- Ys[2, it-1, is]
    
    Xfi  <- cbind(1, f_pi1, w_pi1)
    Xwi  <- cbind(1, f_pi1, w_pi1)
    
    qi <- sample(1:99, 1) 

    yf <- Xfi %*% coeff_f[,qi]
    
    yw <- Xwi %*% coeff_w[,qi]
    
    Ys[1, it, is] <- yf
    Ys[2, it, is] <- yw  
    
  }}

# plot simulated paths
ggplot() +
  geom_line(aes(x = 1:Nt, y = Ys[1, ,1]), color = "blue", alpha = 0.5) +
  geom_line(aes(x = 1:Nt, y = Ys[1, ,2]), color = "blue", alpha = 0.5) +
  geom_line(aes(x = 1:Nt, y = Ys[1, ,3]), color = "blue", alpha = 0.5) +
  labs(y = "Fertilizer Price Change (€/100 kg)", x = "Time Period") +
  theme_minimal()

# Reshape data for ggplot
plot_data <- data.frame(
  time = rep(1:Nt, Ns),
  price = as.vector(Ys[1, , ]),
  simulation = rep(1:Ns, each = Nt)
)

ggplot(plot_data, aes(x = time, y = price, group = simulation)) +
  geom_line(color = "blue", alpha = 0.5) +
  labs(y = "Fertilizer Price Change (€/100 kg)", x = "Time Period") +
  theme_minimal()

### Insurance Contract Price ###########################################################

nitrogen_prices <- Ys[1, 1:Nt, ]  # Extract nitrogen price paths

duration <- 12    # contract duration in months
k1 <- 50          # strike price
amount <- 100     # Notional amount (e.g., tons of nitrogen)

payoffs_call <- array(0, dim=c(Nt, Ns))
for (it in 1:Nt) {
  payoffs_call[it, ] <- pmax(nitrogen_prices[it, ] - k1, 0) * amount
}

# Calculate Present Values 
int <- 0.03  # Interest rate
disc <- (1 + int)**(-(1:Nt/12))  # Monthly discount factors

# Discounted payoffs 
disc_cal_payoffs <- array(0, dim=c(Nt, Ns))
for (it in 1:Nt) {
  disc_cal_payoffs[it, ] <- payoffs_call[it, ] * disc[it]
}

# Average Discounted Payoff
contract_price <- mean(disc_cal_payoffs)

contract_price



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







