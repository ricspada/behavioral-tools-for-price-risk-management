rm(list = ls()) # clean environment

# Load required libraries
library(quantreg)
library(readxl)
library(modelsummary)
library(robustbase)
library(tidyverse)
library(vars)
library(tseries)
library(patchwork)
library(urca)
library(moments)

# Set seed for reproducibility
set.seed(1234)

################################################################################
# 1. NITROGEN-YIELD RELATIONSHIP
################################################################################

### 1.1 Get data from crop simulation model
sim_yield <- read_csv("1_Data/RawData/EPIC2/Riccardo_EPIC.csv") 

sim_yield <- sim_yield |> 
  filter(
    CROP == "WWHT",  # winter wheat only
    SimUID == 52338  # select one grid location (alternative 52339)
    
    # different weather scenarios are simulated only for FTN = 180
    # WTH == "hist"   # historical weather only
  ) |> 
    
  # convert dry matter yield to economic yield (12% moisture content)
  mutate(YLD = YLD_DM / 0.88,
         FTN = round(FTN)) # weird behavior of FTN variable, round to be safe   

# plot the distribution of yields
ggplot(sim_yield, aes(x = YLD)) +
  geom_histogram(binwidth = 0.1) +
  labs(x = "Yield (t/ha)", y = "Frequency") +
  theme_minimal() 

# summary statistics
summary <- sim_yield |>
  group_by(FTN) |>
  summarise(
    mean_YLD = mean(YLD),
    sd_YLD = sd(YLD),
    median_YLD = median(YLD),
    MAD_YLD = mad(YLD),
    min_YLD = min(YLD),
    max_YLD = max(YLD),
    n = n()
  ) |>
  mutate(across(where(is.numeric), ~ round(., 2)))

summary
# write_csv(summary, "3_Outputs/Table_YieldSummary.csv")

# create control variables 
sim_yield <- sim_yield |>
  mutate(
    irr = ifelse(IRR == "ir", 1, 0),
    til_bau = ifelse(TIL == "contill_bau", 1, 0),
    til_min_cons = ifelse(TIL == "mintill_cons", 1, 0),
    til_contill_r00 = ifelse(TIL == "contill_r00", 1, 0),
    til_contill_r30 = ifelse(TIL == "contill_r30", 1, 0),
    til_contill_r60 = ifelse(TIL == "contill_r60", 1, 0),
    til_contill_r90 = ifelse(TIL == "contill_r90", 1, 0)
  )

### Just & Pope production function estimation #######################################

yield_function <- lmrob(
  YLD ~ sqrt(FTN) +
    FTN +
    irr +
    til_min_cons +
    til_contill_r00 +
    til_contill_r30 +
    til_contill_r60 +
    til_contill_r90,
  data = sim_yield,
  method = "MM"
)

summary(yield_function)

# append residuals
sim_yield <- sim_yield |>
  mutate(residuals = yield_function$residuals, abs_residuals = abs(residuals))

variation_function <- lmrob(
  abs_residuals ~ sqrt(FTN) +
    FTN +
    irr +
    til_min_cons +
    til_contill_r00 +
    til_contill_r30 +
    til_contill_r60 +
    til_contill_r90,
  data = sim_yield,
  method = "MM"
)

summary(variation_function)

# table of results
models <- list(
  "Yield function" = yield_function,
  "Variation function" = variation_function
)

modelsummary(
  models,
  output = "3_Outputs/Table_YieldFunction.csv",
  stars = c('*' = 0.1, '**' = 0.05, '***' = 0.01)
)

################################################################################
# 2. NITROGEN PRICES - WHEAT PRICES RELATIONSHIP
################################################################################

### Import data ###############################################

f_data <- read_csv("1_Data/FertilizerPrices.csv")
f_data <- f_data[f_data[[2]] == "CAN", ] # Keep only CAN data

w_data <- read_csv("1_Data/Wheatprices.csv")
w_data <- w_data[w_data[[2]] == "Bread_Wheat", ] # alt. 'Feed_Wheat'

# Check time ranges
range(f_data$Date)
range(w_data$Date)

# Check time intervals (see if regularly spaced)
diff(f_data$Date)
unique(diff(f_data$Date))
diff(w_data$Date)
unique(diff(w_data$Date))

# Merge datasets on common dates
dat <- inner_join(w_data, f_data, by = "Date") |>
  rename(w_p = Avg_Price.x, f_p = Avg_Price.y) |>
  dplyr::select(Date, w_p, f_p) |>
  # Remove pre-2009 obs. (monthly up until that point)
  filter(Date >= as.Date("2009-01-01")) |> 
  # Transform prices to €/kg
  #mutate(w_p = w_p / 100, f_p = f_p / 100) |> 
  # take log prices
  mutate(w_p = log(w_p), f_p = log(f_p))

range(dat$Date)
diff(dat$Date)

# check for time intervals that are not 14 days
unique(diff(dat$Date))
sum(diff(dat$Date) != 14)
sum(diff(dat$Date) == 14)

# we dont have regular 14-day intervals, data is recorder every 1st and 3rd
# monday of the month (sometimes 3 week gaps are introduced)

# plot prices
price_plot <- ggplot(dat, aes(x = Date)) +
  geom_line(
    aes(y = f_p, color = "Calcium Ammonium Nitrate price"),
    linewidth = 0.5
  ) +
  geom_line(aes(y = w_p, color = "Wheat price"), linewidth = 0.5) +
  scale_color_manual(
    values = c(
      "Calcium Ammonium Nitrate price" = "blue",
      "Wheat price" = "orange"
    )
  ) +
  labs(y = "Log of price (€/kg)", x = NULL, color = NULL) +
  theme_minimal() +
  theme(
    legend.position = c(0.1, 0.9),
    legend.justification = c(0, 1),
    legend.background = element_rect(fill = "white", color = "black")
  )

price_plot

#ggsave("3_Outputs/Fig_Prices.pdf", plot = price_plot, width = 6, height = 4)


# check for stationarity
adf.test(dat$f_p, k=24)
adf.test(dat$f_p, k=12)
adf.test(dat$f_p, k=6)
adf.test(dat$f_p, k=4)
adf.test(dat$f_p, k=2)
# reject the null (non-stat) for lags >12 (half year) => stationary 
# fail to reject the null for lags <6 => non-stationary at seasonal frequencies

adf_f <- adf.test(dat$f_p, k=4)

# Dickey-Fuller test may have low power 
# (H0 not rejected, whereas there may not be a unit root)

kpss.test(dat$f_p, null = "Trend")
kpss.test(dat$f_p, null = "Level")
kpss_f <- kpss.test(dat$f_p, null = c("Level", "Trend"))
# reject null of stationarity => time series is non-stationary

pp.test(dat$f_p)
pp_f <- pp.test(dat$f_p)
# fail to reject the null (stationary)=> non-stationarity


adf.test(dat$w_p, k=24)
adf.test(dat$w_p, k=12)
adf.test(dat$w_p, k=6)
adf.test(dat$w_p, k=4)
adf.test(dat$w_p, k=2)
adf_w <- adf.test(dat$w_p, k=4)
# fail to reject the null (non-stat) for every lag => non-stationary

kpss.test(dat$w_p, null = "Trend")
kpss.test(dat$w_p, null = "Level")
kpss_w <- kpss.test(dat$w_p, null = c("Level", "Trend"))
# reject null of stationarity => non-stationary

pp.test(dat$w_p)
pp_w <- pp.test(dat$w_p)
# fail to reject the null of non-stationary => non-stationarity


# time series are non-stationary: take first differences
dif_dat <- dat |>
  mutate(
    w_p = c(NA, diff(w_p)),
    f_p = c(NA, diff(f_p))
  ) |>
  drop_na()

# check for stationarity again
adf.test(dif_dat$f_p, k = 4)
adf.test(dif_dat$w_p, k = 4)
adf_f_dif <- adf.test(dif_dat$f_p, k = 4)
adf_w_dif <- adf.test(dif_dat$w_p, k = 4)

kpss.test(dif_dat$f_p)
kpss.test(dif_dat$w_p)

kpss_f_dif <- kpss.test(dif_dat$f_p, null = c("Level", "Trend"))
kpss_w_dif <- kpss.test(dif_dat$w_p, null = c("Level", "Trend"))

pp.test(dif_dat$f_p)
pp.test(dif_dat$w_p)

pp_f_dif <- pp.test(dif_dat$f_p)
pp_w_dif <- pp.test(dif_dat$w_p)


# plot differenced prices
plot_dif_w <- ggplot(dif_dat, aes(x = Date)) +
  geom_line(aes(y = w_p), linewidth = 0.5, color = "orange") +
  #geom_hline(yintercept = 0, linetype = "dashed") +
  labs(y = "Wheat", x = NULL) +
  theme_minimal()

plot_dif_f <- ggplot(dif_dat, aes(x = Date)) +
  geom_line(aes(y = f_p), linewidth = 0.5, color = "blue") +
  #geom_hline(yintercept = 0, linetype = "dashed") +
  labs(y = "Fertilizer", x = NULL) +
  theme_minimal()

combined_dif_plots <- plot_dif_f / plot_dif_w
combined_dif_plots

# ggsave("3_Outputs/Fig_DifferencedPrices.pdf", plot = combined_dif_plots, width = 8, height = 6)

combined_plots <- price_plot / combined_dif_plots
combined_plots

#ggsave("3_Outputs/Fig_CombinedPrices.pdf", plot = combined_plots, width = 10, height = 8)

# cointegration test
coint_test <- ca.jo(dat[, c("f_p", "w_p")], type = "trace", ecdet = "trend", K = 2)
summary(coint_test)
# weak evidence of cointegration at 10% level

# summary statistics of differenced prices
summary_price_table <- tibble(
  Statistic = c("Mean", "SD", "Min", "Max", "Skew", "Kurtosis", "ADF test", 
                "KPSS test", "PP test"),
  CAN = c(mean(dat$f_p), 
          sd(dat$f_p),
          min(dat$f_p), 
          max(dat$f_p),
          skewness(dat$f_p), 
          kurtosis(dat$f_p),
          adf_f$statistic,
          kpss_f$statistic,
          pp_f$statistic
        ),
  Wheat = c(mean(dat$w_p), 
            sd(dat$w_p),
            min(dat$w_p), 
            max(dat$w_p),
            skewness(dat$w_p), 
            kurtosis(dat$w_p),
            adf_w$statistic,
            kpss_w$statistic,
            pp_w$statistic
          ),
  CAN_dif = c(mean(dif_dat$f_p), 
          sd(dif_dat$f_p),
          min(dif_dat$f_p), 
          max(dif_dat$f_p),
          skewness(dif_dat$f_p), 
          kurtosis(dif_dat$f_p),
          adf_f_dif$statistic,
          kpss_f_dif$statistic,
          pp_f_dif$statistic
        ),
  Wheat_dif = c(mean(dif_dat$w_p), 
            sd(dif_dat$w_p),
            min(dif_dat$w_p), 
            max(dif_dat$w_p),
            skewness(dif_dat$w_p), 
            kurtosis(dif_dat$w_p),
            adf_w_dif$statistic,
            kpss_w_dif$statistic,
            pp_w_dif$statistic
          )
        )

summary_price_table <- summary_price_table %>%
  mutate(across(where(is.numeric), ~ round(., 2)))

summary_price_table

#write_csv(summary_price_table, "3_Outputs/Table_PriceSummary.csv")


### Specify and Estimate VAR ####################################################

# Create regression variables
dif_dat <- dif_dat %>%
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


vardat <- dif_dat %>% dplyr::select(w_p, f_p)

# lag selection
VARselect(vardat, lag.max = 4, type = "const")
# either 1 or 2 lags suggested by info criteria

# estimate VAR
var1 <- VAR(vardat, p = 1)
summary(var1)
serial.test(var1, lags.pt = 12, type = "PT.asymptotic")
# Strong rejection of the null hypothesis - there is significant serial correlation

var2 <- VAR(vardat, p = 2)
summary(var2)
serial.test(var2, lags.pt = 12, type = "PT.asymptotic")

var3 <- VAR(vardat, p = 3)
summary(var3)
serial.test(var3, lags.pt = 12, type = "PT.asymptotic")
# still strong rejection of the null hypothesis - there is significant serial correlation

# test for volatility clustering
arch.test(var1)
arch.test(var2)
# VAR models assume constant variance - ARCH effects present

# create time trends and seasonal dummies
dif_dat <- dif_dat |>
  mutate(
    tt = 1:n(), # general time trend
    tt1 = ifelse(
      Date >= as.Date("2021-09-01"),
      row_number() - which.min(Date < as.Date("2021-09-01")),
      0
    ), # time trend post Sep 2021 (change of data)

    year = as.numeric(format(Date, "%Y")),
    month = as.numeric(format(Date, "%m")),
    day = as.numeric(format(Date, "%d")),

    # quarterly dummies
    Q1 = ifelse(month %in% c(1, 2, 3), 1, 0),
    Q2 = ifelse(month %in% c(4, 5, 6), 1, 0),
    Q3 = ifelse(month %in% c(7, 8, 9), 1, 0),
    Q4 = ifelse(month %in% c(10, 11, 12), 1, 0),

    # covid/ukraine war dummy
    D = ifelse(
      Date >= as.Date("2021-09-01") & Date <= as.Date("2023-09-01"),
      1,
      0
    )
  )


# select lag with exogenous variables
VARselect(vardat, lag.max = 4, type = "const")
VARselect(vardat, lag.max = 4, exogen = cbind(dif_dat[, c("w_ps1", "f_ps1")])) 
VARselect(
  vardat,
  lag.max = 4,
  exogen = cbind(dif_dat[, c("w_ps1", "f_ps1", "tt")])
) 
VARselect(
  vardat,
  lag.max = 4,
  exogen = cbind(dif_dat[, c("w_ps1", "f_ps1", "tt", "tt1")])
) 
VARselect(
  vardat,
  lag.max = 4,
  exogen = cbind(dif_dat[, c("w_ps1", "f_ps1", "tt", "tt1", "Q1", "Q2", "Q3")])
) 
VARselect(
  vardat,
  lag.max = 4,
  exogen = cbind(dif_dat[, c("w_ps1", "f_ps1", "tt", "tt1", "Q1", "Q2", "Q3", "D")])
) 

# best model according to AIC and BIC is p=1/2 with squared terms 
# following JP we stay with the linear model as the squared specification 
# does not produce a stable price simulation 


### Specify and Estimate QVAR ####################################################

# define individual marginal specifications based on information criteria
mf <- f_p ~ f_p1 + w_p1  
mw <- w_p ~ f_p1 + w_p1 

var_f <- lm(mf, data = dif_dat) 
summary(var_f)
var_w <- lm(mw, data = dif_dat)
summary(var_w)

# estimate QVARs for table 
taus_sparse <- c(.1, .3, .5, .7, .9)

qvar_f <- rq(mf, tau = taus_sparse, data = dif_dat)
qvar_w <- rq(mw, tau = taus_sparse, data = dif_dat)
sum_qvar_f <- summary(qvar_f, se = "boot", brmethod = "xy")
sum_qvar_w <- summary(qvar_w, se = "boot", brmethod = "xy")
sum_qvar_f
sum_qvar_w

# Calculate pseudo R^2s for table (both wheat and fertilizer)

psu_r2s_w <- numeric(length(taus_sparse))
psu_r2s_f <- numeric(length(taus_sparse))
rho <- function(u, tau) sum(u * (tau - (u < 0))) # rho calc from JP code

for (i in seq_along(taus_sparse)) {
  tau <- taus_sparse[i]

  # wheat
  fit_w <- rq(mw, tau = tau, data = dif_dat)
  y_w <- dif_dat$w_p
  yhat_w <- predict(fit_w, newdata = dif_dat)
  rho_m_w <- rho(y_w - yhat_w, tau)
  rho0_w  <- rho(y_w - quantile(y_w, probs = tau), tau)
  psu_r2s_w[i] <- 1 - rho_m_w / rho0_w

  # fertilizer
  fit_f <- rq(mf, tau = tau, data = dif_dat)
  y_f <- dif_dat$f_p
  yhat_f <- predict(fit_f, newdata = dif_dat)
  rho_m_f <- rho(y_f - yhat_f, tau)
  rho0_f  <- rho(y_f - quantile(y_f, probs = tau), tau)
  psu_r2s_f[i] <- 1 - rho_m_f / rho0_f
}

round(psu_r2s_f, 3)
round(psu_r2s_w, 3)

# make table
model_list_f <- list("VAR" = var_f)
for (tau in taus_sparse) {
  model_list_f[[paste0("", tau)]] <- rq(mf, tau = tau, data = dif_dat)
}

model_list_w <- list("VAR" = var_w)
for (tau in taus_sparse) {
  model_list_w[[paste0("", tau)]] <- rq(mw, tau = tau, data = dif_dat)
}

modelsummary(
  model_list_f,
  title = "Fertilizer Price",
  output = "3_Outputs/Table_QVAR(1)_FertPrices.csv",
  stars = c('*' = 0.1, '**' = 0.05, '***' = 0.01)
)

modelsummary(
  model_list_w,
  title = "Wheat Price",
  output = "3_Outputs/Table_QVAR(1)_WheatPrices.csv",
  stars = c('*' = 0.1, '**' = 0.05, '***' = 0.01)
)

### full QVAR estimation #################################################

taus_dense <- seq(0.01, 0.99, by = 0.01)
fit_f <- rq(mf, tau = taus_dense, data = dif_dat)
fit_w <- rq(mw, tau = taus_dense, data = dif_dat)

coeff_w <- fit_w$coeff
coeff_f <- fit_f$coeff

# Independent variables
w_p <- dif_dat$w_p
w_p1 <- dif_dat$w_p1

f_p <- dif_dat$f_p
f_p1 <- dif_dat$f_p1

X <- cbind(1, w_p1, f_p1) #col of 1s for b0

# Inverse distributions
y_w <- X %*% coeff_w
y_w <- t(apply(y_w, 1, cummax))

y_f <- X %*% coeff_f
y_f <- t(apply(y_f, 1, cummax))

### Copula estimation #################################################
n = nrow(vardat)

F_w <- array(0, dim = c(n, 1))
F_f <- array(0, dim = c(n, 1))

for (ii in 1:n) {
  F_w[ii] <- taus_dense[min(which(
    min(abs(w_p[ii] - y_w[ii, ])) == abs(w_p[ii] - y_w[ii, ])
  ))]
  F_f[ii] <- taus_dense[min(which(
    min(abs(f_p[ii] - y_f[ii, ])) == abs(f_p[ii] - y_f[ii, ])
  ))]
}

C <- F_w ~ F_f
summary(lm(C))

rC <- rq(C, tau = taus_sparse)
sum_rC <- summary(rC, se = "boot", brmethod = "xy")
sum_rC

rC_full <- rq(C, tau = taus_dense)
coeff_rC <- rC_full$coeff

plot(
  taus_dense,
  coeff_rC[2, ],
  type = "l",
  col = "steelblue",
  lwd = 2,
  xlab = "tau of fertilizer price",
  ylab = "slope"
)
abline(h = 0, lty = 2)

################################################################################
# 3. ADJUSTED CONTRACT PRICING
################################################################################

### Price path simulation ###########################################################

n_dif <- nrow(dif_dat)
n_lev <- nrow(dat)

Nt <- 24 # number of periods to simulate (24 observations twice per month = 1 year)
Ns <- 1000 # number of draws

Ys_changes <- array(0, dim = c(2, Nt, Ns))
Ys_levels <- array(0, dim = c(2, Nt, Ns))

Ys_changes[1, 1, ] <- dif_dat$f_p[n_dif]    # last observed fertilizer change at t=1
Ys_changes[2, 1, ] <- dif_dat$w_p[n_dif]    # last observed wheat change at t=1

# set initial level (only at time 1)
Ys_levels[1, 1, ] <- dat$f_p[n_lev]         # last observed fertilizer level
Ys_levels[2, 1, ] <- dat$w_p[n_lev]         # last observed wheat level

# QVAR Simulation
for (it in 2:Nt) {
  for (is in 1:Ns) {
    f_change1 <- Ys_changes[1, it - 1, is]
    w_change1 <- Ys_changes[2, it - 1, is]

    Xfi <- cbind(1, f_change1, w_change1) # matches mf: intercept, f_p1, w_p1
    Xwi <- cbind(1, f_change1, w_change1) # matches mw: intercept, f_p1, w_p1

    qi <- sample(1:99, 1)

    yf_change <- Xfi %*% coeff_f[, qi]
    yw_change <- Xwi %*% coeff_w[, qi]

    Ys_changes[1, it, is] <- yf_change
    Ys_changes[2, it, is] <- yw_change

    # convert to level using cumulative sum
    Ys_levels[1, it, is] <- dat$f_p[n_lev] + sum(Ys_changes[1, 1:it, is])
    Ys_levels[2, it, is] <- dat$w_p[n_lev] + sum(Ys_changes[2, 1:it, is])
  }
}

# Reshape data for plotting
simulated_prices <- data.frame(
  time = rep(1:Nt, Ns),
  price_f = as.vector(Ys_levels[1, , ]),
  price_w = as.vector(Ys_levels[2, , ]),
  simulation = rep(1:Ns, each = Nt)
)

ggplot(plot_data, aes(x = time, y = price_f, group = simulation)) +
  geom_line(color = "blue", alpha = 0.5) +
  labs(y = "Log of Fertilizer Price", x = "Time Period") +
  theme_minimal()

ggplot(plot_data, aes(x = time, y = price_w, group = simulation)) +
  geom_line(color = "orange", alpha = 0.5) +
  labs(y = "Log of Wheat Price", x = "Time Period") +
  theme_minimal()


ggplot() +
  geom_density(data = simulated_prices, aes(x = price_f), color = "red", size = 1) +
  geom_density(data = dat, aes(x = f_p), color = "blue", size = 1, alpha = 0.5) +
  labs(x = "Fertilizer Price (€/kg)", y = "Density") +
  theme_minimal() +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  ggtitle("Density of Historical vs Simulated Fertilizer Prices") +
  theme(plot.title = element_text(hjust = 0.5))

ggplot() +
  geom_density(data = simulated_prices, aes(x = price_w), color = "red", size = 1) +
  geom_density(data = dat, aes(x = w_p), color = "blue", size = 1, alpha = 0.5) +
  labs(x = "Fertilizer Price (€/kg)", y = "Density") +
  theme_minimal() +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  ggtitle("Density of Historical vs Simulated Fertilizer Prices") +
  theme(plot.title = element_text(hjust = 0.5))


### Insurance Contract Price ###########################################################

nitrogen_prices <- Ys_levels[1, 1:Nt, ] # Extract nitrogen price paths

# Contract parameters
duration <- 12 # contract duration (biweekly periods: 12 = 6 months)
# strike prices
k1 <- mean(dat$f_p) + sd(dat$f_p)
k2 <- mean(dat$f_p) + 1.5 * sd(dat$f_p)
amount <- 100 # Notional amount 
int <- 0.03   # Interest rate
disc <- (1 + int)**(-(1:Nt / 12)) # Monthly discount factors

payoffs_call <- array(0, dim = c(Nt, Ns))
for (it in 1:Nt) {
  payoffs_call[it, ] <- pmax(nitrogen_prices[it, ] - k1, 0) * amount
}

# Discounted payoffs
disc_cal_payoffs <- array(0, dim = c(Nt, Ns))
for (it in 1:Nt) {
  disc_cal_payoffs[it, ] <- payoffs_call[it, ] * disc[it]
}

# Average Discounted Payoff
contract_price <- mean(disc_cal_payoffs)

contract_price

#------- Analytical contract price calc (currently based on simple GBM) ----------
params <- list(
  n <- 10000,       # number of simulated paths
  t <- 10,         # total years
  dt <- 1/12,            # monthly steps
  mu <- 0.0,             # drift
  sigma <- 0.5,          # volatility
  s0 <- 200,               # starting price
  p_cap <- 200,          # price cap
  r <- 0.05             # annual discount rate
) #temporary, estimate based on data etc.

#Expected cost of adjusted contract 
price_ceiling_cost <- function(params){
  with(params, {
    t_vec <- seq(dt, t, by = dt)
    
    #black-scholes
    d1P <- (log(s0 / p_cap) + (mu + 0.5 * sigma^2) * t_vec) / (sigma * sqrt(t_vec))
    d2P <- d1P - sigma * sqrt(t_vec)
    
    cost_t <- s0 * exp(mu * t_vec) * pnorm(d1P) - p_cap * pnorm(d2P)
    cost_PV <- exp(-r * t_vec) * cost_t
    
    sum(cost_PV)
  })
}

# Expected cost of forward contract
forward_cost <- function(params, f) {
  with(params, {
    t_vec <- seq(dt, t, by = dt)
    E_P_t <- s0 * exp(mu * t_vec)
    cost_t <- E_P_t - f
    sum(exp(-r * t_vec) * cost_t)
  })
}

#Finds the cost-equal forward contract (based on adj contract cost)
cost_equal <- function(params) {
  
  target_cost <- price_ceiling_cost(params) #finds cost of adj contract
  
  cost_diff <- function(f) {
    forward_cost(params, f) - target_cost
  }
  
  uniroot(cost_diff, lower = 0, upper = 2000)$root #forward contract price where both contract are cost equivalent
}

f_star <- cost_equal(params)
f_star

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
