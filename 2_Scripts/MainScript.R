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
library(RColorBrewer)

# Set seed for replicability
set.seed(1234)

# ===============================================================================
# 1.0 NITROGEN-YIELD RELATIONSHIP ==============================================
# ===============================================================================

#-------------------------------------------------------------------------------
## 1.1 Get data from crop simulation model -------------------------------------
#-------------------------------------------------------------------------------

sim_yield <- read_csv("1_Data/EPIC.csv")

sim_yield <- sim_yield |>
  filter(
    CROP == "WWHT", # winter wheat only
    SimUID == 52338 # select one grid location (alternative 52339)
  ) |>
  mutate(
    YLD = YLD_DM / 0.88, # convert dry matter yield to economic yield (12% moisture content)
    FTN = round(FTN)
  )

# plot the distribution of yields
ggplot(sim_yield, aes(x = YLD)) +
  geom_histogram(binwidth = 0.1) +
  labs(x = "Yield (t/ha)", y = "Frequency") +
  theme_minimal()

# create control variables
sim_yield <- sim_yield |>
  mutate(
    irr = ifelse(IRR == "ir", 1, 0),
    rf = ifelse(IRR == "rf", 1, 0),
    rot_mono = ifelse(ROT == "MONO", 1, 0),
    rot_1 = ifelse(ROT == "CRS1", 1, 0),
    rot_2 = ifelse(ROT == "CRS2", 1, 0),
    rot_3 = ifelse(ROT == "CRS3", 1, 0),
    rot_4 = ifelse(ROT == "CRS4", 1, 0),
    rot_6 = ifelse(ROT == "CRS6", 1, 0),
    rot_7 = ifelse(ROT == "CRS7", 1, 0),
    wth_hist = ifelse(WTH == "hist", 1, 0),
    cmip6_ipsl_ssp126 = ifelse(WTH == "cmip6_ipsl_ssp126", 1, 0),
    cmip6_ipsl_ssp585 = ifelse(WTH == "cmip6_ipsl_ssp585", 1, 0),
    cmip6_mpi_ssp126 = ifelse(WTH == "cmip6_mpi_ssp126", 1, 0),
    cmip6_mpi_ssp585 = ifelse(WTH == "cmip6_mpi_ssp585", 1, 0),
    til_bau = ifelse(TIL == "contill_bau", 1, 0),
    til_min_cons = ifelse(TIL == "mintill_cons", 1, 0),
    til_contill_r00 = ifelse(TIL == "contill_r00", 1, 0),
    til_contill_r30 = ifelse(TIL == "contill_r30", 1, 0),
    til_contill_r60 = ifelse(TIL == "contill_r60", 1, 0),
    til_contill_r90 = ifelse(TIL == "contill_r90", 1, 0)
  )

# descriptive statistics
sim_yield_summary <- sim_yield |>
  dplyr::select(
    YLD, FTN, irr, rf, rot_mono, rot_1, rot_2, rot_3, rot_4, rot_6, rot_7,
    wth_hist, cmip6_ipsl_ssp126, cmip6_ipsl_ssp585, cmip6_mpi_ssp126, cmip6_mpi_ssp585,
    til_bau, til_min_cons, til_contill_r00, til_contill_r30, til_contill_r60, til_contill_r90
  ) |>
  pivot_longer(everything(), names_to = "variable", values_to = "value") |>
  group_by(variable) |>
  summarise(
    mean = mean(value, na.rm = TRUE),
    sd = sd(value, na.rm = TRUE),
    min = min(value, na.rm = TRUE),
    max = max(value, na.rm = TRUE),
    .groups = "drop"
  ) |>
  mutate(across(where(is.numeric), ~ round(., 2)))

# write_csv(sim_yield_summary, "3_Outputs/Table_YieldDescriptiveStatistics.csv")

#-------------------------------------------------------------------------------
## 1.2 Just & Pope production function estimation ------------------------------
#-------------------------------------------------------------------------------

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

# production function estimation
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

# variation function estimation
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

# export table of results
models <- list(
  "Yield function" = yield_function,
  "Variation function" = variation_function
)

modelsummary(
  models,
  output = "3_Outputs/Table_YieldFunction.csv",
  stars = c("*" = 0.1, "**" = 0.05, "***" = 0.01)
)

# ===============================================================================
# 2.0 Price Dynamics ===========================================================
# ===============================================================================

#-------------------------------------------------------------------------------
## 2.1 Import data -------------------------------------------------------------
#-------------------------------------------------------------------------------

# keep data only for Calcium Ammonium Nitrate (CAN)
f_data <- read_csv("1_Data/FertilizerPrices.csv")
f_data <- f_data[f_data[[2]] == "CAN", ]

# keep data only for Bread wheat
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
  filter(Date >= as.Date("2009-01-01")) |> # Remove pre-2009 obs. (monthly up until that point)
  mutate(w_p = log(w_p), f_p = log(f_p)) # Log prices

range(dat$Date)
diff(dat$Date)

# check for time intervals that are not 14 days
unique(diff(dat$Date))
sum(diff(dat$Date) != 14)
sum(diff(dat$Date) == 14)

# we dont have regular 14-day intervals, data is recorder every 1st and 3rd
# monday of the month (3 week gaps are introduced)

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

# ggsave("3_Outputs/Fig_Prices.pdf", plot = price_plot, width = 6, height = 4)

#-------------------------------------------------------------------------------
## 2.2 Test for stationarity ---------------------------------------------------
#-------------------------------------------------------------------------------

# check for stationarity for CAN prices
adf.test(dat$f_p, k = 24)
adf.test(dat$f_p, k = 12)
adf.test(dat$f_p, k = 6)
adf.test(dat$f_p, k = 4)
adf.test(dat$f_p, k = 2)
# reject the null (non-stat) for lags >12 (half year) => stationary
# fail to reject the null for lags <6 => non-stationary at seasonal frequencies

adf_f <- adf.test(dat$f_p, k = 4)

# Dickey-Fuller test may have low power
# (H0 not rejected, whereas there may not be a unit root)

kpss.test(dat$f_p, null = "Trend")
kpss.test(dat$f_p, null = "Level")
kpss_f <- kpss.test(dat$f_p, null = c("Level", "Trend"))
# reject null of stationarity => time series is non-stationary

pp.test(dat$f_p)
pp_f <- pp.test(dat$f_p)
# fail to reject the null (stationary)=> non-stationarity

# check for stationarity for wheat prices
adf.test(dat$w_p, k = 24)
adf.test(dat$w_p, k = 12)
adf.test(dat$w_p, k = 6)
adf.test(dat$w_p, k = 4)
adf.test(dat$w_p, k = 2)
adf_w <- adf.test(dat$w_p, k = 4)
# fail to reject the null (non-stat) for every lag => non-stationary

kpss.test(dat$w_p, null = "Trend")
kpss.test(dat$w_p, null = "Level")
kpss_w <- kpss.test(dat$w_p, null = c("Level", "Trend"))
# reject null of stationarity => non-stationary

pp.test(dat$w_p)
pp_w <- pp.test(dat$w_p)
# fail to reject the null of non-stationary => non-stationarity

#-------------------------------------------------------------------------------
### 2.2.1 Take first difference if not stationary ------------------------------
#-------------------------------------------------------------------------------

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

#-------------------------------------------------------------------------------
## 2.3 Plot --------------------------------------------------------------------
#-------------------------------------------------------------------------------

plot_dif_w <- ggplot(dif_dat, aes(x = Date)) +
  geom_line(aes(y = w_p), linewidth = 0.5, color = "orange") +
  # geom_hline(yintercept = 0, linetype = "dashed") +
  labs(y = "Wheat", x = NULL) +
  theme_minimal()

plot_dif_f <- ggplot(dif_dat, aes(x = Date)) +
  geom_line(aes(y = f_p), linewidth = 0.5, color = "blue") +
  # geom_hline(yintercept = 0, linetype = "dashed") +
  labs(y = "Fertilizer", x = NULL) +
  theme_minimal()

combined_dif_plots <- plot_dif_f / plot_dif_w
combined_dif_plots

# ggsave("3_Outputs/Fig_DifferencedPrices.pdf", plot = combined_dif_plots, width = 8, height = 6)

combined_plots <- price_plot / combined_dif_plots
combined_plots

# ggsave("3_Outputs/Fig_CombinedPrices.pdf", plot = combined_plots, width = 10, height = 8)

# cointegration test
coint_test <- ca.jo(dat[, c("f_p", "w_p")], type = "trace", ecdet = "trend", K = 2)
summary(coint_test)
# weak evidence of cointegration at 10% level

#-------------------------------------------------------------------------------
## 2.4 Summary statistics of differenced prices --------------------------------
#-------------------------------------------------------------------------------

summary_price_table <- tibble(
  Statistic = c(
    "Mean", "SD", "Min", "Max", "Skew", "Kurtosis", "ADF test",
    "KPSS test", "PP test"
  ),
  CAN = c(
    mean(dat$f_p),
    sd(dat$f_p),
    min(dat$f_p),
    max(dat$f_p),
    skewness(dat$f_p),
    kurtosis(dat$f_p),
    adf_f$statistic,
    kpss_f$statistic,
    pp_f$statistic
  ),
  Wheat = c(
    mean(dat$w_p),
    sd(dat$w_p),
    min(dat$w_p),
    max(dat$w_p),
    skewness(dat$w_p),
    kurtosis(dat$w_p),
    adf_w$statistic,
    kpss_w$statistic,
    pp_w$statistic
  ),
  CAN_dif = c(
    mean(dif_dat$f_p),
    sd(dif_dat$f_p),
    min(dif_dat$f_p),
    max(dif_dat$f_p),
    skewness(dif_dat$f_p),
    kurtosis(dif_dat$f_p),
    adf_f_dif$statistic,
    kpss_f_dif$statistic,
    pp_f_dif$statistic
  ),
  Wheat_dif = c(
    mean(dif_dat$w_p),
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

# write_csv(summary_price_table, "3_Outputs/Table_PriceSummary.csv")

#-------------------------------------------------------------------------------
## 2.5 Create deflated prices for later comparison -----------------------------
#-------------------------------------------------------------------------------

# import PPI data
ppi_output <- read_xlsx("1_Data/RawData/PPI_Cereals_Germany.xlsx", skip = 4)
ppi_input <- read_xlsx("1_Data/RawData/PPI_Fertilizers_Germany.xlsx", skip = 4)

base_year <- 2020

ppi_w <- ppi_output[4, -c(1, 2)] |>
  pivot_longer(
    cols = everything(), # pivot all columns
    names_to = "year",
    values_to = "ppi"
  ) |>
  mutate(
    year = as.numeric(year), # convert year strings to numeric (optional)
    ppi = as.numeric(ppi)
  )

ppi_w <- rbind(ppi_w, data.frame(year = 2025, ppi = tail(ppi_w$ppi, 1)))

ppi_f <- ppi_input[3, -c(1, 2)] |>
  mutate(across(everything(), as.double)) |>
  pivot_longer(
    cols = everything(), # pivot all columns
    names_to = "year",
    values_to = "ppi"
  ) |>
  mutate(
    year = as.numeric(year), # convert year strings to numeric (optional)
    ppi = as.numeric(ppi)
  )

ppi_f <- rbind(ppi_f, data.frame(year = 2025, ppi = tail(ppi_f$ppi, 1)))

ppi_base <- 100

defl_f_p <- f_data |>
  filter(Date >= as.Date("2009-01-01")) |>
  dplyr::select(Date, Avg_Price) |>
  mutate(
    t = row_number(),
    year = as.numeric(format(Date, "%Y")),
    month = as.numeric(format(Date, "%m")),
    year_index = year - min(year) + 1,
    f_p = Avg_Price
  ) |>
  left_join(ppi_f, by = "year") |>
  mutate(real_f_p = f_p * (ppi_base / ppi))

defl_w_p <- w_data |>
  filter(Date >= as.Date("2009-01-01")) |>
  dplyr::select(Date, Avg_Price) |>
  mutate(
    t = row_number(),
    year = as.numeric(format(Date, "%Y")),
    month = as.numeric(format(Date, "%m")),
    year_index = year - min(year) + 1,
    w_p = Avg_Price
  ) |>
  left_join(ppi_w, by = "year") |>
  mutate(real_w_p = w_p * (ppi_base / ppi))

#-------------------------------------------------------------------------------
## 2.6 Test for best QVAR model ------------------------------------------------
#-------------------------------------------------------------------------------

# Create regression variables
dif_dat <- dif_dat |>
  mutate(
    w_p1 = lag(w_p, 1), # wheat lags
    w_p2 = lag(w_p, 2),
    w_p3 = lag(w_p, 3),
    w_p4 = lag(w_p, 4),
    f_p1 = lag(f_p, 1), # fert lags
    f_p2 = lag(f_p, 2),
    f_p3 = lag(f_p, 3),
    f_p4 = lag(f_p, 4),
    w_ps1 = w_p1^2, # squares
    f_ps1 = f_p1^2
  ) |>
  drop_na()

vardat <- dif_dat |> dplyr::select(w_p, f_p)

# lag selection
VARselect(vardat, lag.max = 4, type = "const")
# either 1 or 2 lags suggested by info criteria

# estimate VAR(1)
var1 <- VAR(vardat, p = 1)
summary(var1)
serial.test(var1, lags.pt = 12, type = "PT.asymptotic")
# Strong rejection of the null hypothesis - there is significant serial correlation

# estimate VAR(2)
var2 <- VAR(vardat, p = 2)
summary(var2)
serial.test(var2, lags.pt = 12, type = "PT.asymptotic")

# estimate VAR(3)
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
VARselect(vardat, lag.max = 12, exogen = cbind(dif_dat[, c("w_ps1", "f_ps1")]))
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
# we stay with the linear model as the squared specification does not produce a stable price simulation

#-------------------------------------------------------------------------------
## 2.7 QVAR estimation ---------------------------------------------------------
#-------------------------------------------------------------------------------

# define individual marginal specifications based on information criteria
mf <- f_p ~ f_p1 + w_p1
mw <- w_p ~ f_p1 + w_p1

# var estimation
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
  rho0_w <- rho(y_w - quantile(y_w, probs = tau), tau)
  psu_r2s_w[i] <- 1 - rho_m_w / rho0_w

  # fertilizer
  fit_f <- rq(mf, tau = tau, data = dif_dat)
  y_f <- dif_dat$f_p
  yhat_f <- predict(fit_f, newdata = dif_dat)
  rho_m_f <- rho(y_f - yhat_f, tau)
  rho0_f <- rho(y_f - quantile(y_f, probs = tau), tau)
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
  stars = c("*" = 0.1, "**" = 0.05, "***" = 0.01)
)

modelsummary(
  model_list_w,
  title = "Wheat Price",
  output = "3_Outputs/Table_QVAR(1)_WheatPrices.csv",
  stars = c("*" = 0.1, "**" = 0.05, "***" = 0.01)
)

# full QVAR estimation
taus_dense <- seq(0.01, 0.99, by = 0.01)
fit_f <- rq(mf, tau = taus_dense, data = dif_dat)
fit_w <- rq(mw, tau = taus_dense, data = dif_dat)

# save the coefficients
coeff_w <- fit_w$coeff
coeff_f <- fit_f$coeff

# Get independent variables
w_p <- dif_dat$w_p
w_p1 <- dif_dat$w_p1

f_p <- dif_dat$f_p
f_p1 <- dif_dat$f_p1

X <- cbind(1, w_p1, f_p1) # col of 1s for b0

# Inverse distributions
y_w <- X %*% coeff_w
y_w <- t(apply(y_w, 1, cummax))

y_f <- X %*% coeff_f
y_f <- t(apply(y_f, 1, cummax))

#-------------------------------------------------------------------------------
## 2.8 Copula estimation -------------------------------------------------------
#-------------------------------------------------------------------------------

# get number of rows
n <- nrow(vardat)

# initialize arrays
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

#-------------------------------------------------------------------------------
## 2.9 Dynamic stability analysis ----------------------------------------------
#-------------------------------------------------------------------------------

ntau <- length(taus_sparse)
ncof <- ncol(X)
bootreps <- 1000

rqbs <- array(0, dim = c(2, ntau, bootreps, ncof))
Rev <- array(0, dim = c(bootreps, ntau, ntau))
MRev <- array(0, dim = c(bootreps, ntau, ntau))
res <- array(0, dim = c(ntau, ntau))
qMRev <- array(0, dim = c(3, ntau, ntau))


# bootstrap quantile regressions
for (i in 1:ntau) {
  # Fertilizer
  rqbs[1, i, , ] <- boot.rq(
    y = dif_dat$f_p,
    x = X,
    tau = taus_sparse[i],
    R = bootreps
  )$B

  # Wheat
  rqbs[2, i, , ] <- boot.rq(
    y = dif_dat$w_p,
    x = X,
    tau = taus_sparse[i],
    R = bootreps
  )$B
}

# Projection matrix
for (f in 1:ntau) {
  for (w in 1:ntau) {
    for (r in 1:bootreps) {
      PM <- rbind(
        c(rqbs[1, f, r, 2], rqbs[1, f, r, 3]), # fert eq
        c(rqbs[2, w, r, 2], rqbs[2, w, r, 3]) # wheat eq
      )

      ev <- eigen(PM)$values
      Rev[r, f, w] <- ev[1]
      MRev[r, f, w] <- Mod(ev[1])
    }

    res[f, w] <- mean(Rev[, f, w])
    qMRev[, f, w] <- quantile(MRev[, f, w], c(0.9, 0.95, 0.99))
  }
}

# Test stability at tau = 0.5,0.5
t.test(MRev[, 3, 3], mu = 1)


# ===============================================================================
# 3.0 Price paths simulation ===================================================
# ===============================================================================

set.seed(1234)

n_dif <- nrow(dif_dat)
n_lev <- nrow(dat)

# number of periods to simulate (24 observations twice per month = 1 year)
Nt <- 24

# number of draws
Ns <- 1000

#-------------------------------------------------------------------------------
## 3.1 Simulation --------------------------------------------------------------
#-------------------------------------------------------------------------------

# initialize arrays
Ys_changes <- array(0, dim = c(2, Nt, Ns))
Ys_levels <- array(0, dim = c(2, Nt, Ns))

Ys_changes[1, 1, ] <- dif_dat$f_p[n_dif] # last observed fertilizer change at t=1
Ys_changes[2, 1, ] <- dif_dat$w_p[n_dif] # last observed wheat change at t=1

# set initial level (only at time 1)
Ys_levels[1, 1, ] <- dat$f_p[n_lev] # last observed fertilizer level
Ys_levels[2, 1, ] <- dat$w_p[n_lev] # last observed wheat level

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

#-------------------------------------------------------------------------------
### 3.1.1 Plot simulated prices ------------------------------------------------
#-------------------------------------------------------------------------------

simulated_prices <- data.frame(
  time = rep(1:Nt, Ns),
  price_f = as.vector(Ys_levels[1, , ]),
  price_w = as.vector(Ys_levels[2, , ]),
  price_f_level = exp(as.vector(Ys_levels[1, , ])),
  price_w_level = exp(as.vector(Ys_levels[2, , ])),
  simulation = rep(1:Ns, each = Nt)
)

ggplot(simulated_prices, aes(x = time, y = price_f, group = simulation)) +
  geom_line(color = "blue", alpha = 0.5) +
  labs(y = "Log of Fertilizer Price", x = "Time Period") +
  theme_minimal()

ggplot(simulated_prices, aes(x = time, y = price_w, group = simulation)) +
  geom_line(color = "orange", alpha = 0.5) +
  labs(y = "Log of Wheat Price", x = "Time Period") +
  theme_minimal()

# density distributions vs historical prices
density_plot_f <- ggplot() +
  geom_density(data = simulated_prices, aes(x = price_f_level, fill = "Simulated"), alpha = 0.5) +
  geom_density(data = f_data, aes(x = Avg_Price, fill = "Historical"), alpha = 0.4) +
  scale_fill_manual(
    name = NULL,
    values = c("Simulated" = "#E41A1C", "Historical" = "#377EB8")
  ) +
  labs(
    x = paste0("\nN = ", nrow(simulated_prices)),
    y = "Density\n",
    title = "CAN Price Distribution"
  ) +
  theme_minimal() +
  theme(
    legend.position = c(0.85, 0.85),
    legend.background = element_rect(fill = "white", color = "black"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank()
  )

density_plot_w <- ggplot() +
  geom_density(data = simulated_prices, aes(x = price_w_level, fill = "Simulated"), alpha = 0.5) +
  geom_density(data = w_data, aes(x = Avg_Price, fill = "Historical"), alpha = 0.4) +
  scale_fill_manual(
    name = NULL,
    values = c("Simulated" = "#E41A1C", "Historical" = "#377EB8")
  ) +
  labs(
    x = paste0("\nN = ", nrow(simulated_prices)),
    y = "Density\n",
    title = "Wheat Price Distribution"
  ) +
  theme_minimal() +
  theme(
    legend.position = c(0.85, 0.85),
    legend.background = element_rect(fill = "white", color = "black"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank()
  )

density_plot_f
density_plot_w

# ggsave("3_Outputs/Fig_PriceDensity_CAN.pdf", plot = density_plot_f, width = 10, height = 6)
# ggsave("3_Outputs/Fig_PriceDensity_Wheat.pdf", plot = density_plot_w, width = 10, height = 6)

# density distributions vs deflated historical prices
density_plot_f_deflated <- ggplot() +
  geom_density(data = simulated_prices, aes(x = price_f_level, fill = "Simulated prices"), alpha = 0.5) +
  geom_density(data = defl_f_p, aes(x = real_f_p, fill = "Deflated historical prices"), alpha = 0.4) +
  scale_fill_manual(
    name = NULL,
    values = c("Simulated prices" = "#E41A1C", "Deflated historical prices" = "#377EB8")
  ) +
  labs(
    x = paste0("\n(eur/t), N = ", nrow(simulated_prices)),
    y = "Density\n",
    title = "CAN Price Distributions"
  ) +
  theme_minimal() +
  theme(
    legend.position = "null",
    legend.background = element_rect(fill = "white", color = "black"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  coord_cartesian(ylim = c(0, 0.025))

density_plot_w_deflated <- ggplot() +
  geom_density(data = simulated_prices, aes(x = price_w_level, fill = "Simulated prices"), alpha = 0.5) +
  geom_density(data = defl_w_p, aes(x = real_w_p, fill = "Deflated historical prices"), alpha = 0.4) +
  scale_fill_manual(
    name = NULL,
    values = c("Simulated prices" = "#E41A1C", "Deflated historical prices" = "#377EB8")
  ) +
  labs(
    x = paste0("\n(eur/t), N = ", nrow(simulated_prices)),
    y = NULL,
    title = "Wheat Price Distributions"
  ) +
  theme_minimal() +
  theme(
    legend.position = c(0.7, 0.85),
    legend.background = element_rect(fill = "white", color = "black"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.y = element_blank()
  )

deflated_density_plot <- density_plot_f_deflated | density_plot_w_deflated

# ggsave("3_Outputs/Fig_PriceDensity_deflated.pdf", plot = deflated_density_plot, width = 10, height = 6)


# ===============================================================================
# 4.0 Insurance Contract Pricing ===============================================
# ===============================================================================

# Contract parameters
obs_per_year <- 24
simulated_years <- Nt / obs_per_year
# exercise in early march
exercise_time <- 5 + (0:(simulated_years - 1)) * obs_per_year
int <- 0.03 # Interest rate
T_years <- 0.5 # time to maturity (6 months)
discount_factor <- 1 / (1 + int)^T_years

#-------------------------------------------------------------------------------
## 4.1 pricing based on simulated nitrogen prices ------------------------------
#-------------------------------------------------------------------------------

# get simulated prices
sim_f_p <- array(dim = c(Ns, simulated_years))

for (t in 1:simulated_years) {
  time <- exercise_time[t]
  sim_f_p[, t] <- exp(as.vector(Ys_levels[1, time, ])) / 1000
}
sim_f_p_df <- data.frame(sim_f_p)

# define strike prices
max_prices <- quantile(sim_f_p, probs = c(0.5, 0.6, 0.7, 0.8, 0.9))

premium_array <- numeric(length(max_prices))

for (i in seq_along(max_prices)) {
  strike <- max_prices[i]

  payoff <- pmax(sim_f_p - strike, 0)
  discounted_payoff <- payoff * discount_factor

  premium_array[i] <- mean(discounted_payoff)
}

premium_results <- data.frame(
  strike_price = max_prices,
  option_premium = premium_array
) |>
  mutate(across(everything(), ~ round(., 3)))

# write_csv(premium_results, "3_Outputs/Table_ContractPayoff.csv")

#-------------------------------------------------------------------------------
### 4.1.1 Plot distribution density --------------------------------------------
#-------------------------------------------------------------------------------

strike_vec <- premium_results$strike_price
cb_colors <- brewer.pal(3, "Dark2")[1:2]

# Plot density of original sim_f_p with strike prices
plot_original <- ggplot(data.frame(sim_f_p), aes(x = sim_f_p)) +
  geom_density(alpha = 0.3, linewidth = 0.5, color = cb_colors[2], fill = cb_colors[2]) +
  geom_vline(xintercept = strike_vec, color = cb_colors[1], linetype = "dashed") +
  labs(
    x = "Simulated prices and strike levels\n",
    y = "Density"
  ) +
  theme_minimal() +
  coord_cartesian(ylim = c(0, 80), xlim = c(0.28, 0.43))

# one plot per strike price (capped distribution)
plot_list <- lapply(seq_along(strike_vec), function(i) {
  k <- strike_vec[i]
  sim_f_p_adj <- pmin(sim_f_p, k)

  df <- data.frame(
    value = c(sim_f_p, sim_f_p_adj),
    type = c(
      rep("Original", length(sim_f_p)),
      rep(paste0("Capped (k = ", k, ")"), length(sim_f_p_adj))
    )
  )

  ggplot(df, aes(x = value, color = type, fill = type)) +
    geom_density(alpha = 0.3, linewidth = 0.5) +
    scale_color_manual(values = cb_colors, labels = c("Adjusted contract", "No contract")) +
    scale_fill_manual(values = cb_colors, labels = c("Adjusted contract", "No contract")) +
    labs(
      x = paste("Strike =", k, "\n"),
      y = NULL,
      color = NULL,
      fill = NULL
    ) +
    theme_minimal() +
    theme(legend.position = "top") +
    coord_cartesian(ylim = c(0, 80), xlim = c(0.28, 0.43))
})

plot_list[[2]] <- plot_list[[2]] + theme(legend.position = "none")
plot_list[[3]] <- plot_list[[3]] + labs(y = "Density\n") + theme(legend.position = "none")
plot_list[[4]] <- plot_list[[4]] + labs(x = paste("Strike = 0.375 \n\nPrice (eur/kg)")) + theme(legend.position = "none")
plot_list[[5]] <- plot_list[[5]] + theme(legend.position = "none")

# Combine into 2x3 layout
final_plot <- plot_original + plot_list[[1]] + plot_list[[2]] + plot_list[[3]] +
  plot_list[[4]] + plot_list[[5]]

final_plot

# ggsave("3_Outputs/Fig_PriceDensities_SimulatedStrikes.pdf", plot = final_plot, width = 12, height = 8)

# ===============================================================================
# SIMULATION OF UTILITY OF PROFITS
# ===============================================================================

# ===============================================================================
# 5.0 Broad framing ============================================================
# ===============================================================================

# Parameters
N <- 200 # amount of N applied
CAN <- N / 0.27 # amount of CAN (27% N)

#-------------------------------------------------------------------------------
### 5.0.1 Simulate stochastic yields -------------------------------------------
#-------------------------------------------------------------------------------

# get coeff from production function regression
coeff_0_yield <- yield_function$coefficients[1]
coeff_1_yield <- yield_function$coefficients[2]
coeff_2_yield <- yield_function$coefficients[3]

# generate predicted yield
predicted_yield <- coeff_0_yield + coeff_1_yield * sqrt(N) + coeff_2_yield * N

# get coeff from variation function regression
coeff_0_varf <- variation_function$coefficients[1]
coeff_1_varf <- variation_function$coefficients[2]

# generate predicted variation
predicted_sigma_yield <- coeff_0_varf + coeff_1_varf * sqrt(N)

# derive stochastic error
stochastic_error <- array(dim = c(Ns, simulated_years))

for (t in 1:simulated_years) {
  stochastic_error[, t] <- rnorm(n = Ns, mean = 0, sd = predicted_sigma_yield)
}

# simulation of stochastic yields
simulated_yield <- array(dim = c(Ns, simulated_years))

for (t in 1:simulated_years) {
  simulated_yield[, t] <- predicted_yield + stochastic_error[, t]
}

# get simulated wheat prices (t = 13 -> early July)
output_selling_time <- 13 + (0:(simulated_years - 1)) * obs_per_year

sim_w_p <- array(dim = c(Ns, simulated_years))

# get simulated wheat prices (t = 13 -> early July)
for (t in 1:simulated_years) {
  time <- output_selling_time[t]
  sim_w_p[, t] <- exp(as.vector(Ys_levels[2, time, ]))
}

str(Ys_levels)

#-------------------------------------------------------------------------------
### 5.0.2 Get price simulations and calculate profits --------------------------
#-------------------------------------------------------------------------------

# get forward price
sim_f_p_forward <- rep(tail(f_data$Avg_Price / 1000, 1) * discount_factor, Ns)

# Calculate "no contract" profits
profit_nc <- (sim_w_p * simulated_yield) - (sim_f_p * CAN)

# Calculate "forward" profits
profit_fc <- (sim_w_p * simulated_yield) - (sim_f_p_forward * CAN)

# Calculate "adjusted contract" profits:
# Initialize result data frame
df_simulated_profit <- data.frame(
  nc = profit_nc,
  fc = profit_fc
)

# Loop through each strike price and calculate adjusted contract profits
for (i in 1:nrow(premium_results)) {
  k <- premium_results$strike_price[i]
  contract_cost <- premium_results$option_premium[i]

  sim_f_p_adj_contr <- pmin(sim_f_p, k)

  profit_ac <- (sim_w_p * simulated_yield) - (sim_f_p_adj_contr * CAN) - (contract_cost * CAN)

  # Add result columns for each strike price
  col_name <- paste0("ac_", i)
  df_simulated_profit[[col_name]] <- profit_ac
}

#-------------------------------------------------------------------------------
## 5.1 Summary tables ----------------------------------------------------------
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
### 5.1.1 CAN price paid under different contracts -----------------------------
#-------------------------------------------------------------------------------

base_row_sim_prices <- tibble(
  strike = NA,
  mean = mean(sim_f_p),
  median = median(sim_f_p),
  sd = sd(sim_f_p),
  min = min(sim_f_p),
  max = max(sim_f_p)
)

base_row_forward_prices <- tibble(
  strike = NA,
  mean = mean(sim_f_p_forward),
  median = median(sim_f_p_forward),
  sd = sd(sim_f_p_forward),
  min = min(sim_f_p_forward),
  max = max(sim_f_p_forward)
)

stats_list <- list()

for (i in 1:nrow(premium_results)) {
  k <- premium_results$strike_price[i]

  sim_f_p_adj_contr <- pmin(sim_f_p, k)

  stats_list[[i]] <- tibble(
    strike = k,

    # stats for capped simulated prices
    mean   = mean(sim_f_p_adj_contr),
    median = median(sim_f_p_adj_contr),
    sd     = sd(sim_f_p_adj_contr),
    min    = min(sim_f_p_adj_contr),
    max    = max(sim_f_p_adj_contr)
  )
}

# combine into a single dataframe
price_stats_table <- bind_rows(base_row_sim_prices, base_row_forward_prices, stats_list) |>
  mutate(across(where(is.numeric), ~ round(., 3)))

price_stats_table

# write_csv(price_stats_table, "3_Outputs/Table_ResultsPriceStats.csv")

#-------------------------------------------------------------------------------
### 5.1.2 Total costs under different contracts --------------------------------
#-------------------------------------------------------------------------------

base_row_costs <- tibble(
  strike = NA,
  mean = mean(sim_f_p * CAN),
  median = median(sim_f_p * CAN),
  sd = sd(sim_f_p * CAN),
  min = min(sim_f_p * CAN),
  max = max(sim_f_p * CAN)
)

base_row_forward_costs <- tibble(
  strike = NA,
  mean = mean(sim_f_p_forward * CAN),
  median = median(sim_f_p_forward * CAN),
  sd = sd(sim_f_p_forward * CAN),
  min = min(sim_f_p_forward * CAN),
  max = max(sim_f_p_forward * CAN),
)

stats_list <- list()

for (i in 1:nrow(premium_results)) {
  k <- premium_results$strike_price[i]
  cost <- premium_results$option_premium[i]

  sim_f_p_adj_contr <- pmin(sim_f_p, k)

  stats_list[[i]] <- tibble(
    strike = k,
    mean   = mean((sim_f_p_adj_contr * CAN) + (cost * CAN)),
    median = median((sim_f_p_adj_contr * CAN) + (cost * CAN)),
    sd     = sd((sim_f_p_adj_contr * CAN) + (cost * CAN)),
    min    = min((sim_f_p_adj_contr * CAN) + (cost * CAN)),
    max    = max((sim_f_p_adj_contr * CAN) + (cost * CAN))
  )
}

# combine into a single dataframe
costs_stats_table <- bind_rows(base_row_costs, base_row_forward_costs, stats_list) |>
  mutate(across(where(is.numeric), ~ round(., 2)))

costs_stats_table

# write_csv(costs_stats_table, "3_Outputs/Table_ResultsCostsStats.csv")

#-------------------------------------------------------------------------------
### 5.1.3 Total profits under different contracts ------------------------------
#-------------------------------------------------------------------------------

summary_profits <- df_simulated_profit |>
  pivot_longer(everything(), names_to = "variable", values_to = "value") |>
  group_by(variable) |>
  summarise(
    mean = mean(value),
    median = median(value),
    sd = sd(value),
    min = min(value),
    max = max(value)
  ) |>
  mutate(across(where(is.numeric), ~ round(., 2)))

summary_profits

# write_csv(summary_profits, "3_Outputs/Table_ResultsProfitStats.csv")

#-------------------------------------------------------------------------------
## 5.2 Expected Utility Theory (EUT) values under broad framing ----------------
#-------------------------------------------------------------------------------

# define absolute risk aversion values
delta <- c(0.00001, 0.00005, 0.0001, 0.0005, 0.001, 0.005)

target_cols <- colnames(df_simulated_profit)

# CE and RP with Mean-Variance approach
eut_results_broadframing <- data.frame(
  Variable = character(),
  Delta = numeric(),
  CertaintyEquivalent = numeric(),
  RiskPremium = numeric()
)

for (d in delta) {
  for (c in target_cols) {
    profits <- df_simulated_profit[[c]]

    mean_profit <- mean(profits)
    var_profit <- var(profits)
    RP <- ((d / 2) * var_profit)
    CE <- mean_profit - RP

    eut_results_broadframing <- rbind(eut_results_broadframing, data.frame(
      Variable = c,
      Delta = d,
      CertaintyEquivalent = CE,
      RiskPremium = RP
    ))
  }
}

# make table of results
eut_results_broadframing <- eut_results_broadframing |>
  pivot_wider(
    names_from = Variable,
    values_from = c(CertaintyEquivalent, RiskPremium),
    names_sep = "_"
  ) |>
  mutate(across(where(is.numeric), ~ round(.x, 2)))

# Reorder columns to group CE and RP by contract
contract_columns <- c()
for (contract in target_cols) {
  contract_columns <- c(
    contract_columns,
    paste0("CertaintyEquivalent_", contract),
    paste0("RiskPremium_", contract)
  )
}

eut_results_broadframing <- eut_results_broadframing |>
  dplyr::select(Delta, all_of(contract_columns))

# rename columns
colnames(eut_results_broadframing) <- c(
  "Delta",
  unlist(lapply(target_cols, function(x) c(paste0(x, "_CE"), paste0(x, "_RP"))))
)

eut_results_broadframing

# write_csv(eut_results_broadframing, "3_Outputs/Table_Results_EUT_Broadframing.csv")

#-------------------------------------------------------------------------------
## 5.3 Cumulative Prospect Theory (CPT) under broad framing --------------------
#-------------------------------------------------------------------------------

# risk aversion
a_values <- c(0.1, 0.334, 0.5)

# loss aversion
lambda <- c(0.5, 1.574, 2)

# probability weighting
gamma <- c(0.4, 0.571, 0.7)

# reference point
ref_point <- mean(df_simulated_profit$fc)

# Initialize result data frame
cpt_results_broadframing <- data.frame(
  Variable = character(),
  a = numeric(),
  lambda = numeric(),
  gamma = numeric(),
  CPT_Value = numeric()
)

# value function
v <- function(x, a, l) {
  ifelse(x >= ref_point, (x - ref_point)^a, -l * (-(x - ref_point))^a)
}

# weighting function
weighting_function <- function(p, g) {
  exp(-(-log(p))^g)
}


for (a in a_values) {
  for (l in lambda) {
    for (g in gamma) {
      for (c in target_cols) {
        # get the current profit vector and set the reference point
        profits <- sort(df_simulated_profit[[c]])

        # uniform probabilities
        n <- length(profits)
        probs <- rep(1 / n, n)
        df <- data.frame(profit = profits, prob = probs)

        # Separate into losses and gains and order
        losses <- df |>
          dplyr::filter(profit < ref_point) |>
          dplyr::arrange(profit)
        gains <- df |>
          dplyr::filter(profit >= ref_point) |>
          dplyr::arrange(dplyr::desc(profit))

        # Calculate Decision Weights (dw) for Losses
        if (nrow(losses) > 0) {
          losses <- losses |>
            dplyr::mutate(
              cum_p = cumsum(prob),
              dw = weighting_function(cum_p, g) - weighting_function(cum_p - prob, g)
            )
        }

        # Calculate Decision Weights (dw) for Gains
        if (nrow(gains) > 0) {
          # Get the highest cumulative probability from losses
          cum_p_losses <- if (nrow(losses) > 0) max(losses$cum_p) else 0

          gains <- gains |>
            dplyr::mutate(
              cum_p = cum_p_losses + cumsum(prob),
              dw = weighting_function(cum_p, g) - weighting_function(cum_p - prob, g)
            )
        }

        # Calculate CPT Value
        loss_value <- if (nrow(losses) > 0) {
          sum(v(losses$profit, a, l) * losses$dw)
        } else {
          0
        }

        gain_value <- if (nrow(gains) > 0) {
          sum(v(gains$profit, a, l) * gains$dw)
        } else {
          0
        }

        cpt_val_current <- loss_value + gain_value

        # 7. Store the result
        new_row <- data.frame(
          Variable = c,
          Alpha = a,
          Lambda = l,
          Gamma = g,
          CPT_Value = cpt_val_current
        )

        cpt_results_broadframing <- rbind(cpt_results_broadframing, new_row)
      }
    }
  }
}

#-------------------------------------------------------------------------------
### 5.3.1 display CPT results --------------------------------------------------
#-------------------------------------------------------------------------------

# standard assumptions
alpha_std <- 0.334
lambda_std <- 1.574
gamma_std <- 0.571

# standard scenario
std_row <- cpt_results_broadframing |>
  filter(
    Alpha == alpha_std,
    Lambda == lambda_std,
    Gamma == gamma_std
  ) |>
  dplyr::select(Variable, CPT_Value) |>
  pivot_wider(
    names_from = Variable,
    values_from = CPT_Value
  ) |>
  mutate(Scenario = "Standard Assumptions")

# varying risk aversion
alpha_rows <- cpt_results_broadframing |>
  filter(
    Alpha != alpha_std,
    Lambda == lambda_std,
    Gamma == gamma_std
  ) |>
  dplyr::select(Alpha, Variable, CPT_Value) |>
  pivot_wider(
    names_from = Variable,
    values_from = CPT_Value
  ) |>
  rename(Scenario = Alpha) |>
  mutate(Scenario = paste0("Alpha = ", Scenario))

# varying loss aversion
lambda_rows <- cpt_results_broadframing |>
  filter(
    Alpha == alpha_std,
    Lambda != lambda_std,
    Gamma == gamma_std
  ) |>
  dplyr::select(Lambda, Variable, CPT_Value) |>
  pivot_wider(
    names_from = Variable,
    values_from = CPT_Value
  ) |>
  rename(Scenario = Lambda) |>
  mutate(Scenario = paste0("Lambda = ", Scenario))

# varying prob weighting
gamma_rows <- cpt_results_broadframing |>
  filter(
    Alpha == alpha_std,
    Lambda == lambda_std,
    Gamma != gamma_std
  ) |>
  dplyr::select(Gamma, Variable, CPT_Value) |>
  pivot_wider(
    names_from = Variable,
    values_from = CPT_Value
  ) |>
  rename(Scenario = Gamma) |>
  mutate(Scenario = paste0("Gamma = ", Scenario))

# combine
final_cpt_table <- bind_rows(std_row, alpha_rows, lambda_rows, gamma_rows) |>
  relocate(Scenario) |>
  mutate(across(where(is.numeric), ~ round(., 2)))

final_cpt_table

# write_csv(final_cpt_table, "3_Outputs/Table_Results_CPT_BroadFraming.csv")

# ===============================================================================
# 6.0 Narrow Framing ===========================================================
# ===============================================================================

#-------------------------------------------------------------------------------
## 6.0.1 calculate profits ------------------------------------------------------
#-------------------------------------------------------------------------------

# forward contract
profit_fc_nf <- CAN * (sim_f_p_forward - sim_f_p)

# adjusted contract
# Initialize result data frame
df_profit_narrowframing <- data.frame(
  fc = profit_fc_nf
)

for (i in 1:nrow(premium_results)) {
  k <- premium_results$strike_price[i]
  contract_cost <- premium_results$option_premium[i]

  profit_ac <- CAN * (pmax(sim_f_p - k, 0) - contract_cost)

  # Add column with name ac_1, ac_2, etc.
  col_name <- paste0("ac_", i)
  df_profit_narrowframing[[col_name]] <- profit_ac
}

df_profit_narrowframing

# summary table
summary_profits_narrowframing <- df_profit_narrowframing |>
  pivot_longer(everything(), names_to = "variable", values_to = "value") |>
  group_by(variable) |>
  summarise(
    mean = mean(value),
    median = median(value),
    sd = sd(value),
    min = min(value),
    max = max(value)
  ) |>
  mutate(across(where(is.numeric), ~ round(., 2)))

summary_profits_narrowframing

# write_csv(summary_profits_narrowframing, "3_Outputs/Table_Results_Profits_NarrowFraming.csv")

#-------------------------------------------------------------------------------
## 6.1 Expected Utility Theory (EUT) under narrow framing ----------------------
#-------------------------------------------------------------------------------

delta <- c(0, 0.0005, 0.0007, 0.001, 0.003)

target_cols <- colnames(df_profit_narrowframing)

eut_results_narrowframing <- data.frame(
  Variable = character(),
  Delta = numeric(),
  RiskPremium = numeric(),
  CertaintyEquivalent = numeric()
)

for (d in delta) {
  for (c in target_cols) {
    profits <- df_profit_narrowframing[[c]]

    RP <- ((d / 2) * var(profits))
    CE <- mean(profits) - RP

    temporary_results <- data.frame(
      Variable = c,
      Delta = d,
      RiskPremium = RP,
      CertaintyEquivalent = CE
    )

    eut_results_narrowframing <- rbind(eut_results_narrowframing, temporary_results)
  }
}

# make table of results
final_eut_table_narrowframing <- eut_results_narrowframing |>
  pivot_wider(
    names_from = Variable,
    values_from = c(CertaintyEquivalent, RiskPremium),
    names_sep = "_"
  ) |>
  mutate(across(where(is.numeric), ~ round(.x, 2)))

# Reorder columns to group CE and RP by contract
contract_columns <- c()
for (contract in target_cols) {
  contract_columns <- c(
    contract_columns,
    paste0("CertaintyEquivalent_", contract),
    paste0("RiskPremium_", contract)
  )
}

final_eut_table_narrowframing <- final_eut_table_narrowframing |>
  dplyr::select(Delta, all_of(contract_columns))

# rename columns
colnames(final_eut_table_narrowframing) <- c(
  "Delta",
  unlist(lapply(target_cols, function(x) c(paste0(x, "_CE"), paste0(x, "_RP"))))
)

final_eut_table_narrowframing

# write_csv(final_eut_table_narrowframing, "3_Outputs/Table_Results_EUT_NarrowFraming.csv")

#-------------------------------------------------------------------------------
## 6.2 Cumulative Prospect Theory (CPT) under narrow framing -------------------
#-------------------------------------------------------------------------------

ref_point <- 0

# Initialize result data frame
cpt_results_nf <- data.frame(
  Variable = character(),
  a = numeric(),
  lambda = numeric(),
  gamma = numeric(),
  CPT_Value = numeric()
)

for (a in a_values) {
  for (l in lambda) {
    for (g in gamma) {
      for (c in target_cols) {
        # get the current profit vector and set the reference point
        profits <- sort(df_profit_narrowframing[[c]])

        # uniform probabilities
        n <- length(profits)
        probs <- rep(1 / n, n)
        df <- data.frame(profit = profits, prob = probs)

        # Separate into losses and gains and order
        losses <- df |>
          dplyr::filter(profit < ref_point) |>
          dplyr::arrange(profit)
        gains <- df |>
          dplyr::filter(profit >= ref_point) |>
          dplyr::arrange(dplyr::desc(profit))

        # Calculate Decision Weights (dw) for Losses
        if (nrow(losses) > 0) {
          losses <- losses |>
            dplyr::mutate(
              cum_p = cumsum(prob),
              dw = weighting_function(cum_p, g) - weighting_function(cum_p - prob, g)
            )
        }

        # Calculate Decision Weights (dw) for Gains
        if (nrow(gains) > 0) {
          # Get the highest cumulative probability from losses
          cum_p_losses <- if (nrow(losses) > 0) max(losses$cum_p) else 0

          gains <- gains |>
            dplyr::mutate(
              cum_p = cum_p_losses + cumsum(prob),
              dw = weighting_function(cum_p, g) - weighting_function(cum_p - prob, g)
            )
        }

        # Calculate CPT Value
        loss_value <- if (nrow(losses) > 0) {
          sum(v(losses$profit, a, l) * losses$dw)
        } else {
          0
        }

        gain_value <- if (nrow(gains) > 0) {
          sum(v(gains$profit, a, l) * gains$dw)
        } else {
          0
        }

        cpt_val_current <- loss_value + gain_value

        # 7. Store the result
        new_row <- data.frame(
          Variable = c,
          Alpha = a,
          Lambda = l,
          Gamma = g,
          CPT_Value = cpt_val_current
        )

        cpt_results_nf <- rbind(cpt_results_nf, new_row)
      } # End col_name loop
    } # End gamma loop
  } # End lambda loop
} # End alpha loop

#-------------------------------------------------------------------------------
### 6.2.1 display CPT results --------------------------------------------------
#-------------------------------------------------------------------------------

# standard assumptions
alpha_std <- 0.334
lambda_std <- 1.574
gamma_std <- 0.571

# standard scenario
std_row <- cpt_results_nf |>
  filter(
    Alpha == alpha_std,
    Lambda == lambda_std,
    Gamma == gamma_std
  ) |>
  dplyr::select(Variable, CPT_Value) |>
  pivot_wider(
    names_from = Variable,
    values_from = CPT_Value
  ) |>
  mutate(Scenario = "Standard Assumptions")

# varying risk aversion
alpha_rows <- cpt_results_nf |>
  filter(
    Alpha != alpha_std,
    Lambda == lambda_std,
    Gamma == gamma_std
  ) |>
  dplyr::select(Alpha, Variable, CPT_Value) |>
  pivot_wider(
    names_from = Variable,
    values_from = CPT_Value
  ) |>
  rename(Scenario = Alpha) |>
  mutate(Scenario = paste0("Alpha = ", Scenario))

# varying loss aversion
lambda_rows <- cpt_results_nf |>
  filter(
    Alpha == alpha_std,
    Lambda != lambda_std,
    Gamma == gamma_std
  ) |>
  dplyr::select(Lambda, Variable, CPT_Value) |>
  pivot_wider(
    names_from = Variable,
    values_from = CPT_Value
  ) |>
  rename(Scenario = Lambda) |>
  mutate(Scenario = paste0("Lambda = ", Scenario))

# varying prob weighting
gamma_rows <- cpt_results_nf |>
  filter(
    Alpha == alpha_std,
    Lambda == lambda_std,
    Gamma != gamma_std
  ) |>
  dplyr::select(Gamma, Variable, CPT_Value) |>
  pivot_wider(
    names_from = Variable,
    values_from = CPT_Value
  ) |>
  rename(Scenario = Gamma) |>
  mutate(Scenario = paste0("Gamma = ", Scenario))

# combine
final_cpt_narrowframing <- bind_rows(std_row, alpha_rows, lambda_rows, gamma_rows) |>
  relocate(Scenario) |>
  mutate(across(where(is.numeric), ~ round(., 2)))

final_cpt_narrowframing

# write_csv(final_cpt_narrowframing, "3_Outputs/Table_Results_CPT_NarrowFraming.csv")
