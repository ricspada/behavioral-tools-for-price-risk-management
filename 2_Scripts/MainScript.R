# Load required libraries
library(quantreg)
library(readxl)
library(modelsummary)
library(robustbase)
library(tidyverse)
library(vars)
library(tseries)
library(forecast)
library(patchwork)
library(urca)
library(moments)
library(RColorBrewer)
#library(stargazer)
library(strucchange)

# Toggle to enable/disable writing output files (CSV tables, PDF figures)
SAVE_OUTPUTS <- TRUE

# Set seed for replicability
set.seed(1234)

# ==============================================================================
# 1.0 NITROGEN-YIELD RELATIONSHIP ==============================================
# ==============================================================================

#-------------------------------------------------------------------------------
## 1.1 Get data from crop simulation model -------------------------------------
#-------------------------------------------------------------------------------

sim_yield <- read_csv("1_Data/EPIC.csv")

sim_yield <- sim_yield |>
  filter(
    # winter wheat only
    CROP == "WWHT",
    # grid location of interest (alternative: 52339)
    SimUID == 52338,
    # exlude alternative weather scenario pathways
    WTH == "hist"
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
    YLD,
    FTN,
    irr,
    rf,
    rot_mono,
    rot_1,
    rot_2,
    rot_3,
    rot_4,
    rot_6,
    rot_7,
    wth_hist,
    cmip6_ipsl_ssp126,
    cmip6_ipsl_ssp585,
    cmip6_mpi_ssp126,
    cmip6_mpi_ssp585,
    til_bau,
    til_min_cons,
    til_contill_r00,
    til_contill_r30,
    til_contill_r60,
    til_contill_r90
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

if (SAVE_OUTPUTS) {
  write_csv(sim_yield_summary, "3_Outputs/Table_YieldDescriptiveStatistics.csv")
}

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
# future wheater scenarios only simulate yield for FTN = 180, should we keep
# or exlude them from the estimation? (the benefit is to introduce stochasticity)

if (SAVE_OUTPUTS) {
  write_csv(summary, "3_Outputs/Table_YieldSummary.csv")
}

#-------------------------------------------------------------------------------
## 1.2 Just & Pope production function estimation ------------------------------
#-------------------------------------------------------------------------------

# Weather and rotation dummies are excluded because the downstream simulation
# uses a fixed management scenario (rainfed, conventional tillage); their
# effects are captured in the residual variance.
# The EPIC simulation design varies these factors jointly.
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

residuals <- yield_function$residuals
hist(residuals, breaks = 50)
shapiro.test(sample(residuals, 500))
skewness(residuals)

# Variation function estimation
# Using absolute residuals (not squared) as the dependent variable.
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

if (SAVE_OUTPUTS) {
  modelsummary(
    models,
    output = "3_Outputs/Table_YieldFunction.csv",
    stars = c("*" = 0.1, "**" = 0.05, "***" = 0.01)
  )
}

# ===============================================================================
# 2.0 Price Dynamics ===========================================================
# ===============================================================================

#-------------------------------------------------------------------------------
## 2.1 Import data -------------------------------------------------------------
#-------------------------------------------------------------------------------

# keep data only for Calcium Ammonium Nitrate (CAN)
f_data <- read_csv("1_Data/FertilizerPrices.csv")
f_data <- f_data[f_data[[2]] == "CAN", ]

# keep data only for Bread wheat (alt. 'Feed_Wheat')
w_data <- read_csv("1_Data/WheatPrices.csv")
w_data <- w_data[w_data[[2]] == "Bread_Wheat", ]

# Merge datasets on common dates
# Remove pre-2009 obs. (monthly up until that point)
dat <- inner_join(w_data, f_data, by = "Date") |>
  rename(w_p = Avg_Price.x, f_p = Avg_Price.y) |>
  dplyr::select(Date, w_p, f_p) |>
  filter(Date >= as.Date("2009-01-01")) |>
  mutate(w_p = log(w_p), f_p = log(f_p))

range(dat$Date)
diff(dat$Date)

# check for time intervals that are not 14 days
unique(diff(dat$Date))
sum(diff(dat$Date) != 14)
sum(diff(dat$Date) == 14)

# we dont have regular 14-day intervals, data is recorder every 1st and 3rd
# monday of the month (3 week gaps are introduced)

# plot prices

# consistent colors across plots
colors <- c("darkblue", "darkgreen")

price_plot <- ggplot(dat, aes(x = Date)) +
  geom_line(
    aes(y = f_p, color = "Calcium Ammonium Nitrate price"),
    linewidth = 0.5
  ) +
  geom_line(aes(y = w_p, color = "Wheat price"), linewidth = 0.5) +
  scale_color_manual(
    values = c(
      "Calcium Ammonium Nitrate price" = colors[1],
      "Wheat price" = colors[2]
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

# clear outlier in wheat prices, formally checkin for outliers

# ------------------------------------------------------------------------------
## Outliers --------------------------------------------------------------------
# ------------------------------------------------------------------------------

outliers_f <- tsoutliers(dat$f_p)
outliers_w <- tsoutliers(dat$w_p)

outliers_f$index
outliers_f$replacements

outliers_w$index
outliers_w$replacements

# flag outlier points
dat$outlier_f <- FALSE
dat$outlier_f[outliers_f$index] <- TRUE
dat$outlier_w <- FALSE
dat$outlier_w[outliers_w$index] <- TRUE

# plot identified outliers
outliers_plot <- price_plot +
  geom_point(
    data = subset(dat, outlier_f),
    aes(x = Date, y = f_p),
    color = "red"
  ) +
  geom_point(
    data = subset(dat, outlier_w),
    aes(x = Date, y = w_p),
    color = "red"
  )

outliers_plot

# identified outliers for fertilizer occour after price shocks thus are not removed
# replacing only wheat price outlier
dat$w_p[outliers_w$index] <- outliers_w$replacements

# plot clean time series
clean_price_plot <- ggplot(dat, aes(x = Date)) +
  geom_line(
    aes(y = f_p, color = "Calcium Ammonium Nitrate price"),
    linewidth = 0.5
  ) +
  geom_line(aes(y = w_p, color = "Wheat price"), linewidth = 0.5) +
  scale_color_manual(
    values = c(
      "Calcium Ammonium Nitrate price" = colors[1],
      "Wheat price" = colors[2]
    )
  ) +
  labs(y = "Log of price (€/kg)", x = NULL, color = NULL) +
  theme_minimal() +
  theme(
    legend.position = c(0.1, 0.9),
    legend.justification = c(0, 1),
    legend.background = element_rect(fill = "white", color = "black")
  )

clean_price_plot


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
### 2.2.1 Take first difference ------------------------------------------------
#-------------------------------------------------------------------------------

# taking the first difference as prices are non-stationary
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

plot_dif_f <- ggplot(dif_dat, aes(x = Date)) +
  geom_line(aes(y = f_p), linewidth = 0.5, color = colors[1]) +
  labs(y = "Fertilizer", x = NULL) +
  theme_minimal()

plot_dif_w <- ggplot(dif_dat, aes(x = Date)) +
  geom_line(aes(y = w_p), linewidth = 0.5, color = colors[2]) +
  labs(y = "Wheat", x = NULL) +
  theme_minimal()

combined_dif_plots <- (plot_dif_f +
  theme(axis.text.x = element_blank()) +
  coord_cartesian(xlim = as.Date(c("2009-07-01", "2025-08-01")))) /
  (plot_dif_w +
    coord_cartesian(xlim = as.Date(c("2009-07-01", "2025-08-01"))))
combined_dif_plots

combined_plots <- (clean_price_plot +
  theme(axis.text.x = element_blank()) +
  coord_cartesian(xlim = as.Date(c("2009-07-01", "2025-08-01")))) /
  combined_dif_plots
combined_plots

if (SAVE_OUTPUTS) {
  ggsave(
    "3_Outputs/Fig_CombinedPrices.pdf",
    plot = combined_plots,
    width = 10,
    height = 8
  )
}

# cointegration test
coint_test <- ca.jo(
  dat[, c("f_p", "w_p")],
  type = "trace",
  ecdet = "trend",
  K = 2
)
summary(coint_test)
# weak evidence of cointegration at 10% level

#-------------------------------------------------------------------------------
## 2.4 Summary statistics of differenced prices --------------------------------
#-------------------------------------------------------------------------------

summary_price_table <- tibble(
  Statistic = c(
    "Mean",
    "SD",
    "Min",
    "Max",
    "Skew",
    "Kurtosis",
    "ADF test",
    "KPSS test",
    "PP test"
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

if (SAVE_OUTPUTS) {
  write_csv(summary_price_table, "3_Outputs/Table_PriceSummary.csv")
}

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
var_dat <- dif_dat |>
  mutate(
    w_p1 = lag(w_p, 1), # wheat lags
    w_p2 = lag(w_p, 2),
    w_p3 = lag(w_p, 3),
    w_p4 = lag(w_p, 4),
    f_p1 = lag(f_p, 1), # fert lags
    f_p2 = lag(f_p, 2),
    f_p3 = lag(f_p, 3),
    f_p4 = lag(f_p, 4),
    w_ps1 = w_p1^2, # squared lags
    f_ps1 = f_p1^2
  ) |>
  drop_na()

### 2.6.1 Baseline lag selection -----------------------------------------------

VARselect(var_dat[, c("w_p", "f_p")], lag.max = 12, type = "const")
# all criteria select p = 1

### 2.6.2 VAR(1) baseline diagnostics -----------------------------------------

var1 <- VAR(var_dat[, c("w_p", "f_p")], p = 1)
summary(var1)

serial.test(var1, lags.pt = 12, type = "PT.asymptotic")
serial.test(var1, lags.pt = 12, type = "BG")
# BG test p = 0.04: evidence of serial correlation at 5%

plot(stability(var1))
plot(stability(var1, type = "fluctuation"))
# fluctuation test suggests instability -> check for structural breaks

### 2.6.3 Structural break analysis --------------------------------------------

bp_fp <- breakpoints(f_p ~ f_p1 + w_p1, data = var_dat)
bp_wp <- breakpoints(w_p ~ f_p1 + w_p1, data = var_dat)
summary(bp_fp)
summary(bp_wp)
plot(bp_fp)
plot(bp_wp)
# 2 structural breaks for f_p at obs 198 (~2017-04) and 288 (~2021-01)
# no break for w_p

### 2.6.4 ARCH test ------------------------------------------------------------

arch.test(var1)
# ARCH effects present -> QVAR approach

### 2.6.5 Create exogenous variables -------------------------------------------

var_dat <- var_dat |>
  mutate(
    tt = 1:n(),
    tt1 = ifelse(
      Date >= as.Date("2021-09-01"),
      row_number() - which.min(Date < as.Date("2021-09-01")),
      0
    ),
    year = as.numeric(format(Date, "%Y")),
    month = as.numeric(format(Date, "%m")),
    Q1 = ifelse(month %in% c(1, 2, 3), 1, 0),
    Q2 = ifelse(month %in% c(4, 5, 6), 1, 0),
    Q3 = ifelse(month %in% c(7, 8, 9), 1, 0),
    Q4 = ifelse(month %in% c(10, 11, 12), 1, 0),
    D2 = ifelse(seq_len(nrow(var_dat)) >= 288, 1, 0),
    D_shock = ifelse(seq_len(nrow(var_dat)) %in% 285:315, 1, 0)
  )

### 2.6.6 Models comparison -----------------------------------------

# endogenous variables
vardat_endo <- var_dat[, c("w_p", "f_p")]

# helper to extract information criteria and diagnostics
get_var_ic <- function(mod) {
  aic <- AIC(mod$varresult$w_p) + AIC(mod$varresult$f_p)
  bic <- BIC(mod$varresult$w_p) + BIC(mod$varresult$f_p)
  bg <- serial.test(mod, lags.pt = 12, type = "BG")
  c(AIC = aic, BIC = bic, BG_pval = unname(bg$serial$p.value))
}

# estimate candidate models
var_m1 <- VAR(vardat_endo, p = 1)
var_m2 <- VAR(
  vardat_endo,
  p = 1,
  exogen = as.matrix(var_dat[, c("w_ps1", "f_ps1")])
)
var_m3 <- VAR(
  vardat_endo,
  p = 1,
  exogen = as.matrix(var_dat[, "D2", drop = FALSE])
)
var_m4 <- VAR(
  vardat_endo,
  p = 1,
  exogen = as.matrix(var_dat[, "D_shock", drop = FALSE])
)
var_m5 <- VAR(
  vardat_endo,
  p = 1,
  exogen = as.matrix(var_dat[, c("Q1", "Q2", "Q3")])
)
var_m6 <- VAR(
  vardat_endo,
  p = 1,
  exogen = as.matrix(var_dat[, c("w_ps1", "f_ps1", "D_shock")])
)
var_m7 <- VAR(
  vardat_endo,
  p = 1,
  exogen = as.matrix(var_dat[, c("w_ps1", "f_ps1", "Q1", "Q2", "Q3")])
)
var_m8 <- VAR(vardat_endo, p = 1, exogen = as.matrix(var_dat[, c("tt", "tt1")]))
var_m9 <- VAR(vardat_endo, p = 2)
var_m10 <- VAR(
  vardat_endo,
  p = 2,
  exogen = as.matrix(var_dat[, c("w_ps1", "f_ps1")])
)
var_m11 <- VAR(
  vardat_endo,
  p = 2,
  exogen = as.matrix(var_dat[, "D2", drop = FALSE])
)
var_m12 <- VAR(
  vardat_endo,
  p = 2,
  exogen = as.matrix(var_dat[, c("w_ps1", "f_ps1", "Q1", "Q2", "Q3")])
)

model_list <- list(
  "VAR(1)" = var_m1,
  "VAR(1) + sq" = var_m2,
  "VAR(1) + D2" = var_m3,
  "VAR(1) + D_shock" = var_m4,
  "VAR(1) + Q" = var_m5,
  "VAR(1) + sq + D_shock" = var_m6,
  "VAR(1) + sq + Q" = var_m7,
  "VAR(1) + tt" = var_m8,
  "VAR(2)" = var_m9,
  "VAR(2) + sq" = var_m10,
  "VAR(2) + D2" = var_m11,
  "VAR(2) + sq + Q" = var_m12
)

var_comparison <- as.data.frame(t(sapply(model_list, get_var_ic)))
var_comparison$Model <- rownames(var_comparison)
var_comparison <- var_comparison[, c("Model", "AIC", "BIC", "BG_pval")]
rownames(var_comparison) <- NULL

var_comparison

# Best by BIC: VAR(1) baseline
#
# Squared lags (w_p1^2, f_p1^2) consistently improve fit and remove
# serial correlation
#
# Selected specification for QVAR:
# VAR(1)
# VAR(1) + squared lags

if (SAVE_OUTPUTS) {
  write.csv(
    var_comparison,
    "3_Outputs/Table_VARComparison.csv",
    row.names = FALSE
  )
}

#-------------------------------------------------------------------------------
## 2.7 QVAR estimation ---------------------------------------------------------
#-------------------------------------------------------------------------------

# define individual marginal specifications based on information criteria
mf <- f_p ~ f_p1 + w_p1
mw <- w_p ~ f_p1 + w_p1

# var estimation
var_f <- lm(mf, data = var_dat)
summary(var_f)
var_w <- lm(mw, data = var_dat)
summary(var_w)

# estimate QVARs for table
taus_sparse <- c(.1, .3, .5, .7, .9)

qvar_f <- rq(mf, tau = taus_sparse, data = var_dat)
qvar_w <- rq(mw, tau = taus_sparse, data = var_dat)
sum_qvar_f <- summary(qvar_f, se = "boot", brmethod = "xy")
sum_qvar_w <- summary(qvar_w, se = "boot", brmethod = "xy")
sum_qvar_f
sum_qvar_w

# Calculate pseudo R^2s for table (both wheat and fertilizer)
psu_r2s_w <- numeric(length(taus_sparse))
psu_r2s_f <- numeric(length(taus_sparse))

# rho calculation
rho <- function(u, tau) sum(u * (tau - (u < 0)))

for (i in seq_along(taus_sparse)) {
  tau <- taus_sparse[i]

  # wheat
  fit_w <- rq(mw, tau = tau, data = var_dat)
  y_w <- var_dat$w_p
  yhat_w <- predict(fit_w, newdata = var_dat)
  rho_m_w <- rho(y_w - yhat_w, tau)
  rho0_w <- rho(y_w - quantile(y_w, probs = tau), tau)
  psu_r2s_w[i] <- 1 - rho_m_w / rho0_w

  # fertilizer
  fit_f <- rq(mf, tau = tau, data = var_dat)
  y_f <- var_dat$f_p
  yhat_f <- predict(fit_f, newdata = var_dat)
  rho_m_f <- rho(y_f - yhat_f, tau)
  rho0_f <- rho(y_f - quantile(y_f, probs = tau), tau)
  psu_r2s_f[i] <- 1 - rho_m_f / rho0_f
}

round(psu_r2s_f, 3)
round(psu_r2s_w, 3)

### create QVAR tables ---------------------------------------------------

# Note: substitute QVAR R2 with pseudo R2 derived above

model_list_f <- list("VAR" = var_f)
for (tau in taus_sparse) {
  model_list_f[[paste0("", tau)]] <- rq(mf, tau = tau, data = var_dat)
}

model_list_w <- list("VAR" = var_w)
for (tau in taus_sparse) {
  model_list_w[[paste0("", tau)]] <- rq(mw, tau = tau, data = var_dat)
}

if (SAVE_OUTPUTS) {
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
}

## full QVAR estimation
taus_dense <- seq(0.01, 0.99, by = 0.01)
fit_f <- rq(mf, tau = taus_dense, data = var_dat)
fit_w <- rq(mw, tau = taus_dense, data = var_dat)

# save the coefficients for the copula
coeff_w <- fit_w$coeff
coeff_f <- fit_f$coeff

# Get independent variables
w_p <- var_dat$w_p
w_p1 <- var_dat$w_p1

f_p <- var_dat$f_p
f_p1 <- var_dat$f_p1

X <- cbind(1, f_p1, w_p1)

# Get inverse distributions to build the copula
y_w <- X %*% coeff_w
y_w <- t(apply(y_w, 1, cummax))

y_f <- X %*% coeff_f
y_f <- t(apply(y_f, 1, cummax))

#-------------------------------------------------------------------------------
## 2.8 Copula estimation -------------------------------------------------------
#-------------------------------------------------------------------------------

# get number of rows
n <- nrow(var_dat)

# initialize arrays
F_w <- array(0, dim = c(n, 1))
F_f <- array(0, dim = c(n, 1))

# for each observation, finds the conditional quantile the realized price
# change fall closest to. The result F_w[ii] is a value in [0.01, 0.99]
# used to build the copula.
for (ii in 1:n) {
  F_w[ii] <- taus_dense[min(which(
    min(abs(w_p[ii] - y_w[ii, ])) == abs(w_p[ii] - y_w[ii, ])
  ))]
  F_f[ii] <- taus_dense[min(which(
    min(abs(f_p[ii] - y_f[ii, ])) == abs(f_p[ii] - y_f[ii, ])
  ))]
}

# conditional copula
C <- F_w ~ F_f
summary(lm(C))

# copula table (sparse quantiles with bootstrap SEs)
rC <- rq(C, tau = taus_sparse)
sum_rC <- summary(rC, se = "boot", brmethod = "xy")
sum_rC

# full copula (dense quantiles for simulation)
rC_full <- rq(C, tau = taus_dense)
coeff_rC <- rC_full$coeff

# plot copula coefficients across quantiles
plot(
  taus_dense,
  coeff_rC[2, ],
  type = "l",
  col = "steelblue",
  lwd = 2,
  xlab = "quantile",
  ylab = "coefficient",
  main = "Copula: F_w | F_f"
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
    y = var_dat$f_p,
    x = X,
    tau = taus_sparse[i],
    R = bootreps
  )$B

  # Wheat
  rqbs[2, i, , ] <- boot.rq(
    y = var_dat$w_p,
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


# ==============================================================================
# 3.0 Price paths simulation ===================================================
# ==============================================================================

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

# Copula Simulation
for (it in 2:Nt) {
  for (is in 1:Ns) {
    f_change1 <- Ys_changes[1, it - 1, is]
    w_change1 <- Ys_changes[2, it - 1, is]

    Xfi <- cbind(1, f_change1, w_change1) # matches mf: intercept, f_p1, w_p1
    Xwi <- cbind(1, f_change1, w_change1) # matches mw: intercept, f_p1, w_p1

    # Draw independent marginal quantiles
    qi_f <- sample(1:99, 1)
    qi_w <- sample(1:99, 1)

    yf_change <- Xfi %*% coeff_f[, qi_f]
    yw_change <- Xwi %*% coeff_w[, qi_w]

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
  geom_density(
    data = simulated_prices,
    aes(x = price_f_level, fill = "Simulated"),
    alpha = 0.5
  ) +
  geom_density(
    data = f_data,
    aes(x = Avg_Price, fill = "Historical"),
    alpha = 0.4
  ) +
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
  geom_density(
    data = simulated_prices,
    aes(x = price_w_level, fill = "Simulated"),
    alpha = 0.5
  ) +
  geom_density(
    data = w_data,
    aes(x = Avg_Price, fill = "Historical"),
    alpha = 0.4
  ) +
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

if (SAVE_OUTPUTS) {
  ggsave(
    "3_Outputs/Fig_PriceDensity_CAN.pdf",
    plot = density_plot_f,
    width = 10,
    height = 6
  )
  ggsave(
    "3_Outputs/Fig_PriceDensity_Wheat.pdf",
    plot = density_plot_w,
    width = 10,
    height = 6
  )
}

# density distributions vs deflated historical prices
density_plot_f_deflated <- ggplot() +
  geom_density(
    data = simulated_prices,
    aes(x = price_f_level, fill = "Simulated prices"),
    alpha = 0.5
  ) +
  geom_density(
    data = defl_f_p,
    aes(x = real_f_p, fill = "Deflated historical prices"),
    alpha = 0.4
  ) +
  scale_fill_manual(
    name = NULL,
    values = c(
      "Simulated prices" = "#E41A1C",
      "Deflated historical prices" = "#377EB8"
    )
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
  geom_density(
    data = simulated_prices,
    aes(x = price_w_level, fill = "Simulated prices"),
    alpha = 0.5
  ) +
  geom_density(
    data = defl_w_p,
    aes(x = real_w_p, fill = "Deflated historical prices"),
    alpha = 0.4
  ) +
  scale_fill_manual(
    name = NULL,
    values = c(
      "Simulated prices" = "#E41A1C",
      "Deflated historical prices" = "#377EB8"
    )
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

if (SAVE_OUTPUTS) {
  ggsave(
    "3_Outputs/Fig_PriceDensity_deflated.pdf",
    plot = deflated_density_plot,
    width = 10,
    height = 6
  )
}


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

# get simulated fertilizer prices
sim_f_p <- array(dim = c(Ns, simulated_years))

for (t in 1:simulated_years) {
  time <- exercise_time[t]
  sim_f_p[, t] <- exp(as.vector(Ys_levels[1, time, ]))
}

sim_f_p_df <- data.frame(sim_f_p)

# get simulated wheat prices (t=13, early July)
output_selling_time <- 13 + (0:(simulated_years - 1)) * obs_per_year

sim_w_p <- array(dim = c(Ns, simulated_years))

for (t in 1:simulated_years) {
  sim_w_p[, t] <- exp(as.vector(Ys_levels[2, output_selling_time[t], ]))
}

#-------------------------------------------------------------------------------
## 4.1 Input contract pricing based on simulated prices ------------------------
#-------------------------------------------------------------------------------

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

if (SAVE_OUTPUTS) {
  write_csv(premium_results, "3_Outputs/Table_ContractPayoff.csv")
}

#-------------------------------------------------------------------------------
### 4.1.1 Plot distribution density - fertilizer -------------------------------
#-------------------------------------------------------------------------------

strike_vec <- premium_results$strike_price
cb_colors <- brewer.pal(3, "Dark2")[1:2]

# Plot density of original sim_f_p with strike prices
plot_original <- ggplot(data.frame(sim_f_p), aes(x = sim_f_p)) +
  geom_density(
    alpha = 0.3,
    linewidth = 0.5,
    color = cb_colors[2],
    fill = cb_colors[2]
  ) +
  geom_vline(
    xintercept = strike_vec,
    color = cb_colors[1],
    linetype = "dashed"
  ) +
  labs(
    x = "Simulated prices and strike levels\n",
    y = "Density"
  ) +
  theme_minimal()

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
    scale_color_manual(
      values = cb_colors,
      labels = c("Adjusted contract", "No contract")
    ) +
    scale_fill_manual(
      values = cb_colors,
      labels = c("Adjusted contract", "No contract")
    ) +
    labs(
      x = paste("Strike =", k, "\n"),
      y = NULL,
      color = NULL,
      fill = NULL
    ) +
    theme_minimal() +
    theme(legend.position = "top")
})

plot_list[[2]] <- plot_list[[2]] + theme(legend.position = "none")
plot_list[[3]] <- plot_list[[3]] +
  labs(y = "Density\n") +
  theme(legend.position = "none")
plot_list[[4]] <- plot_list[[4]] +
  labs(x = paste("Strike =", round(strike_vec[4], 3), "\n\nPrice (eur/t)")) +
  theme(legend.position = "none")
plot_list[[5]] <- plot_list[[5]] + theme(legend.position = "none")

# Combine into 2x3 layout
final_plot <- plot_original +
  plot_list[[1]] +
  plot_list[[2]] +
  plot_list[[3]] +
  plot_list[[4]] +
  plot_list[[5]]

final_plot

if (SAVE_OUTPUTS) {
  ggsave(
    "3_Outputs/Fig_PriceDensities_SimulatedStrikes.pdf",
    plot = final_plot,
    width = 12,
    height = 8
  )
}

#-------------------------------------------------------------------------------
## 4.2 Output Contract Pricing based on simulated wheat prices -----------------
#-------------------------------------------------------------------------------

# Wheat is sold at t=13 (early July), approximately 12 months from the autumn
# reference at which contracts are initiated. The put option matures at delivery.
T_years_w <- 1.0
discount_factor_w <- 1 / (1 + int)^T_years_w

# Strike prices: lower quantiles provide downside (floor) protection on revenue.
# Symmetric to the fertilizer call's 0.5-0.9 range: same 40-pct span, lower tail.
min_prices <- quantile(sim_w_p, probs = c(0.1, 0.2, 0.3, 0.4, 0.5))

premium_array_w <- numeric(length(min_prices))

for (i in seq_along(min_prices)) {
  strike_w <- min_prices[i]

  # Put option payoff: max(k_w - P_w, 0) per tonne
  payoff_w <- pmax(strike_w - sim_w_p, 0)
  premium_array_w[i] <- mean(payoff_w * discount_factor_w)
}

premium_results_w <- data.frame(
  strike_price = min_prices,
  option_premium = premium_array_w
) |>
  mutate(across(everything(), ~ round(., 3)))

premium_results_w

if (SAVE_OUTPUTS) {
  write_csv(premium_results_w, "3_Outputs/Table_ContractPayoff_Wheat.csv")
}

#-------------------------------------------------------------------------------
### 4.2.1 Plot distribution density - wheat ------------------------------------
#-------------------------------------------------------------------------------

strike_vec_w <- premium_results_w$strike_price

# Plot density of sim_w_p with strike price levels
plot_original_w <- ggplot(data.frame(sim_w_p), aes(x = sim_w_p)) +
  geom_density(
    alpha = 0.3,
    linewidth = 0.5,
    color = cb_colors[2],
    fill = cb_colors[2]
  ) +
  geom_vline(
    xintercept = strike_vec_w,
    color = cb_colors[1],
    linetype = "dashed"
  ) +
  labs(
    x = "Simulated wheat prices and strike levels\n",
    y = "Density"
  ) +
  theme_minimal() +
  coord_cartesian(ylim = c(0, 0.06))

# One plot per strike price (floored distribution)
plot_list_w <- lapply(seq_along(strike_vec_w), function(i) {
  k_w <- strike_vec_w[i]
  sim_w_p_adj <- pmax(sim_w_p, k_w)

  df_w <- data.frame(
    value = c(sim_w_p, sim_w_p_adj),
    type = c(
      rep("Original", length(sim_w_p)),
      rep(paste0("Floored (k = ", k_w, ")"), length(sim_w_p_adj))
    )
  )

  ggplot(df_w, aes(x = value, color = type, fill = type)) +
    geom_density(alpha = 0.3, linewidth = 0.5) +
    scale_color_manual(
      values = cb_colors,
      labels = c("Min-price contract", "No contract")
    ) +
    scale_fill_manual(
      values = cb_colors,
      labels = c("Min-price contract", "No contract")
    ) +
    labs(
      x = paste("Strike =", k_w, "\n"),
      y = NULL,
      color = NULL,
      fill = NULL
    ) +
    theme_minimal() +
    theme(legend.position = "top") +
    coord_cartesian(ylim = c(0, 0.06))
})

plot_list_w[[2]] <- plot_list_w[[2]] + theme(legend.position = "none")
plot_list_w[[3]] <- plot_list_w[[3]] +
  labs(y = "Density\n") +
  theme(legend.position = "none")
plot_list_w[[4]] <- plot_list_w[[4]] +
  labs(x = paste("Strike =", round(strike_vec_w[4], 3), "\n\nPrice (eur/t)")) +
  theme(legend.position = "none")
plot_list_w[[5]] <- plot_list_w[[5]] + theme(legend.position = "none")

# Combine into 2x3 layout
final_plot_w <- plot_original_w +
  plot_list_w[[1]] +
  plot_list_w[[2]] +
  plot_list_w[[3]] +
  plot_list_w[[4]] +
  plot_list_w[[5]]

final_plot_w

if (SAVE_OUTPUTS) {
  ggsave(
    "3_Outputs/Fig_PriceDensities_SimulatedStrikes_Wheat.pdf",
    plot = final_plot_w,
    width = 12,
    height = 8
  )
}

# ==============================================================================
# 5.0 Simulation of Utility of profits -----------------------------------------
# ==============================================================================

# ------------------------------------------------------------------------------
## 5.1 Functions for EUT, CPT, and Ambiguity calculations ----------------------
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
## 5.1.1 EUT -------------------------------------------------------------------
# ------------------------------------------------------------------------------

# Compute EUT certainty equivalents using Mean-Variance approximation
compute_eut <- function(profit_df, delta_values) {
  target_cols <- colnames(profit_df)
  results_list <- vector("list", length(delta_values) * length(target_cols))
  idx <- 0L
  for (d in delta_values) {
    for (col_name in target_cols) {
      idx <- idx + 1L
      profits <- profit_df[[col_name]]
      mean_profit <- mean(profits)
      var_profit <- var(profits)
      RP <- (d / 2) * var_profit
      CE <- mean_profit - RP
      results_list[[idx]] <- data.frame(
        Variable = col_name,
        Delta = d,
        CertaintyEquivalent = CE,
        RiskPremium = RP
      )
    }
  }
  bind_rows(results_list)
}

# ------------------------------------------------------------------------------
## 5.1.2 CPT -------------------------------------------------------------------
# ------------------------------------------------------------------------------

# Value function (Tversky & Kahneman 1992)
v <- function(x, a, l, ref) {
  ifelse(x >= ref, (x - ref)^a, -l * (-(x - ref))^a)
}

# Probability weighting function (Prelec 1998)
weighting_function <- function(p, g) {
  exp(-(-log(p))^g)
}

# Compute CPT values for all parameter combinations
compute_cpt <- function(
  profit_df,
  ref_point,
  a_values,
  lambda_values,
  gamma_values
) {
  target_cols <- colnames(profit_df)
  n_combos <- length(a_values) *
    length(lambda_values) *
    length(gamma_values) *
    length(target_cols)
  results_list <- vector("list", n_combos)
  idx <- 0L

  for (a in a_values) {
    for (l in lambda_values) {
      for (g in gamma_values) {
        for (col_name in target_cols) {
          idx <- idx + 1L
          profits <- sort(profit_df[[col_name]])

          n <- length(profits)
          probs <- rep(1 / n, n)
          df <- data.frame(profit = profits, prob = probs)

          losses <- df |>
            dplyr::filter(profit < ref_point) |>
            dplyr::arrange(profit)
          gains <- df |>
            dplyr::filter(profit >= ref_point) |>
            dplyr::arrange(dplyr::desc(profit))

          # Decision weights for losses
          if (nrow(losses) > 0) {
            losses <- losses |>
              dplyr::mutate(
                cum_p = cumsum(prob),
                dw = weighting_function(cum_p, g) -
                  weighting_function(cum_p - prob, g)
              )
          }

          # Decision weights for gains (weighted separately per TK1992)
          if (nrow(gains) > 0) {
            gains <- gains |>
              dplyr::mutate(
                cum_p = cumsum(prob),
                dw = weighting_function(cum_p, g) -
                  weighting_function(cum_p - prob, g)
              )
          }

          loss_value <- if (nrow(losses) > 0) {
            sum(v(losses$profit, a, l, ref_point) * losses$dw)
          } else {
            0
          }

          gain_value <- if (nrow(gains) > 0) {
            sum(v(gains$profit, a, l, ref_point) * gains$dw)
          } else {
            0
          }

          results_list[[idx]] <- data.frame(
            Variable = col_name,
            Alpha = a,
            Lambda = l,
            Gamma = g,
            CPT_Value = loss_value + gain_value
          )
        }
      }
    }
  }
  bind_rows(results_list)
}

# ------------------------------------------------------------------------------
## 5.1.3 Ambiguity -------------------------------------------------------------
# ------------------------------------------------------------------------------

# Compute alpha-EU certainty equivalents under ambiguity.
# p_intervals refer to Pr(low), i.e. the probability of the adverse state
# (profits < ref_point). Intervals should be centered on the empirical
# frequency of that state.
compute_alpha_eu <- function(
  profit_df,
  ref_point,
  delta_values,
  alpha_values,
  p_intervals
) {
  target_cols <- colnames(profit_df)
  results_list <- list()
  idx <- 0L

  for (d in delta_values) {
    for (col_name in target_cols) {
      profits <- profit_df[[col_name]]

      high <- profits[profits >= ref_point]
      low <- profits[profits < ref_point]

      # Expected CARA utility within each state
      eu_high <- if (length(high) > 0) mean(-exp(-d * high)) else NA
      eu_low <- if (length(low) > 0) mean(-exp(-d * low)) else NA

      # If all observations fall in one state, E^P[U] is constant in P
      # and alpha-EU collapses to standard EU
      if (is.na(eu_high) || is.na(eu_low)) {
        eu_single <- if (!is.na(eu_high)) eu_high else eu_low
        ce_single <- -log(-eu_single) / d

        for (a in alpha_values) {
          for (pi in p_intervals) {
            idx <- idx + 1L
            results_list[[idx]] <- data.frame(
              Variable = col_name,
              Delta = d,
              Alpha = a,
              P_interval = paste0("[", pi[1], ", ", pi[2], "]"),
              AlphaEU_CE = ce_single
            )
          }
        }
        next
      }

      for (a in alpha_values) {
        for (pi in p_intervals) {
          p_min <- pi[1]
          p_max <- pi[2]

          # p refers to Pr(low): EU is decreasing in p, so
          # p_max gives min EU (worst case), p_min gives max EU (best case)
          eu_at_pmin <- p_min * eu_low + (1 - p_min) * eu_high
          eu_at_pmax <- p_max * eu_low + (1 - p_max) * eu_high

          # alpha weights the pessimistic (min) prior
          alpha_eu_val <- a * eu_at_pmax + (1 - a) * eu_at_pmin

          ce <- -log(-alpha_eu_val) / d

          idx <- idx + 1L
          results_list[[idx]] <- data.frame(
            Variable = col_name,
            Delta = d,
            Alpha = a,
            P_interval = paste0("[", pi[1], ", ", pi[2], "]"),
            AlphaEU_CE = ce
          )
        }
      }
    }
  }
  bind_rows(results_list)
}

# ===============================================================================
## 5.2 Broad framing ============================================================
# ===============================================================================

# Parameters

# kg of N/ha applied
N <- 200
# tonnes of CAN/ha (contains 27% N)
CAN <- N / 0.27 / 1000

#-------------------------------------------------------------------------------
### 5.2.1 Simulate stochastic yields -------------------------------------------
#-------------------------------------------------------------------------------

# Get coeff from production function regression
# Only intercept + nitrogen coefficients are used, which implicitly predicts
# for rainfed + conventional BAU tillage (the reference categories, where
# irr = 0 and all tillage dummies = 0).
coeff_0_yield <- yield_function$coefficients[1]
coeff_1_yield <- yield_function$coefficients[2]
coeff_2_yield <- yield_function$coefficients[3]

# generate predicted yield
predicted_yield <- coeff_0_yield + coeff_1_yield * sqrt(N) + coeff_2_yield * N

# Get coeff from variation function regression
# Same reference scenario: rainfed + conventional BAU tillage.
coeff_0_varf <- variation_function$coefficients[1]
coeff_1_varf <- variation_function$coefficients[2]

# generate predicted variation
predicted_sigma_yield <- coeff_0_varf + coeff_1_varf * sqrt(N)

# derive stochastic error
stochastic_error <- array(dim = c(Ns, simulated_years))

for (t in 1:simulated_years) {
  # NOTE: Assumes normally distributed yield errors. If residuals from the
  # yield function are significantly skewed (e.g., left-skew from drought),
  # consider using a skewed distribution or adding a third-moment function.
  stochastic_error[, t] <- rnorm(n = Ns, mean = 0, sd = predicted_sigma_yield)
}

# simulation of stochastic yields
simulated_yield <- array(dim = c(Ns, simulated_years))

for (t in 1:simulated_years) {
  simulated_yield[, t] <- predicted_yield + stochastic_error[, t]
}

# plot simulated yields
df_sim_yield <- data.frame(yield = simulated_yield[, 1])

sim_yield_plot <- ggplot(df_sim_yield, aes(x = yield)) +
  geom_histogram() +
  theme_minimal() +
  labs(x = "Simulated yield (t/ha)", y = "Density")

sim_yield_plot

#-------------------------------------------------------------------------------
### 5.2.2 Get price simulations and calculate profits --------------------------
#-------------------------------------------------------------------------------

# get forward price
sim_f_p_forward <- rep(tail(f_data$Avg_Price, 1) * discount_factor, Ns)

# wheat forward price: last observed wheat price discounted to harvest
# w_data$Avg_Price is in eur/t (no /1000 conversion needed unlike fertilizer)
sim_w_p_forward <- rep(tail(w_data$Avg_Price, 1) * discount_factor_w, Ns)

# Calculate "no contract" profits
profit_nc <- (sim_w_p * simulated_yield) - (sim_f_p * CAN)

# Calculate "forward" profits
profit_fc <- (sim_w_p * simulated_yield) - (sim_f_p_forward * CAN)
profit_wfc <- (sim_w_p_forward * simulated_yield) - (sim_f_p * CAN)
profit_both_fc <- (sim_w_p_forward * simulated_yield) - (sim_f_p_forward * CAN)

# Calculate "adjusted contract" profits:

# Initialize result data frame
df_simulated_profit <- data.frame(
  nc = profit_nc,
  f_fc = profit_fc,
  w_fc = profit_wfc,
  both_fc = profit_both_fc
)

# Loop through each strike price and calculate adjusted contract profits
for (i in 1:nrow(premium_results)) {
  k <- premium_results$strike_price[i]
  contract_cost <- premium_results$option_premium[i]

  sim_f_p_adj_contr <- pmin(sim_f_p, k)

  profit_ac <- (sim_w_p * simulated_yield) -
    (sim_f_p_adj_contr * CAN) -
    (contract_cost * CAN)

  # Add result columns for each strike price
  col_name <- paste0("fc_", i)
  df_simulated_profit[[col_name]] <- profit_ac
}

# Wheat minimum-price contracts (wc_j): spot fertilizer + floored wheat revenue ----
for (j in 1:nrow(premium_results_w)) {
  k_w <- premium_results_w$strike_price[j]
  c_w <- premium_results_w$option_premium[j]

  profit_wc <- (pmax(sim_w_p, k_w) * simulated_yield) -
    (c_w * simulated_yield) - # put premium (eur/t) x yield (t/ha) = eur/ha
    (sim_f_p * CAN) # fertilizer at spot price

  col_name <- paste0("wc_", j)
  df_simulated_profit[[col_name]] <- profit_wc
}

# Combined fertilizer max-price + wheat min-price contracts (both_i, diagonal) ----
# Strike i for fertilizer is paired with strike i for wheat (rank-matched)
for (i in 1:nrow(premium_results)) {
  k_f <- premium_results$strike_price[i]
  c_f <- premium_results$option_premium[i]
  k_w <- premium_results_w$strike_price[i]
  c_w <- premium_results_w$option_premium[i]

  profit_both <- (pmax(sim_w_p, k_w) * simulated_yield) -
    (c_w * simulated_yield) - # wheat put premium (eur/ha)
    (pmin(sim_f_p, k_f) * CAN) - # fertilizer capped at k_f
    (c_f * CAN) # fertilizer call premium (eur/ha)

  col_name <- paste0("both_", i)
  df_simulated_profit[[col_name]] <- profit_both
}

#-------------------------------------------------------------------------------
## 5.3 Summary tables ----------------------------------------------------------
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
### 5.3.1 CAN price paid under different contracts -----------------------------
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
    mean = mean(sim_f_p_adj_contr),
    median = median(sim_f_p_adj_contr),
    sd = sd(sim_f_p_adj_contr),
    min = min(sim_f_p_adj_contr),
    max = max(sim_f_p_adj_contr)
  )
}

# combine into a single dataframe
price_stats_table <- bind_rows(
  base_row_sim_prices,
  base_row_forward_prices,
  stats_list
) |>
  mutate(across(where(is.numeric), ~ round(., 3)))

price_stats_table

if (SAVE_OUTPUTS) {
  write_csv(price_stats_table, "3_Outputs/Table_ResultsPriceStats.csv")
}

#-------------------------------------------------------------------------------
### 5.3.2 Total costs under different contracts --------------------------------
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
    mean = mean((sim_f_p_adj_contr * CAN) + (cost * CAN)),
    median = median((sim_f_p_adj_contr * CAN) + (cost * CAN)),
    sd = sd((sim_f_p_adj_contr * CAN) + (cost * CAN)),
    min = min((sim_f_p_adj_contr * CAN) + (cost * CAN)),
    max = max((sim_f_p_adj_contr * CAN) + (cost * CAN))
  )
}

# combine into a single dataframe
costs_stats_table <- bind_rows(
  base_row_costs,
  base_row_forward_costs,
  stats_list
) |>
  mutate(across(where(is.numeric), ~ round(., 2)))

costs_stats_table

if (SAVE_OUTPUTS) {
  write_csv(costs_stats_table, "3_Outputs/Table_ResultsCostsStats.csv")
}

#-------------------------------------------------------------------------------
### 5.3.3 Wheat price received under different contracts -----------------------
#-------------------------------------------------------------------------------

base_row_wheat_prices <- tibble(
  strike = NA,
  mean = mean(sim_w_p),
  median = median(sim_w_p),
  sd = sd(sim_w_p),
  min = min(sim_w_p),
  max = max(sim_w_p)
)

base_row_wheat_forward_prices <- tibble(
  strike = NA,
  mean = mean(sim_w_p_forward),
  median = median(sim_w_p_forward),
  sd = sd(sim_w_p_forward),
  min = min(sim_w_p_forward),
  max = max(sim_w_p_forward)
)

stats_list <- list()

for (i in 1:nrow(premium_results_w)) {
  k_w <- premium_results_w$strike_price[i]

  sim_w_p_adj_contr <- pmax(sim_w_p, k_w)

  stats_list[[i]] <- tibble(
    strike = k_w,
    mean = mean(sim_w_p_adj_contr),
    median = median(sim_w_p_adj_contr),
    sd = sd(sim_w_p_adj_contr),
    min = min(sim_w_p_adj_contr),
    max = max(sim_w_p_adj_contr)
  )
}

wheat_price_stats_table <- bind_rows(
  base_row_wheat_prices,
  base_row_wheat_forward_prices,
  stats_list
) |>
  mutate(across(where(is.numeric), ~ round(., 3)))

wheat_price_stats_table

if (SAVE_OUTPUTS) {
  write_csv(
    wheat_price_stats_table,
    "3_Outputs/Table_ResultsWheatPriceStats.csv"
  )
}

#-------------------------------------------------------------------------------
### 5.3.4 Total wheat revenue under different contracts ------------------------
#-------------------------------------------------------------------------------

base_row_wheat_revenue <- tibble(
  strike = NA,
  mean = mean(sim_w_p * simulated_yield),
  median = median(sim_w_p * simulated_yield),
  sd = sd(sim_w_p * simulated_yield),
  min = min(sim_w_p * simulated_yield),
  max = max(sim_w_p * simulated_yield)
)

base_row_wheat_forward_revenue <- tibble(
  strike = NA,
  mean = mean(sim_w_p_forward * simulated_yield),
  median = median(sim_w_p_forward * simulated_yield),
  sd = sd(sim_w_p_forward * simulated_yield),
  min = min(sim_w_p_forward * simulated_yield),
  max = max(sim_w_p_forward * simulated_yield)
)

stats_list <- list()

for (i in 1:nrow(premium_results_w)) {
  k_w <- premium_results_w$strike_price[i]
  c_w <- premium_results_w$option_premium[i]

  wheat_revenue <- pmax(sim_w_p, k_w) * simulated_yield - c_w * simulated_yield

  stats_list[[i]] <- tibble(
    strike = k_w,
    mean = mean(wheat_revenue),
    median = median(wheat_revenue),
    sd = sd(wheat_revenue),
    min = min(wheat_revenue),
    max = max(wheat_revenue)
  )
}

wheat_revenue_stats_table <- bind_rows(
  base_row_wheat_revenue,
  base_row_wheat_forward_revenue,
  stats_list
) |>
  mutate(across(where(is.numeric), ~ round(., 2)))

wheat_revenue_stats_table

if (SAVE_OUTPUTS) {
  write_csv(
    wheat_revenue_stats_table,
    "3_Outputs/Table_ResultsWheatRevenueStats.csv"
  )
}

#-------------------------------------------------------------------------------
### 5.3.5 Total profits under different contracts ------------------------------
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

if (SAVE_OUTPUTS) {
  write_csv(summary_profits, "3_Outputs/Table_ResultsProfitStats.csv")
}

#-------------------------------------------------------------------------------
## 5.4 Expected Utility Theory (EUT) values under broad framing ----------------
#-------------------------------------------------------------------------------

# define absolute risk aversion values
delta <- c(0.00001, 0.00005, 0.0001, 0.0005, 0.001, 0.005)

target_cols <- colnames(df_simulated_profit)

# CE and RP with Mean-Variance approximation: CE = E[π] - (δ/2)·Var[π]
# This is a second-order Taylor approximation of expected utility, valid when
# the distribution of profits is approximately normal or risk aversion is small.
eut_results_broadframing <- compute_eut(df_simulated_profit, delta)

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

if (SAVE_OUTPUTS) {
  write_csv(
    eut_results_broadframing,
    "3_Outputs/Table_Results_EUT_Broadframing.csv"
  )
}

#-------------------------------------------------------------------------------
## 5.5 Cumulative Prospect Theory (CPT) under broad framing --------------------
#-------------------------------------------------------------------------------

# risk aversion
a_values <- c(0.1, 0.334, 0.5)

# loss aversion
lambda <- c(0.5, 1.574, 2)

# probability weighting
gamma <- c(0.4, 0.571, 0.7)

# Reference point: expected profit under the status quo hedging strategy
# (forward contract/buy and hold strategy). This represents the farmer's natural
# benchmark for evaluating alternative risk management contracts.
ref_point_broad <- mean(df_simulated_profit$both_fc)

cpt_results_broadframing <- compute_cpt(
  df_simulated_profit,
  ref_point_broad,
  a_values,
  lambda,
  gamma
)

#-------------------------------------------------------------------------------
### 5.5.1 display CPT results --------------------------------------------------
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

if (SAVE_OUTPUTS) {
  write_csv(final_cpt_table, "3_Outputs/Table_Results_CPT_BroadFraming.csv")
}

#-------------------------------------------------------------------------------
## 5.6 Alpha-EU under broad framing --------------------------------------------
#-------------------------------------------------------------------------------

# Ambiguity aversion parameter (0 = full ambiguity loving, 1 = full aversion)
alpha_amb <- c(0.1, 0.5, 0.9)

# Risk aversion (CARA)
delta_amb <- c(0.00001, 0.00005, 0.0001, 0.0005, 0.001, 0.005)

# Empirical probability of the adverse state (profits < ref point)
p_hat <- mean(df_simulated_profit$both_fc < ref_point_broad)

# Belief bounds on Pr(low) centered on the empirical frequency
# Three levels of ambiguity: high (±0.20), medium (±0.10), low (±0.05)
p_intervals <- list(
  c(max(p_hat - 0.20, 0), min(p_hat + 0.20, 1)),
  c(max(p_hat - 0.10, 0), min(p_hat + 0.10, 1)),
  c(max(p_hat - 0.05, 0), min(p_hat + 0.05, 1))
)

alpha_eu_results <- compute_alpha_eu(
  df_simulated_profit,
  ref_point_broad,
  delta_amb,
  alpha_amb,
  p_intervals
)

# Pivot to wide format: one row per (Delta, Alpha, P_interval), columns = scenarios
alpha_eu_table <- alpha_eu_results |>
  pivot_wider(
    names_from = Variable,
    values_from = AlphaEU_CE,
    names_prefix = "CE_"
  ) |>
  mutate(across(where(is.numeric), ~ round(., 2)))

alpha_eu_table

if (SAVE_OUTPUTS) {
  write_csv(
    alpha_eu_table,
    "3_Outputs/Table_Results_AlphaEU_BroadFraming.csv"
  )
}

# ===============================================================================
# 6.0 Narrow Framing ===========================================================
# ===============================================================================

#-------------------------------------------------------------------------------
## 6.0.1 calculate profits ------------------------------------------------------
#-------------------------------------------------------------------------------

# forward contracts
profit_fc_nf <- CAN * (sim_f_p_forward - sim_f_p)
profit_wfc_nf <- predicted_yield * (sim_w_p_forward - sim_w_p)
profit_both_fc_nf <- CAN *
  (sim_f_p_forward - sim_f_p) +
  predicted_yield * (sim_w_p_forward - sim_w_p)

# adjusted contract
# Initialize result data frame
df_profit_narrowframing <- data.frame(
  f_fc = profit_fc_nf,
  w_fc = profit_wfc_nf,
  both_fc = profit_both_fc_nf
)

for (i in 1:nrow(premium_results)) {
  k <- premium_results$strike_price[i]
  contract_cost <- premium_results$option_premium[i]

  profit_ac <- CAN * (pmax(sim_f_p - k, 0) - contract_cost)

  # Add column with name fc_1, fc_2, etc.
  col_name <- paste0("fc_", i)
  df_profit_narrowframing[[col_name]] <- profit_ac
}

# Wheat put option: narrow framing payoff (per ha) ----
# Uses predicted_yield (fixed) to isolate price risk,
# consistent with CAN being fixed on the fertilizer side.
for (j in 1:nrow(premium_results_w)) {
  k_w <- premium_results_w$strike_price[j]
  c_w <- premium_results_w$option_premium[j]

  profit_wc_nf <- predicted_yield * (pmax(k_w - sim_w_p, 0) - c_w)

  col_name <- paste0("wc_", j)
  df_profit_narrowframing[[col_name]] <- profit_wc_nf
}

# Combined narrow framing: diagonal — fertilizer call i + wheat put i ----
for (i in 1:nrow(premium_results)) {
  k_f <- premium_results$strike_price[i]
  c_f <- premium_results$option_premium[i]
  k_w <- premium_results_w$strike_price[i]
  c_w <- premium_results_w$option_premium[i]

  profit_both_nf <- CAN *
    (pmax(sim_f_p - k_f, 0) - c_f) +
    predicted_yield * (pmax(k_w - sim_w_p, 0) - c_w)

  col_name <- paste0("both_", i)
  df_profit_narrowframing[[col_name]] <- profit_both_nf
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

if (SAVE_OUTPUTS) {
  write_csv(
    summary_profits_narrowframing,
    "3_Outputs/Table_Results_Profits_NarrowFraming.csv"
  )
}

#-------------------------------------------------------------------------------
## 6.1 Expected Utility Theory (EUT) under narrow framing ----------------------
#-------------------------------------------------------------------------------

delta <- c(0, 0.0005, 0.0007, 0.001, 0.003)

target_cols <- colnames(df_profit_narrowframing)

# CE and RP with Mean-Variance approximation: CE = E[π] - (δ/2)·Var[π]
eut_results_narrowframing <- compute_eut(df_profit_narrowframing, delta)

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

if (SAVE_OUTPUTS) {
  write_csv(
    final_eut_table_narrowframing,
    "3_Outputs/Table_Results_EUT_NarrowFraming.csv"
  )
}

#-------------------------------------------------------------------------------
## 6.2 Cumulative Prospect Theory (CPT) under narrow framing -------------------
#-------------------------------------------------------------------------------

ref_point_narrow <- 0

cpt_results_nf <- compute_cpt(
  df_profit_narrowframing,
  ref_point_narrow,
  a_values,
  lambda,
  gamma
)

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
final_cpt_narrowframing <- bind_rows(
  std_row,
  alpha_rows,
  lambda_rows,
  gamma_rows
) |>
  relocate(Scenario) |>
  mutate(across(where(is.numeric), ~ round(., 2)))

final_cpt_narrowframing

if (SAVE_OUTPUTS) {
  write_csv(
    final_cpt_narrowframing,
    "3_Outputs/Table_Results_CPT_NarrowFraming.csv"
  )
}

#-------------------------------------------------------------------------------
## 6.3 Alpha-EU under narrow framing -------------------------------------------
#-------------------------------------------------------------------------------

# Empirical probability of the adverse state (payoff < 0) under forward contract
p_hat_nf <- mean(df_profit_narrowframing$f_fc < ref_point_narrow)

# Belief bounds on Pr(low) centered on empirical frequency
p_intervals_nf <- list(
  c(max(p_hat_nf - 0.20, 0), min(p_hat_nf + 0.20, 1)),
  c(max(p_hat_nf - 0.10, 0), min(p_hat_nf + 0.10, 1)),
  c(max(p_hat_nf - 0.05, 0), min(p_hat_nf + 0.05, 1))
)

alpha_eu_results_nf <- compute_alpha_eu(
  df_profit_narrowframing,
  ref_point_narrow,
  delta_amb,
  alpha_amb,
  p_intervals_nf
)

# Pivot to wide format: one row per (Delta, Alpha, P_interval), columns = scenarios
alpha_eu_table_nf <- alpha_eu_results_nf |>
  pivot_wider(
    names_from = Variable,
    values_from = AlphaEU_CE,
    names_prefix = "CE_"
  ) |>
  mutate(across(where(is.numeric), ~ round(., 2)))

alpha_eu_table_nf

if (SAVE_OUTPUTS) {
  write_csv(
    alpha_eu_table_nf,
    "3_Outputs/Table_Results_AlphaEU_NarrowFraming.csv"
  )
}

# ===============================================================================
# 7.0 MULTI-FARM SIMULATION ====================================================
# ===============================================================================

# Link to location-specific yields from EPIC when we get the data. For now farm
# heterogeneity comes from varying nitrogen application rates (N).
n_farms <- 500
N_mean <- 200
N_sd <- 30
N_min <- 50
N_max <- 350

set.seed(5678)
N_farms <- rnorm(n_farms, mean = N_mean, sd = N_sd)
N_farms <- pmax(pmin(N_farms, N_max), N_min)

# -------------------------------------------------------------------------------
## 7.0.1 Vectorized farm parameters (vectors of length n_farms) -----------------
# -------------------------------------------------------------------------------

CAN_farms <- N_farms / 0.27 / 1000
pred_yield_farms <- coeff_0_yield +
  coeff_1_yield * sqrt(N_farms) +
  coeff_2_yield * N_farms
pred_sigma_farms <- coeff_0_varf + coeff_1_varf * sqrt(N_farms)

# -------------------------------------------------------------------------------
## 7.0.2 Vectorized yield simulation (Ns x n_farms matrix) ----------------------
# -------------------------------------------------------------------------------

# Each column is one farm's Ns (Number of simulations) yield draws
sim_yield_mat <- matrix(rnorm(Ns * n_farms), nrow = Ns, ncol = n_farms)
# each column gets farm-specific mean and sd
sim_yield_mat <- sweep(sim_yield_mat, 2, pred_sigma_farms, "*")
sim_yield_mat <- sweep(sim_yield_mat, 2, pred_yield_farms, "+")

# -------------------------------------------------------------------------------
## 7.0.3 Vectorized profit matrices (Ns x n_farms each) -------------------------
# -------------------------------------------------------------------------------

# Price vectors (drop the single-year matrix dimension to plain vectors)
sim_f_p_vec <- drop(sim_f_p)
sim_w_p_vec <- drop(sim_w_p)
sim_f_p_fwd_vec <- drop(sim_f_p_forward)

# Helper: outer products give Ns x n_farms matrices
fert_cost_spot <- outer(sim_f_p_vec, CAN_farms) # sim_f_p * CAN_f
fert_cost_fwd <- outer(sim_f_p_fwd_vec, CAN_farms) # sim_f_p_forward * CAN_f

# Wheat revenue: sim_w_p_vec (Ns) * sim_yield_mat (Ns x n_farms)
wheat_rev <- sim_w_p_vec * sim_yield_mat

# --- Broad framing profit matrices ---
n_strikes <- nrow(premium_results)
n_strikes_w <- nrow(premium_results_w)
# Total scenarios: nc, f_fc, w_fc, both_fc, fc_1..5, wc_1..5, both_1..5 = 19
sim_w_p_fwd_vec <- rep(sim_w_p_forward[1], Ns)
broad_names <- c(
  "nc",
  "f_fc",
  "w_fc",
  "both_fc",
  paste0("fc_", 1:n_strikes),
  paste0("wc_", 1:n_strikes_w),
  paste0("both_", 1:n_strikes)
)

# Store as a named list of Ns x n_farms matrices
wheat_rev_fwd <- sim_w_p_fwd_vec * sim_yield_mat
broad_mats <- list()
broad_mats[["nc"]] <- wheat_rev - fert_cost_spot
broad_mats[["f_fc"]] <- wheat_rev - fert_cost_fwd
broad_mats[["w_fc"]] <- wheat_rev_fwd - fert_cost_spot
broad_mats[["both_fc"]] <- wheat_rev_fwd - fert_cost_fwd

for (i in 1:n_strikes) {
  k <- premium_results$strike_price[i]
  cost <- premium_results$option_premium[i]
  fert_capped <- outer(pmin(sim_f_p_vec, k), CAN_farms)
  prem_cost <- outer(rep(cost, Ns), CAN_farms)
  broad_mats[[paste0("fc_", i)]] <- wheat_rev - fert_capped - prem_cost
}

for (j in 1:n_strikes_w) {
  k_w <- premium_results_w$strike_price[j]
  c_w <- premium_results_w$option_premium[j]
  wheat_floored <- pmax(sim_w_p_vec, k_w) * sim_yield_mat
  wheat_prem <- c_w * sim_yield_mat
  broad_mats[[paste0("wc_", j)]] <- wheat_floored - wheat_prem - fert_cost_spot
}

for (i in 1:n_strikes) {
  k_f <- premium_results$strike_price[i]
  c_f <- premium_results$option_premium[i]
  k_w <- premium_results_w$strike_price[i]
  c_w <- premium_results_w$option_premium[i]
  broad_mats[[paste0("both_", i)]] <-
    pmax(sim_w_p_vec, k_w) *
    sim_yield_mat -
    c_w * sim_yield_mat -
    outer(pmin(sim_f_p_vec, k_f), CAN_farms) -
    outer(rep(c_f, Ns), CAN_farms)
}

# --- Narrow framing profit matrices ---
narrow_names <- c(
  "f_fc",
  "w_fc",
  "both_fc",
  paste0("fc_", 1:n_strikes),
  paste0("wc_", 1:n_strikes_w),
  paste0("both_", 1:n_strikes)
)

# Forward contract payoffs
fwd_minus_spot <- sim_f_p_fwd_vec - sim_f_p_vec
wfwd_minus_spot <- sim_w_p_fwd_vec - sim_w_p_vec

narrow_mats <- list()
narrow_mats[["f_fc"]] <- outer(fwd_minus_spot, CAN_farms)
narrow_mats[["w_fc"]] <- outer(wfwd_minus_spot, pred_yield_farms)
narrow_mats[["both_fc"]] <- outer(fwd_minus_spot, CAN_farms) +
  outer(wfwd_minus_spot, pred_yield_farms)

for (i in 1:n_strikes) {
  k <- premium_results$strike_price[i]
  cost <- premium_results$option_premium[i]
  narrow_mats[[paste0("fc_", i)]] <-
    outer(pmax(sim_f_p_vec - k, 0) - cost, CAN_farms)
}

for (j in 1:n_strikes_w) {
  k_w <- premium_results_w$strike_price[j]
  c_w <- premium_results_w$option_premium[j]
  narrow_mats[[paste0("wc_", j)]] <-
    outer(pmax(k_w - sim_w_p_vec, 0) - c_w, pred_yield_farms)
}

for (i in 1:n_strikes) {
  k_f <- premium_results$strike_price[i]
  c_f <- premium_results$option_premium[i]
  k_w <- premium_results_w$strike_price[i]
  c_w <- premium_results_w$option_premium[i]
  narrow_mats[[paste0("both_", i)]] <-
    outer(pmax(sim_f_p_vec - k_f, 0) - c_f, CAN_farms) +
    outer(pmax(k_w - sim_w_p_vec, 0) - c_w, pred_yield_farms)
}

# ------------------------------------------------------------------------------
## 7.0.4 Vectorized EUT  -------------------------------------------------------
# ------------------------------------------------------------------------------

# EUT CE = mean - (delta/2) * var, computed directly from profit matrices.

compute_eut_vectorized <- function(
  profit_mats,
  scenario_names,
  delta_values,
  N_farms,
  n_farms
) {
  # profit_mats: named list of Ns x n_farms matrices
  # Returns a data frame with columns: Variable, Delta, CertaintyEquivalent,
  #   RiskPremium, farm_id, N
  col_means <- lapply(profit_mats[scenario_names], colMeans) # list of n_farms vectors
  col_vars <- lapply(profit_mats[scenario_names], function(m) apply(m, 2, var))

  results <- vector("list", length(delta_values) * length(scenario_names))
  idx <- 0L
  for (d in delta_values) {
    for (s in scenario_names) {
      idx <- idx + 1L
      mu <- col_means[[s]]
      vr <- col_vars[[s]]
      rp <- (d / 2) * vr
      ce <- mu - rp
      results[[idx]] <- data.frame(
        Variable = s,
        Delta = d,
        CertaintyEquivalent = ce,
        RiskPremium = rp,
        farm_id = 1:n_farms,
        N = N_farms
      )
    }
  }
  bind_rows(results)
}

delta_bf <- c(0.00001, 0.00005, 0.0001, 0.0005, 0.001, 0.005)
delta_nf <- c(0, 0.0005, 0.0007, 0.001, 0.003)

eut_broad_all <- compute_eut_vectorized(
  broad_mats,
  broad_names,
  delta_bf,
  N_farms,
  n_farms
)
eut_narrow_all <- compute_eut_vectorized(
  narrow_mats,
  narrow_names,
  delta_nf,
  N_farms,
  n_farms
)

# -------------------------------------------------------------------------------
## 7.0.5 Fast CPT for multi-farm ----------------------------------------
# -------------------------------------------------------------------------------

# CPT value for a single sorted numeric vector (base R only)
cpt_single <- function(profits_sorted, ref, a, l, g) {
  n <- length(profits_sorted)
  prob <- 1 / n

  # Split into losses (< ref) and gains (>= ref)
  idx_loss <- profits_sorted < ref
  n_loss <- sum(idx_loss)
  n_gain <- n - n_loss

  loss_val <- 0
  if (n_loss > 0) {
    # losses already sorted ascending; cumulative prob from worst
    cum_p <- (1:n_loss) * prob
    prev_p <- cum_p - prob
    dw <- exp(-(-log(cum_p))^g) - ifelse(prev_p > 0, exp(-(-log(prev_p))^g), 0)
    x_loss <- profits_sorted[idx_loss]
    loss_val <- sum(-l * (-(x_loss - ref))^a * dw)
  }

  gain_val <- 0
  if (n_gain > 0) {
    # gains: reverse order (best first) for decumulative weighting
    x_gain <- rev(profits_sorted[!idx_loss])
    cum_p <- (1:n_gain) * prob
    prev_p <- cum_p - prob
    dw <- exp(-(-log(cum_p))^g) - ifelse(prev_p > 0, exp(-(-log(prev_p))^g), 0)
    gain_val <- sum((x_gain - ref)^a * dw)
  }

  loss_val + gain_val
}

# Compute CPT across all farms for one scenario using its Ns x n_farms matrix
compute_cpt_multifarm <- function(
  profit_mats,
  scenario_names,
  ref_vecs,
  a_values,
  lambda_values,
  gamma_values,
  N_farms,
  n_farms
) {
  # ref_vecs: named list of length-n_farms vectors (one ref point per farm per scenario set)
  #           OR a single numeric vector of length n_farms applied to all scenarios
  n_combos <- length(a_values) *
    length(lambda_values) *
    length(gamma_values) *
    length(scenario_names)
  results <- vector("list", n_combos)
  idx <- 0L

  for (a in a_values) {
    for (l in lambda_values) {
      for (g in gamma_values) {
        for (s in scenario_names) {
          idx <- idx + 1L
          mat <- profit_mats[[s]] # Ns x n_farms
          cpt_vals <- numeric(n_farms)
          for (f in 1:n_farms) {
            profits_sorted <- sort(mat[, f])
            cpt_vals[f] <- cpt_single(profits_sorted, ref_vecs[f], a, l, g)
          }
          results[[idx]] <- data.frame(
            Variable = s,
            Alpha = a,
            Lambda = l,
            Gamma = g,
            CPT_Value = cpt_vals,
            farm_id = 1:n_farms,
            N = N_farms
          )
        }
      }
    }
  }
  bind_rows(results)
}

# -------------------------------------------------------------------------------
## 7.0.6 Alpha-EU for multi-farm -------------------------------------------
# -------------------------------------------------------------------------------

compute_alpha_eu_multifarm <- function(
  profit_mats,
  scenario_names,
  ref_vecs,
  delta_values,
  alpha_values,
  p_interval_offsets,
  N_farms,
  n_farms
) {
  # p_interval_offsets: vector of half-widths, e.g. c(0.20, 0.10, 0.05)
  # ref_vecs: numeric vector length n_farms
  results <- list()
  idx <- 0L

  for (d in delta_values) {
    for (s in scenario_names) {
      mat <- profit_mats[[s]] # Ns x n_farms

      # Pre-compute per-farm EU for high and low states
      eu_high <- numeric(n_farms)
      eu_low <- numeric(n_farms)
      n_high <- integer(n_farms)
      n_low <- integer(n_farms)

      for (f in 1:n_farms) {
        x <- mat[, f]
        ref <- ref_vecs[f]
        hi <- x[x >= ref]
        lo <- x[x < ref]
        n_high[f] <- length(hi)
        n_low[f] <- length(lo)
        eu_high[f] <- if (length(hi) > 0) mean(-exp(-d * hi)) else NA
        eu_low[f] <- if (length(lo) > 0) mean(-exp(-d * lo)) else NA
      }

      for (a in alpha_values) {
        for (hw in p_interval_offsets) {
          ce_vec <- numeric(n_farms)

          for (f in 1:n_farms) {
            if (is.na(eu_high[f]) || is.na(eu_low[f])) {
              eu_s <- if (!is.na(eu_high[f])) eu_high[f] else eu_low[f]
              ce_vec[f] <- -log(-eu_s) / d
            } else {
              p_hat_f <- n_low[f] / (n_low[f] + n_high[f])
              p_min <- max(p_hat_f - hw, 0)
              p_max <- min(p_hat_f + hw, 1)
              eu_at_pmin <- p_min * eu_low[f] + (1 - p_min) * eu_high[f]
              eu_at_pmax <- p_max * eu_low[f] + (1 - p_max) * eu_high[f]
              alpha_eu_val <- a * eu_at_pmax + (1 - a) * eu_at_pmin
              ce_vec[f] <- -log(-alpha_eu_val) / d
            }
          }

          p_label <- paste0("[", hw, "]")
          idx <- idx + 1L
          results[[idx]] <- data.frame(
            Variable = s,
            Delta = d,
            Alpha = a,
            P_interval = p_label,
            AlphaEU_CE = ce_vec,
            farm_id = 1:n_farms,
            N = N_farms
          )
        }
      }
    }
  }
  bind_rows(results)
}

# -------------------------------------------------------------------------------
## 7.0.7 Run CPT and Alpha-EU across all farms ----------------------------------
# -------------------------------------------------------------------------------

# Broad framing reference points: per-farm mean of fc profits
ref_broad_farms <- colMeans(broad_mats[["both_fc"]])

cpt_broad_all <- compute_cpt_multifarm(
  broad_mats,
  broad_names,
  ref_broad_farms,
  a_values,
  lambda,
  gamma,
  N_farms,
  n_farms
)

alphaeu_broad_all <- compute_alpha_eu_multifarm(
  broad_mats,
  broad_names,
  ref_broad_farms,
  delta_amb,
  alpha_amb,
  c(0.20, 0.10, 0.05),
  N_farms,
  n_farms
)

# Narrow framing reference point: 0 for all farms
ref_narrow_farms <- rep(0, n_farms)

cpt_narrow_all <- compute_cpt_multifarm(
  narrow_mats,
  narrow_names,
  ref_narrow_farms,
  a_values,
  lambda,
  gamma,
  N_farms,
  n_farms
)

alphaeu_narrow_all <- compute_alpha_eu_multifarm(
  narrow_mats,
  narrow_names,
  ref_narrow_farms,
  delta_amb,
  alpha_amb,
  c(0.20, 0.10, 0.05),
  N_farms,
  n_farms
)

# -------------------------------------------------------------------------------
## 7.2 Summary statistics across farms -------------------------------------------
# -------------------------------------------------------------------------------

eut_broad_summary <- eut_broad_all |>
  group_by(Delta, Variable) |>
  summarise(
    mean_CE = mean(CertaintyEquivalent),
    sd_CE = sd(CertaintyEquivalent),
    q05_CE = quantile(CertaintyEquivalent, 0.05),
    q25_CE = quantile(CertaintyEquivalent, 0.25),
    median_CE = median(CertaintyEquivalent),
    q75_CE = quantile(CertaintyEquivalent, 0.75),
    q95_CE = quantile(CertaintyEquivalent, 0.95),
    .groups = "drop"
  ) |>
  mutate(across(where(is.numeric), ~ round(., 2)))

cpt_broad_summary <- cpt_broad_all |>
  group_by(Alpha, Lambda, Gamma, Variable) |>
  summarise(
    mean_CPT = mean(CPT_Value),
    sd_CPT = sd(CPT_Value),
    q05_CPT = quantile(CPT_Value, 0.05),
    q25_CPT = quantile(CPT_Value, 0.25),
    median_CPT = median(CPT_Value),
    q75_CPT = quantile(CPT_Value, 0.75),
    q95_CPT = quantile(CPT_Value, 0.95),
    .groups = "drop"
  ) |>
  mutate(across(where(is.numeric), ~ round(., 2)))

alphaeu_broad_summary <- alphaeu_broad_all |>
  group_by(Delta, Alpha, P_interval, Variable) |>
  summarise(
    mean_CE = mean(AlphaEU_CE),
    sd_CE = sd(AlphaEU_CE),
    q05_CE = quantile(AlphaEU_CE, 0.05),
    q25_CE = quantile(AlphaEU_CE, 0.25),
    median_CE = median(AlphaEU_CE),
    q75_CE = quantile(AlphaEU_CE, 0.75),
    q95_CE = quantile(AlphaEU_CE, 0.95),
    .groups = "drop"
  ) |>
  mutate(across(where(is.numeric), ~ round(., 2)))

eut_narrow_summary <- eut_narrow_all |>
  group_by(Delta, Variable) |>
  summarise(
    mean_CE = mean(CertaintyEquivalent),
    sd_CE = sd(CertaintyEquivalent),
    q05_CE = quantile(CertaintyEquivalent, 0.05),
    q25_CE = quantile(CertaintyEquivalent, 0.25),
    median_CE = median(CertaintyEquivalent),
    q75_CE = quantile(CertaintyEquivalent, 0.75),
    q95_CE = quantile(CertaintyEquivalent, 0.95),
    .groups = "drop"
  ) |>
  mutate(across(where(is.numeric), ~ round(., 2)))

cpt_narrow_summary <- cpt_narrow_all |>
  group_by(Alpha, Lambda, Gamma, Variable) |>
  summarise(
    mean_CPT = mean(CPT_Value),
    sd_CPT = sd(CPT_Value),
    q05_CPT = quantile(CPT_Value, 0.05),
    q25_CPT = quantile(CPT_Value, 0.25),
    median_CPT = median(CPT_Value),
    q75_CPT = quantile(CPT_Value, 0.75),
    q95_CPT = quantile(CPT_Value, 0.95),
    .groups = "drop"
  ) |>
  mutate(across(where(is.numeric), ~ round(., 2)))

alphaeu_narrow_summary <- alphaeu_narrow_all |>
  group_by(Delta, Alpha, P_interval, Variable) |>
  summarise(
    mean_CE = mean(AlphaEU_CE),
    sd_CE = sd(AlphaEU_CE),
    q05_CE = quantile(AlphaEU_CE, 0.05),
    q25_CE = quantile(AlphaEU_CE, 0.25),
    median_CE = median(AlphaEU_CE),
    q75_CE = quantile(AlphaEU_CE, 0.75),
    q95_CE = quantile(AlphaEU_CE, 0.95),
    .groups = "drop"
  ) |>
  mutate(across(where(is.numeric), ~ round(., 2)))

# -------------------------------------------------------------------------------
## 7.3 Save outputs --------------------------------------------------------------
# -------------------------------------------------------------------------------

if (SAVE_OUTPUTS) {
  write_csv(eut_broad_summary, "3_Outputs/Table_MultiFarm_EUT_BroadFraming.csv")
  write_csv(cpt_broad_summary, "3_Outputs/Table_MultiFarm_CPT_BroadFraming.csv")
  write_csv(
    alphaeu_broad_summary,
    "3_Outputs/Table_MultiFarm_AlphaEU_BroadFraming.csv"
  )
  write_csv(
    eut_narrow_summary,
    "3_Outputs/Table_MultiFarm_EUT_NarrowFraming.csv"
  )
  write_csv(
    cpt_narrow_summary,
    "3_Outputs/Table_MultiFarm_CPT_NarrowFraming.csv"
  )
  write_csv(
    alphaeu_narrow_summary,
    "3_Outputs/Table_MultiFarm_AlphaEU_NarrowFraming.csv"
  )
}

# -------------------------------------------------------------------------------
## 7.4 Distribution plots -------------------------------------------------------
# -------------------------------------------------------------------------------
library(ggridges)
# EUT certainty equivalents under broad framing — one delta value
plot_eut_broad_dist <- eut_broad_all |>
  filter(Delta == 0.0001) |>
  ggplot(aes(x = CertaintyEquivalent, y = Variable)) +
  geom_density_ridges_gradient(alpha = 0.4) +
  labs(
    title = "EUT Certainty Equivalents across farms (broad framing, delta = 0.0001)",
    x = "Certainty Equivalent (eur/ha)",
    y = "Density"
  ) +
  theme_minimal()

plot_eut_broad_dist

# CPT values under broad framing — standard assumptions
plot_cpt_broad_dist <- cpt_broad_all |>
  filter(Alpha == 0.334, Lambda == 1.574, Gamma == 0.571) |>
  ggplot(aes(x = CPT_Value, y = Variable)) +
  geom_density_ridges(alpha = 0.4) +
  labs(
    title = "CPT Values across farms (broad framing, standard assumptions)",
    x = "CPT Value",
    y = "Density"
  ) +
  theme_minimal()

plot_cpt_broad_dist

# EUT certainty equivalents under narrow framing — one delta value
plot_eut_narrow_dist <- eut_narrow_all |>
  filter(Delta == 0.001) |>
  ggplot(aes(x = CertaintyEquivalent, y = Variable)) +
  geom_density_ridges(alpha = 0.4) +
  labs(
    title = "EUT Certainty Equivalents across farms (narrow framing, delta = 0.001)",
    x = "Certainty Equivalent (eur/ha)",
    y = "Density"
  ) +
  theme_ridges()

plot_eut_narrow_dist

# CPT values under narrow framing — standard assumptions
plot_cpt_narrow_dist <- cpt_narrow_all |>
  filter(Alpha == 0.334, Lambda == 1.574, Gamma == 0.571) |>
  ggplot(aes(x = CPT_Value, y = Variable)) +
  geom_density_ridges(alpha = 0.4) +
  labs(
    title = "CPT Values across farms (narrow framing, standard assumptions)",
    x = "CPT Value",
    y = "Density"
  ) +
  theme_minimal()

plot_cpt_narrow_dist

if (SAVE_OUTPUTS) {
  ggsave(
    "3_Outputs/Fig_MultiFarm_EUT_BroadFraming.pdf",
    plot = plot_eut_broad_dist,
    width = 14,
    height = 10
  )
  ggsave(
    "3_Outputs/Fig_MultiFarm_CPT_BroadFraming.pdf",
    plot = plot_cpt_broad_dist,
    width = 14,
    height = 10
  )
  ggsave(
    "3_Outputs/Fig_MultiFarm_EUT_NarrowFraming.pdf",
    plot = plot_eut_narrow_dist,
    width = 14,
    height = 10
  )
  ggsave(
    "3_Outputs/Fig_MultiFarm_CPT_NarrowFraming.pdf",
    plot = plot_cpt_narrow_dist,
    width = 14,
    height = 10
  )
}
