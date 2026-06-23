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
library(ggridges)
library(fixest)

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

sim_yield <- read_csv("1_Data/Germany_EPIC.csv.gz")

# keep winter wheat only, and exclude alternative weather scenario pathways
# convert dry matter yield to economic yield (12% moisture content)
sim_yield <- sim_yield |>
  filter(
    CROP == "WWHT",
    WTH == "hist"
  ) |>
  mutate(
    YLD = YLD_DM / 0.88,
    FTN = round(FTN)
  )

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

# full descriptive statistics
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

selected_summary <- summary |>
  filter(FTN %in% c(10, 50, 100, 180, 250))

if (SAVE_OUTPUTS) {
  write_csv(selected_summary, "3_Outputs/Table_YieldSummary.csv")
}

#-------------------------------------------------------------------------------
## 1.2 Just & Pope production function estimation ------------------------------
#-------------------------------------------------------------------------------

# Weather and rotation dummies are excluded
# Fixed-effects model: SimUID absorbs location-level heterogeneity
yield_function <- feols(
  YLD ~ sqrt(FTN) +
    FTN +
    irr +
    til_min_cons +
    til_contill_r00 +
    til_contill_r30 +
    til_contill_r60 +
    til_contill_r90 |
    SimUID,
  data = sim_yield,
  vcov = "hetero"
)

summary(yield_function)

# append residuals
sim_yield$residuals <- residuals(yield_function)
sim_yield$abs_residuals <- abs(sim_yield$residuals)

# Variation function estimation
variation_function <- feols(
  abs_residuals ~ sqrt(FTN) +
    FTN +
    irr +
    til_min_cons +
    til_contill_r00 +
    til_contill_r30 +
    til_contill_r60 +
    til_contill_r90 |
    SimUID,
  data = sim_yield,
  vcov = "hetero"
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

# plot prices

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
  labs(y = bquote(Delta[1] ~ "CAN"), x = NULL) +
  theme_minimal()

plot_dif_w <- ggplot(dif_dat, aes(x = Date)) +
  geom_line(aes(y = w_p), linewidth = 0.5, color = colors[2]) +
  labs(y = bquote(Delta[1] ~ "Wheat"), x = NULL) +
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


#-------------------------------------------------------------------------------
## 2.4 Summary statistics of prices --------------------------------------------
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

# estimate models
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
var_comparison$AIC <- round(as.numeric(var_comparison$AIC), 2)
var_comparison$BIC <- round(as.numeric(var_comparison$BIC), 2)
var_comparison$BG_pval <- round(as.numeric(var_comparison$BG_pval), 3)

var_comparison

# Best by BIC: VAR(1) + squared lags + D_shock
#
# Selecting VAR(1) + D_shock
# as QVAR altready capture non-linear dynamics

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
mf <- f_p ~ f_p1 + w_p1 + D_shock
mw <- w_p ~ f_p1 + w_p1 + D_shock

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
D_shock_hist <- var_dat$D_shock

X <- cbind(1, f_p1, w_p1, D_shock_hist)

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

# zero initial change (mean of stationary differences ≈ 0)
Ys_changes[1, 1, ] <- 0
Ys_changes[2, 1, ] <- 0

# set initial level to sample mean (geometric mean of raw prices)
f_p_start <- mean(dat$f_p)
w_p_start <- mean(dat$w_p)
Ys_levels[1, 1, ] <- f_p_start
Ys_levels[2, 1, ] <- w_p_start

# Simulation based on marginal distributions
for (it in 2:Nt) {
  for (is in 1:Ns) {
    f_change1 <- Ys_changes[1, it - 1, is]
    w_change1 <- Ys_changes[2, it - 1, is]

    # D_shock = 0 post-shock
    Xfi <- cbind(1, f_change1, w_change1, 0)
    Xwi <- cbind(1, f_change1, w_change1, 0)

    # Draw fertilizer marginal quantile, then wheat conditional on it via copula
    qi_f <- sample(1:99, 1)
    u_f <- taus_dense[qi_f]
    v <- sample(1:99, 1)
    F_w_pred <- coeff_rC[1, v] + coeff_rC[2, v] * u_f
    F_w_pred <- max(taus_dense[1], min(taus_dense[99], F_w_pred))
    qi_w <- which.min(abs(taus_dense - F_w_pred))

    yf_change <- Xfi %*% coeff_f[, qi_f]
    yw_change <- Xwi %*% coeff_w[, qi_w]

    Ys_changes[1, it, is] <- yf_change
    Ys_changes[2, it, is] <- yw_change

    # convert to level using cumulative sum
    Ys_levels[1, it, is] <- f_p_start + sum(Ys_changes[1, 1:it, is])
    Ys_levels[2, it, is] <- w_p_start + sum(Ys_changes[2, 1:it, is])
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

# Historical prices in level form (same space as simulation output)
dat_level <- dat |>
  mutate(f_p_level = exp(f_p), w_p_level = exp(w_p))

# Future date axis: 24 biweekly periods from last historical date
last_date <- max(dat$Date)
sim_dates <- seq(last_date, by = "2 weeks", length.out = Nt + 1)[-1]

# Per-period quantile summary of simulated paths
sim_fan <- simulated_prices |>
  group_by(time) |>
  summarise(
    q10_f = quantile(price_f_level, 0.10),
    q25_f = quantile(price_f_level, 0.25),
    q50_f = quantile(price_f_level, 0.50),
    q75_f = quantile(price_f_level, 0.75),
    q90_f = quantile(price_f_level, 0.90),
    q10_w = quantile(price_w_level, 0.10),
    q25_w = quantile(price_w_level, 0.25),
    q50_w = quantile(price_w_level, 0.50),
    q75_w = quantile(price_w_level, 0.75),
    q90_w = quantile(price_w_level, 0.90),
    .groups = "drop"
  ) |>
  mutate(Date = sim_dates)

# Shared y-axis limits per commodity
ylim_f <- range(c(dat_level$f_p_level, simulated_prices$price_f_level))
ylim_w <- range(c(dat_level$w_p_level, simulated_prices$price_w_level))

# Fertilizer plot

p_ts_f <- ggplot() +
  geom_ribbon(
    data = sim_fan,
    aes(x = Date, ymin = q10_f, ymax = q90_f),
    fill = "#E41A1C",
    alpha = 0.15
  ) +
  geom_ribbon(
    data = sim_fan,
    aes(x = Date, ymin = q25_f, ymax = q75_f),
    fill = "#E41A1C",
    alpha = 0.25
  ) +
  geom_line(
    data = sim_fan,
    aes(x = Date, y = q50_f),
    color = "#E41A1C",
    linewidth = 0.5
  ) +
  geom_line(
    data = dat_level,
    aes(x = Date, y = f_p_level),
    color = "darkblue",
    linewidth = 0.5
  ) +
  coord_cartesian(ylim = ylim_f) +
  labs(x = NULL, y = "CAN price (€/t)") +
  theme_minimal() +
  theme(panel.grid.minor = element_blank())

p_dens_f <- ggplot() +
  geom_density(
    data = simulated_prices,
    aes(x = price_f_level, fill = "Simulated"),
    alpha = 0.5
  ) +
  geom_density(
    data = dat_level,
    aes(x = f_p_level, fill = "Historical"),
    alpha = 0.4
  ) +
  scale_fill_manual(
    name = NULL,
    values = c("Simulated" = "#E41A1C", "Historical" = "darkblue")
  ) +
  coord_flip(xlim = ylim_f) +
  labs(x = "Density", y = NULL) +
  theme_minimal() +
  theme(
    panel.grid.minor = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    legend.position = c(0.8, 0.85),
    legend.background = element_rect(fill = "white", color = "black")
  )

fig_f <- p_ts_f | p_dens_f

# Wheat plot

p_ts_w <- ggplot() +
  geom_ribbon(
    data = sim_fan,
    aes(x = Date, ymin = q10_w, ymax = q90_w),
    fill = "#E41A1C",
    alpha = 0.15
  ) +
  geom_ribbon(
    data = sim_fan,
    aes(x = Date, ymin = q25_w, ymax = q75_w),
    fill = "#E41A1C",
    alpha = 0.25
  ) +
  geom_line(
    data = sim_fan,
    aes(x = Date, y = q50_w),
    color = "#E41A1C",
    linewidth = 0.5
  ) +
  geom_line(
    data = dat_level,
    aes(x = Date, y = w_p_level),
    color = "darkgreen",
    linewidth = 0.5
  ) +
  coord_cartesian(ylim = ylim_w) +
  labs(x = NULL, y = "Wheat price (€/t)") +
  theme_minimal() +
  theme(panel.grid.minor = element_blank())

p_dens_w <- ggplot() +
  geom_density(
    data = simulated_prices,
    aes(x = price_w_level, fill = "Simulated"),
    alpha = 0.5
  ) +
  geom_density(
    data = dat_level,
    aes(x = w_p_level, fill = "Historical"),
    alpha = 0.4
  ) +
  scale_fill_manual(
    name = NULL,
    values = c("Simulated" = "#E41A1C", "Historical" = "darkgreen")
  ) +
  coord_flip(xlim = ylim_w) +
  labs(x = "Density", y = NULL) +
  theme_minimal() +
  theme(
    panel.grid.minor = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    legend.position = c(0.8, 0.85),
    legend.background = element_rect(fill = "white", color = "black")
  )

fig_w <- p_ts_w | p_dens_w

plot_sim_prices <- fig_f / fig_w

plot_sim_prices

if (SAVE_OUTPUTS) {
  ggsave(
    "3_Outputs/Fig_SimulationFan.pdf",
    plot_sim_prices,
    width = 10,
    height = 4
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
# 5.0 MULTI-FARM SIMULATION ----------------------------------------------------
# ==============================================================================

# ------------------------------------------------------------------------------
## 5.0.0 Parameters needed by multi-farm simulation ----------------------------
# ------------------------------------------------------------------------------

# Get nitrogen coefficients from production function.
fe_yield <- fixef(yield_function)$SimUID
fe_varf <- fixef(variation_function)$SimUID

coeff_1_yield <- coef(yield_function)["sqrt(FTN)"]
coeff_2_yield <- coef(yield_function)["FTN"]

coeff_1_varf <- coef(variation_function)["sqrt(FTN)"]
coeff_2_varf <- coef(variation_function)["FTN"]

# Forward prices: risk-neutral expectation E[S_T] from the simulation
sim_f_p_forward <- rep(mean(sim_f_p), Ns)
sim_w_p_forward <- rep(mean(sim_w_p), Ns)

n_farms <- length(unique(sim_yield$SimUID))

# Joint bootstrap: sample (SimUID, FTN) pairs from EPIC rows
set.seed(5678)
boot_rows <- sim_yield[sample(nrow(sim_yield), n_farms, replace = TRUE), ]
N_farms <- boot_rows$FTN
fe_yield_farms <- fe_yield[as.character(boot_rows$SimUID)]
fe_varf_farms <- fe_varf[as.character(boot_rows$SimUID)]

# convert from N to tonnes of CAN (assuming 27% N content)
CAN_farms <- N_farms / 0.27 / 1000

# -------------------------------------------------------------------------------
## 5.0.1 Vectorized farm parameters (vectors of length n_farms) -----------------
# -------------------------------------------------------------------------------

# Implicitly predicts for rainfed + conventional BAU tillage
# (irr = 0, all tillage dummies = 0).
pred_yield_farms <- fe_yield_farms +
  coeff_1_yield * sqrt(N_farms) +
  coeff_2_yield * N_farms
pred_sigma_farms <- fe_varf_farms +
  coeff_1_varf * sqrt(N_farms) +
  coeff_2_varf * N_farms

# -------------------------------------------------------------------------------
## 5.0.2 Vectorized yield simulation (Ns x n_farms matrix) ----------------------
# -------------------------------------------------------------------------------

# Each column is one farm's Ns (Number of simulations) yield draws
sim_yield_mat <- matrix(rnorm(Ns * n_farms), nrow = Ns, ncol = n_farms)
# each column gets farm-specific mean and sd
sim_yield_mat <- sweep(sim_yield_mat, 2, pred_sigma_farms, "*")
sim_yield_mat <- sweep(sim_yield_mat, 2, pred_yield_farms, "+")

# -------------------------------------------------------------------------------
## 5.0.3 Vectorized profit matrices (Ns x n_farms each) -------------------------
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
## 5.0.4 Profits summary statistics --------------------------------
# ------------------------------------------------------------------------------

# Broad framing: for each scenario, pool all (Ns x n_farms) profit values
broad_profit_summary <- map_dfr(broad_names, function(s) {
  x <- as.vector(broad_mats[[s]])
  tibble(
    scenario = s,
    mean = mean(x),
    median = median(x),
    sd = sd(x),
    min = min(x),
    max = max(x),
    skewness = moments::skewness(x),
    q05 = quantile(x, 0.05),
    q25 = quantile(x, 0.25),
    q75 = quantile(x, 0.75),
    q95 = quantile(x, 0.95)
  )
}) |>
  mutate(across(where(is.numeric), ~ round(., 2)))

# Narrow framing
narrow_profit_summary <- map_dfr(narrow_names, function(s) {
  x <- as.vector(narrow_mats[[s]])
  tibble(
    scenario = s,
    mean = mean(x),
    median = median(x),
    sd = sd(x),
    min = min(x),
    max = max(x),
    skewness = moments::skewness(x),
    q05 = quantile(x, 0.05),
    q25 = quantile(x, 0.25),
    q75 = quantile(x, 0.75),
    q95 = quantile(x, 0.95)
  )
}) |>
  mutate(across(where(is.numeric), ~ round(., 2)))

if (SAVE_OUTPUTS) {
  write_csv(
    broad_profit_summary,
    "3_Outputs/Table_MultiFarm_ProfitStats_BroadFraming.csv"
  )
  write_csv(
    narrow_profit_summary,
    "3_Outputs/Table_MultiFarm_ProfitStats_NarrowFraming.csv"
  )
}

# -------------------------------------------------------------------------------
## 5.4d Across-farm distribution of per-farm profit statistics -------------------
# -------------------------------------------------------------------------------

farm_stat_df <- function(mats, names, framing) {
  bind_rows(lapply(names, function(s) {
    m <- mats[[s]]
    bind_rows(
      data.frame(value = apply(m, 2, min), statistic = "Min", scenario = s),
      data.frame(value = colMeans(m), statistic = "Mean", scenario = s),
      data.frame(value = apply(m, 2, max), statistic = "Max", scenario = s)
    )
  })) |>
    mutate(
      # compute group from character scenario before reassigning to a factor
      group = factor(group_scenario(scenario), levels = group_order),
      scenario = factor(scenario, levels = rev(broad_names)),
      statistic = factor(statistic, levels = c("Min", "Mean", "Max")),
      framing = framing
    )
}

df_stat_broad <- farm_stat_df(broad_mats, broad_names, "Broad framing")
df_stat_narrow <- farm_stat_df(narrow_mats, narrow_names, "Narrow framing")

mk_stat_plot <- function(df, xlab) {
  ggplot(df, aes(x = value, y = scenario, fill = group)) +
    geom_density_ridges(alpha = 0.4, scale = 1.1) +
    facet_wrap(~statistic, scales = "free_x") +
    scale_fill_brewer(palette = "Set2", drop = FALSE) +
    labs(x = xlab, y = NULL, fill = NULL) +
    theme_minimal()
}

plot_stat_broad <- mk_stat_plot(df_stat_broad, "Profit (€/ha)") +
  labs(title = "Broad framing")
plot_stat_narrow <- mk_stat_plot(df_stat_narrow, "Contract payoff (€/ha)") +
  labs(title = "Narrow framing")

plot_stat_broad
plot_stat_narrow

# -------------------------------------------------------------------------------
## 5.4c Single-farm profit distributions ----------------------------------------
# -------------------------------------------------------------------------------

farm_idx <- 3
narrow_names <- setdiff(broad_names, "nc")

# Contract-type grouping for colour coding
group_order <- c(
  "No contract",
  "Forward",
  "Fert. call",
  "Wheat put",
  "Both options"
)
group_scenario <- function(s) {
  case_when(
    s == "nc" ~ "No contract",
    s %in% c("f_fc", "w_fc", "both_fc") ~ "Forward",
    grepl("^fc_", s) ~ "Fert. call",
    grepl("^wc_", s) ~ "Wheat put",
    grepl("^both_", s) ~ "Both options"
  )
}

# Stack the single farm's column across scenarios into a long data frame
single_farm_df <- function(mats, names, framing) {
  bind_rows(lapply(names, function(s) {
    data.frame(
      value = mats[[s]][, farm_idx],
      scenario = s,
      group = group_scenario(s),
      framing = framing
    )
  })) |>
    mutate(
      # reverse levels so the first scenario sits at the top of the ridge
      scenario = factor(scenario, levels = rev(broad_names)),
      group = factor(group, levels = group_order)
    )
}

df_broad <- single_farm_df(broad_mats, broad_names, "Broad framing")
df_narrow <- single_farm_df(narrow_mats, narrow_names, "Narrow framing")

plot_broad_sf <- ggplot(df_broad, aes(x = value, y = scenario, fill = group)) +
  geom_density_ridges(alpha = 0.4, scale = 1.1) +
  scale_fill_brewer(palette = "Set2", drop = FALSE) +
  labs(title = "Broad framing", x = "Profit (€/ha)", y = NULL, fill = NULL) +
  theme_minimal()

plot_narrow_sf <- ggplot(
  df_narrow,
  aes(x = value, y = scenario, fill = group)
) +
  geom_density_ridges(alpha = 0.4, scale = 1.1) +
  scale_fill_brewer(palette = "Set2", drop = FALSE) +
  labs(
    title = "Narrow framing",
    x = "Contract payoff (€/ha)",
    y = NULL,
    fill = NULL
  ) +
  theme_minimal()

plot_broad_sf
plot_narrow_sf


# parameters
delta_bf <- c(0.0005, 0.001, 0.002, 0.003, 0.004)
delta_nf <- delta_bf

a_values <- c(0.1, 0.334, 0.5)
lambda <- c(0.5, 1.574, 2)
gamma <- c(0.4, 0.571, 0.7)

alpha_amb <- c(0.1, 0.5, 0.9)

# Ambiguity half-widths: probability ± around each farm's own p_hat.
# hw=0.3 → [p-0.3, p+0.3] (high ambiguity), hw=0.2 → medium, hw=0.1 → low
amb_halfwidths <- c(0.3, 0.2, 0.1)

# Broad framing reference points: per-farm mean profit under no contract (status quo)
ref_broad_farms <- colMeans(broad_mats[["nc"]])

# Narrow framing reference point: 0 for all farms
ref_narrow_farms <- rep(0, n_farms)


# ------------------------------------------------------------------------------
## 5.0.4 Compute EUT  ----------------------------------------------------------
# ------------------------------------------------------------------------------

# EUT CE = mean - (delta/2) * var

compute_eut_vectorized <- function(
  profit_mats,
  scenario_names,
  delta_values,
  N_farms,
  n_farms
) {
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


# ------------------------------------------------------------------------------
## 5.0.5 Compute CPT  ----------------------------------------------------------
# ------------------------------------------------------------------------------

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
## 5.0.6 Alpha-EU  -------------------------------------------------------------
# -------------------------------------------------------------------------------

compute_alpha_eu_multifarm <- function(
  profit_mats,
  scenario_names,
  ref_vecs,
  delta_values,
  alpha_values,
  amb_halfwidths,
  N_farms,
  n_farms
) {
  # ref_vecs: numeric vector length n_farms
  results <- list()
  idx <- 0L

  for (d in delta_values) {
    for (s in scenario_names) {
      mat <- profit_mats[[s]] # Ns x n_farms
      Ns_obs <- nrow(mat)

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

      # Per-farm estimated loss probability
      p_hat_all <- n_low / Ns_obs # vector length n_farms

      for (a in alpha_values) {
        for (hw in amb_halfwidths) {
          # Farm-specific ambiguity interval: p_hat[f] ± hw, clipped to [0,1].
          # hw is a probability half-width (not a quantile index), representing
          # the agent's subjective uncertainty about the true loss probability.
          p_min_vec <- pmax(0, p_hat_all - hw)
          p_max_vec <- pmin(1, p_hat_all + hw)
          ce_vec <- numeric(n_farms)

          for (f in 1:n_farms) {
            if (is.na(eu_high[f]) || is.na(eu_low[f])) {
              eu_s <- if (!is.na(eu_high[f])) eu_high[f] else eu_low[f]
              ce_vec[f] <- -log(-eu_s) / d
            } else {
              eu_at_pmin <- p_min_vec[f] *
                eu_low[f] +
                (1 - p_min_vec[f]) * eu_high[f]
              eu_at_pmax <- p_max_vec[f] *
                eu_low[f] +
                (1 - p_max_vec[f]) * eu_high[f]
              alpha_eu_val <- a * eu_at_pmax + (1 - a) * eu_at_pmin
              ce_vec[f] <- -log(-alpha_eu_val) / d
            }
          }

          p_label <- paste0("hw=", hw)
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

# ------------------------------------------------------------------------------
## 5.0.7 Run EUT, CPT, and Alpha-EU across all farms ---------------------------
# ------------------------------------------------------------------------------

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
  delta_bf,
  alpha_amb,
  amb_halfwidths,
  N_farms,
  n_farms
)

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
  delta_nf,
  alpha_amb,
  amb_halfwidths,
  N_farms,
  n_farms
)

# -------------------------------------------------------------------------------
## 5.2 Summary statistics across farms -------------------------------------------
# -------------------------------------------------------------------------------

eut_broad_summary <- eut_broad_all |>
  group_by(Delta, Variable) |>
  summarise(
    mean_CE = mean(CertaintyEquivalent),
    # RP: risk premium in €/ha (mean profit - CE)
    RP = mean(RiskPremium),
    sd_CE = sd(CertaintyEquivalent),
    q05_CE = quantile(CertaintyEquivalent, 0.05),
    q25_CE = quantile(CertaintyEquivalent, 0.25),
    median_CE = median(CertaintyEquivalent),
    q75_CE = quantile(CertaintyEquivalent, 0.75),
    q95_CE = quantile(CertaintyEquivalent, 0.95),
    .groups = "drop"
  ) |>
  mutate(across(where(is.numeric) & !any_of("Delta"), ~ round(., 2)))

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
  mutate(across(where(is.numeric) & !any_of("Delta"), ~ round(., 2)))

eut_narrow_summary <- eut_narrow_all |>
  group_by(Delta, Variable) |>
  summarise(
    mean_CE = mean(CertaintyEquivalent),
    RP = mean(RiskPremium),
    sd_CE = sd(CertaintyEquivalent),
    q05_CE = quantile(CertaintyEquivalent, 0.05),
    q25_CE = quantile(CertaintyEquivalent, 0.25),
    median_CE = median(CertaintyEquivalent),
    q75_CE = quantile(CertaintyEquivalent, 0.75),
    q95_CE = quantile(CertaintyEquivalent, 0.95),
    .groups = "drop"
  ) |>
  mutate(across(where(is.numeric) & !any_of("Delta"), ~ round(., 2)))

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
  mutate(across(where(is.numeric) & !any_of("Delta"), ~ round(., 2)))

# -------------------------------------------------------------------------------
## 5.3 Save outputs --------------------------------------------------------------
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
## 5.4 Distribution plots -------------------------------------------------------
# -------------------------------------------------------------------------------

# EUT certainty equivalents under broad framing
plot_eut_broad_dist <- eut_broad_all |>
  filter(Delta == delta_bf[3]) |>
  ggplot(aes(x = CertaintyEquivalent, y = Variable)) +
  geom_density_ridges(alpha = 0.4) +
  labs(
    title = "EUT",
    x = NULL,
    y = NULL
  ) +
  theme_minimal()


# CPT values under broad framing — standard assumptions
plot_cpt_broad_dist <- cpt_broad_all |>
  filter(Alpha == 0.334, Lambda == 1.574, Gamma == 0.571) |>
  ggplot(aes(x = CPT_Value, y = Variable)) +
  geom_density_ridges(alpha = 0.4) +
  labs(
    title = "CPT",
    x = NULL,
    y = NULL
  ) +
  theme_minimal()

# EUT certainty equivalents under narrow framing
plot_eut_narrow_dist <- eut_narrow_all |>
  filter(Delta == delta_nf[3]) |>
  ggplot(aes(x = CertaintyEquivalent, y = Variable)) +
  geom_density_ridges(alpha = 0.4) +
  labs(
    x = NULL,
    y = NULL
  ) +
  theme_ridges()


# CPT values under narrow framing — standard assumptions
plot_cpt_narrow_dist <- cpt_narrow_all |>
  filter(Alpha == 0.334, Lambda == 1.574, Gamma == 0.571) |>
  ggplot(aes(x = CPT_Value, y = Variable)) +
  geom_density_ridges(alpha = 0.4) +
  labs(
    x = NULL,
    y = NULL
  ) +
  theme_minimal()


# Alpha-EU certainty equivalents under broad framing
plot_alphaeu_broad_dist <- alphaeu_broad_all |>
  filter(Delta == delta_bf[3], Alpha == 0.9, P_interval == "hw=0.2") |>
  ggplot(aes(x = AlphaEU_CE, y = Variable)) +
  geom_density_ridges(alpha = 0.4) +
  labs(
    title = "Alpha-EU",
    x = NULL,
    y = NULL
  ) +
  theme_minimal()


# Alpha-EU certainty equivalents under narrow framing
plot_alphaeu_narrow_dist <- alphaeu_narrow_all |>
  filter(Delta == delta_nf[3], Alpha == 0.9, P_interval == "hw=0.2") |>
  ggplot(aes(x = AlphaEU_CE, y = Variable)) +
  geom_density_ridges(alpha = 0.4) +
  labs(
    x = NULL,
    y = NULL
  ) +
  theme_ridges()

utility_plots <- ((plot_eut_broad_dist + labs(y = "Broad Framing")) +
  plot_cpt_broad_dist +
  plot_alphaeu_broad_dist) /
  ((plot_eut_narrow_dist + labs(y = "Narrow Framing")) +
    plot_cpt_narrow_dist +
    plot_alphaeu_narrow_dist)

utility_plots

if (SAVE_OUTPUTS) {
  ggsave(
    "3_Outputs/Fig_UtilityPlots.pdf",
    plot = utility_plots,
    width = 14,
    height = 10
  )
}

# ==============================================================================
# 5.4b SENSITIVITY ANALYSIS PLOTS ==============================================
# ==============================================================================

# EUT — broad framing: one facet per delta, scenarios on y-axis
plot_eut_broad_sens <- eut_broad_summary |>
  mutate(Delta = factor(Delta)) |>
  ggplot(aes(x = mean_CE, y = Variable)) +
  geom_point() +
  facet_wrap(~Delta, labeller = label_both, scales = "free_x") +
  labs(title = "EUT - Broad Framing", x = "Mean CE (€)", y = NULL) +
  theme_minimal()

# EUT — narrow framing
plot_eut_narrow_sens <- eut_narrow_summary |>
  mutate(Delta = factor(Delta)) |>
  ggplot(aes(x = mean_CE, y = Variable)) +
  geom_point() +
  facet_wrap(~Delta, labeller = label_both, scales = "free_x") +
  labs(title = "EUT - Narrow Framing", x = "Mean CE (€)", y = NULL) +
  theme_minimal()

# CPT — broad framing: facet_grid(Lambda ~ Gamma), colour by Alpha
plot_cpt_broad_sens <- cpt_broad_summary |>
  mutate(
    Alpha = factor(Alpha),
    Lambda = factor(Lambda),
    Gamma = factor(Gamma)
  ) |>
  ggplot(aes(x = mean_CPT, y = Variable, colour = Alpha, group = Alpha)) +
  geom_point(position = position_dodge(width = 0.6)) +
  facet_grid(Lambda ~ Gamma, labeller = label_both) +
  labs(
    title = "CPT - Broad Framing",
    x = "Mean CPT Value",
    y = NULL,
    colour = "α"
  ) +
  theme_minimal()

# CPT — narrow framing
plot_cpt_narrow_sens <- cpt_narrow_summary |>
  mutate(
    Alpha = factor(Alpha),
    Lambda = factor(Lambda),
    Gamma = factor(Gamma)
  ) |>
  ggplot(aes(x = mean_CPT, y = Variable, colour = Alpha, group = Alpha)) +
  geom_point(position = position_dodge(width = 0.6)) +
  facet_grid(Lambda ~ Gamma, labeller = label_both) +
  labs(
    title = "CPT - Narrow Framing",
    x = "Mean CPT Value",
    y = NULL,
    colour = "α"
  ) +
  theme_minimal()

# Alpha-EU — broad framing: facet_grid(Alpha ~ P_interval), colour by Delta
plot_alphaeu_broad_sens <- alphaeu_broad_summary |>
  mutate(Alpha = factor(Alpha), Delta = factor(Delta)) |>
  ggplot(aes(x = mean_CE, y = Variable, colour = Delta, group = Delta)) +
  geom_point(position = position_dodge(width = 0.6)) +
  facet_grid(Alpha ~ P_interval, labeller = label_both) +
  scale_colour_brewer(palette = "Blues", direction = 1) +
  labs(
    title = "Alpha-EU - Broad Framing",
    x = "Mean CE (€)",
    y = NULL,
    colour = "δ"
  ) +
  theme_minimal()

# Alpha-EU — narrow framing
plot_alphaeu_narrow_sens <- alphaeu_narrow_summary |>
  mutate(Alpha = factor(Alpha), Delta = factor(Delta)) |>
  ggplot(aes(x = mean_CE, y = Variable, colour = Delta, group = Delta)) +
  geom_point(position = position_dodge(width = 0.6)) +
  facet_grid(Alpha ~ P_interval, labeller = label_both) +
  scale_colour_brewer(palette = "Blues", direction = 1) +
  labs(
    title = "Alpha-EU - Narrow Framing",
    x = "Mean CE (€)",
    y = NULL,
    colour = "δ"
  ) +
  theme_minimal()

plot_eut_broad_sens
plot_eut_narrow_sens
plot_cpt_broad_sens
plot_cpt_narrow_sens
plot_alphaeu_broad_sens
plot_alphaeu_narrow_sens

sens_plots <- (plot_eut_broad_sens +
  plot_cpt_broad_sens +
  plot_alphaeu_broad_sens) /
  (plot_eut_narrow_sens + plot_cpt_narrow_sens + plot_alphaeu_narrow_sens)

if (SAVE_OUTPUTS) {
  ggsave(
    "3_Outputs/Fig_Sensitivity_Plots.pdf",
    plot = sens_plots,
    width = 14,
    height = 10
  )
}

# ==============================================================================
# 5.5 FORMAL STATISTICAL TESTS =================================================
# ==============================================================================

# One representative parameter slice per model used throughout this section.
DELTA_EUT_B <- delta_bf[3]
DELTA_AEU_B <- delta_bf[3]
ALPHA_AEU <- 0.5
P_AEU <- "hw=0.3"


# -------------------------------------------------------------------------------
## paired Wilcoxon signed-rank (each contract vs "nc") ----
# -------------------------------------------------------------------------------

nc_eut <- eut_broad_all |>
  filter(Variable == "nc", Delta == DELTA_EUT_B) |>
  arrange(farm_id) |>
  pull(CertaintyEquivalent)

nc_aeu <- alphaeu_broad_all |>
  filter(
    Variable == "nc",
    Delta == DELTA_AEU_B,
    Alpha == ALPHA_AEU,
    P_interval == P_AEU
  ) |>
  arrange(farm_id) |>
  pull(AlphaEU_CE)

wilcox_broad <- map_dfr(setdiff(broad_names, "nc"), function(s) {
  eut_ce <- eut_broad_all |>
    filter(Variable == s, Delta == DELTA_EUT_B) |>
    arrange(farm_id) |>
    pull(CertaintyEquivalent)
  aeu_ce <- alphaeu_broad_all |>
    filter(
      Variable == s,
      Delta == DELTA_AEU_B,
      Alpha == ALPHA_AEU,
      P_interval == P_AEU
    ) |>
    arrange(farm_id) |>
    pull(AlphaEU_CE)

  te <- wilcox.test(eut_ce, nc_eut, paired = TRUE, alternative = "greater")
  ta <- wilcox.test(aeu_ce, nc_aeu, paired = TRUE, alternative = "greater")

  tibble(
    scenario = s,
    p_EUT = te$p.value,
    p_AlphaEU = ta$p.value
  )
})

wilcox_broad

# -------------------------------------------------------------------------------
## 5.5.2 Broad framing: paired Wilcoxon: option vs forward baseline -----------
# -------------------------------------------------------------------------------

# alternative = "two.sided": direction not predetermined EUT agents may prefer
# forwards (greater variance reduction); behavioural agents may prefer options.

# Forward baselines: one CE vector per farm per model
fwd_baselines <- list(
  fc = list(
    eut = eut_broad_all |>
      filter(Variable == "f_fc", Delta == DELTA_EUT_B) |>
      arrange(farm_id) |>
      pull(CertaintyEquivalent),
    aeu = alphaeu_broad_all |>
      filter(
        Variable == "f_fc",
        Delta == DELTA_AEU_B,
        Alpha == ALPHA_AEU,
        P_interval == P_AEU
      ) |>
      arrange(farm_id) |>
      pull(AlphaEU_CE)
  ),
  wc = list(
    eut = eut_broad_all |>
      filter(Variable == "w_fc", Delta == DELTA_EUT_B) |>
      arrange(farm_id) |>
      pull(CertaintyEquivalent),
    aeu = alphaeu_broad_all |>
      filter(
        Variable == "w_fc",
        Delta == DELTA_AEU_B,
        Alpha == ALPHA_AEU,
        P_interval == P_AEU
      ) |>
      arrange(farm_id) |>
      pull(AlphaEU_CE)
  ),
  both = list(
    eut = eut_broad_all |>
      filter(Variable == "both_fc", Delta == DELTA_EUT_B) |>
      arrange(farm_id) |>
      pull(CertaintyEquivalent),
    aeu = alphaeu_broad_all |>
      filter(
        Variable == "both_fc",
        Delta == DELTA_AEU_B,
        Alpha == ALPHA_AEU,
        P_interval == P_AEU
      ) |>
      arrange(farm_id) |>
      pull(AlphaEU_CE)
  )
)

wilcox_vs_fwd <- map_dfr(c("fc", "wc", "both"), function(prefix) {
  map_dfr(1:5, function(i) {
    s <- paste0(prefix, "_", i)

    eut_ce <- eut_broad_all |>
      filter(Variable == s, Delta == DELTA_EUT_B) |>
      arrange(farm_id) |>
      pull(CertaintyEquivalent)
    aeu_ce <- alphaeu_broad_all |>
      filter(
        Variable == s,
        Delta == DELTA_AEU_B,
        Alpha == ALPHA_AEU,
        P_interval == P_AEU
      ) |>
      arrange(farm_id) |>
      pull(AlphaEU_CE)

    te <- wilcox.test(
      eut_ce,
      fwd_baselines[[prefix]]$eut,
      paired = TRUE,
      alternative = "greater"
    )
    ta <- wilcox.test(
      aeu_ce,
      fwd_baselines[[prefix]]$aeu,
      paired = TRUE,
      alternative = "greater"
    )

    tibble(
      scenario = s,
      p_EUT = te$p.value,
      p_AlphaEU = ta$p.value
    )
  })
})

wilcox_vs_fwd


# -------------------------------------------------------------------------------
## 5.5.5 Save test results -------------------------------------------------------
# -------------------------------------------------------------------------------

if (SAVE_OUTPUTS) {
  write_csv(wilcox_broad, "3_Outputs/Table_WilcoxBroad.csv")
  write_csv(wilcox_vs_fwd, "3_Outputs/Table_WilcoxVsForward.csv")
}
