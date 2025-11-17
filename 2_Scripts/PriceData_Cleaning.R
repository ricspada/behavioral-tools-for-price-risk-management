rm(list=ls())

library(readxl)
library(tidyverse)
library(RColorBrewer)

#####################################################################
# Fertilizer Price Data #############################################
#####################################################################

# Get the files and sheets names
fert_old_file <- "1_Data/RawData/PriceData/Düngemittelabfrage Altdaten - KAB - ab 2007 bis August 2021.xlsx"
fert_new_file <- "1_Data/RawData/PriceData/Düngemittelabfrage Neudaten - KAB - ab September 2021 bis September 2025.xlsx"

fert_old_sheets <- excel_sheets(fert_old_file)
fert_new_sheets <- excel_sheets(fert_new_file)

# Setting variables names 
eng_colnames <- c("Year", "Date", "Product", "Min_Price", "Max_Price", "Avg_Price")

# process old file, adding avg column
old_dfs <- list() # initialize list

for (sheet in fert_old_sheets) {
  df <- read_excel(fert_old_file, sheet = sheet, skip = 2)
  df <- df[, 1:5]   # keep only first 5 columns
  colnames(df) <- eng_colnames # renaming columns
  df$Avg_Price <- (df$Min_Price + df$Max_Price) / 2 # calculate avg values
  old_dfs[[sheet]] <- df
}

# process new file
new_dfs <- list()

for (sheet in fert_new_sheets) {
  df <- read_excel(fert_new_file, sheet = sheet, skip = 2)
  df <- df[, 1:6]   # keep only first 6 columns
  colnames(df) <- eng_colnames
  new_dfs[[sheet]] <- df
}

#change name, modify date format
old_dfs <- old_dfs %>%
  set_names(~ case_when(
    . == "KAS ab März 2007" ~ "CAN", # Calcium Ammonium Nitrate (CAN), 27% N
    . == "AHL ab März 2007" ~ "UAN", # Urea Ammonium Nitrate (UAN), 28% N
    . == "ASS ab Mai 2007"  ~ "ASN", # Ammonium Sulphate Nitrate (ASN),26% N + 13% S
    . == "SSA ab August 2009" ~ "AS", # Ammonium Sulphate (AS),21% N + 24% S
    . == "Harnstoff gesch. ab Febr. 2011" ~ "urea_P", # Urea 46% N, granulated, protected (Urea_P), 46% N
    . == "Harnstoff gekörnt ab März 2007" ~ "urea", # Urea 46% N, granulated, (Urea), 46% N
    . == "DAP ab März 2007" ~ "DAP", # Diammonium phosphate, (DAP) 18% N + 46% P2O5
    . == "TSP ab März 2007" ~ "TSP", # Triple superphosphate (TSP) 45 % P2O5
    . == "Kornkali ab März 2007" ~ "Kornkali" # Kornkali + Mg, 38 % K2O + 6 % MgO
  )) %>%
  map(~ .x %>% mutate(Date = as.Date(Date, format = "%d.%m.%Y"))) %>% 
  imap(~ .x %>% mutate(Product = .y))

new_dfs <- new_dfs %>%
  set_names(~ case_when(
    . == "KAS" ~ "CAN", # Calcium Ammonium Nitrate (CAN), 27% N
    . == "AHL" ~ "UAN", # Urea Ammonium Nitrate (UAN), 28% N
    . == "ASS"  ~ "ASN", # Ammonium Sulphate Nitrate (ASN),26% N + 13% S
    . == "SSA" ~ "AS", # Ammonium Sulphate (AS),21% N + 24% S
    . == "Harnstoff gek. gesch." ~ "urea_P", # Urea 46% N, granulated, protected (Urea_P), 46% N
    . == "Harnstoff gek." ~ "urea", # Urea 46% N, granulated, (Urea), 46% N
    . == "DAP" ~ "DAP", # Diammonium phosphate, (DAP) 18% N + 46% P2O5
    . == "TSP" ~ "TSP", # Triple superphosphate (TSP) 45 % P2O5
    . == "Kornkali" ~ "Kornkali" # Kornkali + Mg, 38 % K2O + 6 % MgO
  )) %>%
  map(~ .x %>% mutate(Date = as.Date(Date, format = "%d.%m.%Y"))) %>% 
  imap(~ .x %>% mutate(Product = .y))

# combine in a single df all the data
fert_df <- map2_dfr(old_dfs, new_dfs, rbind)

#drop the Year column (redundant data)
fert_df$Year <- NULL

### Plot fert data ########################################################

# Create a named vector for titles
product_titles <- c(
  "UAN" = "Urea Ammonium Nitrate (UAN)",
  "CAN" = "Calcium Ammonium Nitrate (CAN)",
  "ASN" = "Ammonium Sulphate Nitrate (ASN)",
  "AS" = "Ammonium Sulphate (AS)",
  "urea_P" = "Urea protected (Urea_P)",
  "urea" = "Urea (Urea)",
  "DAP" = "Diammonium phosphate (DAP)",
  "TSP" = "Triple superphosphate (TSP)",
  "Kornkali" = "Kornkali")

ggplot(fert_df, aes(x = Date, color = Product, fill = Product)) +
  geom_ribbon(aes(ymin = Min_Price, ymax = Max_Price), alpha = 0.2, color = NA) +
  geom_line(aes(y = Avg_Price), linewidth = 0.8) +
  facet_wrap(~ Product, scales = "fixed", ncol = 3,
             labeller = labeller(Product = product_titles)) +
  labs(y = "Price (€/100 kg)",
       x = "Date") +
  theme_minimal() +
  theme(legend.position = "none")

# Save fertilizers price  data 
# write.csv(fert_df, "1_Data/FertilizerPrices.csv", row.names = FALSE)

###############################################################################
### Crop data #################################################################
###############################################################################

# Get the files and sheets names
crop_old_file <- "1_Data/RawData/PriceData/Getreideabfrage Altdaten - KAB - ab 2000 bis August 2021.xlsx"                
crop_new_file <- "1_Data/RawData/PriceData/Getreideabfrage Neudaten - KAB - ab September 2021 bis September 2025.xlsx" 

crop_old_sheets <- excel_sheets(crop_old_file)
crop_new_sheets <- excel_sheets(crop_old_file)

# process old file, adding avg column
old_dfs <- list() # emptying list
for (sheet in crop_old_sheets) {
  df <- read_excel(crop_old_file, sheet = sheet, skip = 2)
  df <- df[, 1:5]   # keep only first 5 columns
  colnames(df) <- eng_colnames # renaming columns
  df$Avg_Price <- (df$Min_Price + df$Max_Price) / 2 # calculate avg values
  old_dfs[[sheet]] <- df
}

# process new file
new_dfs <- list()
for (sheet in crop_new_sheets) {
  df <- read_excel(crop_new_file, sheet = sheet, skip = 2)
  df <- df[, 1:6]   # keep only first 6 columns
  colnames(df) <- eng_colnames
  new_dfs[[sheet]] <- df
}

#change name, modify date format
old_dfs <- old_dfs %>%
  set_names(~ case_when(
    . == "Futterweizen" ~ "Feed_Wheat", 
    . == "Brotweizen" ~ "Bread_Wheat")) %>%
  map(~ .x %>% mutate(Date = as.Date(Date, format = "%d.%m.%Y"))) %>% 
  imap(~ .x %>% mutate(Product = .y))

new_dfs <- new_dfs %>%
  set_names(~ case_when(
    . == "Futterweizen" ~ "Feed_Wheat", 
    . == "Brotweizen" ~ "Bread_Wheat")) %>%
  map(~ .x %>% mutate(Date = as.Date(Date, format = "%d.%m.%Y"))) %>% 
  imap(~ .x %>% mutate(Product = .y))

# combine in a single df all the data
crop_df <- map2_dfr(old_dfs, new_dfs, rbind)

#drop the Year column (redundant data)
crop_df$Year <- NULL

# save crop price data
# write.csv(crop_df, "1_Data/WheatPrices.csv", row.names = FALSE)


### Plot Crop data ############################################################

ggplot(crop_df, aes(x = Date, color = Product, fill = Product)) +
  geom_ribbon(aes(ymin = Min_Price, ymax = Max_Price), alpha = 0.2, color = NA) +
  geom_line(aes(y = Avg_Price), linewidth = 0.8) +
  scale_color_brewer(palette = "Dark2") +
  scale_fill_brewer(palette = "Set2") +
  labs(y = "Price (€/100 kg)",
       x = "Date") +
  theme_minimal() +
  theme(legend.position = "bottom")


###############################################################################
### Resample to regular time intervals ########################################
###############################################################################

f_data <- fert_df[fert_df[[2]] == "CAN", ] # Keep only CAN data
w_data <- crop_df[crop_df[[2]] == "Bread_Wheat", ]

# Check time ranges
range(f_data$Date)
range(w_data$Date)

# Check time intervals (see if regularly spaced)
diff(f_data$Date) 
diff(w_data$Date)  


f_df <- f_data |> 
  rename(f_p = Avg_Price) |> 
  dplyr::select(Date, f_p) |>
  filter(Date >= as.Date("2009-01-01"))           # keep observations from 2009 onwards

w_df <- w_data |> 
  rename(w_p = Avg_Price) |> 
  dplyr::select(Date, w_p) |>
  filter(Date >= as.Date("2009-01-01"))           # keep observations from 2009 onwards


# Find optimal time interval that minimizes distance to observations ##############################################################

find_optimal_interval <- function(dates, search_range = c(14, 21)) {
  
  # Convert dates to days since first observation
  start_date <- dates[1]
  days_from_start <- as.numeric(difftime(dates, start_date, units = "days"))
  
  # Objective function: total squared distance to nearest grid point
  objective <- function(interval) {
    total_error <- 0
    for (day in days_from_start) {
      # Find nearest multiple of interval
      nearest_multiple <- round(day / interval) * interval
      error <- abs(day - nearest_multiple)
      total_error <- total_error + error^2
    }
    return(total_error)
  }
  
  # Optimize to find best interval
  result <- optimize(
    f = objective,
    interval = search_range,
    maximum = FALSE
  )
  
  optimal_interval <- round(result$minimum)
  
  # Generate dates based on optimal interval
  max_days <- max(days_from_start)
  num_periods <- floor(max_days / optimal_interval) + 2 # extra periods to cover range
  new_times <- start_date + (0:(num_periods - 1)) * optimal_interval
  
  list(
    optimal_interval_days = optimal_interval,
    new_times = new_times
  )
}


result <- find_optimal_interval(f_df$Date)

result$optimal_interval_days
result$new_times


length(result$new_times)
length(f_df$Date)

diff(result$new_times)


# Map CAN observations to new dates based on biweekly proximity ###############################################################

#' @param df Data frame with Date, f_p, and w_p columns
#' @param new_times Vector of grid dates
#' @param biweek_threshold Number of days to consider "same bi-week" (+/- 7 days)
#' @return Data frame with grid_date, f_p (NA if no matching observation)
map_f_to_grid <- function(df, new_times, biweek_threshold = 7, interpolate_na = TRUE) {
  
  # Initialize result data frame
  result <- data.frame(
    grid_date = new_times,
    f_p = NA_real_
  )
  
  # For each new date, find if there's a matching observation
  for (i in seq_along(new_times)) {
    grid_date <- new_times[i]
    
    # Calculate days difference from this grid date to all observations
    days_diff <- abs(as.numeric(difftime(df$Date, grid_date, units = "days")))
    
    # Find observations within the biweek threshold
    within_week <- which(days_diff <= biweek_threshold)
    
    # If there are matching observations, take the closest one
    if (length(within_week) > 0) {
      closest_idx <- within_week[which.min(days_diff[within_week])]
      result$f_p[i] <- df$f_p[closest_idx]
    }
    # Otherwise, f_p and w_p remain NA (already initialized)
  }
 
  # Interpolate NA values in f_p 
  if (interpolate_na && any(is.na(result$f_p))) {
    # Convert dates to numeric for interpolation
    df_numeric <- as.numeric(df$Date)
    grid_numeric <- as.numeric(result$grid_date)
    
    # Identify which values need interpolation
    na_indices <- which(is.na(result$f_p))
    
    # Interpolate f_p for NA values using linear interpolation
    interpolated_f_p <- approx(
      x = df_numeric,
      y = df$f_p,
      xout = grid_numeric[na_indices],
      rule = 2  # Use nearest value for points outside range
    )$y
    
    # Fill in the interpolated values
    result$f_p[na_indices] <- interpolated_f_p
  }

  return(result)
}


# Map wheat observations to new dates based on weekly proximity ###############################################################


#' @param df Data frame with Date, f_p, and w_p columns
#' @param new_times Vector of grid dates
#' @param week_threshold Number of days to consider "same week" (+/- 3.5 days)
#' @return Data frame with grid_date, f_p (NA if no matching observation)
map_w_to_grid <- function(df, new_times, week_threshold = 3.5) {
  
  # Initialize result data frame
  result <- data.frame(
    grid_date = new_times,
    w_p = NA_real_
  )
  
  # For each new date, find if there's a matching observation
  for (i in seq_along(new_times)) {
    grid_date <- new_times[i]
    
    # Calculate days difference from this grid date to all observations
    days_diff <- abs(as.numeric(difftime(df$Date, grid_date, units = "days")))
    
    # Find observations within the week threshold
    within_week <- which(days_diff <= week_threshold)
    
    # If there are matching observations, take the closest one
    if (length(within_week) > 0) {
      closest_idx <- within_week[which.min(days_diff[within_week])]
      result$w_p[i] <- df$w_p[closest_idx]
    }
    # Otherwise, f_p and w_p remain NA (already initialized)
  }
  
  return(result)
}


# Map data to new times ###############################################################
dat_f <- map_f_to_grid(f_df, result$new_times, biweek_threshold = 7)
dat_w <- map_w_to_grid(w_df, result$new_times, week_threshold = 3.5)

# combine datasets
combined_df <- merge(dat_f, dat_w, by = "grid_date")

# Save combined data
# write.csv(combined_df, "1_Data/Adjusted_Prices.csv", row.names = FALSE)
