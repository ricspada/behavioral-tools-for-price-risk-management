rm(list=ls())

library(readxl)
library(tidyverse)
library(RColorBrewer)

#####################################################################
### Fertilizer Data #################################################
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

###########################################################################
### Plot fert data ########################################################
###########################################################################

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
write.csv(fert_df, "1_Data/FertilizerPrices.csv", row.names = FALSE)

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
write.csv(crop_df, "1_Data/WheatPrices.csv", row.names = FALSE)

###############################################################################
### Plot Crop data ############################################################
###############################################################################

ggplot(crop_df, aes(x = Date, color = Product, fill = Product)) +
  geom_ribbon(aes(ymin = Min_Price, ymax = Max_Price), alpha = 0.2, color = NA) +
  geom_line(aes(y = Avg_Price), linewidth = 0.8) +
  scale_color_brewer(palette = "Dark2") +
  scale_fill_brewer(palette = "Set2") +
  labs(y = "Price (€/100 kg)",
       x = "Date") +
  theme_minimal() +
  theme(legend.position = "bottom")




