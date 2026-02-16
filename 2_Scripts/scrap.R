# contract costs based on historical nitrogen prices ###########################################################
hist_f_p <- f_data |> 
  filter(Date >= as.Date("2009-01-01")) |> 
  dplyr::select(Date, Avg_Price) |> 
  mutate(t = row_number(),
         year = as.numeric(format(Date, "%Y")),
         month = as.numeric(format(Date, "%m")),
         year_index = year - min(year) + 1,
         f_p = Avg_Price/1000 # convert to €/kg
         )
         
# define strike prices
max_prices <- quantile(hist_f_p$f_p, probs = c(0.5, 0.6, 0.7, 0.8, 0.9))

# flag exercise times
hist_f_p$exercise <- hist_f_p$t == exercise_time + obs_per_year * (hist_f_p$year_index - 1)

# Initialize results data frame
premium_results <- data.frame(
  strike_price = numeric(),
  option_premium = numeric()
)

# Loop through each strike price
for (i in 1:length(max_prices)) {
  strike <- max_prices[i]
  
  # Calculate payoffs for this strike price
  hist_f_p$payoff <- ifelse(
    hist_f_p$exercise,
    pmax(0, hist_f_p$f_p - strike),
    0
  )
  
  # Apply discount factor
  hist_f_p$discounted_payoff <- hist_f_p$payoff * discount_factor
  
  # Calculate option premium
  premium <- mean(hist_f_p[hist_f_p$exercise, ]$discounted_payoff)
  
  # Store results
  premium_results <- rbind(
    premium_results,
    data.frame(strike_price = strike, option_premium = premium)
  )
}

premium_results <- premium_results |>
  mutate(across(where(is.numeric), ~ round(.x, 3)))

premium_results

#write_csv(premium_results, "3_Outputs/Table_ContractPayoff_HistoricalPrices.csv")

## density plot of prices #########################################
strike_vec <- premium_results$strike_price

sim_f_p <- exp(as.vector(Ys_levels[1, exercise_time, ])) / 1000

# Plot density of original sim_f_p with strike prices 
plot_original <- ggplot(data.frame(sim_f_p), aes(x = sim_f_p)) +
  geom_density(alpha = 0.3, linewidth = 0.5) +
  geom_vline(xintercept = strike_vec, color = "red", linetype = "dashed") +
  labs(
    title = "Simulated prices with strike levels",
    x = NULL,
    y = "Density"
  ) +
  theme_minimal()

# one plot per strike price (capped distribution)
plot_list <- lapply(seq_along(strike_vec), function(i) {
  k <- strike_vec[i]
  sim_f_p_adj <- pmin(sim_f_p, k)
  
  df <- data.frame(
    value = c(sim_f_p, sim_f_p_adj),
    type  = c(rep("Original", length(sim_f_p)),
              rep(paste0("Capped (k = ", k, ")"), length(sim_f_p_adj)))
  )
  
  ggplot(df, aes(x = value, color = type, fill = type)) +
    geom_density(alpha = 0.3, linewidth = 0.5) +
    labs(
      title = paste("Strike price =", k),
      x = NULL,
      y = NULL
    ) +
    theme_minimal() +
    theme(legend.position = "none")
})

plot_list[[3]] <- plot_list[[3]] + labs(
    y = "Density"
  ) 

plot_list[[4]] <- plot_list[[4]] + labs(
    x = "Price"
  ) 

# Combine into 2x3 layout 
final_plot <- plot_original + plot_list[[1]] + plot_list[[2]] + plot_list[[3]] +
              plot_list[[4]] + plot_list[[5]] 

final_plot

#ggsave("3_Outputs/Fig_PriceDensities_SimulatedStrikes.pdf", plot = final_plot, width = 12, height = 8)


# pricing based on deflated historical nitrogen prices ###########################################################

ppi_f_base <- ppi_f |> filter(year == base_year) |> pull(ppi)

defl_f_p <- hist_f_p |> 
  left_join(ppi_f, by = "year") |> 
  mutate(real_f_p = f_p * (ppi_f_base / ppi)) 


# define strike prices
max_prices <- quantile(defl_f_p$real_f_p, probs = c(0.5, 0.6, 0.7, 0.8, 0.9))

# flag exercise times
hist_f_p$exercise <- hist_f_p$t == exercise_time + obs_per_year * (hist_f_p$year_index - 1)

# Initialize results data frame
premium_results <- data.frame(
  strike_price = numeric(),
  option_premium = numeric()
)

# Loop through each strike price
for (i in 1:length(max_prices)) {
  strike <- max_prices[i]
  
  # Calculate payoffs for this strike price
  hist_f_p$payoff <- ifelse(
    hist_f_p$exercise,
    pmax(0, hist_f_p$f_p - strike),
    0
  )
  
  # Apply discount factor
  hist_f_p$discounted_payoff <- hist_f_p$payoff * discount_factor
  
  # Calculate option premium
  premium <- mean(hist_f_p[hist_f_p$exercise, ]$discounted_payoff)
  
  # Store results
  premium_results <- rbind(
    premium_results,
    data.frame(strike_price = strike, option_premium = premium)
  )
}

premium_results <- premium_results |>
  mutate(across(where(is.numeric), ~ round(.x, 3)))

premium_results

#write_csv(premium_results, "3_Outputs/Table_ContractPayoff_HistoricalPrices.csv")

## density plot of prices #########################################
strike_vec <- premium_results$strike_price

# Plot density of original sim_f_p with strike prices 
plot_original <- ggplot(data.frame(sim_f_p), aes(x = sim_f_p)) +
  geom_density(alpha = 0.3, linewidth = 0.5) +
  geom_vline(xintercept = strike_vec, color = "red", linetype = "dashed") +
  labs(
    title = "Simulated prices with strike levels",
    x = NULL,
    y = "Density"
  ) +
  theme_minimal()

# one plot per strike price (capped distribution)
plot_list <- lapply(seq_along(strike_vec), function(i) {
  k <- strike_vec[i]
  sim_f_p_adj <- pmin(sim_f_p, k)
  
  df <- data.frame(
    value = c(sim_f_p, sim_f_p_adj),
    type  = c(rep("Original", length(sim_f_p)),
              rep(paste0("Capped (k = ", k, ")"), length(sim_f_p_adj)))
  )
  
  ggplot(df, aes(x = value, color = type, fill = type)) +
    geom_density(alpha = 0.3, linewidth = 0.5) +
    labs(
      title = paste("Strike price =", k),
      x = NULL,
      y = NULL
    ) +
    theme_minimal() +
    theme(legend.position = "none")
})

plot_list[[3]] <- plot_list[[3]] + labs(
    y = "Density"
  ) 

plot_list[[4]] <- plot_list[[4]] + labs(
    x = "Price"
  ) 

# Combine into 2x3 layout 
final_plot <- plot_original + plot_list[[1]] + plot_list[[2]] + plot_list[[3]] +
              plot_list[[4]] + plot_list[[5]] 

final_plot

#ggsave("3_Outputs/Fig_PriceDensities_DeflatedPrices.pdf", plot = final_plot, width = 12, height = 8)



# EUT calculation broad framing ###########################################################

# risk aversion values
delta <- c(0, 0.1, 0.229, 0.5, 1, 2)        

target_cols <- colnames(df_simulated_profit)

# Power utility function (to derive expected utilities)
utility_function <- function(delta, profit){
  if(delta == 1){
    log(profit)
  } else {
    ((profit^(1 - delta)) / (1 - delta))
  }
}

# Inverse utility function (to derive certainty equivalents)
inverse_utility_function <- function(delta, eu){
  if(delta == 1){
    exp(eu)
  } else {
    (eu*(1-delta))^(1/(1-delta))
  }
}

## calculate utilities ######################################################################### 

# initialize results array
utility <- array(dim=c(length(delta),nrow(df_simulated_profit), length(target_cols)))
dimnames(utility)[[3]] <- target_cols

# Utility 
for (d in 1:length(delta)){
  for (i in 1:nrow(df_simulated_profit)){
    for (c in 1:ncol(df_simulated_profit)){
      utility[d,i,c] <- utility_function(delta = delta[d], profit = df_simulated_profit[i,c])
    }
  }
}

utility[5, ,]

### calculate expected utility, CE, and risk premium
summary_EU <- array(dim = c(length(delta), nrow(df_simulated_profit), length(target_cols), 4))
dimnames(summary_EU)[[1]] <- delta
dimnames(summary_EU)[[3]] <- target_cols
dimnames(summary_EU)[[4]] <- c("EU","Expected Profits", "Certainty Equivalent", "Risk Premium")

for (d in 1:length(delta)){
  for (i in 1:nrow(df_simulated_profit)){
    for (c in 1:length(target_cols)){
    
      # Expected Utility 
      summary_EU[d,i,c,1] <- mean(as.numeric(utility[d,i,c]), na.rm=T)
      
      # Expected profits
      summary_EU[d,i,c,2] <- mean(as.numeric(df_simulated_profit[i,c]), na.rm=T) 
      
      # Certainty Equivalent
      summary_EU[d,i,c,3] <- inverse_utility_function(delta = delta[d], eu = as.numeric(summary_EU[d,i,c,1])) 
      
      # Risk Premium
      summary_EU[d,i,c,4] <- summary_EU[d,i,c,2] - summary_EU[d,i,c,3] 
    }
  }
}

summary_EU[5, , ,1]
summary_EU[5, , ,2]
summary_EU[5, , ,3]
summary_EU[5, , ,4]
summary(summary_EU[5, , ,4])


### testing for differences in utilty values #####################################################

contract_vs_contract_EU <- array(dim=c(length(delta),length(target_cols),length(target_cols)))
dimnames(contract_vs_contract_EU)[[1]] <- delta
dimnames(contract_vs_contract_EU)[[2]] <- target_cols
dimnames(contract_vs_contract_EU)[[3]] <- target_cols

# Calculating the p values
for (d in 2:length(delta)){
  for (c in 1:length(target_cols)){
    for (cc in 1:length(target_cols)){
      
      # Calculate the p-value of the test statistic
      contract_vs_contract_EU[d,c,cc] <- wilcox.test(utility[d,i,c], utility[d,i,cc], paired=T,alternative="l")$p.value
    }
  }
}

contract_vs_contract_EU[5,,]

### testing for differences in risk premium #####################################################

# testing against no contract  
nc_vs_contracts_RP <- array(dim = c(length(delta), length(target_cols), 1))
dimnames(nc_vs_contracts_RP)[[1]] <- delta
dimnames(nc_vs_contracts_RP)[[2]] <- target_cols
dimnames(nc_vs_contracts_RP)[[3]] <- ("fc")

# Calculating the p-values
for (d in 1:length(delta)){
  for (c in 1:length(target_cols)){
    nc_vs_contracts_RP[d,c,1] <- wilcox.test(summary_EU[d, ,c,4], summary_EU[d, ,1,4], paired=T, alternative="l")$p.value
  }
}

# testing the contracts against each other
contract_vs_contract_RP <- array(dim=c(length(delta),length(target_cols),length(target_cols)))
dimnames(contract_vs_contract_RP)[[1]] <- delta
dimnames(contract_vs_contract_RP)[[2]] <- target_cols
dimnames(contract_vs_contract_RP)[[3]] <- target_cols

# Calculating the p values
for (d in 2:length(delta)){
  for (c in 1:length(target_cols)){
    for (cc in 1:length(target_cols)){
      
      # Calculate the p-value of the test statistic
      contract_vs_contract_RP[d,c,cc] <- wilcox.test(summary_EU[d, ,c,4], summary_EU[d, ,cc,4], paired=T,alternative="l")$p.value
    }
  }
}

#### calculating absolute average changes in risk premiums #####################################

# contracts against no contract
average_RP_change <- as.data.frame(matrix(NA, nrow=length(target_cols), ncol=length(delta)))
row.names(average_RP_change) <- target_cols
colnames(average_RP_change) <- delta

for (d in 1:length(delta)){
  for (c in 1:length(target_cols)){
    temp1 <- summary_EU[d, ,c,4] - summary_EU[d, ,1,4]
    average_RP_change[c,d] <- mean(as.numeric(temp1))
    rm(temp1)
  }
}

average_RP_change


# contracts against each other
average_RP_change_contracts <- array(dim = c(length(delta), length(target_cols), length(target_cols)))
dimnames(average_RP_change_contracts)[[1]] <- delta
dimnames(average_RP_change_contracts)[[2]] <- target_cols
dimnames(average_RP_change_contracts)[[3]] <- target_cols

for (d in 1:length(delta)){
  for (c in 1:length(target_cols)){
    for (cc in 1:length(target_cols)){
      temp1 <- summary_EU[d, ,c,4] - summary_EU[d, ,cc,4]
      average_RP_change_contracts [d,c,cc] <- mean(as.numeric(temp1))
      rm(temp1)
    }
  }
}


average_RP_change_contracts[2,,]
average_RP_change_contracts[3,,]
average_RP_change_contracts[4,,]
average_RP_change_contracts[5,,]

#### calculating relative average changes in risk premiums ##################################################

#### changes in CE
#### average changes in the certainty eqivalent
#### relative changes in the certainty eqivalent

#### lower partial moments of first and second moment

# EUT calculations narrow framing ############################################

# Power utility function (to derive expected utilities)
u_cara <- function(x, delta) {
  return(-(1/delta) * exp(-delta * x))
}

u_cara_inv <- function(eu, delta) {
  return(-(1/delta) * log(-delta * eu))
}

## calculate utilities ######################################################################### 

# initialize results array
utility_narrow_framing <- array(dim=c(length(delta),nrow(df_profit_narrowframing), length(target_cols)))
dimnames(utility_narrow_framing)[[3]] <- target_cols

# Utility 
for (d in 2:length(delta)){
  for (i in 1:nrow(df_profit_narrowframing)){
    for (c in 1:ncol(df_profit_narrowframing)){
      utility_narrow_framing[d,i,c] <- u_cara(x = df_profit_narrowframing[i,c], delta = delta[d])
    }
  }
}

utility_narrow_framing[5, ,]
round(utility_narrow_framing[5,,],2)

### calculate expected utility, CE, and risk premium
summary_EU_narrow <- array(dim = c(length(delta), nrow(df_profit_narrowframing), length(target_cols), 4))
dimnames(summary_EU_narrow)[[1]] <- delta
dimnames(summary_EU_narrow)[[3]] <- target_cols
dimnames(summary_EU_narrow)[[4]] <- c("EU","Expected Profits", "Certainty Equivalent", "Risk Premium")

for (d in 1:length(delta)){
  for (i in 1:nrow(df_simulated_profit)){
    for (c in 1:length(target_cols)){
    
      # Expected Utility 
      summary_EU_narrow[d,i,c,1] <- mean(as.numeric(utility_narrow_framing[d,i,c]), na.rm=T)
      
      # Expected profits
      summary_EU_narrow[d,i,c,2] <- mean(as.numeric(df_profit_narrowframing[i,c]), na.rm=T) 
      
      # Certainty Equivalent
      summary_EU_narrow[d,i,c,3] <- u_cara_inv(delta = delta[d], eu = as.numeric(summary_EU_narrow[d,i,c,1])) 
      
      # Risk Premium
      summary_EU_narrow[d,i,c,4] <- summary_EU_narrow[d,i,c,2] - summary_EU_narrow[d,i,c,3] 
    }
  }
}

summary_EU_narrow[5, , ,1]
summary_EU_narrow[5, , ,2]
summary_EU_narrow[5, , ,3]
summary_EU_narrow[5, , ,4]

round(summary_EU_narrow[5, , ,4])

### testing for differences in risk premium #####################################################

# testing against no contract  
nc_vs_contracts_RP_narrow <- array(dim = c(length(delta), length(target_cols), 1))
dimnames(nc_vs_contracts_RP_narrow)[[1]] <- delta
dimnames(nc_vs_contracts_RP_narrow)[[2]] <- target_cols
dimnames(nc_vs_contracts_RP_narrow)[[3]] <- ("fc")

# Calculating the p-values
for (d in 1:length(delta)){
  for (c in 1:length(target_cols)){
    nc_vs_contracts_RP_narrow[d,c,1] <- wilcox.test(summary_EU[d, ,c,4], summary_EU[d, ,1,4], paired=T, alternative="l")$p.value
  }
}

# testing the contracts against each other
contract_vs_contract_RP <- array(dim=c(length(delta),length(target_cols),length(target_cols)))
dimnames(contract_vs_contract_RP)[[1]] <- delta
dimnames(contract_vs_contract_RP)[[2]] <- target_cols
dimnames(contract_vs_contract_RP)[[3]] <- target_cols

# Calculating the p values
for (d in 2:length(delta)){
  for (c in 1:length(target_cols)){
    for (cc in 1:length(target_cols)){
      
      # Calculate the p-value of the test statistic
      contract_vs_contract_RP[d,c,cc] <- wilcox.test(summary_EU[d, ,c,4], summary_EU[d, ,cc,4], paired=T,alternative="l")$p.value
    }
  }
}

contract_vs_contract_RP[5,,]

#### calculating absolute average changes in risk premiums #####################################

# contracts against no contract
average_RP_change <- as.data.frame(matrix(NA, nrow=length(target_cols), ncol=length(delta)))
row.names(average_RP_change) <- target_cols
colnames(average_RP_change) <- delta

for (d in 1:length(delta)){
  for (c in 1:length(target_cols)){
    temp1 <- summary_EU[d, ,c,4] - summary_EU[d, ,1,4]
    average_RP_change[c,d] <- mean(as.numeric(temp1))
    rm(temp1)
  }
}

average_RP_change <- round(average_RP_change,2)

# contracts against each other
average_RP_change_contracts <- array(dim = c(length(delta), length(target_cols), length(target_cols)))
dimnames(average_RP_change_contracts)[[1]] <- delta
dimnames(average_RP_change_contracts)[[2]] <- target_cols
dimnames(average_RP_change_contracts)[[3]] <- target_cols

for (d in 1:length(delta)){
  for (c in 1:length(target_cols)){
    for (cc in 1:length(target_cols)){
      temp1 <- summary_EU[d, ,c,4] - summary_EU[d, ,cc,4]
      average_RP_change_contracts [d,c,cc] <- mean(as.numeric(temp1))
      rm(temp1)
    }
  }
}


delta_0 <- round(average_RP_change_contracts[1,,],2)
delta_1 <- round(average_RP_change_contracts[2,,],2)
delta_2 <- round(average_RP_change_contracts[3,,],2)
delta_3 <- round(average_RP_change_contracts[4,,],2)
delta_4 <- round(average_RP_change_contracts[5,,],2)



