# =====================================================================================================================
# =====================================================================================================================
# 
# R Script for the publication: "The optimal drought index for designing weather index insurance"
#
# Authors:  Janic Bucheli, Tobias Dalhaus, Robert Finger
# Paper:    The optimal drought index for designing weather index insurance
# Journal:  European Review of Agricultural Economics
# 
# Last modification: May 19, 2020
# 
# Please report any bugs to <jbucheli@ethz.ch>.
#
# The script uses detrended yield data, farm-individual phenology observations and farm-individual meteorological data
# as input data.
#
# All definitions and dataset can be defined in section 1.0. From there on, do not change the code.
#
# =====================================================================================================================
# =====================================================================================================================

# =====================================================================================================================
# 1.0 Packages and data
# =====================================================================================================================

# Load packages if necessary
#install.packages("lubridate")
#install.packages("SPEI")
#install.packages("quantreg")

library(lubridate)
library(SPEI)
library(quantreg)

# Define quantile of interest (for quantile estimator and derivation of tick size)
quantile_of_interest <- 0.3

# Define premium load (e.g. 0.1 for 10%)
premium_load <- 0.0

# Define wheat price:
wheat_price     <- 15.80

# Load the yield, phenology and meteorological data
# Detrended wheat yields
# Each row represents a farm and each column a year
farm_yield_detr <- read.csv("detrended_farm_wheat_yield.csv", row.names=1)
colnames(farm_yield_detr) <- seq(from=1995,to=2015, by=1)

# Phenology from network of growth phase reporters from year 1992
# A colum with Farm, year, date of stem elongation (dat_stem) and date of milk ripeness (dat_milk) 
pheno <- read.csv("phenology_wheat_farm.csv", row.names= 1)
# We remove the year 1991 because of unreliable reports
pheno <- pheno[-which(pheno[,1]==1991),]
# Read dates as dates
pheno <- transform(pheno, dat_stem = as.Date(as.character(dat_stem), "%Y%m%d"))
pheno <- transform(pheno, dat_ear = as.Date(as.character(dat_ear), "%Y%m%d"))
pheno <- transform(pheno, dat_milk = as.Date(as.character(dat_milk), "%Y%m%d"))

# Get farm ID from phenology dataset
farm_id <- unique(pheno$Farm)

# Meteorological data:each row is a date specified in first column. Subsequent columns represent farms. 
# Precipitation at farms
PP <- read.csv("PP_farm.csv", row.names= 1)
colnames(PP) <- c("Date",farm_id)
PP[,1] <- as.Date(PP[,1])

# Potential evapotranspiration at farms
PET <- read.csv("PET_farm.csv", row.names= 1)
colnames(PET) <- c("Date",farm_id)
PET[,1] <- as.Date(PET[,1])
# Dividing the PET by 10 to convert into unit mm
PET[,2:ncol(PET)] <- PET[,2:ncol(PET)] / (10)

# Actual evapotranspiration at farms
AET <- read.csv("AET_farm.csv", row.names= 1)
colnames(AET) <- c("Date",farm_id)
AET[,1] <- as.Date(AET[,1])
# Dividing the AET by 10 to convert into unit mm
AET[,2:ncol(AET)] <- AET[,2:ncol(AET)] / (10)

# Soil moisture at farms
SM <- read.csv("SM_farm.csv", row.names= 1)
colnames(SM) <- c("Date",farm_id)
SM[,1] <- as.Date(SM[,1])

# =====================================================================================================================
# 2.0 Calculations of drought indices
# =====================================================================================================================

# ---------------------------------------------------------------------------------------------------------------------
# 2.1 Defining the flexible index measurement period from stem elongation to milk ripeness and substitutional dates
# ---------------------------------------------------------------------------------------------------------------------

# First, we remove typos indicated by a negative length of the index measurement period (here only one removal)
pheno$riskdays <- as.numeric(pheno$dat_milk - pheno$dat_stem)
pheno[which(pheno$riskdays < 0),4:6] <- NA

# Second, we calculate substitutional dates when data is not available.
substitute_dates <- as.data.frame(matrix(NA,nrow=length(unique(pheno$year)), ncol=3))
substitute_dates[,1] <- sort(unique(pheno$year))
colnames(substitute_dates) <- c("year","start.pheno","end.pheno")

# Calculation of the substitutional dates (i.e. average dates) 
for(t in 1:length(unique(pheno$year))){
  sub_pheno.final <- subset(pheno, year==substitute_dates[t,1])
  substitute_dates[t,2] <- mean.Date(sub_pheno.final$dat_stem, na.rm=T)
  substitute_dates[t,3] <- mean.Date(sub_pheno.final$dat_milk, na.rm=T)
}

# Correct formate
class.date <- class(pheno$dat_stem)
class(substitute_dates$start.pheno) <- class.date
class(substitute_dates$end.pheno) <- class.date
format(substitute_dates$start.pheno, format="%Y-%m-%d")
format(substitute_dates$end.pheno, format="%Y-%m-%d")

# Tidy up global environment
rm(sub_pheno.final)

# ---------------------------------------------------------------------------------------------------------------------
# 2.2 Index calculation
# ---------------------------------------------------------------------------------------------------------------------
#
# We calculate the idnices for each farm and year in accordance to the flexible index measurement period
# If no phenology data is available, we substitute with the average date for the particular year
# Phenology is either not available because there is no reliable report available or there is an obvious typo
# Note that the index measurement period ends 1 day before enetering the growth phase milk ripeness.

# We define a unique farm-year identifier to check the need for substitution in the following calculations
pheno$key <- as.numeric(paste(pheno$year,pheno$Farm,sep=""))

#----------------------------------------------------------------------------------------------------------------------
# 2.2.1 Cumulative precipitation index (CPI)
#----------------------------------------------------------------------------------------------------------------------

# Empty df with row=farm and col=years
CPI_df <- data.frame(matrix(NA,nrow=length(farm_id), ncol= length(unique(pheno[,1]))))
row.names(CPI_df) <- farm_id
colnames(CPI_df) <- unique(pheno[,1])

for (i in 1:length(farm_id)){
  for (t in 1:length(unique(pheno[,1]))){
    
    # Creation of unique identifier for farm i and year t
    obs <- as.numeric(paste(unique(pheno[,1])[t],farm_id[i],sep=""))
    
    # Is there an observation for farm i in year t? If no then substitute.
    if(obs %in% pheno$key == FALSE)
    {
      start.pheno <- (substitute_dates[which(substitute_dates$year== unique(pheno[,1])[t]),2])
      end.pheno <- (substitute_dates[which(substitute_dates$year== unique(pheno[,1])[t]),3]-1)
      
      
      # Have we a negative risk period and thus incerted NA? If yes, then substitute.
    }else if (is.na(pheno[which(pheno[,1]==unique(pheno[,1])[t] & pheno[,8]==farm_id[i]), which(colnames(pheno)=="dat_stem")]))
    {
      start.pheno <- (substitute_dates[which(substitute_dates$year== unique(pheno[,1])[t]),2])
      end.pheno <- (substitute_dates[which(substitute_dates$year== unique(pheno[,1])[t]),3]-1)
      
      
      # If the observation for farm i and year t is available and not an outlier, take this value 
    }else
    {
      start.pheno <- pheno[which(pheno[,1]==unique(pheno[,1])[t] & pheno[,8]==farm_id[i]), which(colnames(pheno)=="dat_stem")]
      end.pheno <- (pheno[which(pheno[,1]==unique(pheno[,1])[t] & pheno[,8]==farm_id[i]), which(colnames(pheno)=="dat_milk")]-1)
    }
    
    # Make a subset of the risk period for farm i in year t (each column after the date is a farm!)
    sub_pp <- subset(PP, PP[,1] >= start.pheno & PP[,1] <= end.pheno)
    
    # We fill the pp_df with the cumulative precipitation for each farm i and year t
    CPI_df[i,t] <- sum(sub_pp[,which(row.names(CPI_df)[i]==colnames(sub_pp))])
    
    # We remove the subsets we do not need anymore
    rm(sub_pp)
  }
}

#----------------------------------------------------------------------------------------------------------------------
# 2.2.2 Standardized precipitation index (SPI)
#----------------------------------------------------------------------------------------------------------------------
#
# The SPEI-package currently calculates monthly SPI values at different time scales.It interpretes the records as monthly
# values of the cumulative precipitation index starting in January and then each following record represents the consecutive months.
# We have a flexible period of index measurement (not monthly) and set the time scale to 1. Hence, we can insert the CPI values at  
# every 12th position of a vector and the SPEI-package subsequently calculates the SPEI-values for the flexible index measurement 
# periods correctly. This approach only works when the time scale is equal to 1! Values have to be ordered in columns. 

# DF with nrow = (12 months) * (years of available data)
# This is just a vector for CPI values at every 12th location of a column.
input_SPI <- data.frame(matrix(NA,nrow=ncol(CPI_df)*12, ncol=nrow(CPI_df)))
colnames(input_SPI) <- seq(from=1, to= ncol(input_SPI), by=1)

# Here we assign the CPI-values to every 12th location in a column.
for(i in 1:ncol(input_SPI)){
  input_SPI[1,i] <- CPI_df[i,1]
  
  for(t in 1:(ncol(CPI_df)-1)){
    input_SPI[12*t+1,i] <- CPI_df[i,t+1]
  }
}

# The columns in input_SPI cannot contain NA-values. We replace them by 0.
for (i in 1:ncol(input_SPI)){
  for (t in 1:nrow(input_SPI)){
    
    if(is.na(input_SPI[t,i])){input_SPI[t,i] <- 0}
    else {}
  }
}

# Calculation of SPI values for time scale 1.Default distribution is Gamma.
SPI_values <- spi(input_SPI,1)

# Extraction of the correct values (every 12th location for each farm)
SPI <- SPI_values[[2]][seq(from=1,to=(nrow(input_SPI)-11),by=12),]

# Transposition to consistent df.
SPI_df <- as.data.frame(t(SPI))
colnames(SPI_df) <- colnames(CPI_df)

# Tidying up the global environment
rm(input_SPI,SPI,SPI_values)

#----------------------------------------------------------------------------------------------------------------------
# 2.2.3 Standardized precipitation evapotranspiration index (SPEI)
#----------------------------------------------------------------------------------------------------------------------
#
# We use cumulative precipitation and cumulative potential evapotranspiration to calculate the climatic water balance 
# (precipitation minus potential evapotranspiration) than is subsequently standardized to the SPEI.

# Calculation of cumulative potential evapotranspiration (analog to cumulative precipitation)

# Empty df with row=farm and col=years
CPET_df <- data.frame(matrix(NA,nrow=length(farm_id), ncol= length(unique(pheno[,1]))))
row.names(CPET_df) <- farm_id
colnames(CPET_df) <- unique(pheno[,1])


# Calculation of the CPET
for (i in 1:length(farm_id)){
  for (t in 1:length(unique(pheno[,1]))){
    
    # Creation of unique identifier for farm i and year t
    obs <- as.numeric(paste(unique(pheno[,1])[t],farm_id[i],sep=""))
    
    # Is there an observation for farm i in year t? If no then substitute.
    if(obs %in% pheno$key == FALSE)
    {
      start.pheno <- (substitute_dates[which(substitute_dates$year== unique(pheno[,1])[t]),2])
      end.pheno <- (substitute_dates[which(substitute_dates$year== unique(pheno[,1])[t]),3]-1)
      
      
      # Have we a negative risk period and thus incerted NA? If yes, then substitute.
    }else if (is.na(pheno[which(pheno[,1]==unique(pheno[,1])[t] & pheno[,8]==farm_id[i]), which(colnames(pheno)=="dat_stem")]))
    {
      
      start.pheno <- (substitute_dates[which(substitute_dates$year== unique(pheno[,1])[t]),2])
      end.pheno <- (substitute_dates[which(substitute_dates$year== unique(pheno[,1])[t]),3]-1)
      
      
      # If the observation for farm i and year t is available and not an outlier, take this value 
    }else
    {
      
      start.pheno <- pheno[which(pheno[,1]==unique(pheno[,1])[t] & pheno[,8]==farm_id[i]), which(colnames(pheno)=="dat_stem")]
      end.pheno <- (pheno[which(pheno[,1]==unique(pheno[,1])[t] & pheno[,8]==farm_id[i]), which(colnames(pheno)=="dat_milk")]-1)
    }
    
    # Make a subset of the risk period for farm i in year t (each column after the date is a farm!)
    sub_pet <- subset(PET, PET[,1] >= start.pheno & PET[,1] <= end.pheno)
    
    # We fill the pp_df with the cumulative precipitation for each farm i and year t
    
    CPET_df[i,t] <- sum(sub_pet[,which(row.names(CPET_df)[i]==colnames(sub_pet))])
    
    # We remove the subsets we do not need anymore
    rm(sub_pet)
  }
}

# Calculation of the cumulative climatic water balance

WBI_df <- data.frame(matrix(NA,nrow=length(farm_id), ncol= length(unique(pheno[,1]))))
row.names(WBI_df) <- farm_id
colnames(WBI_df) <- unique(pheno[,1])
WBI_df <- CPI_df - CPET_df

# The SPEI-package currently calculates monthly SPEI values at different time scales.It interpretes the records as monthly
# values of the climatic water balance starting in January and then each following record represents the consecutive months.
# We have a flexible period of index measurement (not monthly) and set the time scale to 1. Hence, we can insert the climatic  
# water balance values at every 12th position of a vector and the SPEI-package subsequently calculates the SPEI-values for 
# the flexible index measurement periods correctly. This approach only works when the time scale is equal to 1! Values 
# have to be ordered in columns.

# DF with nrow = (12 months) * (years of available data)
# This is just a vector for CPI values at every 12th location of a column.
input_SPEI <- data.frame(matrix(NA,nrow=ncol(WBI_df)*12, ncol=nrow(WBI_df)))
colnames(input_SPEI) <- seq(from=1, to= ncol(input_SPEI), by=1)

# Here we assign the climatic water balance values to every 12th location in a column.
for(i in 1:ncol(input_SPEI)){
  input_SPEI[1,i] <- WBI_df[i,1]
  
  for(t in 1:(ncol(CPI_df)-1)){
    input_SPEI[12*t+1,i] <- WBI_df[i,t+1]
  }
}

# The columns in input_SPEI cannot contain NA-values. We replace them by 0.
for (i in 1:ncol(input_SPEI)){
  for (t in 1:nrow(input_SPEI)){
    
    if(is.na(input_SPEI[t,i])){input_SPEI[t,i] <- 0}
    else {}
  }
}

# Calculation of SPI values for time scale 1.Default distribution is log-logistic.
SPEI_values <- spei(input_SPEI,1)

# Extraction of the correct values (every 12th location for each farm)
SPEI <- SPEI_values[[2]][seq(from=1,to=(nrow(input_SPEI)-11),by=12),]

# Transposition to consistent df.
SPEI_df <- as.data.frame(t(SPEI))
colnames(SPEI_df) <- colnames(WBI_df)

# Tidying up the global environment
rm(input_SPEI,SPEI,SPEI_values, WBI_df)

#----------------------------------------------------------------------------------------------------------------------
# 2.2.4 Soil moisture index (SMI)
#----------------------------------------------------------------------------------------------------------------------
#
# Representing the average of plant available field capacity

# Empty df with row=farm and col=years
SMI_df <- data.frame(matrix(NA,nrow=length(farm_id), ncol= length(unique(pheno[,1]))))
row.names(SMI_df) <- farm_id
colnames(SMI_df) <- unique(pheno[,1])


# Calculation of the average plant available field capacity
for (i in 1:length(farm_id)){
  for (t in 1:length(unique(pheno[,1]))){
    
    # Creation of unique identifier for farm i and year t
    obs <- as.numeric(paste(unique(pheno[,1])[t],farm_id[i],sep=""))
    
    # Is there an observation for farm i in year t? If no then substitute.
    if(obs %in% pheno$key == FALSE)
    {
      
      start.pheno <- (substitute_dates[which(substitute_dates$year== unique(pheno[,1])[t]),2])
      end.pheno <- (substitute_dates[which(substitute_dates$year== unique(pheno[,1])[t]),3]-1)
      
      
      # Have we a negative risk period and thus incerted NA? If yes, then substitute.
    }else if (is.na(pheno[which(pheno[,1]==unique(pheno[,1])[t] & pheno[,8]==farm_id[i]), which(colnames(pheno)=="dat_stem")]))
    {
      start.pheno <- (substitute_dates[which(substitute_dates$year== unique(pheno[,1])[t]),2])
      end.pheno <- (substitute_dates[which(substitute_dates$year== unique(pheno[,1])[t]),3]-1)
      
      
      # If the observation for farm i and year t is available and not an outlier, take this value 
    }else
    {
      
      start.pheno <- pheno[which(pheno[,1]==unique(pheno[,1])[t] & pheno[,8]==farm_id[i]), which(colnames(pheno)=="dat_stem")]
      end.pheno <- (pheno[which(pheno[,1]==unique(pheno[,1])[t] & pheno[,8]==farm_id[i]), which(colnames(pheno)=="dat_milk")]-1)
    }
    
    # Make a subset of the risk period for farm i in year t (each column after the date is a farm!)
    sub_smi <- subset(SM, SM[,1] >= start.pheno & SM[,1] <= end.pheno)
    
    # We fill the pp_df with the cumulative precipitation for each farm i and year t
    
    SMI_df[i,t] <- mean(sub_smi[,which(row.names(SMI_df)[i]==colnames(sub_smi))])
    
    #We remove the subsets we do not need anymore
    rm(sub_smi)
  }
}

#----------------------------------------------------------------------------------------------------------------------
# 2.2.4 Evaporative stress index
#----------------------------------------------------------------------------------------------------------------------
#
# We use cumulative evapotranspiration calculated in 2.2.3 and cumulative actual evapotranspiration to calculate the ratio of  
# actual to potential evapotranspiration that is subsequently standardized. Note that the standardization is not based on a
# probabilistic approach.

# Calculation of cumulative actual evapotranspiration (analog to cumulative precipitation or potential evapotranspiration)


# Empty df with row=farm and col=years
CAET_df <- data.frame(matrix(NA,nrow=length(farm_id), ncol= length(unique(pheno[,1]))))
row.names(CAET_df) <- farm_id
colnames(CAET_df) <- unique(pheno[,1])


# Calculation of the CPET
for (i in 1:length(farm_id)){
  for (t in 1:length(unique(pheno[,1]))){
    
    # Creation of unique identifier for farm i and year t
    obs <- as.numeric(paste(unique(pheno[,1])[t],farm_id[i],sep=""))
    
    # Is there an observation for farm i in year t? If no then substitute.
    if(obs %in% pheno$key == FALSE)
    {
      
      start.pheno <- (substitute_dates[which(substitute_dates$year== unique(pheno[,1])[t]),2])
      end.pheno <- (substitute_dates[which(substitute_dates$year== unique(pheno[,1])[t]),3]-1)
      
      
      # Have we a negative risk period and thus incerted NA? If yes, then substitute.
    }else if (is.na(pheno[which(pheno[,1]==unique(pheno[,1])[t] & pheno[,8]==farm_id[i]), which(colnames(pheno)=="dat_stem")]))
    {
      
      start.pheno <- (substitute_dates[which(substitute_dates$year== unique(pheno[,1])[t]),2])
      end.pheno <- (substitute_dates[which(substitute_dates$year== unique(pheno[,1])[t]),3]-1)
      
      
      # If the observation for farm i and year t is available and not an outlier, take this value 
    }else
    {
      
      start.pheno <- pheno[which(pheno[,1]==unique(pheno[,1])[t] & pheno[,8]==farm_id[i]), which(colnames(pheno)=="dat_stem")]
      end.pheno <- (pheno[which(pheno[,1]==unique(pheno[,1])[t] & pheno[,8]==farm_id[i]), which(colnames(pheno)=="dat_milk")]-1)
    }
    
    # Make a subset of the risk period for farm i in year t (each column after the date is a farm!)
    sub_aet <- subset(AET, AET[,1] >= start.pheno & AET[,1] <= end.pheno)
    
    # We fill the pp_df with the cumulative precipitation for each farm i and year t
    
    CAET_df[i,t] <- sum(sub_aet[,which(row.names(CAET_df)[i]==colnames(sub_aet))])
    
  }
  
  #We remove the subsets we do not need anymore
  rm(sub_aet)
}

# Empty df for water requirement stress index (WRSI = ratio of actual to potential evapot.) with row=farm and col=years
WRSI_df <- data.frame(matrix(NA,nrow=length(farm_id), ncol= length(unique(pheno[,1]))))
row.names(WRSI_df) <- farm_id
colnames(WRSI_df) <- unique(pheno[,1])

WRSI_df <- CAET_df/CPET_df * 100

# Standardized anomalies of WRSI = ESI
ESI_df <- data.frame(matrix(NA,nrow=length(farm_id), ncol= length(unique(pheno[,1]))))
row.names(ESI_df) <- farm_id
colnames(ESI_df) <- unique(pheno[,1])

for (i in 1:nrow(farm_yield_detr)){
  ESI_df [i,] <- (WRSI_df [i,] - mean(as.numeric(WRSI_df[i,]))) / sd(as.numeric(WRSI_df[i,]))
}


# We tidy up the global environment after index calculations
rm (WRSI_df)

# =====================================================================================================================
# 3.0 Contract specifics: Tick size, strike level and premium
# =====================================================================================================================
#
# Here, we define all farm-individual contract specifics. 

# We design contracts for the following indices:
indices     <- c("CPI","SPI","SPEI","SMI","ESI")

# We create an empty array. Each level of the array represents an index
contract_specifics                <- array(dim=c(length(indices),nrow(farm_yield_detr),ncol=4))
dimnames(contract_specifics)[[1]] <- indices 
dimnames(contract_specifics)[[3]] <- c("Intercept QR","tick size", "strike", "premium")

# Array for historical payouts. Each level of the array represents an index.
historical_payouts                <- array(dim=c(length(indices),nrow(farm_yield_detr),ncol(farm_yield_detr)))
dimnames(historical_payouts)[[1]] <- indices 
dimnames(historical_payouts)[[3]] <- colnames(farm_yield_detr)

# Years with yield entries
yield_years <- colnames(farm_yield_detr)

# Calculation of contract specifics
for (k in 1:length(indices)){
  
  # Load the input index k and subset from years 1995-2015. Originally we would have data  from 1992 to 2016 but only yield data from 1995 to 2015.
  temp_index <- get(paste(indices[k],"df",sep="_"))
  sub_temp_index <- temp_index[,yield_years]
  
  for (i in 1:nrow(farm_yield_detr)){
    
    # ------------------------------------------------------------------------------------------------------------------
    # 3.1 Tick size 
    # ------------------------------------------------------------------------------------------------------------------
    
    # Run quantile regression and save intercept & tick size in array contract_specifics (first/ second column for each index)
    temp_reg <- rq(as.numeric(farm_yield_detr[i,]) ~ as.numeric(sub_temp_index[i,]), tau=quantile_of_interest)
    contract_specifics[k,i,1] <- temp_reg$coefficients[1]
    contract_specifics[k,i,2] <- temp_reg$coefficients[2]
    
    # ------------------------------------------------------------------------------------------------------------------
    # 3.2 Strike level 
    # ------------------------------------------------------------------------------------------------------------------
    
    contract_specifics[k,i,3] <- ((quantile(as.numeric(farm_yield_detr[i,]), na.rm=T, quantile_of_interest, type=1)-(temp_reg$coefficients[1])) / temp_reg$coefficients[2])
    
    # ------------------------------------------------------------------------------------------------------------------
    # 3.3 Premium 
    # ------------------------------------------------------------------------------------------------------------------
    
    # Create a sub_subset for pricing -> only calculate payouts in years the yield is reported (section 2.3)
    available_data <- which(!is.na(farm_yield_detr[i,])==TRUE)
    sub_temp_index_pricing <- sub_temp_index[i,available_data]
    
    # Calculation of historical payouts of farm i
    payout_vector <- vector(length = ncol(sub_temp_index_pricing))
    names(payout_vector) <- colnames(sub_temp_index_pricing)
    
    for (l in 1:length(payout_vector)){
      
      payout_vector[l] <- as.numeric(wheat_price)*as.numeric(contract_specifics[k,i,2])*
        max((as.numeric(contract_specifics[k,i,3]) - as.numeric(sub_temp_index_pricing[l])), 0)
    }
    
    # Farmers with a negative payout do not show a drought risk exposure from the empirical data.
    # Hence, they will get a premium of 0 Euros, i.e. should not buy WII for drought risk management.
    payout_vector[payout_vector < 0] <- 0
    
    contract_specifics[k,i,4] <- mean(payout_vector)
    
    for (l in 1:length(payout_vector)){
      historical_payouts[k,i,which(dimnames(historical_payouts)[[3]]==names(payout_vector[l]))] <- payout_vector[l]
    }
  }}

### Add a loading:

contract_specifics[,,4] <- contract_specifics [,,4] * (1+premium_load)

# =====================================================================================================================
# 4. Empirical risk analysis
# =====================================================================================================================

# ---------------------------------------------------------------------------------------------------------------------
# 4.1 Revenues without and with insurance
# ---------------------------------------------------------------------------------------------------------------------

# Revenue without insurance
revenue_crop <- as.data.frame(matrix(NA,nrow=nrow(farm_yield_detr), ncol=ncol(farm_yield_detr)))
colnames(revenue_crop) <- colnames(farm_yield_detr)

for (i in 1:nrow(revenue_crop)){
  revenue_crop[i,] <- farm_yield_detr[i,]*wheat_price
}

# Revenues with insurance
revenue_insured_crop    <- array(dim=c(length(indices),nrow(farm_yield_detr),ncol(farm_yield_detr)))
dimnames(revenue_insured_crop)[[3]] <- yield_years

for (k in 1:length(indices)){
  for (i in 1:nrow(farm_yield_detr)){
    
    revenue_insured_crop[k,i,] <- (as.numeric(revenue_crop[i,]) + as.numeric(historical_payouts[k,i,]) - as.numeric(contract_specifics[k,i,4]))
  }
}

# ---------------------------------------------------------------------------------------------------------------------
# 4.2 Risk analysis I: Expected utility model
# ---------------------------------------------------------------------------------------------------------------------

# Coefficients of constant relative risk aversion
alpha <- c(0, 0.5, 1, 2, 3, 4)

# Power utility function (to derive expected utilities)
utility_function <- function(alpha,revenue){
  
  if(alpha == 1){
    log(revenue)
  } else {
    ((revenue^(1 - alpha)) / (1 - alpha))
  }
}

# Inverse utility function (to derive certainty equivalents)
inverse_utility_function <- function(alpha,eu){
  
  if(alpha == 1){
    exp(eu)
  } else {
    (eu*(1-alpha))^(1/(1-alpha))
  }
}

# ---------------------------------------------------------------------------------------------------------------------
# 4.2.1 Calculating utilities
# ---------------------------------------------------------------------------------------------------------------------


# We create an empty array for each index. Each level of the array represent a CRRA-coefficient
utility_uninsured <- array(dim=c(length(alpha),nrow(farm_yield_detr),ncol(farm_yield_detr)))
utility_CPI <- array(dim=c(length(alpha),nrow(farm_yield_detr),ncol(farm_yield_detr)))
utility_SPI <- array(dim=c(length(alpha),nrow(farm_yield_detr),ncol(farm_yield_detr)))
utility_SPEI <- array(dim=c(length(alpha),nrow(farm_yield_detr),ncol(farm_yield_detr)))
utility_SMI <- array(dim=c(length(alpha),nrow(farm_yield_detr),ncol(farm_yield_detr)))
utility_ESI <- array(dim=c(length(alpha),nrow(farm_yield_detr),ncol(farm_yield_detr)))

# Utility of being uninsured
for (a in 1:length(alpha)){
  for (i in 1:nrow(farm_yield_detr)){
    for (t in 1:ncol(farm_yield_detr)){
      utility_uninsured[a,i,t] <- utility_function(alpha = alpha[a], revenue = revenue_crop[i,t])
    }
  }
}

# Utility for each index 
for (k in 1:length(indices)){
  # Get the correct utility array
  temp <- get(paste("utility",indices[k],sep="_"))
  
  for (a in 1:length(alpha)){
    for (i in 1:nrow(farm_yield_detr)){
      for (t in 1:ncol(farm_yield_detr)){
        temp[a,i,t] <- utility_function(alpha = alpha[a], revenue = revenue_insured_crop[k,i,t]) 
      }
    }
  }
  assign(paste("utility",indices[k],sep="_"),temp, envir = .GlobalEnv)
}

# ---------------------------------------------------------------------------------------------------------------------
# 4.2.2 Calculating expected utility, certainty equivalent and risk premium
# ---------------------------------------------------------------------------------------------------------------------

# For being uninsured
summary_EU_uninsured <- array(dim=c(length(alpha),nrow(farm_yield_detr),4))
dimnames(summary_EU_uninsured)[[1]] <- alpha
dimnames(summary_EU_uninsured)[[3]] <- c("EU","Expected Revenue", "Certainty Equivalent", "Risk Premium")

for (a in 1:length(alpha)){
  for (i in 1:nrow(farm_yield_detr)){
    
    # Expected Utility 
    summary_EU_uninsured [a,i,1] <- mean(as.numeric(utility_uninsured[a,i,]), na.rm=T)
    
    # Expected Revenue
    summary_EU_uninsured [a,i,2] <- mean(as.numeric(revenue_crop[i,]), na.rm=T) 
    
    # Certainty Equivalent
    summary_EU_uninsured [a,i,3] <- inverse_utility_function(alpha = alpha[a], eu = as.numeric(summary_EU_uninsured[a,i,1])) 
    
    # Risk Premium
    summary_EU_uninsured [a,i,4] <- summary_EU_uninsured [a,i,2] - summary_EU_uninsured [a,i,3] 
  }
}

# For each drought index
for (k in 1:length(indices)){
  
  temp_utility <- get(paste("utility", indices[k], sep="_"))
  temp_revenue <- revenue_insured_crop [k,,]
  temp_risk <- array(dim=c(length(alpha),nrow(farm_yield_detr),4)) 
  
  for (a in 1:length(alpha)){
    for (i in 1:nrow(farm_yield_detr)){
      
      # Expected Utility    
      temp_risk[a,i,1] <- mean(as.numeric(temp_utility[a,i,]), na.rm=T)
      
      # Expected Revenue
      temp_risk[a,i,2] <- mean(as.numeric(temp_revenue[i,]), na.rm=T)
      
      # Certainty Equivalent
      temp_risk[a,i,3] <- inverse_utility_function(alpha = alpha[a], eu = as.numeric(temp_risk[a,i,1])) 
      
      # Risk Premium
      temp_risk[a,i,4] <- temp_risk [a,i,2] - temp_risk [a,i,3] 
    }
  }
  dimnames(temp_risk)[[1]] <- alpha
  dimnames(temp_risk)[[3]] <- c("EU","Expected Revenue", "Certainty Equivalent", "Risk Premium")
  assign(paste("summary_EU",indices[k],sep="_"),temp_risk, envir = .GlobalEnv)
}

# ---------------------------------------------------------------------------------------------------------------------
# 4.2.3 DF for risk premium and identification of BEST
# ---------------------------------------------------------------------------------------------------------------------

# For an overview and further data handling we put the risk premiums for each CRRA in a df
report_RP <- array(dim=c(length(alpha), nrow(farm_yield_detr), ncol=length(indices)+1))
dimnames(report_RP)[[3]] <- c("uninsured","CPI", "SPI", "SPEI", "SMI", "ESI")

for (a in 1:length(alpha)){
  for (i in 1:nrow(farm_yield_detr)){
    
    # Uninsured
    report_RP[a,i,1] <- summary_EU_uninsured[a,i,4]
    
    # Insured with index k
    for (k in 1:length(indices)){
      report_RP[a,i,k+1] <- get(paste("summary_EU", indices[k],sep="_"))[a,i,4] 
    }
  }
}

# Identification of BEST (each farm most risk-reducing index)
best_index <- as.data.frame(matrix(NA,nrow=nrow(farm_yield_detr), ncol=length(alpha))) 
colnames(best_index) <- alpha

# If several indices have the same RP, we choose the simplest one.
for (a in 1:length(alpha)){
  for (i in 1:nrow(farm_yield_detr)){
    temp1 <- which(report_RP[a,i,] == min(report_RP[a,i,])) 
    best_index [i,a] <- as.character(names(temp1[1]))
    rm(temp1)
  }
}

# Assign the RP value of the best drought index
best_index_RP <- as.data.frame(matrix(NA,nrow=nrow(farm_yield_detr), ncol=length(alpha)))
colnames(best_index_RP) <- alpha
for (a in 1:length(alpha)){
  for (i in 1:nrow(farm_yield_detr)){
    best_index_RP[i,a] <- report_RP[a,i,best_index[i,a]]
  }
}

# ---------------------------------------------------------------------------------------------------------------------
# 4.2.4 Testing for differences in risk premium
# ---------------------------------------------------------------------------------------------------------------------

# Testing the indices against being uninsured
indices_vs_uninsured_RP <- array(dim=c(length(alpha),length(indices),1))
dimnames(indices_vs_uninsured_RP)[[1]] <- alpha
dimnames(indices_vs_uninsured_RP)[[2]] <- indices
dimnames(indices_vs_uninsured_RP)[[3]] <- ("uninsured")

# Calculating the p-values
for (a in 1:length(alpha)){
  for (k in 1:length(indices)){
    
    indices_vs_uninsured_RP[a,k,1] <- wilcox.test(report_RP[a,,k+1], report_RP[a,,1], paired=T, alternative="l")$p.value
  }
}

# Testing the indices against each other
indices_vs_indices_RP <- array(dim=c(length(alpha),length(indices),length(indices)))
dimnames(indices_vs_indices_RP)[[1]] <- alpha
dimnames(indices_vs_indices_RP)[[2]] <- indices
dimnames(indices_vs_indices_RP)[[3]] <- indices

# Calculating the p values
for (a in 2:length(alpha)){
  for (k in 1:length(indices)){
    for (kk in 1:length(indices)){
      
      # Calculate the p-value of the test statistic
      indices_vs_indices_RP[a,k,kk] <- wilcox.test(report_RP[a,,k+1], report_RP[a,,kk+1], paired=T,alternative="l")$p.value
    }
  }
}

# Testing BEST against being uninsured and the indices
testing_best_RP <- as.data.frame(matrix(NA, nrow=length(alpha), ncol= (length(indices)+1)))
colnames(testing_best_RP) <- c("uninsured", "CPI", "SPI", "SPEI","SMI" ,"ESI")
row.names(testing_best_RP)<- alpha

for (a in 1:length(alpha)){
  for (k in 1:(length(indices)+1))
    testing_best_RP[a,k] <- wilcox.test(best_index_RP[,a], report_RP[a,,k], paired=T, alternative="l")$p.value 
}

# ---------------------------------------------------------------------------------------------------------------------
# 4.2.5 Calculating absolute average changes in the risk premium
# ---------------------------------------------------------------------------------------------------------------------

# Indices against uninsured
average_RP_change_uninsured <- as.data.frame(matrix(NA, nrow=length(indices), ncol=length(alpha)))
row.names(average_RP_change_uninsured) <- indices
colnames(average_RP_change_uninsured) <- alpha

for (a in 1:length(alpha)){
  for (k in 1:length(indices)){
    temp1 <- report_RP[a,,k+1] - report_RP[a,,1]
    average_RP_change_uninsured [k,a] <- mean(as.numeric(temp1))
    rm(temp1)
  }
}

# Indices against each other
average_RP_change_indices <- array(dim = c(length(alpha), length(indices), length(indices)))
dimnames(average_RP_change_indices)[[1]] <- alpha
dimnames(average_RP_change_indices)[[2]] <- indices
dimnames(average_RP_change_indices)[[3]] <- indices

for (a in 1:length(alpha)){
  for (k in 1:length(indices)){
    for (kk in 1:length(indices)){
      temp1 <- report_RP[a,,k+1] - report_RP[a,,kk+1]
      average_RP_change_indices [a,k,kk] <- mean(as.numeric(temp1))
      rm(temp1)
    }
  }
}

# BEST against uninsured
average_RP_change_BEST_uninsured <- vector(length= length(alpha))
names(average_RP_change_BEST_uninsured) <- alpha

for (a in 1:length(alpha)){
  average_RP_change_BEST_uninsured [a] <- mean(as.numeric(best_index_RP[,a] - report_RP[a,,1]))
}

# BEST against the indices

# BEST against indices
average_RP_change_BEST_indices <- as.data.frame(matrix(NA,nrow=length(alpha), ncol=length(indices)))
row.names(average_RP_change_BEST_indices) <- alpha
colnames(average_RP_change_BEST_indices) <- indices

for (a in 1:length(alpha)){
  for (k in 1:length(indices)){
    average_RP_change_BEST_indices [a,k] <- mean(best_index_RP[,a] - report_RP[a,,k+1]) 
  }
}

# ---------------------------------------------------------------------------------------------------------------------
# 4.2.5 Calculating relative average changes in the risk premium
# ---------------------------------------------------------------------------------------------------------------------

# Indices against uninsured
rel_average_RP_change_uninsured <- as.data.frame(matrix(NA, nrow=length(indices), ncol=length(alpha)))
row.names(rel_average_RP_change_uninsured) <- indices
colnames(rel_average_RP_change_uninsured) <- alpha

for (a in 1:length(alpha)){
  for (k in 1:length(indices)){
    temp1 <- (report_RP[a,,k+1] - report_RP[a,,1]) / report_RP[a,,1] 
    rel_average_RP_change_uninsured [k,a] <- mean(as.numeric(temp1))
    rm(temp1)
  }
}

# Indices against each other
rel_average_RP_change_indices <- array(dim = c(length(alpha), length(indices), length(indices)))
dimnames(rel_average_RP_change_indices)[[1]] <- alpha
dimnames(rel_average_RP_change_indices)[[2]] <- indices
dimnames(rel_average_RP_change_indices)[[3]] <- indices

for (a in 1:length(alpha)){
  for (k in 1:length(indices)){
    for (kk in 1:length(indices)){
      temp1 <- (report_RP[a,,k+1] - report_RP[a,,kk+1]) / report_RP[a,,kk+1]
      rel_average_RP_change_indices [a,k,kk] <- mean(as.numeric(temp1))
      rm(temp1)
    }
  }
}

# BEST against uninsured
rel_average_RP_change_BEST_uninsured <- vector(length= length(alpha))
names(rel_average_RP_change_BEST_uninsured) <- alpha

for (a in 1:length(alpha)){
  rel_average_RP_change_BEST_uninsured [a] <- mean(as.numeric(best_index_RP[,a] - report_RP[a,,1]) / report_RP[a,,1])
}

# BEST against the indices

# BEST against indices
rel_average_RP_change_BEST_indices <- as.data.frame(matrix(NA,nrow=length(alpha), ncol=length(indices)))
row.names(rel_average_RP_change_BEST_indices) <- alpha
colnames(rel_average_RP_change_BEST_indices) <- indices

for (a in 1:length(alpha)){
  for (k in 1:length(indices)){
    rel_average_RP_change_BEST_indices [a,k] <- mean(as.numeric(best_index_RP[,a] - report_RP[a,,k+1]) / report_RP[a,,k+1])
  }
}

# ---------------------------------------------------------------------------------------------------------------------
# 4.2.7 Changes in certainty equivalent (use this risk measure when actuarially fair premium is loaded)
# ---------------------------------------------------------------------------------------------------------------------

# For an overview and further data handling we put the risk premiums for each CRRA in a df
report_CE <- array(dim=c(length(alpha), nrow(farm_yield_detr), ncol=length(indices)+1))
dimnames(report_CE)[[3]] <- c("uninsured","CPI", "SPI", "SPEI", "SMI", "ESI")

for (a in 1:length(alpha)){
  for (i in 1:nrow(farm_yield_detr)){
    
    # Uninsured
    report_CE[a,i,1] <- summary_EU_uninsured[a,i,3]
    
    # Insured with index k
    for (k in 1:length(indices)){
      report_CE[a,i,k+1] <- get(paste("summary_EU", indices[k],sep="_"))[a,i,3] 
    }
  }
}

# Identification of BEST (each farm most risk-reducing index)
best_index_CE <- as.data.frame(matrix(NA,nrow=nrow(farm_yield_detr), ncol=length(alpha))) 
colnames(best_index_CE) <- alpha

# If several indices have the same CE, we choose the simplest one.
for (a in 1:length(alpha)){
  for (i in 1:nrow(farm_yield_detr)){
    temp1 <- which(report_CE[a,i,] == max(report_CE[a,i,])) 
    best_index_CE [i,a] <- as.character(names(temp1[1]))
    rm(temp1)
  }
}

# Assign the CE value of the best drought index
best_index_CE_value <- as.data.frame(matrix(NA,nrow=nrow(farm_yield_detr), ncol=length(alpha)))
colnames(best_index_CE_value) <- alpha
for (a in 1:length(alpha)){
  for (i in 1:nrow(farm_yield_detr)){
    best_index_CE_value[i,a] <- report_CE[a,i,best_index_CE[i,a]]
  }
}

# Testing the indices against being uninsured
indices_vs_uninsured_CE <- array(dim=c(length(alpha),length(indices),1))
dimnames(indices_vs_uninsured_CE)[[1]] <- alpha
dimnames(indices_vs_uninsured_CE)[[2]] <- indices
dimnames(indices_vs_uninsured_CE)[[3]] <- ("uninsured")

# Calculating the p-values
for (a in 1:length(alpha)){
  for (k in 1:length(indices)){
    
    indices_vs_uninsured_CE[a,k,1] <- wilcox.test(report_CE[a,,k+1], report_CE[a,,1], paired=T, alternative="g")$p.value
  }
}

# Testing the indices against each other
indices_vs_indices_CE <- array(dim=c(length(alpha),length(indices),length(indices)))
dimnames(indices_vs_indices_CE)[[1]] <- alpha
dimnames(indices_vs_indices_CE)[[2]] <- indices
dimnames(indices_vs_indices_CE)[[3]] <- indices

# Calculating the p values
for (a in 2:length(alpha)){
  for (k in 1:length(indices)){
    for (kk in 1:length(indices)){
      
      # Calculate the p-value of the test statistic
      indices_vs_indices_CE[a,k,kk] <- wilcox.test(report_CE[a,,k+1], report_CE[a,,kk+1], paired=T,alternative="g")$p.value
    }
  }
}

# Testing BEST against being uninsured and the indices
testing_best_CE <- as.data.frame(matrix(NA, nrow=length(alpha), ncol= (length(indices)+1)))
colnames(testing_best_CE) <- c("uninsured", "CPI", "SPI", "SPEI","SMI" ,"ESI")
row.names(testing_best_CE)<- alpha

for (a in 1:length(alpha)){
  for (k in 1:(length(indices)+1))
    testing_best_CE[a,k] <- wilcox.test(best_index_CE_value[,a], report_CE[a,,k], paired=T, alternative="g")$p.value 
}

# ---------------------------------------------------------------------------------------------------------------------
# 4.2.8 Calculating average changes in the certainty eqivalent
# ---------------------------------------------------------------------------------------------------------------------

# Indices against uninsured
average_CE_change_uninsured <- as.data.frame(matrix(NA, nrow=length(indices), ncol=length(alpha)))
row.names(average_CE_change_uninsured) <- indices
colnames(average_CE_change_uninsured) <- alpha

for (a in 1:length(alpha)){
  for (k in 1:length(indices)){
    temp1 <- report_CE[a,,k+1] - report_CE[a,,1]
    average_CE_change_uninsured [k,a] <- mean(as.numeric(temp1))
    rm(temp1)
  }
}

# Indices against each other
average_CE_change_indices <- array(dim = c(length(alpha), length(indices), length(indices)))
dimnames(average_CE_change_indices)[[1]] <- alpha
dimnames(average_CE_change_indices)[[2]] <- indices
dimnames(average_CE_change_indices)[[3]] <- indices

for (a in 1:length(alpha)){
  for (k in 1:length(indices)){
    for (kk in 1:length(indices)){
      temp1 <- report_CE[a,,k+1] - report_CE[a,,kk+1]
      average_CE_change_indices [a,k,kk] <- mean(as.numeric(temp1))
      rm(temp1)
    }
  }
}

# BEST against uninsured
average_CE_change_BEST_uninsured <- vector(length= length(alpha))
names(average_CE_change_BEST_uninsured) <- alpha

for (a in 1:length(alpha)){
  average_CE_change_BEST_uninsured [a] <- mean(as.numeric(best_index_CE_value[,a] - report_CE[a,,1]))
}

# BEST against the indices

# BEST against indices
average_CE_change_BEST_indices <- as.data.frame(matrix(NA,nrow=length(alpha), ncol=length(indices)))
row.names(average_CE_change_BEST_indices) <- alpha
colnames(average_CE_change_BEST_indices) <- indices

for (a in 1:length(alpha)){
  for (k in 1:length(indices)){
    average_CE_change_BEST_indices [a,k] <- mean(best_index_CE_value[,a] - report_CE[a,,k+1]) 
  }
}

# ---------------------------------------------------------------------------------------------------------------------
# 4.2.9 Calculating relative changes in certainty equivalent
# ---------------------------------------------------------------------------------------------------------------------

# Indices against uninsured
rel_average_CE_change_uninsured <- as.data.frame(matrix(NA, nrow=length(indices), ncol=length(alpha)))
row.names(rel_average_CE_change_uninsured) <- indices
colnames(rel_average_CE_change_uninsured) <- alpha

for (a in 1:length(alpha)){
  for (k in 1:length(indices)){
    temp1 <- (report_CE[a,,k+1] - report_CE[a,,1]) / report_CE[a,,1] 
    rel_average_CE_change_uninsured [k,a] <- mean(as.numeric(temp1))
    rm(temp1)
  }
}

# Indices against each other
rel_average_CE_change_indices <- array(dim = c(length(alpha), length(indices), length(indices)))
dimnames(rel_average_CE_change_indices)[[1]] <- alpha
dimnames(rel_average_CE_change_indices)[[2]] <- indices
dimnames(rel_average_CE_change_indices)[[3]] <- indices

for (a in 1:length(alpha)){
  for (k in 1:length(indices)){
    for (kk in 1:length(indices)){
      temp1 <- (report_CE[a,,k+1] - report_CE[a,,kk+1]) / report_CE[a,,kk+1]
      rel_average_CE_change_indices [a,k,kk] <- mean(as.numeric(temp1))
      rm(temp1)
    }
  }
}

# BEST against uninsured
rel_average_CE_change_BEST_uninsured <- vector(length= length(alpha))
names(rel_average_CE_change_BEST_uninsured) <- alpha

for (a in 1:length(alpha)){
  rel_average_CE_change_BEST_uninsured [a] <- mean(as.numeric(best_index_CE_value[,a] - report_CE[a,,1]) / report_CE[a,,1])
}

# BEST against the indices

# BEST against indices
rel_average_CE_change_BEST_indices <- as.data.frame(matrix(NA,nrow=length(alpha), ncol=length(indices)))
row.names(rel_average_CE_change_BEST_indices) <- alpha
colnames(rel_average_CE_change_BEST_indices) <- indices

for (a in 1:length(alpha)){
  for (k in 1:length(indices)){
    rel_average_CE_change_BEST_indices [a,k] <- mean(as.numeric(best_index_CE_value[,a] - report_CE[a,,k+1]) / report_CE[a,,k+1])
  }
}

# ---------------------------------------------------------------------------------------------------------------------
# 4.3 Lower partial moments of first and second moment
# ---------------------------------------------------------------------------------------------------------------------

# DF for LPM 1 & 2
LPM1 <- as.data.frame(matrix(NA,nrow=nrow(farm_yield_detr), ncol=length(indices)+1))
LPM2 <- as.data.frame(matrix(NA,nrow=nrow(farm_yield_detr), ncol=length(indices)+1))
colnames(LPM1) <- c("uninsured", "CPI", "SPI", "SPEI","SMI", "ESI")
colnames(LPM2) <- c("uninsured", "CPI", "SPI", "SPEI","SMI", "ESI")

# For being uninsured
for (i in 1:nrow(farm_yield_detr)){
  
  temp1 <- revenue_crop[i,!is.na(revenue_crop[i,])]
  expected_revenue <- quantile(as.numeric(temp1), 0.3, type=1)
  
  shortfall_uninsured <-expected_revenue - temp1
  shortfall_uninsured <-shortfall_uninsured[shortfall_uninsured>0]
  LPM1 [i,1] <- mean(as.numeric(shortfall_uninsured))
  
  variance_uninsured <- (shortfall_uninsured)^2
  LPM2 [i,1] <- mean(variance_uninsured)
  
  rm(temp1)
}

# For being insured
for (k in 1:length(indices)){
  for (i in 1:nrow(farm_yield_detr)){
    
    expected_revenue <- quantile(as.numeric(revenue_crop[i,!is.na(revenue_crop[i,])]), 0.3, type=1)
    temp1 <- revenue_insured_crop[k,i,!is.na(revenue_insured_crop[k,i,])]
    
    temp2 <- expected_revenue - temp1
    temp2 <- temp2[temp2>0]
    LPM1[i,k+1] <- mean(as.numeric(temp2))
    
    temp3 <- (temp2)^2
    LPM2[i,k+1] <- mean(as.numeric(temp3))
  }
}

# ---------------------------------------------------------------------------------------------------------------------
# 4.3.1 Identification of BEST
# ---------------------------------------------------------------------------------------------------------------------

# For LPM1
best_index_LPM1 <- vector(length=nrow(farm_yield_detr))

# If several indexes with same LPM1 then the simplest one is taken
for (i in 1:nrow(farm_yield_detr)){
  temp1 <- which(LPM1[i,] == min(LPM1[i,])) 
  best_index_LPM1 [i] <- as.character(colnames(LPM1)[temp1[1]])
  rm(temp1)
}

# Assign value of BEST for LPM1
best_index_LPM1_value <- vector(length=nrow(farm_yield_detr))
for (i in 1:nrow(farm_yield_detr)){
  best_index_LPM1_value [i] <- LPM1[i,best_index_LPM1[i]]
}


# For LPM2
best_index_LPM2 <- vector(length=nrow(farm_yield_detr))

# If several indexes with same LPM2 then the simplest one is taken
for (i in 1:nrow(farm_yield_detr)){
  
  temp1 <- which(LPM2[i,] == min(LPM2[i,])) 
  best_index_LPM2 [i] <- as.character(colnames(LPM2)[temp1[1]])
  rm(temp1)
}

# Assign value of BEST for LPM2
best_index_LPM2_value <- vector(length=nrow(farm_yield_detr))

for (i in 1:nrow(farm_yield_detr)){
  best_index_LPM2_value [i] <- LPM2[i,best_index_LPM2[i]]
}

# ---------------------------------------------------------------------------------------------------------------------
# 4.3.2 Testing for differences in the lower partial moments
# ---------------------------------------------------------------------------------------------------------------------

# LPM1: Indices vs. uninsured
indices_vs_uninsured_LPM1 <- as.data.frame(matrix(0,nrow=1, ncol=length(indices)))
colnames(indices_vs_uninsured_LPM1) <- indices 

for (k in 1:length(indices)){
  indices_vs_uninsured_LPM1[1,k] <- wilcox.test(LPM1[,k+1], LPM1[,1], paired=T, alternative="l")$p.value
}

# LPM2: Indices vs. uninsured
indices_vs_uninsured_LPM2 <- as.data.frame(matrix(0,nrow=1, ncol=length(indices)))
colnames(indices_vs_uninsured_LPM2) <- indices 

for (k in 1:length(indices)){
  indices_vs_uninsured_LPM2[1,k] <- wilcox.test(LPM2[,k+1], LPM2[,1], paired=T, alternative="l")$p.value
}


# LPM1: Indices vs. indices
indices_vs_indices_LPM1 <- as.data.frame(matrix(0,nrow=length(indices), ncol=length(indices)))
row.names(indices_vs_indices_LPM1) <- indices
colnames(indices_vs_indices_LPM1) <- indices    

for (k in 1:length(indices)){
  for (kk in 1:length(indices)){
    indices_vs_indices_LPM1[k,kk] <- wilcox.test(LPM1[,k+1], LPM1[,kk+1], paired=T,alternative="l")$p.value
  }
}

# LPM2: Indices vs. indices
indices_vs_indices_LPM2 <- as.data.frame(matrix(0,nrow=length(indices), ncol=length(indices)))
row.names(indices_vs_indices_LPM2) <- indices
colnames(indices_vs_indices_LPM2) <- indices    

for (k in 1:length(indices)){
  for (kk in 1:length(indices)){
    indices_vs_indices_LPM2[k,kk] <- wilcox.test(LPM2[,k+1], LPM2[,kk+1], paired=T,alternative="l")$p.value
  }
}

# BEST for LPM1
testing_best_LPM1 <- vector(length= length(indices)+1)
names(testing_best_LPM1) <- c("uninsured", "CPI", "SPI", "SPEI","SMI", "ESI")

for (k in 1:(length(indices)+1)){
  testing_best_LPM1[k] <- wilcox.test(best_index_LPM1_value, LPM1[,k], paired=T, alternative="l")$p.value 
}

# BEST for LPM2
testing_best_LPM2 <- vector(length= length(indices)+1)
names(testing_best_LPM2) <- c("uninsured", "CPI", "SPI", "SPEI","SMI","ESI")

for (k in 1:(length(indices)+1)){
  testing_best_LPM2[k] <- wilcox.test(best_index_LPM2_value, LPM2[,k], paired=T, alternative="l")$p.value 
}

# ---------------------------------------------------------------------------------------------------------------------
# 4.3.3 Average absolute changes in the lower partial moments
# ---------------------------------------------------------------------------------------------------------------------

# Indices vs. uninsured for LPM 1
average_LPM1_change_indices_uninsured <- as.data.frame(matrix(0,nrow=1, ncol=length(indices)))
colnames(average_LPM1_change_indices_uninsured) <- indices 

for (k in 1:length(indices)){
  average_LPM1_change_indices_uninsured[1,k] <- mean(LPM1[,k+1] - LPM1[,1])
}

# Indices against each other LPM 1
average_LPM1_change_indices <- as.data.frame(matrix(0,nrow=length(indices), ncol=length(indices)))
row.names(average_LPM1_change_indices) <- indices
colnames(average_LPM1_change_indices) <- indices    

for (k in 1:length(indices)){
  for (kk in 1:length(indices)){
    # Calculate the p-value of the test statistic
    average_LPM1_change_indices[k,kk] <- mean (LPM1[,k+1] - LPM1[,kk+1])
  }
}


# Indices vs. uninsured for LPM1
average_LPM2_change_indices_uninsured <- as.data.frame(matrix(0,nrow=1, ncol=length(indices)))
colnames(average_LPM2_change_indices_uninsured) <- indices 

for (k in 1:length(indices)){
  average_LPM2_change_indices_uninsured[1,k] <- mean(LPM2[,k+1] - LPM2[,1])
}

# Indices vs. indices LPM2
average_LPM2_change_indices <- as.data.frame(matrix(0,nrow=length(indices), ncol=length(indices)))
row.names(average_LPM2_change_indices) <- indices
colnames(average_LPM2_change_indices) <- indices    

for (k in 1:length(indices)){
  for (kk in 1:length(indices)){
    
    # Calculate the p-value of the test statistic
    average_LPM2_change_indices[k,kk] <- mean (LPM2[,k+1] - LPM2[,kk+1])
  }
}


# BEST against uninsured for LPM 1 & 2
average_LPM1_change_BEST_uninsured <- mean(best_index_LPM1_value - LPM1[,1])
average_LPM2_change_BEST_uninsured <- mean(best_index_LPM2_value - LPM2[,1])

# BEST against indices for LPM 1 & 2

average_LPM1_change_BEST_indices <- vector(length = length(indices))
names(average_LPM1_change_BEST_indices) <- indices
for (k in 1:length(indices)){
  average_LPM1_change_BEST_indices [k] <- mean(best_index_LPM1_value - LPM1 [,k+1])
}

average_LPM2_change_BEST_indices <- vector(length = length(indices))
names(average_LPM2_change_BEST_indices) <- indices
for (k in 1:length(indices)){
  average_LPM2_change_BEST_indices [k] <- mean(best_index_LPM2_value - LPM2 [,k+1])
}

# --------------------------------------------------------------------------------------------------------------------- 
# 4.3.4 Average relative changes in lower partial moments
# ---------------------------------------------------------------------------------------------------------------------

# Indices vs. uninsured for LPM 1
rel_average_LPM1_change_indices_uninsured <- as.data.frame(matrix(0,nrow=1, ncol=length(indices)))
colnames(rel_average_LPM1_change_indices_uninsured) <- indices 

for (k in 1:length(indices)){
  rel_average_LPM1_change_indices_uninsured[1,k] <- mean((LPM1[,k+1] - LPM1[,1]) / LPM1 [,1])
}

# Indices against each other LPM 1
rel_average_LPM1_change_indices <- as.data.frame(matrix(0,nrow=length(indices), ncol=length(indices)))
row.names(rel_average_LPM1_change_indices) <- indices
colnames(rel_average_LPM1_change_indices) <- indices    

for (k in 1:length(indices)){
  for (kk in 1:length(indices)){
    # Calculate the p-value of the test statistic
    rel_average_LPM1_change_indices[k,kk] <- mean ((LPM1[,k+1] - LPM1[,kk+1]) / LPM1[,kk+1] )
  }
}


# Indices vs. uninsured for LPM2
rel_average_LPM2_change_indices_uninsured <- as.data.frame(matrix(0,nrow=1, ncol=length(indices)))
colnames(rel_average_LPM2_change_indices_uninsured) <- indices 

for (k in 1:length(indices)){
  rel_average_LPM2_change_indices_uninsured[1,k] <- mean((LPM2[,k+1] - LPM2[,1]) / LPM2[,1])
}

# Indices vs. indices LPM2
rel_average_LPM2_change_indices <- as.data.frame(matrix(0,nrow=length(indices), ncol=length(indices)))
row.names(rel_average_LPM2_change_indices) <- indices
colnames(rel_average_LPM2_change_indices) <- indices    

for (k in 1:length(indices)){
  for (kk in 1:length(indices)){
    
    # Calculate the p-value of the test statistic
    rel_average_LPM2_change_indices[k,kk] <- mean ((LPM2[,k+1] - LPM2[,kk+1]) / LPM2[,kk+1])
  }
}


# BEST against uninsured for LPM 1 & 2
rel_average_LPM1_change_BEST_uninsured <- mean((best_index_LPM1_value - LPM1[,1]) / LPM1[,1])
rel_average_LPM2_change_BEST_uninsured <- mean((best_index_LPM2_value - LPM2[,1]) / LPM2[,1])

# BEST against indices for LPM 1 & 2

rel_average_LPM1_change_BEST_indices <- vector(length = length(indices))
names(rel_average_LPM1_change_BEST_indices) <- indices
for (k in 1:length(indices)){
  rel_average_LPM1_change_BEST_indices [k] <- mean((best_index_LPM1_value - LPM1 [,k+1]) / LPM1 [,k+1])
}

rel_average_LPM2_change_BEST_indices <- vector(length = length(indices))
names(rel_average_LPM2_change_BEST_indices) <- indices
for (k in 1:length(indices)){
  rel_average_LPM2_change_BEST_indices [k] <- mean((best_index_LPM2_value - LPM2[,k+1]) / LPM2[,k+1])
}

# =====================================================================================================================
# =====================================================================================================================
#
# End of the script 
#
# =====================================================================================================================
# =====================================================================================================================

# =====================================================================================================================
# =====================================================================================================================
#
# Appendix: Cross-validation
#
# See online supplementary for explanations.
#
# =====================================================================================================================
# =====================================================================================================================

library(reshape2)

# ---------------------------------------------------------------------------------------------------------------------
# A.1 Contract calibration
# ---------------------------------------------------------------------------------------------------------------------

quantile_of_interest <- 0.3

# We create an empty array. Each level of the array represents an index
contract_specifics_crossval <- array(dim=c(length(indices),nrow(farm_yield_detr),ncol=4))
dimnames(contract_specifics_crossval)[[1]] <- indices 
dimnames(contract_specifics_crossval)[[3]] <- c("Intercept QR","tick size", "strike", "premium")

# Years with yield entries
yield_years <- colnames(farm_yield_detr)

# Merging yield and weather data for pooled quantile regression

melted_yield  <- melt(as.matrix(farm_yield_detr))
melted_CPI    <- melt(as.matrix(CPI_df[,yield_years]))
melted_SPI    <- melt(as.matrix(SPI_df[,yield_years]))
melted_SPEI   <- melt(as.matrix(SPEI_df[,yield_years]))
melted_SMI    <- melt(as.matrix(SMI_df[,yield_years]))
melted_ESI    <- melt(as.matrix(ESI_df[,yield_years]))

# Array for historical payouts. Each level of the array represents an index.
historical_payouts_crossval                <- array(dim=c(length(indices),nrow(farm_yield_detr),ncol(farm_yield_detr)))
dimnames(historical_payouts_crossval)[[1]] <- indices 
dimnames(historical_payouts_crossval)[[3]] <- colnames(farm_yield_detr)

# Start of cross-validation contract calibration

for (i in 1:nrow(farm_yield_detr)){
  for (k in 1:length(indices)){
    
    # Tick size
    
    temp_index <- get(paste("melted",indices[k],sep="_"))
    temp_reg   <- rq(melted_yield[melted_yield[,1] != i ,3] ~ temp_index [temp_index[,1] != i,3], tau=quantile_of_interest)
    contract_specifics_crossval[k,i,1] <- temp_reg$coefficients[1]
    contract_specifics_crossval[k,i,2] <- temp_reg$coefficients[2]
    
    # Strike level
    
    contract_specifics_crossval[k,i,3] <- ((quantile(melted_yield[melted_yield[,1] != i,3], na.rm=T, quantile_of_interest, type=1)-(temp_reg$coefficients[1])) / temp_reg$coefficients[2])
    
    
    # Premium
    
    available_data <- which(!is.na(farm_yield_detr[i,])==TRUE)
    temp1 <- seq(1995,2015,1)
    take_years <- as.character(temp1[available_data]) 
    rm(temp1)
    temp_index_pricing <-  get(paste(indices[k],"df",sep="_"))[i,take_years]
    temp_payouts <- vector(length=length(temp_index_pricing))
    names(temp_payouts) <- colnames(temp_index_pricing)
    
    for (l in 1:length(temp_index_pricing)){
    temp_payouts [l] <- as.numeric(wheat_price)*as.numeric(contract_specifics_crossval[k,i,2])*
      max((as.numeric(contract_specifics_crossval[k,i,3]) - as.numeric(temp_index_pricing[l])), 0)
    }
    
    contract_specifics_crossval[k,i,4] <- mean(temp_payouts)
    
    for (l in 1:length(temp_payouts)){
      historical_payouts_crossval[k,i,which(dimnames(historical_payouts_crossval)[[3]]==names(temp_payouts)[l])] <- temp_payouts[l]
    }
    
    rm(temp_payouts, temp_index_pricing)
    
  }
}

# ------------------------------------------------------------------------------------------------------------------
# A.2 Calculation of wealth
# ------------------------------------------------------------------------------------------------------------------

# Revenue without insurance
revenue_crop <- as.data.frame(matrix(NA,nrow=nrow(farm_yield_detr), ncol=ncol(farm_yield_detr)))
colnames(revenue_crop) <- colnames(farm_yield_detr)

for (i in 1:nrow(revenue_crop)){
  revenue_crop[i,] <- farm_yield_detr[i,]*wheat_price
}

# Revenues with insurance
revenue_insured_crop_crossval    <- array(dim=c(length(indices),nrow(farm_yield_detr),ncol(farm_yield_detr)))
dimnames(revenue_insured_crop_crossval)[[3]] <- yield_years

for (k in 1:length(indices)){
  for (i in 1:nrow(farm_yield_detr)){
    
    revenue_insured_crop_crossval[k,i,] <- (as.numeric(revenue_crop[i,]) + as.numeric(historical_payouts_crossval[k,i,]) - as.numeric(contract_specifics_crossval[k,i,4]))
  }
}

# ------------------------------------------------------------------------------------------------------------------
# A.3 Expected utility calculations
# ------------------------------------------------------------------------------------------------------------------

# We create an empty array for each index. Each level of the array represent a CRRA-coefficient
utility_uninsured_crossval <- array(dim=c(length(alpha),nrow(farm_yield_detr),ncol(farm_yield_detr)))
utility_CPI_crossval       <- array(dim=c(length(alpha),nrow(farm_yield_detr),ncol(farm_yield_detr)))
utility_SPI_crossval       <- array(dim=c(length(alpha),nrow(farm_yield_detr),ncol(farm_yield_detr)))
utility_SPEI_crossval      <- array(dim=c(length(alpha),nrow(farm_yield_detr),ncol(farm_yield_detr)))
utility_SMI_crossval       <- array(dim=c(length(alpha),nrow(farm_yield_detr),ncol(farm_yield_detr)))
utility_ESI_crossval       <- array(dim=c(length(alpha),nrow(farm_yield_detr),ncol(farm_yield_detr)))

# Utility of being uninsured
for (a in 1:length(alpha)){
  for (i in 1:nrow(farm_yield_detr)){
    for (t in 1:ncol(farm_yield_detr)){
      utility_uninsured_crossval[a,i,t] <- utility_function(alpha = alpha[a], revenue = revenue_crop[i,t])
    }
  }
}

# Utility for each index 
for (k in 1:length(indices)){
  # Get the correct utility array
  temp <- get(paste("utility",indices[k],"crossval",sep="_"))
  
  for (a in 1:length(alpha)){
    for (i in 1:nrow(farm_yield_detr)){
      for (t in 1:ncol(farm_yield_detr)){
        temp[a,i,t] <- utility_function(alpha = alpha[a], revenue = revenue_insured_crop_crossval[k,i,t]) 
      }
    }
  }
  assign(paste("utility",indices[k],"crossval",sep="_"),temp, envir = .GlobalEnv)
  rm(temp)
}

# ---------------------------------------------------------------------------------------------------------------------
# A.4 Calculating expected utility, certainty equivalent and risk premium
# ---------------------------------------------------------------------------------------------------------------------

# For being uninsured
summary_EU_uninsured_crossval <- array(dim=c(length(alpha),nrow(farm_yield_detr),4))
dimnames(summary_EU_uninsured_crossval)[[1]] <- alpha
dimnames(summary_EU_uninsured_crossval)[[3]] <- c("EU","Expected Revenue", "Certainty Equivalent", "Risk Premium")

for (a in 1:length(alpha)){
  for (i in 1:nrow(farm_yield_detr)){
    
    # Expected Utility 
    summary_EU_uninsured_crossval [a,i,1] <- mean(as.numeric(utility_uninsured_crossval[a,i,]), na.rm=T)
    
    # Expected Revenue
    summary_EU_uninsured_crossval [a,i,2] <- mean(as.numeric(revenue_crop[i,]), na.rm=T) 
    
    # Certainty Equivalent
    summary_EU_uninsured_crossval [a,i,3] <- inverse_utility_function(alpha = alpha[a], eu = as.numeric(summary_EU_uninsured_crossval[a,i,1])) 
    
    # Risk Premium
    summary_EU_uninsured_crossval [a,i,4] <- summary_EU_uninsured_crossval [a,i,2] - summary_EU_uninsured_crossval [a,i,3] 
  }
}

# For each drought index
for (k in 1:length(indices)){
  
  temp_utility_crossval <- get(paste("utility", indices[k],"crossval", sep="_"))
  temp_revenue_crossval <- revenue_insured_crop_crossval [k,,]
  temp_risk_crossval <- array(dim=c(length(alpha),nrow(farm_yield_detr),4)) 
  
  for (a in 1:length(alpha)){
    for (i in 1:nrow(farm_yield_detr)){
      
      # Expected Utility    
      temp_risk_crossval[a,i,1] <- mean(as.numeric(temp_utility_crossval[a,i,]), na.rm=T)
      
      # Expected Revenue
      temp_risk_crossval[a,i,2] <- mean(as.numeric(temp_revenue_crossval[i,]), na.rm=T)
      
      # Certainty Equivalent
      temp_risk_crossval[a,i,3] <- inverse_utility_function(alpha = alpha[a], eu = as.numeric(temp_risk_crossval[a,i,1])) 
      
      # Risk Premium
      temp_risk_crossval[a,i,4] <- temp_risk_crossval [a,i,2] - temp_risk_crossval [a,i,3] 
    }
  }
  dimnames(temp_risk_crossval)[[1]] <- alpha
  dimnames(temp_risk_crossval)[[3]] <- c("EU","Expected Revenue", "Certainty Equivalent", "Risk Premium")
  assign(paste("summary_EU",indices[k],"crossval",sep="_"),temp_risk_crossval, envir = .GlobalEnv)
  rm(temp_risk_crossval)
}

# For an overview and further data handling, we put the risk premiums for each level of risk aversion in a df
report_RP_crossval <- array(dim=c(length(alpha), nrow(farm_yield_detr), ncol=length(indices)+1))
dimnames(report_RP_crossval)[[3]] <- c("uninsured","CPI", "SPI", "SPEI", "SMI", "ESI")

for (a in 1:length(alpha)){
  for (i in 1:nrow(farm_yield_detr)){
    
    # Uninsured
    report_RP_crossval[a,i,1] <- summary_EU_uninsured_crossval[a,i,4]
    
    # Insured with index k
    for (k in 1:length(indices)){
      report_RP_crossval[a,i,k+1] <- get(paste("summary_EU", indices[k],"crossval",sep="_"))[a,i,4] 
    }
  }
}

# ---------------------------------------------------------------------------------------------------------------------
# A.5 Identification of BEST
# ---------------------------------------------------------------------------------------------------------------------

best_index_crossval <- as.data.frame(matrix(NA,nrow=nrow(farm_yield_detr), ncol=length(alpha))) 
colnames(best_index_crossval) <- alpha

# If several indices have the same RP, we choose the simplest one.
for (a in 1:length(alpha)){
  for (i in 1:nrow(farm_yield_detr)){
    temp1 <- which(report_RP_crossval[a,i,] == min(report_RP_crossval[a,i,])) 
    best_index_crossval [i,a] <- as.character(names(temp1[1]))
    rm(temp1)
  }
}

# Assign the RP value of the best drought index
best_index_RP_crossval <- as.data.frame(matrix(NA,nrow=nrow(farm_yield_detr), ncol=length(alpha)))
colnames(best_index_RP_crossval) <- alpha
for (a in 1:length(alpha)){
  for (i in 1:nrow(farm_yield_detr)){
    best_index_RP_crossval[i,a] <- report_RP_crossval[a,i,best_index_crossval[i,a]]
  }
}

# ---------------------------------------------------------------------------------------------------------------------
# A.6 Changes in the risk premium
# ---------------------------------------------------------------------------------------------------------------------

# Indices against uninsured
average_RP_change_uninsured_crossval <- as.data.frame(matrix(NA, nrow=length(indices), ncol=length(alpha)))
row.names(average_RP_change_uninsured_crossval) <- indices
colnames(average_RP_change_uninsured_crossval)  <- alpha

for (a in 1:length(alpha)){
  for (k in 1:length(indices)){
    temp1 <- report_RP_crossval[a,,k+1] - report_RP_crossval[a,,1]
    average_RP_change_uninsured_crossval [k,a] <- mean(as.numeric(temp1))
    rm(temp1)
  }
}

# Indices against each other
average_RP_change_indices_crossval <- array(dim = c(length(alpha), length(indices), length(indices)))
dimnames(average_RP_change_indices_crossval)[[1]] <- alpha
dimnames(average_RP_change_indices_crossval)[[2]] <- indices
dimnames(average_RP_change_indices_crossval)[[3]] <- indices

for (a in 1:length(alpha)){
  for (k in 1:length(indices)){
    for (kk in 1:length(indices)){
      temp1 <- report_RP_crossval[a,,k+1] - report_RP_crossval[a,,kk+1]
      average_RP_change_indices_crossval [a,k,kk] <- mean(as.numeric(temp1))
      rm(temp1)
    }
  }
}

# BEST against uninsured
average_RP_change_BEST_uninsured_crossval <- vector(length= length(alpha))
names(average_RP_change_BEST_uninsured_crossval) <- alpha

for (a in 1:length(alpha)){
  average_RP_change_BEST_uninsured_crossval [a] <- mean(as.numeric(best_index_RP_crossval[,a] - report_RP_crossval[a,,1]))
}

# BEST against the indices

# BEST against indices
average_RP_change_BEST_indices_crossval <- as.data.frame(matrix(NA,nrow=length(alpha), ncol=length(indices)))
row.names(average_RP_change_BEST_indices_crossval) <- alpha
colnames(average_RP_change_BEST_indices_crossval) <- indices

for (a in 1:length(alpha)){
  for (k in 1:length(indices)){
    average_RP_change_BEST_indices_crossval [a,k] <- mean(best_index_RP_crossval[,a] - report_RP_crossval[a,,k+1]) 
  }
}

# ---------------------------------------------------------------------------------------------------------------------
# Relative changes in RP
# ---------------------------------------------------------------------------------------------------------------------

# Indices against uninsured
rel_average_RP_change_uninsured_crossval <- as.data.frame(matrix(NA, nrow=length(indices), ncol=length(alpha)))
row.names(rel_average_RP_change_uninsured_crossval) <- indices
colnames(rel_average_RP_change_uninsured_crossval) <- alpha

for (a in 1:length(alpha)){
  for (k in 1:length(indices)){
    temp1 <- (report_RP_crossval[a,,k+1] - report_RP_crossval[a,,1]) / report_RP_crossval[a,,1] 
    rel_average_RP_change_uninsured_crossval [k,a] <- mean(as.numeric(temp1))
    rm(temp1)
  }
}

# Indices against each other
rel_average_RP_change_indices_crossval                <- array(dim = c(length(alpha), length(indices), length(indices)))
dimnames(rel_average_RP_change_indices_crossval)[[1]] <- alpha
dimnames(rel_average_RP_change_indices_crossval)[[2]] <- indices
dimnames(rel_average_RP_change_indices_crossval)[[3]] <- indices

for (a in 1:length(alpha)){
  for (k in 1:length(indices)){
    for (kk in 1:length(indices)){
      temp1 <- (report_RP_crossval[a,,k+1] - report_RP_crossval[a,,kk+1]) / report_RP_crossval[a,,kk+1]
      rel_average_RP_change_indices_crossval [a,k,kk] <- mean(as.numeric(temp1))
      rm(temp1)
    }
  }
}

# BEST against uninsured
rel_average_RP_change_BEST_uninsured_crossval <- vector(length= length(alpha))
names(rel_average_RP_change_BEST_uninsured_crossval) <- alpha

for (a in 1:length(alpha)){
  rel_average_RP_change_BEST_uninsured_crossval [a] <- mean(as.numeric(best_index_RP_crossval[,a] - report_RP_crossval[a,,1]) / report_RP_crossval[a,,1])
}

# BEST against the indices

# BEST against indices
rel_average_RP_change_BEST_indices_crossval            <- as.data.frame(matrix(NA,nrow=length(alpha), ncol=length(indices)))
row.names(rel_average_RP_change_BEST_indices_crossval) <- alpha
colnames(rel_average_RP_change_BEST_indices_crossval)  <- indices

for (a in 1:length(alpha)){
  for (k in 1:length(indices)){
    rel_average_RP_change_BEST_indices_crossval [a,k]  <- mean(as.numeric(best_index_RP_crossval[,a] - report_RP_crossval[a,,k+1]) / report_RP_crossval[a,,k+1])
  }
}

# ---------------------------------------------------------------------------------------------------------------------
# Significant differences in RP
# ---------------------------------------------------------------------------------------------------------------------

# Testing the indices against being uninsured
indices_vs_uninsured_RP_crossval <- array(dim=c(length(alpha),length(indices),1))
dimnames(indices_vs_uninsured_RP_crossval)[[1]] <- alpha
dimnames(indices_vs_uninsured_RP_crossval)[[2]] <- indices
dimnames(indices_vs_uninsured_RP_crossval)[[3]] <- ("uninsured")

# Calculating the p-values
for (a in 1:length(alpha)){
  for (k in 1:length(indices)){
    
    indices_vs_uninsured_RP_crossval[a,k,1] <- wilcox.test(report_RP_crossval[a,,k+1], report_RP_crossval[a,,1], paired=T, alternative="l")$p.value
  }
}

# Testing the indices against each other
indices_vs_indices_RP_crossval <- array(dim=c(length(alpha),length(indices),length(indices)))
dimnames(indices_vs_indices_RP_crossval)[[1]] <- alpha
dimnames(indices_vs_indices_RP_crossval)[[2]] <- indices
dimnames(indices_vs_indices_RP_crossval)[[3]] <- indices

# Calculating the p values
for (a in 2:length(alpha)){
  for (k in 1:length(indices)){
    for (kk in 1:length(indices)){
      
      # Calculate the p-value of the test statistic
      indices_vs_indices_RP_crossval[a,k,kk] <- wilcox.test(report_RP_crossval[a,,k+1], report_RP_crossval[a,,kk+1], paired=T,alternative="l")$p.value
    }
  }
}

# Testing BEST against being uninsured and the indices
testing_best_RP_crossval <- as.data.frame(matrix(NA, nrow=length(alpha), ncol= (length(indices)+1)))
colnames(testing_best_RP_crossval) <- c("uninsured", "CPI", "SPI", "SPEI","SMI" ,"ESI")
row.names(testing_best_RP_crossval)<- alpha

for (a in 1:length(alpha)){
  for (k in 1:(length(indices)+1))
    testing_best_RP_crossval[a,k] <- wilcox.test(best_index_RP_crossval[,a], report_RP_crossval[a,,k], paired=T, alternative="l")$p.value 
}


# =====================================================================================================================
# =====================================================================================================================
# End of the cross-validation
# =====================================================================================================================
# =====================================================================================================================
