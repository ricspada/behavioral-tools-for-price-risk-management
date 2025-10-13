#set working directory as folder location

#setwd()

library(lubridate)
library(readxl)

#load and transform data ----
data <- read_excel("CMO-Historical-Data-Monthly.xlsx", sheet = "Monthly Prices")
data1 <- read_excel("CMO-Historical-Data-Monthly.xlsx", sheet = "Monthly Indices")

gas <- as.numeric(unlist(data[380:788 , 9 ]))
dap <- as.numeric(unlist(data[380:788 , 59 ]))
food <- as.numeric(unlist(data1[383:791 , 7 ]))
time <- data[380:788 , 1 ]
time <- ym(time$`World Bank Commodity Price Data (The Pink Sheet)`)

yr0 <- year(time)
mm  <- month(time)
yr <-  yr0 + (mm-1)/12

dat <- cbind(gas, food, dap, yr)
dat <- na.omit(dat)
write.csv(dat, 'data.csv')

