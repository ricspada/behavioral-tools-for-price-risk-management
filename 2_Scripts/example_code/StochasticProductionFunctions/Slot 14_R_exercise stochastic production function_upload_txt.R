# Risk Analysis and Risk Management in Agriculture: Theory and Applications 2025. Robert Finger
# Empirical Application  - Estimate stochastic production functions
# Goal: you estimate stochastic production functions and visualize and interpret findings 
# We use the following dataset: Stephanie Schmid & Robert Finger (2008). Dataset on simulated maize and winter wheat yields in Switzerland under current and future climate conditions. https://www.research-collection.ethz.ch/handle/20.500.11850/344460 


## Step 1 - read in data from ETH research Collection 

raw_data<-read.table("https://www.research-collection.ethz.ch/bitstream/handle/20.500.11850/344460/maize_base.txt", sep=";", header=TRUE, na.string="NA")

#Note: you may also import data directly in R Studio. 

# We restrict our analysis to a subset of this data based on soil type 1. 
# We initially also limit our first analysis to one input, namely Nitrogen (i.e. Water, W=0)
maize_data_nitrogen_only <-subset(raw_data, raw_data$Dboden==1 & raw_data$Water==0)

##Step 2 - plotting  and exploring the data 

plot(Yield~Nitrogen,data=maize_data_nitrogen_only, main = "Yield Response to Nitrogen raw data",    xlab = "Nitrogen (kg/ha)",    ylab = "Yield (tons/ha)")

# what is mean, standard deviation and skewness of yields?
# hint: to calculate skewness, please load the library: library(e1071)


#Step 3 - conducting regression analysis

# We use a OLS regression to estimate the production function that is the effect of inputs on the expected yield level
# We use nitrogen in levels and  square roots. This allows to reflect non-linearity

reg1 <- lm(Yield ~ Nitrogen + sqrt(Nitrogen), data=maize_data_nitrogen_only)
summary(reg1)

#Question: what is your interpretation of these coefficients? 


#plot the regression line from estimated coefficients for nitrogen 

x<-seq(from=0, to=350, by=1)
predictedyieldnitrogen<- coef(reg1)[1]+coef(reg1)[2]*x+coef(reg1)[3]*x^0.5
plot(Yield~Nitrogen,data=maize_data_nitrogen_only, main = "Yield and Fitted Production Function",    xlab = "Nitrogen (kg/ha)",    ylab = "Yield (tons/ha)")
points(x,predictedyieldnitrogen , type='l', col="red")
legend("topleft", legend = "Fitted Production Function", col = "red", lty = 1, bty = "n")

#step 4 - investigate risk properties graphically
#squared residuals reflect the estimated yield variance
# we now give heteroscedasticity an economic meaning 

plot(maize_data_nitrogen_only$Nitrogen, residuals(reg1)^2, main = "Estimated Yield Variance vs Nitrogen",   xlab = "Nitrogen (kg/ha)", ylab = "Squared Residuals (Proxy for Variance)")

##Question: what can you say about risk properties on yield variance of nitrogen graphically?



#Step 5 - investigate risk properties empirically

#estimate variance function - note we proceed only with square root terms for simplicity. Other specifications could be tested
reg2 <- lm((residuals(reg1))^2 ~ sqrt(Nitrogen), data = maize_data_nitrogen_only)
summary(reg2)
#Question: what is your interpretation of these coefficients?

##plot the regression line from estimated coefficients for yield variance for nitrogen 

x<-seq(from=0, to=350, by=1)
predictedvariancenitrogen<- coef(reg2)[1]+coef(reg2)[2]*x^0.5
plot(maize_data_nitrogen_only$Nitrogen, residuals(reg1)^2, main = "Fitted Yield Variance Function",  xlab = "Nitrogen (kg/ha)", ylab = "Estimated Variance")
points(x,predictedvariancenitrogen , type='l', col="red")



# Step 6 - find optimal levels of inputs use

#define the goal function - certainty equivalent maximization 
# We assume an absolute Arrow Pratt coefficient of risk aversion of 0.002 (the absolute Arrow Pratt coefficient of risk aversion reflect the relative Arrow Pratt coefficient of risk aversion divided by an expected wealth proxy) --> e.g. 2/1000=0.002
#Price of Maize is  300 Fr/ton; the price of nitrogen is 1Fr/kg; 

#optimal nitrogen use - using R function optimize 

CENitrogen <- function(N) {
  price_maize <- 300
  price_nitrogen <- 1
  risk_aversion <- 0.002
  expected_yield <- coef(reg1)[1] + coef(reg1)[2]*N + coef(reg1)[3]*sqrt(N)
  yield_var <- coef(reg2)[1] + coef(reg2)[2]*sqrt(N)
  CE <- price_maize * expected_yield - price_nitrogen * N - 0.5 * risk_aversion * price_maize^2 * yield_var
  return(CE)
}

optimize(CENitrogen,c(1,500),maximum=T)

# Questions you shall address 

#layer 1

# 1 How large are the differences in optimal Nitrogen use between a risk neutral decision maker and a risk averse decision maker with an absolute Arrow Pratt coefficient of risk aversion of 0.002?
# 2 What input levels would a more risk averse decision maker with an absolute Arrow Pratt coefficient of risk aversion of 0.004 chose?
# 3 How do differences in optimal input use between a risk neutral decision maker and a risk averse decision maker with an absolute Arrow Pratt coefficient of risk aversion of 0.002 change if the price for maize is 500FR/t instead of 300Fr/t?


#layer 2

# 4 what is the role of crop yield skewness effects? 

#Note: cubic residuals reflect the yield skewness
# Visual presentation 
plot(maize_data_nitrogen_only$Nitrogen, (reg1$residuals)^3, main = "Estimated Yield Skewness vs Nitrogen",  xlab = "Nitrogen (kg/ha)", ylab = "Cubic Residuals (Proxy for Skewness)")

##Question: what can you say about risk properties on yield skewness graphically?

#Estimate skewness function
reg3 <- lm((residuals(reg1))^3 ~ sqrt(Nitrogen), data = maize_data_nitrogen_only)
summary(reg3)
#Question: what is your interpretation of these coefficients? 
#next: these could be integrated in optimization problems?
#e.g. for a downside coefficient of risk aversion of -2/1000^2=0.000002 (see Finger 2013) 

# see also: Finger, R. (2013). Expanding risk consideration in integrated models - the role of downside risk aversion in irrigation decisions. Environmental Modelling & Software 43: 169-172 


#layer 3
# what about using water instead of nitrogen? --> start with a dataset with only water, dat<-subset(x, x$Dboden==1 & x$Nitrogen==0)
