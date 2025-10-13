rm(list = ls())

#setwd()

## load packages ----
library(moments)
library(vars)
library(quantreg)# VAR models
library(tseries)
library(forecast)

#load and transform data ----
dat <- read.csv("data.csv")
yr  <- dat[,5]
dat <- dat[,2:4]
#transform data
#transform Dollar/mmBTU in EUR/GJ (multiplying by mmBTU/GJ and EUR/Dollar (05-04-2024))
dat[ ,1] <- dat[ ,1] * 1.055056 * 1.0837
#transform Dollar/mt in EUR/mt
dat[ ,3] <- dat[ ,3] * 1.0837

nn <- nrow(dat)-1 #number of observations
nL <- 4 # trim for lags
n <- nn-4

#Create variables and lags

pg  <- dat[nL    :nn    ,1]#gas
pg1 <- dat[(nL-1):(nn-1),1]
pg2 <- dat[(nL-2):(nn-2),1]
pg3 <- dat[(nL-3):(nn-3),1]
pd  <- dat[nL    :nn    ,3]#DAP
pd1 <- dat[(nL-1):(nn-1),3]
pd2 <- dat[(nL-2):(nn-2),3]
pd3 <- dat[(nL-3):(nn-3),3]
pf  <- diff(dat[nL    :(nn+1)   ,2])#Food price index
pf1 <- diff(dat[(nL-1):(nn),2])
pf2 <- diff(dat[(nL-2):(nn-1),2])
pf3 <- diff(dat[(nL-3):(nn-2),2])
tt  <- 1:nrow(dat[nL :nn,])
yr  <- yr[nL     :      nn]
pgs1 <- pg1**2
pds1 <- pd1**2
pfs1 <- pf1**2


#adf.test(pd, k = 3)
#adf.test(pg, k = 3)


pp.test(pd)
pp.test(pg)
pp.test(pf)


vardat <- cbind(pg, pd, pf)

#check length
nrow(vardat)
length(pgs1)


# VAR Specification ----
#select var order based on SBIC

a <- VARselect(vardat, lag.max = 10, type = "const")
summary(a)
summary(nvar <- VAR(vardat, p=2))
#add exogenous variables
VARselect(vardat, lag.max =4, exogen = cbind(pgs1))
VARselect(vardat, lag.max =4, exogen = cbind(pds1))
VARselect(vardat, lag.max =4, exogen = cbind(pgs1, pds1))
VARselect(vardat, lag.max =4, exogen = cbind(pgs1, pds1, pfs1))

# Specify and Estimate Q/VAR ----

#define individual marginal specifications
mg <- pg ~ pg1 + pd1 + pf1 + pg2 + pd2 + pf2
md <- pd ~ pg1 + pd1 + pf1 + pg2 + pd2 + pf2
mf <- pf ~ pg1 + pd1 + pf1 + pg2 + pd2 + pf2

#select quantiles for tables
aa <- c(.1, .3, .5, .7, .9)
nqa <- length(aa)

#Caclulate Pseude-R^2 for QVAR
fits1 <-rq(mg, tau = aa)
fits0 <-rq(pg~1, tau = aa)
rho <- function(u,tau=.5)u*(tau - (u < 0))
R1g <- 1 - fits1$rho/fits0$rho
fits1 <-rq(md, tau = aa)
fits0 <-rq(pd~1, tau = aa)
R1f <- 1 - fits1$rho/fits0$rho
fits1 <-rq(mf, tau = aa)
fits0 <-rq(pf~1, tau = aa)
R1o <- 1 - fits1$rho/fits0$rho

R1g
R1f
R1o

#estimate QVAR for integer quantiles
nqb <- 99
ab <- c(1:nqb)/100
qgb <- rq(mg, tau = ab)
qdb <- rq(md, tau = ab)
qfb <- rq(mf, tau = ab)
bgb <- qgb$coeff
bdb <- qdb$coeff
bfb <- qfb$coeff

#create dataset of all independent variables
Xg <- cbind(1, pg1, pd1, pf1, pg2, pd2, pf2)
Xd <- cbind(1, pg1, pd1, pf1, pg2, pd2, pf2)
Xf <- cbind(1, pg1, pd1, pf1, pg2, pd2, pf2)

#create inverse distribution 
yg <- Xg %*% bgb
yg <- t(apply(yg, 1, cummax))
yg1 <- yg
ygcc <- cbind(yg1, yg1[,nqb])
yd <- Xd %*% bdb
yd <- t(apply(yd, 1, cummax))
yd1 <- yd
ydcc <- cbind(yd1, yd1[,nqb])
yf <- Xf %*% bfb
yf <- t(apply(yf, 1, cummax))
yf1 <- yf

# Obtain the copula 
Fg <- array(0, dim=c(n,1))
Fd <- array(0, dim=c(n,1))
Ff <- array(0, dim=c(n,1))

for (ii in 1:n) {
  Fg[ii] <- ab[min(which(min(abs(pg[ii] - yg[ii,])) == abs(pg[ii] - yg[ii,])))]
  Fd[ii] <- ab[min(which(min(abs(pd[ii] - yd[ii,])) == abs(pd[ii] - yd[ii,])))]
  Ff[ii] <- ab[min(which(min(abs(pf[ii] - yf[ii,])) == abs(pf[ii] - yf[ii,])))]
}

## Estimate both conditional copulas to Food price index

#Copula between gas price and FPI marginals
CQb <- Fg ~ Ff
summary(lm(CQb))
rCQb <- rq(CQb , tau=aa)
summary(rCQb, se = "boot", brmethod="xy")
rCQbb <- rq(CQb , tau=ab)
bCQbb <- rCQbb$coeff

#Copula between DAP price and FPI marginals
CQa <- Fd ~ Ff
summary(lm(CQa))
rCQa <- rq(CQa , tau=aa)
summary(rCQa, se = "boot", brmethod="xy")
rCQab <- rq(CQa , tau=ab)
bCQab <- rCQab$coeff

# Simulation ----
Nt <- 183 # number of periods to simulate
Ns <- 1000 # number of draws
Ys <- array(0, dim=c(3, Nt, Ns))

set.seed(132)
# set initial conditions
pg[395:397]
Ys[1, ,] <- pg[n] 
Ys[2, ,] <- pd[n] 
Ys[3, ,] <- pf[n] 

for (it in (3:Nt)) {
  for (is in (1:Ns)) {
    
    pgi1 <- Ys[1, it-1, is]
    pgi2 <- Ys[1, it-2, is]
    pdi1 <- Ys[2, it-1, is]
    pdi2 <- Ys[2, it-2, is]
    pfi1 <- Ys[3, it-1, is]
    pfi2 <- Ys[3, it-2, is]
    
    Xgi  <- cbind(1, pgi1, pdi1, pfi1, pgi2, pdi2, pfi2)
    Xdi  <- cbind(1, pgi1, pdi1, pfi1, pgi2, pdi2, pfi2)
    Xfi  <- cbind(1, pgi1, pdi1, pfi1, pgi2, pdi2, pfi2)
    
    qif  <- floor(100*runif(1, 0.01, 1))
    Fif  <- qif/100
    yf   <- Xfi %*% bfb[,qif]
    XQ   <- cbind(1, Fif)
    
    Qid  <- XQ %*% bCQab[,floor(100*runif(1,.01,1))]  
    qid0 <- ab[min(which(min(abs(Qid - ab)) == (abs(Qid - ab))))]
    qid  <- qid0*100
    yd   <- Xdi %*% bdb[,qid]
    
    Qig  <- XQ %*% bCQbb[,floor(100*runif(1,.01,1))]  
    qig0 <- ab[min(which(min(abs(Qig - ab)) == (abs(Qig - ab))))]
    qig  <- qig0*100
    yg   <- Xgi %*% bgb[,qig]
    
    Ys[1, it, is] <- yg 
    Ys[2, it, is] <- yd
    Ys[3, it, is] <- yf
    
  }}

#set parameters from manuscript appendix
q <- 7.480177 #produced dap in mt
x <- 296.4973 #gas input in GJ
Z <- 4285.687 #Other costs
I <- 1038000
sa <- (83333.33/1000) *245 #savings of disposal costs
int <- 0.03 #Interest rate
Ntt <- 180

R <- array(0, dim=c( Ntt, Ns))
Rd <- array(0, dim=c( Ntt, Ns))
npv <- array(0, dim = c(Ns))
disc <- (1+ int)**(-(1:Ntt/12))#

R[, ] <- Ys[2,1:Ntt,] * q - Ys[1,1:Ntt,] * x - Z #+ sa
for (it in (1:Ntt)) {
  Rd[it,] <-  R[it,] * disc[it]
}

npv[] <- colSums(Rd) - I

# Ornstein-Uhlenbeck ----

dt <- 1

pgd <- pg-pg1
pdd <- pd-pd1

arlmg <- lm(pgd ~ pg1)
arlmd <- lm(pdd ~ pd1)


mug <- (-log(1+as.numeric(arlmg$coefficients[2])))
mud <- (-log(1+as.numeric(arlmd$coefficients[2])))

thetag <- as.numeric(-arlmg$coefficients[1])/(as.numeric(arlmg$coefficients[2]))
thetad <- as.numeric(-arlmd$coefficients[1])/(as.numeric(arlmd$coefficients[2]))

epsilong <- arlmg$residuals
sigmag   <- sd(epsilong)*(log(1+as.numeric(arlmg$coefficients[2]))/((1+as.numeric(arlmg$coefficients[2]))**2-1))**0.5

epsilond <- arlmd$residuals
sigmad   <- sd(epsilond)*(log(1+as.numeric(arlmd$coefficients[2]))/((1+as.numeric(arlmd$coefficients[2]))**2-1))**0.5

rho <- cor(pg, pd)

price <- cbind(pg, pd)

T <- 180                        # Investment horizon in years
num_simulations <- 1000        # Number of simulation paths

initial_values <- as.numeric(price[nrow(price),])

# Simulate revenue paths using GBM no correlation
simulate_price_paths_no <- function(initial_values, mug, thetag, sigmag,
                                    mud, thetad, sigmad, T, num_simulations) {
  
  dt <- 1
  paths <- array(0, c(2, T + 1, num_simulations))
  paths[1, ,1] <- initial_values[1]
  paths[2, ,1] <- initial_values[2]
  set.seed(132)
  for (i in 1:num_simulations) {
    for (t in 2:(T + 1)) {
      dWg <- rnorm(1, mean = 0, sd = sqrt(dt))
      dWd <- rnorm(1, mean = 0, sd = sqrt(dt))
      
      paths[1, t,i] <- paths[1, t-1, i] + mug * (thetag - paths[1, t-1, i]) + sigmag * dWg
      paths[2, t,i] <- paths[2, t-1, i] + mud * (thetad - paths[2, t-1, i]) + sigmad * dWd
    }
  }
  
  return(paths)
}

price_paths_no <- simulate_price_paths_no(initial_values, mug, thetag, sigmag,
                                          mud, thetad, sigmad, T, num_simulations)

sa <- (83333.33/1000) *245

R_no <- array(0, dim=c( Ntt, Ns))
Rd_no <- array(0, dim=c( Ntt, Ns))
npv_no <- array(0, dim = c(Ns))

R_no[, ] <- price_paths_no[2,1:Ntt,] * q - price_paths_no[1, 1:Ntt,] * x - Z #+ sa
for (it in (1:Ntt)) {
  Rd_no[it,] <- R_no[it,] * disc[it]
}

npv_no[] <- colSums(Rd_no[])  - I

# Dynamic Stability Analysis ----
naa <- length(aa)
ncof <- ncol(Xg)
REP <- 1000
rqbs <- array(0, dim=c(3, naa, REP, ncof))
Rev <- array(0, dim= c(REP, naa, naa, naa))
MRev <- array(0, dim= c(REP, naa, naa, naa))
res <- array(0, dim= c( naa, naa, naa))
qMRev <- array(0, dim= c(3, naa, naa, naa))

set.seed(132)

for (ia in (1:naa)){
  rqbs[1,ia,,] <- boot.rq(y= pg, x= cbind(1, pg1 , pd1 , pf1 , pg2 , pd2 , pf2), R= REP, tau= aa[ia])$B
  rqbs[2,ia,,] <- boot.rq(y= pd, x= cbind(1, pg1 , pd1 , pf1 , pg2 , pd2 , pf2), R= REP, tau= aa[ia])$B
  rqbs[3,ia,,] <- boot.rq(y= pf, x= cbind(1, pg1 , pd1 , pf1 , pg2 , pd2 , pf2), R= REP, tau= aa[ia])$B
}

#PM Projection Matrix
for (iag in (1:naa)) {
  for (iad in (1:naa)) {
    for (iaf in (1:naa)) {
      for (r in (1:REP))   {
        PM <- rbind(cbind(rqbs[1, iag, r, 2], rqbs[1, iag, r, 3],
                          rqbs[1, iag, r, 4], rqbs[1, iag, r, 5],
                          rqbs[1, iag, r, 6], rqbs[1, iag, r, 7]),
                    cbind(rqbs[2, iad, r, 2], rqbs[2, iad, r, 3],
                          rqbs[2, iad, r, 4], rqbs[2, iad, r, 5],
                          rqbs[2, iad, r, 6], rqbs[2, iad, r, 7]),
                    cbind(rqbs[3, iaf, r, 2], rqbs[3, iaf, r, 3],
                          rqbs[3, iaf, r, 4], rqbs[3, iaf, r, 5],
                          rqbs[3, iaf, r, 6], rqbs[3, iaf, r, 7]),
                    cbind(diag(3), array(0, dim=c(3,3))) )
        Rev[r, iag, iad, iaf]  <- eigen(PM)$values[1]
        MRev[r, iag, iad, iaf] <- Mod(eigen(PM)$values[1])
      }
      res[iag, iad, iaf]<- mean(Rev[, iag, iad, iaf])
      qMRev[, iag, iad, iaf] <- quantile(MRev[, iag, iad, iaf],c(0.9,0.95,0.99))
    } } }
#extract estimates and standard error with varying [,,,]
t.test(MRev[,3,3,3], mu = 1)

# Visualization ----

pdf('plots/figure_1.pdf',
    family = "serif",
    width = 8, height = 6)

# Adjust margins to allow space for the right-side y-axis label
par(mar = c(5, 5, 4, 5) + 0.1)

# Set up the primary y-axis for gas price (pg)
plot(yr, pg, lwd = 1, xlab = "Year", cex.lab=1.5, cex.axis=1.5,
     ylab = "Gas", ylim = c(min(pg), max(pg) + 50), type = "l")

# Add secondary y-axis for diammonium phosphate price (pd)
par(new = TRUE)
plot(yr, pd, lwd = 1, xaxt = "n", yaxt = "n", xlab = "", ylab = "", 
     ylim = c(min(pd), max(pd) + 350), type = "l", lty = 5)

# Add the right-side y-axis for pd
axis(4, cex.axis=1.5)
mtext("Diammonium Phosphate", side = 4, line = 3, cex=1.5)

# Add legend
legend("topleft", c("Gas price", "Diammonium Phosphate price"), cex=1.8, 
       lwd=2, lty=c(1,5))

dev.off()

pdf(file="plots/figure_2a.pdf",family = "serif", width = 8, height = 5)
plot(ecdf(npv/1e6),
     xlab = "",xlim =c(-2.3, -0.9), axes = F, ann = F, lty = 1  )
box()
dd1 <- ecdf(npv_no/1e6)
lines(dd1, lty=5, col= "gray70")
legend("topleft", c("QVAR","OUP"), cex=1.8, lwd=2.5,  col= c("black","gray70"),
       lty= c(1,1))
pts <- seq(-2.3, -1.1, 0.2)
axis(1, at= pts, labels = paste(pts, "MM", sep = ""))
dev.off()

lim <- rbind(npv, npv_no)
pdf(file="plots/figure_2b.pdf",family = "serif", width = 8, height = 5)
plot(density(npv/1e6, na.rm=T),
     xlab = "",xlim =c(-2.3, -0.9), axes = F, ann = F
     , type = "l" )
box()
dd1 <- density(npv_no/1e6, na.rm=T)
lines(dd1, lty=5)
legend("topleft", c("QVAR","OUP"), cex=1.8, lwd=2.5, 
       lty= c(1,5))
pts <- seq(-2.3, -1.1, 0.2)
axis(1, at= pts, labels = paste(pts, "MM", sep = ""))
dev.off() 



# Necessary disposal saving for min, mean, and sd-tradeoff
min(npv)
mint <- int/12
ml <- min(npv) * (((1+mint)**T *mint) / ((1+mint)**T-1))
be_sa <- ml / 83.33333 

mean(npv)
mint <- int/12
ml <- mean(npv) * (((1+mint)**T *mint) / ((1+mint)**T-1))
be_sa <- ml / 83.33333 


RP <- array(0, c(3,6))
RP[1,] <- seq(0,0.0005,0.0001)
length(seq(0,0.0005,0.0001))
RP[2,] <- 1/2 * RP[1,] * var(npv) - 1/6 * RP[1,]**2 * skewness(npv)
RP[3,] <- mean(npv) - (1/2 * RP[1,] * var(npv) - 1/6 * RP[1,]**2 * skewness(npv))
pts <- seq(-1.5, -2.75, -0.25)


pdf('plots/figure_3.pdf',
    family = "serif",
    width = 8, height = 6)
plot(RP[1,], RP[3,]/1e6, type = "b", pch = 16, ylab = "Certainty Equivalent", xlab = expression("Risk aversion" ~ r[a]), main = "", ylim = range(RP[3,]/1e6), axes = F)
axis(1)
box()
axis(2, at= pts, labels = paste(pts, "MM", sep = " "))
dev.off()

