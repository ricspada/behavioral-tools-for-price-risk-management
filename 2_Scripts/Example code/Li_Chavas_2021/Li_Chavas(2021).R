## This code is used to replicate the main results of Li and Chavas (2021)

setwd("C:/Users/lijian/Desktop")

library(quantreg)
library(moments)
library(rlist)
library(stargazer)
library(Qtools)

dat1 <- read.csv("data_soybean.csv", header = TRUE )
dat1[,c(4,5)] <- log(dat1[,c(4,5)])

nn <- nrow(dat1)
n  <- nn-2
t1 <- dat1[,1]

nL  <- 3
yr0 <- dat1[nL:nn,1]
ww  <- dat1[nL:nn,3]
yr  <- yr0 + (ww-1)/52
tt  <- yr-1980
tt1 <- 0*(tt<18)+(tt-18)*(tt>=18)
tt2 <- 0*(tt<25)+(tt-25)*(tt>=25)
tt <- tt/10
mm <- dat1[nL:nn,2]
Q1 <- (mm==1) + (mm==2) + (mm==3)
Q2 <- (mm==4) + (mm==5) + (mm==6)
Q3 <- (mm==7) + (mm==8) + (mm==9)
Q4 <- (mm==10) + (mm==11) + (mm==12)

ps  <- dat1[nL    :nn    ,4]
ps1 <- dat1[(nL-1):(nn-1),4]
ps2 <- dat1[(nL-2):(nn-2),4]
pf  <- dat1[nL    :nn    ,5]
pf1 <- dat1[(nL-1):(nn-1),5]
pf2 <- dat1[(nL-2):(nn-2),5]
ma  <- dat1[nL    :nn    ,6]
ma1 <- dat1[(nL-1):(nn-1),6]
maf1  <- ma*pf1
ma1f1 <- ma1*pf1
ma1s1 <- ma*ps1
maf2 <- ma*pf2
ps1s <- (ps1-median(ps))^2
ps2s <- (ps2-median(ps))^2
pf1s <- (pf1-median(pf))^2
pf2s <- (pf2-median(pf))^2

pdf('figure1.pdf',
    width = 10, height = 6)
plot(yr, ps, lwd = 1, xlab = "year", cex.lab=1.5,cex.axis=1.5, ylab = "log price", ylim = c(min(ps), max(pf)), col="red", type = "l")
lines(yr, pf, lwd =1, col="blue" ,lty=5)
legend("topleft", c("Spot price", "Futures price"), cex=1.8, lwd=2, 
       col=c("red", "blue"), lty=c(1,5))
dev.off() 

ms  <- ps ~ ps1 + pf1 + ps2 + pf2 + ps1s + ps2s + ma1 + ma1f1 + Q1+Q2+Q3  +tt +tt1 +tt2 
mf  <- pf ~ ps1 + pf1 + ps2 + pf2 + pf1s + pf2s + ma1 + ma1f1 + Q1+Q2+Q3  +tt +tt1 +tt2 

summary(lm(ms))
summary(lm(mf))

aa <- c(.1, .3, .5, .7, .9)
nqa <- length(aa)
qsa <- rq(ms, tau = aa)
qfa <- rq(mf, tau = aa)
sum1 <- summary(qsa, se = "boot", brmethod="xy" )
sum2 <- summary(qfa, se = "boot", brmethod="xy" )
sum1
sum2
# goodness of fit
fits1 <-rq(ms, tau = aa)
fits0 <-rq(ps~1, tau = aa)
rho <- function(u,tau=.5)u*(tau - (u < 0))
R1s <- 1 - fits1$rho/fits0$rho
R1s
fitf1 <-rq(mf, tau = aa)
fitf0 <-rq(pf~1, tau = aa)
R1f <- 1 - fitf1$rho/fitf0$rho
R1f

bsa <- qsa$coeff
bfa <- qfa$coeff

nqb <- 99
ab <- c(1:nqb)/100
qsb <- rq(ms, tau = ab)
qfb <- rq(mf, tau = ab)
bsb <- qsb$coeff
bfb <- qfb$coeff

# evaluating roots of [d Y_t / d Y_t-1] 
ns <- 3
Rcn1 <- array(0, dim=c(ns, ns, nqb, nqb ) )
Rcn2 <- Rcn1

Qps <- quantile(ps, c(.1, .5, .9)) - median(ps)
Qpf <- quantile(pf, c(.1, .5, .9)) - median(pf)
Qma <- quantile(ma, c(.1, .5, .9)) 
Qq <- c(10, 30, 50, 70, 90)

# scenarios: scp = price scenarios, scm = ma scenarios, iq1 = sp quantiles, iq2 = fp quantiles
for (scp in (1:ns)) {
  # scp = price scenarios : 1=low, 2=med, 3=high prices for (ps, pf)
  for (scm in (1:ns)) {
    # scm = ma scenarios : 1=low, 2=med, 3=high ma
    for (iq1 in 1:nqb) {
      for (iq2 in 1:nqb) {
        Ac <- rbind(cbind(bsb[2,iq1]+2*bsb[6,iq1]*Qps[scp], bsb[3,iq1]+bsb[9, iq1]*Qma[scm], bsb[4,iq1]+2*bsb[7,iq1]*Qps[scp], bsb[5,iq1]),
                    cbind(bfb[2,iq2], bfb[3,iq2]+2*bfb[6,iq2]*Qpf[scp]+bfb[9, iq2]*Qma[scm], bfb[4,iq2], bfb[5,iq2]+2*bfb[7,iq2]*Qpf[scp]),
                    cbind(diag(2), 0*diag(2)) )
        Rcn1[scp, scm, iq1, iq2 ] <- eigen(Ac)$values[1]
        Rcn2[scp, scm, iq1, iq2 ] <- eigen(Ac)$values[2]
      } }
  } }


# SAMPLE SIMULATION for ps and pf
t90 <- min(which(t1 >= 1990))
t00 <- min(which(t1 >= 2000))
t10 <- min(which(t1 >= 2010))
t15 <- min(which(t1 >= 2015))
ts <- c(t90, t00, t10)
Xs <- cbind(1,  ps1, pf1, ps2, pf2, ps1s, ps2s, ma1, ma1f1, Q1, Q2, Q3,  tt, tt1, tt2 )
Xf <- cbind(1,  ps1, pf1, ps2, pf2, pf1s, pf2s, ma1, ma1f1, Q1, Q2, Q3,  tt, tt1, tt2  )
ys <- Xs %*% bsb
ys <- t(apply(ys, 1, cummax))
ys1 <- ys
yscc <- cbind(ys1, ys1[,nqb])
yf <- Xf %*% bfb
yf <- t(apply(yf, 1, cummax))
yf1 <- yf
yfcc <- cbind(yf1, yf1[,nqb])

pq <- rep(1/(nqb+1), nqb+1)
Myf <- rowSums(pq*yfcc)
Vyf <- rowSums(pq*(yfcc-Myf)^2)
Syf <- rowSums(pq*(yfcc-Myf)^3)
Kyf <- rowSums(pq*(yfcc-Myf)^4)
MMf <- cbind(Myf, sqrt(Vyf), Syf/(Vyf^1.5), Kyf/(Vyf^2))
apply(MMf,2,mean)

pdf('figure_b3.pdf')
par(mar=c(5.1,4.1,4.1,4.1))
plot(   yr, MMf[,1], lwd = 3, xlab = "year", ylab = "mean / skewness / kurtosis", 
        cex.axis=1.5, cex.lab=1.5,cex.main=2, ylim = c(-2,12), col="red", type = "l",main = "futures price")
lines( yr, MMf[,3], lwd = 3, col="green" ,lty=3)
lines( yr, MMf[,4], lwd = 3, col="brown",lty=4)
legend("topleft", c("mean", "stand. dev.","Skewness","Kurtosis"), cex=1.5, lwd=2.5,
       col=c("red","blue","green","brown"), lty=c(1,2,3,4))
par(new = TRUE)
plot(   yr, MMf[,2], lwd = 1.5, ylim = c(-0.05,0.3), xaxt = "n", yaxt = "n", cex.axis=1.5,
        col="blue", type = "l", ylab = "", xlab = "", lty=2)
abline(a=0,b=0)
axis(side=4)
mtext("stand. dev.", cex=1.5, side=4, line=3)
dev.off()

min5 <- min(ys1[t90,], ys1[t00,], ys1[t10,])
max5 <- max(yf1[t90,], yf1[t00,], yf1[t10,])
dif <- max((max(yf1[t90,])- min(ys1[t90,])),(max(yf1[t00,])- min(ys1[t00,])),(max(yf1[t10,])- min(ys1[t10,])) )
pdf('figure_b1.pdf')
plot(   ys1[t10,], ab, lwd = 3, xlab = "price", ylab = "quantile", xlim = c(min(ys1[t10,]), min(ys1[t10,])+dif), 
        col="red", type = "l", main="2010", cex.lab=1.8, cex.axis=1.5,cex.main=2)
lines( yf1[t10,], ab, lwd = 3, col="blue",lty=2 )
legend("bottomright", c( "pf","ps"), cex=2, lwd=2.5,lty=c(2,1),
       col=c("blue","red"))
dev.off()

qq1 <- c(10,25, 50, 75,90)
ry1 <- ys1[,qq1[3]]/ys1[,qq1[3]]
ry2 <- ys1[,qq1[1]]/ys1[,qq1[3]]
ry3 <- ys1[,qq1[2]]/ys1[,qq1[3]]
ry4 <- ys1[,qq1[4]]/ys1[,qq1[3]]
ry5 <- ys1[,qq1[5]]/ys1[,qq1[3]]

pdf('figure_b2.pdf')
plot( yr, ry1, xlim = c(min(yr), max(yr)), ylim = c(0.94,1.06), xlab = "year", ylab = "Relative ratio",
      lwd=2, col="red", type="l", main="Spot price", cex.lab=1.5, cex.axis=1.5,cex.main=2)
lines(yr, ry2, lwd=1, col="green", lty = 2)
lines(yr, ry3, lwd=2, col="blue", lty = 3)
lines(yr, ry4, lwd=2, col="blue", lty = 4)
lines(yr, ry5, lwd=1, col="green", lty = 5)
legend("bottomleft", 
       c("median","0.1 quantile","0.25 quantile","0.75 quantile","0.9 quantile"), 
       col=c("red","green","blue","blue","green"), cex=1.3, lwd = 2, lty=c(1,2,3,4,5))
dev.off()


# Obtain the copula 
Fs <- array(0, dim=c(n,1))
Ff <- array(0, dim=c(n,1))
for (ii in 1:n) {
  Fs[ii] <- ab[min(which(min(abs(ps[ii] - ys[ii,])) == abs(ps[ii] - ys[ii,])))]
  Ff[ii] <- ab[min(which(min(abs(pf[ii] - yf[ii,])) == abs(pf[ii] - yf[ii,])))]
}

## Estimate conditional copula: Fs (Ff, ma) 
Fs2 <- (Fs-.5)^2
Fsma <- Fs*ma

CQ <- Ff ~ Fs + Fs2 + ma
summary(lm(CQ ))
rCQ <- rq(CQ , tau=aa)
summary(rCQ, se = "boot", brmethod="xy")
rCQb <- rq(CQ , tau=ab)
bCQb <- rCQb$coeff

# FORWARD SIMULATION for ps and pf
nts <- length(ts)
Qma <- quantile(ma, c(0, .25, .5, .75, 0.95)) 
nma <- length(Qma)
Nt <- 150
Ns <- 50
MOps <- array(0, dim=c(nts,nma, 2,Nt))
MOpf <- MOps
MObas <- array(0, dim=c(nts,nma, Nt))
QRbas <- array(0, dim=c(nts,nma, Nt, Ns))
Ys <- array(0, dim=c(2, nts, nma, Nt, Ns,3))

# set initial conditions : 1=t90, 2=t00, 3=t10 for (ps, pf)
for (its in (1:nts)) {
  Ys[1, its, , , ,] <- ps[ts[its]] 
  Ys[2, its, , , ,] <- pf[ts[its]]
}

Ys[2, , , , , ] <- pf[ts]
dc <- 1.3
Ys[1, , , 2, ,2] <- Ys[1, , , 2, ,2]*dc
Ys[2, , , 2, ,3] <- Ys[2, , , 2, ,3]*dc

set.seed(12345)
for (its in (1:nts)) {
  # its = initial conditions : 1=t90, 2=t00, 3=t10 for (ps, pf)
  for (ima in (1:nma)) {
    # ima = ma scenarios : 1=low, 2=med, 3=high ma
    for (it in (3:Nt)) {
      for (is in (1:Ns)) {
        psi1 <- Ys[1,its,ima,it-1,is,]
        psi2 <- Ys[1,its,ima,it-2,is,]
        pfi1 <- Ys[2,its,ima,it-1,is,]
        pfi2 <- Ys[2,its,ima,it-2,is,]
        Xsi  <- rbind(cbind(1,  psi1[1], pfi1[1], psi2[1], pfi2[1], (psi1[1]-median(ps))^2,  (psi2[1]-median(ps))^2, Qma[ima], Qma[ima]*pfi1[1], 
                            Q1[its], Q2[its], Q3[its],  tt[its], tt1[its], tt2[its] ),
                      cbind(1,  psi1[2], pfi1[2], psi2[2], pfi2[2], (psi1[2]-median(ps))^2,  (psi2[2]-median(ps))^2, Qma[ima], Qma[ima]*pfi1[2], 
                            Q1[its], Q2[its], Q3[its],  tt[its], tt1[its], tt2[its] ),
                      cbind(1,  psi1[3], pfi1[3], psi2[3], pfi2[3], (psi1[3]-median(ps))^2,  (psi2[3]-median(ps))^2, Qma[ima], Qma[ima]*pfi1[3], 
                            Q1[its], Q2[its], Q3[its],  tt[its], tt1[its], tt2[its] ))
        Xfi  <- rbind(cbind(1,  psi1[1], pfi1[1], psi2[1], pfi2[1], (pfi1[1]-median(pf))^2,  (pfi2[1]-median(pf))^2, Qma[ima], Qma[ima]*pfi1[1], 
                            Q1[its], Q2[its], Q3[its],  tt[its], tt1[its], tt2[its] ),
                      cbind(1,  psi1[2], pfi1[2], psi2[2], pfi2[2], (pfi1[2]-median(pf))^2,  (pfi2[2]-median(pf))^2, Qma[ima], Qma[ima]*pfi1[2], 
                            Q1[its], Q2[its], Q3[its],  tt[its], tt1[its], tt2[its] ),
                      cbind(1,  psi1[3], pfi1[3], psi2[3], pfi2[3], (pfi1[3]-median(pf))^2,  (pfi2[3]-median(pf))^2, Qma[ima], Qma[ima]*pfi1[3], 
                            Q1[its], Q2[its], Q3[its],  tt[its], tt1[its], tt2[its] ))
        qis <- floor(100*runif(1,.01,1))
        Fis <- qis/100            
        ys <- Xsi %*% bsb[,qis]
        XQ <- cbind(1, Fis, (Fis-.5)^2, Qma[ima] ) 
        Qif <- XQ %*% bCQb[,floor(100*runif(1,.01,1))]  
        qif0 <- ab[min(which(min(abs(Qif - ab)) == (abs(Qif - ab))))]
        qif <- qif0*100
        yf <- Xfi %*% bfb[,qif]
        Ys[1, its, ima, it, is,] <- ys  
        Ys[2, its, ima, it, is,] <- yf
        
        Acs <- rbind( cbind( bsb[2,qis]+2*bsb[6,qis]*(psi1[1]-median(ps)), bsb[3,qis]+bsb[9, qis]*Qma[ima], 
                             bsb[4,qis]+2*bsb[7,qis]*(psi2[1]-median(ps)), bsb[5,qis] ),
                      cbind( bfb[2,qif], bfb[3,qif]+2*bfb[6,qif]*(pfi1[1]-median(pf))+bfb[9, qif]*Qma[ima], 
                             bfb[4,qif], bfb[5,qif]+2*bfb[7,qif]*(pfi2[1]-median(pf)) ),
                      cbind( diag(2), 0*diag(2) ) )
        
      } } 
    Mps <- apply(Ys[1, its, ima, , , 1], 1, mean)
    Sps <- apply(Ys[1, its, ima, , , 1], 1, sd)
    Mpf <- apply(Ys[2, its, ima, , , 1], 1, mean)
    Spf <- apply(Ys[2, its, ima, , , 1], 1, sd)
    
    MOps[its,ima,,] <- rbind(Mps, Sps)
    MOpf[its,ima,,] <- rbind(Mpf, Spf)
    MObas[its,ima,] <- Mpf-Mps 
    QRbas[its,ima,,] <- Ys[2, its, ima, , , 1]-Ys[1, its, ima, , , 1]
  } }
CQRbas <-apply (QRbas[3,1,,], 1, quantile, c(0.1,0.5,0.9))
pdf('figure_3a.pdf')
plot (CQRbas[2,], type="l", xlab= "time", ylab="basis",  ylim=c(-0.5*max(CQRbas), 1.2*max(CQRbas)),col="red"
      , cex.lab=1.5, cex.axis=1.5,cex.main=2, main="ma_low")
lines(CQRbas[1,],col="green",lty=6)
lines(CQRbas[3,],col="blue",lty=5)
abline(a=0,b=0)
legend("bottomleft", c("q=0.1", "q=0.5","q=0.9"), cex=1.8, lwd=2.5, 
       col=c("green","red","blue"), lty=c(6,1,5))
dev.off() 


IRpsps <- apply(Ys[1, 1, , 2:Nt, , 2], c(1,2),  quantile, c(0.1, 0.5, 0.9))-apply(Ys[1, 1, ,2:Nt , , 1], c(1,2),  quantile, c(0.1, 0.5, 0.9))
IRpspf <- apply(Ys[2, 1, , 2:Nt, , 2], c(1,2),  quantile, c(0.1, 0.5, 0.9))-apply(Ys[2, 1, ,2:Nt , , 1], c(1,2), quantile, c(0.1, 0.5, 0.9))
IRpfpf <- apply(Ys[2, 1, , 2:Nt, , 3], c(1,2),  quantile, c(0.1, 0.5, 0.9))-apply(Ys[2, 1, ,2:Nt , , 1], c(1,2),  quantile, c(0.1, 0.5, 0.9))
IRpfps <- apply(Ys[1, 1, , 2:Nt, , 3], c(1,2),  quantile, c(0.1, 0.5, 0.9))-apply(Ys[1, 1, ,2:Nt , , 1], c(1,2),  quantile, c(0.1, 0.5, 0.9))
max1 <- max(IRpsps, IRpspf,IRpfpf, IRpfps)
min1 <- min(IRpsps, IRpspf,IRpfpf, IRpfps)

pdf('figure_2a.pdf')
plot (IRpsps[2,3,],lwd=2, type="l", ylim=c(-.2,0.8),col="blue",main="Impulse Response of ps on ps and pf",
      xlab = "time", ylab = "impulse response", cex.lab=1.5, cex.axis=1.5,cex.main=1.5 )
abline(a=0, b=0)
lines(IRpsps[1,3,],lwd=2, col="blue",lty=2)
lines(IRpsps[3,3,],lwd=2, col="blue",lty=3)
lines(IRpspf[2,3,],lwd=2,  col="green",lty=1)
lines(IRpspf[1,3,],lwd=2,  col="green",lty=2)
lines(IRpspf[3,3,],lwd=2,  col="green",lty=3)
legend("topright", c("IR_psps, q=0.5", "IR_psps, q=0.1", "IR_psps, q=0.9", "IR_pspf, q=0.5", "IR_pspf, q=0.1", "IR_pspf, q=0.9"), cex=1.5,
       lwd=c(2,2,2,2,2,2), col=c("blue","blue","blue","green","green","green"), lty=c(1,2,3,1,2,3))
dev.off() 


