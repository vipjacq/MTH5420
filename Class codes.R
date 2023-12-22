# December 9 class ----
# Set directory to where bradfield.txt is located or use Import dataset.
brad <- scan("bradfield.txt") 
brad2 <- brad[!is.na(brad)] # remove NA in brad and change to brad2
pacf(brad2) # partial autocorrelation plot of brad2
plot(brad[1:87671], brad[2:87672]) # plot of brad against its immediate past

# Modification of the cluster10 function
# This function newcluster outputs the declustered observation given any value of kappa.
newcluster <- function(dataset, threshold, kappa) { 
  x=list() 
  z=list() 
  j=1 
  { 
    for(i in (kappa+1):length(dataset)) { 
      if(all(dataset[i-kappa]>threshold&dataset[i-((kappa-1):0)]<=threshold)) { 
        x=max(dataset[j:i]) 
        ifelse(i !=length(dataset), j<-i+1, NA) 
        z=c(z,x)}}} 
  return(as.numeric(z))
}

library("ismev")
cluster.peaks <- newcluster(brad2, 50, 10) # Cluster peaks of the brad2 data with threshold=50 and clustering parameter kappa=10
A <- gpd.fit(cluster.peaks, 50) 
A$rate <- length(cluster.peaks)/length(brad2)

# This function returns the exact return level given the parameters indicated. M is set at 365.25*24
return.level <- function(threshold, scale, shape, rate, M) {
  threshold+(scale/shape)*((threshold*M*rate)^(shape)-1)
}
ret <- return.level(threshold=50, scale=A$mle[1], shape=A$mle[2], rate=A$rate, M=365.25*24) # The value here is 101.5816.

gpd.prof(A, 50, npy=365.25*24, xlow = 92, xup = 130) 
abline(v=ret, lty=2)

# December 16 class ----
# Import venice.txt first
attach(venice)
head(venice)
colnames(venice)[1] <- c("Year")
colnames(venice)[2] <- c("Annual Maxima")
venice.anmax <- venice$`Annual Maxima`
plot(venice.anmax, type="l")
ti <- as.matrix(seq(1, 51, 1))
Y <- lm(venice.anmax~ti)
abline(Y, lty=2, col='red')

venice.gev <- gev.fit(venice.anmax, ydat=ti, mul=1) # nonstationary GEV
venice.gev.st <- gev.fit(venice.anmax) # stationary GEV
mu.t <- venice.gev$mle[1]+venice.gev$mle[2]*ti
D <- 2*(-venice.gev$nllh-(-venice.gev.st$nllh))

# reject H0 if D > crit(chi)
# Nonstationary model is better than the stationary model

ti2 <- cbind(seq(1,51,1), seq(1,51,1)^2)
venice.gev.quad <- gev.fit(venice.anmax, ydat=ti2, mul=c(1,2), show=F) # nonstationary GEV (quadratic form) 
2*(-venice.gev.quad$nllh-(-venice.gev$nllh)) # value of D

# Do not reject H0. Quadratic model did not improve the model fit.

gev.diag(venice.gev)
gev.prof(venice.gev, m=100, xlow=-10, xup=90, conf=0.95, nint=100)

# December 21 class ----
library(ismev)
data("fremantle")
head(fremantle)
par(mfrow = c(1,2))
plot(fremantle$SeaLevel, xlab='SOI', ylab='Sea level (m)')
plot(fremantle$SOI, type='l', xlab='Year', ylab='Sea level (m)')
par(mfrow = c(1,1))

covar <- cbind(fremantle$SOI, seq(1, 86, 1)) # 2x86 covariance matrix
A <- gev.fit(fremantle$SeaLevel) # stationary model
B <- gev.fit(fremantle[,2], ydat=covar, mul=1, show=F) # SOI covariate
C <- gev.fit(fremantle[,2], ydat=covar, mul=2, show=F) # time covariate
D <- gev.fit(fremantle[,2], ydat=covar, mul=c(1,2), show=F) # both SOI and time

2*(-B$nllh-(-A$nllh)) # D=7.289023, significant improvement over the stationary model
2*(-C$nllh-(-A$nllh)) # D=12.44618, dependence on time gives a more significant improvement than allowing for a dependence on SOI
2*(-D$nllh-(-C$nllh)) # D=8.071955, we should include both time and SOI as covariates

# Diagnostic plots
gev.diag(A)
gev.diag(B)
gev.diag(C)
gev.diag(D)

data("rain")
ti <- matrix(seq(1, 17531, 1), ncol=1, nrow=17531)
f <- gpd.fit(rain, threshold=30, ydat=ti, sigl=1, siglink=exp, show=F)
f$mle # [1] 1.8039221883 0.0000196265 0.1977440505, sigma(t)=exp{1.804 + 0.00002t}
E <- gpd.fit(rain, threshold=30) # stationary model for rain
2*(-f$nllh-(-E$nllh)) # D=0.9841494, there is no evidence of a (linear) trend in the log–scale parameter.
