# 1st quarter ----

wassaw <- c(8.5, 8.9, 9.1, 8.9, 8.4, 9.7, 9.1, 9.6, 8.7, 9.3,
             9.6, 9.3, 8.7, 9.0, 8.8, 8.9, 8.9, 12.2, 7.8, 7.7,
             8.3, 8.1, 7.3, 6.8, 6.7, 7.3, 7.6, 8.2, 8.6, 9.8,
             9.5, 7.4, 7.3, 10.2, 10.3, 10.4, 8.8, 9.7, 10.0, 10.8,
             11.1, 12.7, 11.5, 11.8, 12.6, 13.0, 10.5, 10.5, 10.0, 9.4)

#log likelihood
gev.loglik <- function(theta, dataset) {
  mu <- theta[1]
  sigma <- theta[2]
  xi <- theta[3]
  m <- min((1 + (xi * (dataset - mu) / sigma)))
  if (m < 0.00001) return(1e6)
  if (sigma < 0.00001) return(1e6)
  if (xi == 0) {
    loglik <- -length(dataset) * log(sigma) - sum((dataset - mu) /   sigma) - 
      sum(exp(-((dataset - mu) / sigma)))
  } else {
    loglik <- -length(dataset) * log(sigma) - (1 / xi + 1) * sum(log(1 + (xi * (dataset - mu) / sigma))) -
      sum((1 + (xi * (dataset - mu) / sigma))^(-1 / xi))
  }
  return(-loglik)
}


# Initial parameter values
theta.init <- c(mean(wassaw), sd(wassaw), 0.01)

# Optimization using nlminb
opt <- nlminb(start = theta.init, objective = gev.loglik, dataset = wassaw)
opt$par

# Load the extRemes package
library(extRemes)
fit <- fevd(x=wassaw, type="GEV", method="MLE")
fit
sqrt(0.043872195)
fit$results
inf.mat<-fit$results$hessian #information matrix
varcov<-solve(inf.mat) #inverse of information matrix
diag(varcov) #variance of estimated parameters
sqrt(diag(varcov)) #standard error of MLE estimates


fit2 <- fevd(x=wassaw, type="GEV", method="Lmoments")
fit2
fit2$results
ci(fit2, type="parameter")

n<-length(wassaw)
emp<-(1:n)/(n+1) #empirical cumulative probabilities
library(extRemes)
#cumulative probabilities based of GEV
mod<-pevd(sort(wassaw), loc = 8.7113, scale = 1.3115, shape = -0.1084, type = "GEV") 
plot(emp,mod,main="Probability plot",xlab="Empirical",
     ylab="Model")
#adding a straight line
abline(0, 1, col = "red", lwd=2) 
ci(fit, type="parameter")

#quantiles based of GEV
mod.q <- qevd(emp,loc=8.7113, scale=1.3115, shape=0.1084, type="GEV")
#plot model quantile vs ordered obs.
plot(sort(wassaw), mod.q, main="QQ plot", xlab ="Ordered obs.", ylab ="Model")
#adding a straight line
abline(0, 1, col = "red", lwd =2)

fit <- fevd(wassaw, type="GEV", method="MLE") #MLE
plot(fit, "probprob") #prob. plot
plot(fit, "qq") #Q Q plot
plot(fit, "hist", ylim =c(0,0.4)) #histogram against GEV density
fit2 <- fevd(wassaw, type="GEV", method="MLE") #MLE
plot(fit2,"probprob") #prob. plot
plot(fit, "qq") #Q Q plot
plot(fit, "hist", ylim =c(0,0.4)) #histogram against GEV density


fit <- fevd(x=wassaw, type="GEV", method="MLE")
return.level(fit,return.period=c(10, 100))
rlevd(c(10, 100, 1000), loc = 8.7112817, scale = 1.3114838, shape = -0.1084465)
inf.mat <- fit$results$hessian
varcov <- solve(inf.mat) #inverse of information matrix
par <- fit$results$par
par[1]; par[2]; par[3]

y10 <- -log(1-(1/10))
del <- matrix(ncol=1,nrow=3)
del[1,1] <- 1
del[2,1] <- -((par[3])^(-1))*(1-(y10^(-par[3])))
del[3, 1] <- ((par[2])*((par[3])^(-2))*(1-((y10)^(-par[3]))))-((par[2])*((par[3])^(-1))*((y10)^(-(par[3])))*log(y10))
del.transpose = t(del)
var.rl <- del.transpose%*%varcov%*%del #variance of RL
var.rl
SE10 <- as.numeric(sqrt(var.rl));SE10 #standard error of RL

y100 <- -log(1-(1/100))
del100 <- matrix(ncol=1,nrow=3)
del100[1,1] <- 1
del100[2,1] <- -((par[3])^(-1))*(1-(y10^(-par[3])))
del100[3, 1] <- ((par[2])*((par[3])^(-2))*(1-((y100)^(-par[3]))))-((par[2])*((par[3])^(-1))*((y100)^(-(par[3])))*log(y100))
del.transpose100 = t(del100)
var.rl100 <- del.transpose100%*%varcov%*%del100 #variance of RL
var.rl100
SE100 <- as.numeric(sqrt(var.rl100));SE100 #standard error of RL


rl10<-rlevd(10, loc = 8.7112817, scale = 1.3114838, shape = -0.1084465)
low10<-rl10-1.96*SE10
upp10<-rl10+1.96*SE10
cat(low10,upp10)

fit <- fevd(x=wassaw, type="GEV", method="MLE")
ci(fit,return.period = 10, method = "normal") #extRemes package

# 100 year return level
rl100 <- rlevd(100, loc=8.7112817, scale=1.3114838, shape=-0.1084465)
low100 <- rl100-1.96*SE100
upp100 <- rl100+1.96*SE100
cat(low100,upp100)

fit <- fevd(x=wassaw, type="GEV", method="MLE")
ci(fit,return.period = 100, method = "normal") #extRemes package

#plot RL
rp<-seq(5,100,5)
length(rp)
plot(fit, type = "rl", rperiods=rp) #extRemes package

rl<-rlevd(rp, loc = 8.7112817, scale = 1.3114838, shape = -0.1084465)
ci<-ci(fit, type = "return.level", return.period = rp)
ci[1:20]  #lower 95% CI
ci[21:40] #RL
ci[41:60] #upper 95% CI
plot(rp,rl,xlab="Return period", ylab="Return level", xlim=c(5,100), ylim=c(9, 16))
lines(rp, ci[1:20] , type="l", col="red")
lines(rp, ci[41:60] , type="l", col="red")

lake<-c(333, 213, 790, 343, 351, 521, 307, 305, 352, 277, 319, 319, 
        339, 262, 285, 297, 327, 620, 350, 545, 258)
fit3 <- fevd(x=lake, type="GEV", method="MLE")
#delta/normal
ci(fit3, type = "return.level", return.period = c(100,1000))
#bootstrap
ci(fit3, type = "return.level", return.period = c(100,1000), method="boot")
#profile likelihood
ci(fit3, type = "return.level", method = "proflik", return.period = 100,
   xrange = c(600,1800), verbose=T)
ci(fit3, type = "return.level", method = "proflik", return.period = 1000,
   xrange = c(500,10400), verbose=T)

# Extract parameter estimates
mu <- fit$results$par[1]
sigma <- fit$results$par[2]
xi <- fit$results$par[3]

# Specify a return level (replace this with the return level you're interested in)
z <- 11.33

# Calculate the return period for the given return level
return_period <- 1 / (1 - exp(-((1 + xi * ((z - mu) / sigma)) ^ (-1 / xi))))
as.numeric(return_period)

(rt <- 1/(1-pevd(z, loc=mu, scale=sigma, shape=xi, type="GEV")))
fit5 <- fevd(x=wassaw, type="Gumbel", method="MLE")
ks.test(wassaw,pevd, loc=8.636136, scale=1.274520)
return.level(fit5, return.period=c(10,100), do.ci=T)

#block minima
data(Denmint)
head(Denmint)
Denmint$Minneg=-Denmint$Min
minDenmint<-blockmaxxer(Denmint, blocks=Denmint$Year, which="Minneg")
head(minDenmint)
hist(minDenmint$Minneg)
fit4 <- fevd(x=minDenmint$Minneg, type="GEV", method="MLE")
fit4
fit7 <- fevd(x=minDenmint$Minneg, type="GEV", method="MLE")
fit7
ks.test(minDenmint$Minneg, pevd, loc=8.6151698, scale=6.1164944, shape=-0.1716314)

plot(fit4,"probprob")
plot(fit4,"qq")
plot(fit4,"hist",ylim=c(0,0.1),lwd=1,col="grey")


F_gpd <- function(x, u, xi, sigma) {
  if (xi != 0) {
    return(1 - (1 + xi * (x - u) / sigma)^(-1 / xi))
  } else {
    return(1 - exp(-(x - u) / sigma))
  }
}

F_gpd(2, u=1, sigma=0.5, xi=0.8)
F_gpd(10, u=1, sigma=0.5, xi=0.8)
F_gpd(5:10, u=1, sigma=0.5, xi=0.8)

library(extRemes)
pevd(2, scale=0.5, shape=0.8, threshold=1, type="GP")
pevd(5:10, scale=0.5, shape=0.8, threshold=1, type="GP")
devd(2:4, scale=0.5, shape=0.8, threshold=1, type="GP") #density
qevd(c(0.1,0.9), scale=0.5, shape=0.8, threshold=1, type="GP")
revd(10, scale=0.5, shape=0.8, threshold=1, type="GP")

library(ismev)
data(rain)
plot(rain, type='l')

u <- seq(0, max(rain), 0.1)
mrlplot(rain, nint=1000)

threshrange.plot(rain, r=c(0,40), nint=15)


# 2nd quarter ----
#return level
#Obs above the threshold
threshold <- 30 #Set the threshold
rainext<-rain[rain>30]
rainext

# Optimize the parameters using nlminb
start_params <- c(sigma = sd(rainext), xi = 0.1) 
result <- nlminb(start = start_params, objective = loglik_gpd, data = rainext, threshold = threshold)

# The estimated parameters
estimated_params <- result$par
estimated_params

#load extRemes package
library(extRemes)
fit<-fevd(rainext, threshold = 30, method = "MLE", type="GP")
fit
fit$results$par
fit$results$hessian
varcov<-solve(fit$results$hessian)
sqrt(diag(varcov))

#95%C1
sd<-sqrt(diag(varcov))
par1<-fit$results$par[1]
lowsc<-par1-1.96*sd[1]
upsc<-par1+1.96*sd[1]
cat(lowsc, upsc)

par2<-fit$results$par[2]
lowsh<-par2-1.96*sd[2]
upsh<-par2+1.96*sd[2]
cat(lowsh, upsh)

#using extRemes package
fit<-fevd(rainext, threshold = 30, method = "MLE", type="GP")
ci(fit,type="parameter")

#L-moments
library(extRemes)
fit2<-fevd(rainext, threshold = 30, method = "Lmoments", type="GP")
fit2
ci(fit2,type="parameter")

x<-sort(rainext)
x
p<-pevd(x, scale=7.440252, shape=0.184498, threshold=30, type="GP")
p
emp<-(1:152)/153
plot(emp, p)
abline(0, 1, col = "red", lwd=2)

fit<-fevd(rainext, threshold = 30, method = "MLE", type="GP")
plot(fit, "prob")
fit2<-fevd(rainext, threshold = 30, method = "Lmoments", type="GP")
plot(fit2, "prob")

emp<-(1:152)/153
q<-qevd(emp, scale=7.440252, shape=0.184498, threshold=30, type="GP")
q
x<-sort(rainext)
x
plot(x, q)
abline(0, 1, col = "red", lwd=2)


library(extRemes)
lam<-length(rainext)/length(rain)
#MLE
rlevd(50, scale = 7.440252, shape = 0.184498, threshold = 30, npy = 365.25, rate = lam, type="GP")
#L-moment
rlevd(50, scale = 7.299019, shape = 0.1965159, threshold = 30, npy = 365.25, rate = lam, type="GP")

#std error of return level
varlam<-(lam*(1-lam))/length(rain)
fit<-fevd(rainext, threshold = 30, method = "MLE", type="GP")
varcov<-solve(fit$results$hessian)
newvarcov<-matrix(0, ncol=3, nrow=3)
newvarcov[1,1]<-varlam
newvarcov[2:3,2:3]<-varcov
newvarcov
sqrt(diag(newvarcov))

# Define your variables
r <- 50      # The return period
ny <- 365.25  # The number of years
lam <- length(rainext)/length(rain) # Threshold exceedance rate
sigma <- as.numeric(fit$results$par[1])  # Scale parameter
xi <- as.numeric(fit$results$par[2])   # Shape parameter

# Calculate each element of the gradient
del_1 <- (sigma*(ny*r*lam)^xi)/lam
del_2 <- ((ny*r*lam)^xi-1)/xi
del_3 <-(sigma*((ny*r*lam)^xi*log(ny*r*lam)*xi-(ny*r*lam)^xi+1))/xi^2

# Create a matrix with these values
del <- matrix(c(del_1, del_2, del_3), ncol=1)
del

#Variance of return level
varrl<-t(del)%*%newvarcov%*%del
varrl
sqrt(varrl)


rl<-rlevd(50, scale = 7.440252, shape = 0.184498, threshold = 30, npy = 365.25, rate = lam, type="GP")
sd<-sqrt(varrl)
low<-rl-1.96*sd 
upp<-rl+1.96*sd 
cat(low, upp)


# Chapter 4 ----

data<-read.csv(file.choose(), header=T)
data$API
api <- as.numeric(data$API)
api <- na.omit(api)
pacf(api)
acf(api)

n <- length(api)
plot(api[1:(n-1)],api[2:n],xlab="API(t)",ylab="API(t+1)",main="API")
cor(api[1:(n-1)], api[2:n])

library(extRemes)
mrlplot(api, nint=1000)
plot(api, xlab = "No.of observation", ylab ="API", cex= 1.25, cex.lab = 1.25, col = "darkblue", bg ="lightblue", pch = 21, xaxt = "n", type="l")
abline(h=100, col="red", lwd=1.5) #add threshold line
index<-which(api>100) #position
points(index,api[index], col = "darkred", cex =1) #data above threshold

brad <- scan("bradfield.txt")
brad2 <- brad[!is.na(brad)]
pacf(brad2)
plot(brad[1:87671], brad[2:87672])

cluster10 <- function(dataset,threshold) { 
  x=list() 
  z=list() 
  j=1 
  { 
    for(i in (11):length(dataset)) { 
      if(dataset[i-10]>threshold 
         & dataset[i-9]<=threshold & dataset[i-8]<=threshold 
         & dataset[i-7]<=threshold & dataset[i-6]<=threshold 
         & dataset[i-5]<=threshold & dataset[i-4]<=threshold 
         & dataset[i-3]<=threshold & dataset[i-2]<=threshold 
         & dataset[i-1]<=threshold 
         & dataset[i]<=threshold) { 
            x=max(dataset[j:i]) 
            ifelse(i !=length(dataset), j<-i+1, NA) 
            z=c(z,x)}}} 
  return(z)
}

all(as.numeric(cluster10(brad2, 50))==newcluster(brad2, 50, 10))

cluster4 <- function(dataset,threshold,kappa) { 
  x=list() 
  z=list()
  kappa.int=kappa+1
  j=1 
  { 
    for(i in kappa.int:length(dataset)) {
      Seq <- k:(kappa.int+(k-1))
      if(dataset[i-kappa]>threshold) {
        
      }
         # & dataset[i-3]<=threshold & dataset[i-2]<=threshold 
         # & dataset[i-1]<=threshold 
         # & dataset[i]<=threshold) { 
        x=max(dataset[j:i]) 
        ifelse(i !=length(dataset), j<-i+1, NA) 
        z=c(z,x)}}} 
  return(z)
}

brad2[1]>50&(brad2[2]<=50)

as.numeric(cluster10(brad2, 50))

library("ismev")
cluster.peaks <- as.numeric(cluster10(brad2,50))
A <- gpd.fit(cluster.peaks, 50)
A$rate <- length(cluster.peaks)/length(brad2)

plot(decluster(brad2,50))
plot(decluster(brad2, 50, r=10))

gpd.prof(A, 50, npy=365.25*24, xlow = 92, xup = 130)
abline(v=101.533, lty=2)
data(wooster)

cluster10(-wooster, 10)


cluster2 = function(dataset, threshold) {
  x = list()
  z = list()
  j = 1
  {
    for(i in (3):length(dataset)) {
      if(dataset[i-2]>threshold & dataset[i-1]<=threshold & dataset[i] <= threshold) {
        x = max(dataset[j:i])
        ifelse(i!=length(dataset), j <- i+1, NA)
        z = c(z,x)}
    }
  }
  return(z)
}


cluster2(-wooster, 10)
cluster2(-wooster, 5)
cluster2(-wooster, 20) # no values

cluster4(-wooster, 10)
cluster4(-wooster, 20)

gpd.fit(as.numeric(cluster2(-wooster, 5)), 5)
gpd.fit(as.numeric(cluster2(-wooster, 20)), 20)
gpd.fit(as.numeric(cluster2(-wooster, 10)), 10)
gpd.fit(as.numeric(cluster4(-wooster, 10)), 10)
gpd.fit(as.numeric(cluster4(-wooster, 20)), 20)

wooster[wooster<0]
plot(wooster)

B <- gpd.fit(as.numeric(cluster4(-wooster, 10)), 10)
length(as.numeric(cluster4(-wooster, 10)))/length(wooster)


# Run GEV given a kappa
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

newcluster(brad2, 50, 20)

cluster.peaks <- newcluster(brad2, threshold = 50, kappa = 10)
library("ismev")
A <- gpd.fit(cluster.peaks, 50)
A$rate <- length(clusterpeaks)/length(brad2)

return.level <- function(threshold, scale, shape, rate) {
  threshold+(scale/shape)*((threshold*(365.25*24)*rate)^(shape)-1)
}
return.level(50, A$mle[1], A$mle[2], A$rate)
gpd.prof(A, 50, npy=365.25*24, xlow=92, xup=140)
abline(v=return.level(50, A$mle[1], A$mle[2], A$rate), lty=2)


data("wooster")
newcluster(-wooster, 10, 2)

data("Tphap")
plot(Tphap)
head(Tphap)

extremalindex(Tphap$MaxT, threshold=105)
pacf(Tphap$MaxT) # significant autocorrelation at lag 1
y <- decluster(Tphap$MaxT, threshold=105, r=2)
plot(y)
gpd.fit(y, threshold=105)

mrlplot(Tphap$MinT)
extremalindex(Tphap$MinT, threshold = 73)
decluster(Tphap$MinT, threshold=77, r=1)
gpd.fit(decluster(Tphap$MinT, threshold=77, r=1), threshold = 77)

attach(venice)
head(venice)
colnames(venice)[1] <- c("Year")
colnames(venice)[2] <- c("Annual Maxima")
anmax <- venice$`Annual Maxima`
plot(anmax, type="l")
ti <- matrix(ncol=1, nrow=51)
ti[,1] <- seq(1, 51, 1)
Y <- lm(anmax~seq(1, 51, 1))
abline(Y)
venice.gev <- gev.fit(anmax, ydat=ti, mul=1) # nonstationary GEV
venice.gev.st <- gev.fit(anmax)
mu.t <- venice.gev$mle[1]+venice.gev$mle[2]*ti[,1]
D <- 2*(-venice.gev$nllh-(-venice.gev.st$nllh))
# reject H0 if D > crit(chi)
# Nonstationary model is better than the stationary model

ti2 <- cbind(seq(1,51,1), seq(1,51,1)^2)
venice.gev.quad <- gev.fit(anmax, ydat=ti2, mul=c(1,2), show=F)
2*(-venice.gev.quad$nllh-(-venice.gev$nllh))

gev.diag(venice.gev)
gev.prof(venice.gev, m=100, xlow=-10, xup=90, conf=0.95, nint=100)
erlevd(venice.gev)
