# Probit, Logit comparison

# working directory
setwd("~/Documents/ProbLogAnalysis")

# ************** QQ Plot for Probit/Logit *******************
# create vector for plotting
v<-seq(from = -3, to = 3, by = 0.02)
# vector for fitting
fv<-seq(from = -0.5, to = 0.5, by = 0.02)

# probit cdf
pcdf<-function(x){
  pnorm(x, mean = 0, sd = 1, lower.tail = TRUE, log.p = FALSE)
  
}

# logit cdf
lcdf<-function(x){
  1 / (1 + exp(-x))
  
}

pv<-pcdf(v)
lv<-lcdf(v)

# fit highly linear region
vfit<-lm(lcdf(fv) ~ pcdf(fv))

# plotting
plot(pv, lv, main="QQ Plot for Probit vs Logit", xlab = "Probit Quantile", ylab = "Logit Quantile")
abline(vfit, col="red")

# ***************** Beta Comparison Chart *******************

# define p to beta for probit
pbeta<-function(x){
  qnorm(x)
  
}

# define p to beta for logit
lbeta<-function(x){
  log(x/(1-x))
  
}

# probabilities from evenly spaced beta
p<-pv
fp<-pcdf(fv)

pp<-pbeta(p)
lp<-lbeta(p)

# fit highly linear region
pfit<-lm(lbeta(fp) ~ pbeta(fp))

# fit all
pafit<-lm(lp ~ pp)

# learning about beta distribution

# plotting
plot(pp, lp, main="Beta-Beta Plot for Probit vs Logit", xlab = "Probit Beta", ylab = "Logit Beta")
abline(pfit, col="red")
abline(pafit, col="green")

# *************** Analysis on AJ Data ****************

# import the effects and fequencies
myData = read.csv("effects.txt", sep=" ", header=FALSE)

# set beta0 to try to get p0 around .005
beta0 = -20.5

# calculate mu and sigma^2 for aj and non-aj in logit
ajmu = sum(2 * myData[2] * myData[4]) + beta0
ajsigma = sqrt(sum(myData[2] * (1 - myData[2]) * myData[4]^2))

najmu = sum(2 * myData[3] * myData[4]) + beta0
najsigma = sqrt(sum(myData[3] * (1 - myData[3]) * myData[4]^2))

# convert mu and sigma to probit model
ajmup = pbeta(lcdf(ajmu))
ajsigmap = pbeta(lcdf(ajsigma))

najmup = pbeta(lcdf(najmu))
najsigmap = pbeta(lcdf(najsigma))

# start generating f(p) for both
probs<-seq(0, 0.05, 0.00005)

logitf<-function(x, sigma, mu){
  1 / (sigma * x * (1 - x)) * dnorm(1 / sigma * (lbeta(x) - mu))
  
}

probitf<-function(x, sigma, mu){
  1 / (sigma * dnorm(pbeta(x))) * dnorm(1 / sigma * (pbeta(x) - mu))
  
}

# area under curve funciton (trapezoid)
auc<-function(x, y){
  sum((x[-1] - x[-length(x)]) * (y[-1] + y[-length(y)]) / 2)
  
}

# calculate f(p)
ajlf = logitf(probs, ajsigma, ajmu)
najlf = logitf(probs, najsigma, najmu)

ajpf = probitf(probs, ajsigmap, ajmup)
najpf = probitf(probs, najsigmap, najmup)

# fix NaN issue (sometimes first element is NaN)
probs = probs[-1]
ajlf = ajlf[-1]
najlf = najlf[-1]
ajpf = ajpf[-1]
najpf = najpf[-1]

# plot f(p) for logit
plot(probs, najlf, type="l", main="f(p) for Logit Model", xlab = "p", ylab = "f(p)")
abline(v = 0.005465162)
lines(probs, ajlf, col="red")
abline(v = 0.008242534, col="red")

# plot f(p) for probit
#plot(probs, najpf, type="l", main="f(p) for Probit Model", xlab = "p", ylab = "f(p)")
#abline(v = 0.007158909)
#lines(probs, ajpf, col="red")
#abline(v = 0.009774161, col="r")