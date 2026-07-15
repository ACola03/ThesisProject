library(checkPlotR)
library(dplyr)
library(purrr)
library(reshape2)
library(shellpipes)
library(patchwork)
source("trianglePlot_Poisson.R")

# ----- random 

# dat <- generatePoisson(lambdas = c(5,10,100), dat = NULL, numSims = 1e4, numReps = 1, testv = "poisson.test")
# dat <- fuzzPoisson(dat, testv = "poisson.test", fuzz.type = "unsupervised", use.fuzz = TRUE)
# checkplotStats(dat, 0.01, 3)
# checkplotWrapper(dat, 1e4, 1, "poisson.test", 0.05, 3)
# the pdf plots are for smaller sample size (?)
# maybe convert to MEAN error now

# ----- viewing beta changes

epsilon <- seq(-0.1,0.1,0.025)
numSims <- 1e4

for (e in epsilon){
  beta.data <- data.frame(lambda = 0, p=rbeta(numSims, 1+e, 1))
  title <- paste0(sprintf("Beta(%.0f, %.3f)", 1, 1+e))
  checkplotWrapper(beta.data, numSims = numSims, 1, title, 0.025)
}

# ----- pinning 0 and 1 improves the statistic

z <- rbeta(numSims, 1, 1.25)
pv <- sort(c(0,z,1))
d <- (numSims+1)*diff(pv)
var(d)
mean(d)

z <- rbeta(numSims, 1, 1.25)
pv <- sort(z)
d <- (numSims+1)*diff(pv)
var(d)
mean(d)

# ----- the statistic distributions under the small perturbations

empDist.beta <- function(numSims, numReps, epsilon){
  
  stat.space <- c()
  stat.location <- c()
  
  for (x in 1:numReps){
    beta.dat <- data.frame(lambda = 0, p=rbeta(numSims, 1+epsilon, 1))
    stats <- checkplotStats(dat = beta.dat, binwidth = 0.01, 0)
    stat.space <- c(stat.space, stats$var2[[1]])
    stat.location <- c(stat.location, stats$var3[[1]])
  }
  
  return(list("stat.space" = stat.space, "stat.location" = stat.location))
}


for (e in seq(-0.1,0,0.025)){
  betaStats <- empDist.beta(1e4, 1e4, e)
  
  space <- betaStats$stat.space
  location <- betaStats$stat.location
  
  min <- min(space)
  max <- max(space)
  title <- paste0(sprintf("Beta(%.0f, %.3f): Space Range: (%.3f, %.3f) ", 1, 1+e, min, max))
  hist(space, main = title)
  
  min <- min(location)
  max <- max(location)
  title <- paste0(sprintf("Beta(%.0f, %.3f): Location Range: (%.3f, %.3f) ", 1, 1+e, min, max))
  hist(location, main = title)
}











