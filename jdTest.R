#library(shellpipes) ## to load the package
#startGraphics()

library(checkPlotR)
library(dplyr)
library(purrr)
library(reshape2)
library(shellpipes)

source("trianglePlot_Poisson.R")

trianglePlot(lambda = c(1,3,5), plot = "both", testv = "poisson.test", fuzz.x = FALSE, add.checkplot = TRUE)

trianglePlot(lambda = c(1,3,5), plot = "both", testv = "wald.intercept", fuzz.x = FALSE, add.checkplot = TRUE)
trianglePlot(lambda = c(1,3,5), plot = "both",  numReps = 10, testv = "wald.intercept", fuzz.x = FALSE, add.checkplot = TRUE)
trianglePlot(lambda = c(1,3,5), plot = "both",  numReps = 100, testv = "wald.intercept", fuzz.x = FALSE, add.checkplot = TRUE)

dat <- rpois(1e4, lambda = 3)
trianglePlot(lambda = c(1,3,5), dat = dat, plot = "both", testv = "poisson.test", fuzz.x = FALSE, add.checkplot = TRUE)







