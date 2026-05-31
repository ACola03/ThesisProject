#library(shellpipes) ## to load the package
#startGraphics()

library(checkPlotR)
library(dplyr)
library(purrr)
library(reshape2)
library(shellpipes)

source("trianglePlot_Poisson.R")

## trianglePlot(lambda = c(1,3,5), plot = "both", testv = "poisson.test", fuzz.x = TRUE, add.checkplot = TRUE)
trianglePlot(lambda = c(1,3,5), plot = "both", testv = "poisson.test", fuzz.x = FALSE, add.checkplot = TRUE)
trianglePlot(lambda = c(1,3, 5), plot = "one", testv = "wald.intercept", fuzz.x = FALSE, add.checkplot = TRUE)







