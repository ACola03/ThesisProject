library(checkPlotR)
library(dplyr)
library(purrr)
library(reshape2)
library(shellpipes)
library(patchwork)
source("trianglePlot_Poisson.R")

play <- function(numSims){
  unif.dat <- data.frame(lambda = 1, p=runif(numSims))
  checkplotWrapper(unif.dat, numSims, 1, "poisson.test", 0.01)
}

play(1e4)

