library(checkPlotR)
library(dplyr)
library(purrr)
library(reshape2)
library(shellpipes)
library(patchwork)
source("trianglePlot_Poisson.R")

set.seed(19)

## Please fix case_match (see warning)
play <- function(numSims){
  unif.dat <- data.frame(lambda = 0, p=runif(numSims))
  checkplotWrapper(unif.dat, numSims, 1, "I hate supervisors!", 0.01)
}

play(1e6)

