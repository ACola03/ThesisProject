library(checkPlotR)
library(dplyr)
library(purrr)
library(reshape2)
library(shellpipes)

source("trianglePlot_Poisson.R")

# Notes:
# i) focus on non-fuzzed cases
# ii) the checkplots now represent the data shown in the triangles

# Example 1: Discrete Exact Tails
trianglePlot(lambda = c(1,3,5), plot = "both", testv = "poisson.test", fuzz.x = FALSE, add.checkplot = TRUE)

# ---

# Example 2: Wald Approximation - across increasing 'sample size' for 1e4 experiments
trianglePlot(lambda = c(1,3,5), plot = "two", testv = "wald.intercept", fuzz.x = FALSE, add.checkplot = TRUE)
trianglePlot(lambda = c(1,3,5), plot = "two",  numReps = 10, testv = "wald.intercept", fuzz.x = FALSE, add.checkplot = TRUE)
trianglePlot(lambda = c(1,3,5), plot = "two",  numReps = 100, testv = "wald.intercept", fuzz.x = FALSE, add.checkplot = TRUE)

# ---

# Example 3: Misspecified Null - truth is lambda 3 
# dat <- rpois(1e4, lambda = 3)
# trianglePlot(lambda = c(1,3,5), dat = dat, plot = "both", testv = "poisson.test", fuzz.x = FALSE, add.checkplot = TRUE)

# ---

# Example 4: Intro to Wald Fuzz
wald.data <- generatePoisson(lambdas = c(1,2,3,10), testv = "wald.intercept")

# Non-Fuzzed 
checkPlot(wald.data, facets = 4) + facet_grid(~lambda)

# Fuzzed: excluding 0 count
wald.fuzzed <- wald.fuzz(wald.data, filter.zero = TRUE)
wald.fuzzed$p <- wald.fuzzed$rp
wald.fuzzed %>%
  select(pois.mean, p, rp, lambda, p.lower, p.upper) %>%
  split(.$lambda) %>%
  purrr::map(checkPlot) 

# Fuzzed: including 0 count
wald.fuzzed <- wald.fuzz(wald.data, filter.zero = FALSE)
wald.fuzzed$p <- wald.fuzzed$rp
checkPlot(wald.fuzzed, facets = 4) + facet_grid(~lambda) 



