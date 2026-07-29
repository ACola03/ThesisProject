# ----- required
library(tidyr)
library(ggplot2)
library(dplyr)
library(patchwork)
source("pianoPowerF.R")

# ----- u.shape

shapes1 <- c(1.0, 0.90, 0.85, 0.80)
shapes2 <- c(1.0, 0.90, 0.85, 0.80)

zpa.ushape <- zhang.powerAnalysis(1e4, c(50, 100, 250, 500, 1000, 2000), 
                                  shapes1, shapes2, c(0, 0.01, 0.025, 0.05))

powerPlot(zpa.ushape, 2000, 0.05, TRUE)

# ----- bell.shape

shapes1 <- c(1.0, 1.2, 1.15, 1.1)
shapes2 <- c(1.0, 1.2, 1.15, 1.1)

zpa.bell <- zhang.powerAnalysis(1e4, c(50, 100, 250, 500, 1000, 2000), 
                                shapes1, shapes2, c(0, 0.01, 0.025, 0.05))

powerPlot(zpa.bell, 2000, 0.05, TRUE)

# ----- dushoff.shape

shapes1 <- c(1.0, 1.25, 2, 0.75)
shapes2 <- c(1.0, 0.75, 2, 0.75)

zpa.dush <- zhang.powerAnalysis(1e4, c(50, 100, 250, 500, 1000, 2000), 
                                shapes1, shapes2, c(0, 0.01, 0.025, 0.05))

powerPlot(zpa.dush, 2000, 0, TRUE) 
powerPlot(zpa.dush, 2000, 0.01, TRUE) 
powerPlot(zpa.dush, 2000, 0.05, TRUE) 

# ----- statistic distribution from violin plot

shapes1 <- c(1.0, 1.25, 2, 0.75)
shapes2 <- c(1.0, 0.75, 2, 0.75)
N <- 1e4; n <- 2000

zks.dist <- zhang.dist(N, n, shapes1, shapes2, "KS", TRUE)
zad.dist <- zhang.dist(N, n, shapes1, shapes2, "AD", TRUE)
zcvm.dist <- zhang.dist(N, n, shapes1, shapes2, "CVM", TRUE)

p3 <- statDistPlot(zks.dist, "ZKS", N, n)
p1 <- statDistPlot(zad.dist, "ZAD", N, n)
p2 <- statDistPlot(zcvm.dist, "ZCVM", N, n)

p1 / p2 / p3

