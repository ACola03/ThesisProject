library(ggplot2)
library(dplyr)
source("statCurvesF.R")

# ===== Evaluate Across Sample Sizes: 50, 100, 250, 500, 1000, 2000

# Grid of lambda (including 1)
lseq <- sort(c(exp(seq(-3, 0.5, length.out=35)),1))

# ----- n = 25
sym.25 <- sim_ss_over_lambda(lseq, 25, 1e4, "n")
skew.25 <- sim_ss_over_lambda(lseq, 25, 1e4, "skew")

# ----- n = 50
sym.50 <- sim_ss_over_lambda(lseq, 50, 1e4, "n")
skew.50 <- sim_ss_over_lambda(lseq, 50, 1e4, "skew")

# ----- n = 100
sym.100 <- sim_ss_over_lambda(lseq, 100, 1e4, "n")
skew.100 <- sim_ss_over_lambda(lseq, 100, 1e4, "skew")

# ----- n = 250
sym.250 <- sim_ss_over_lambda(lseq, 250, 1e4, "n")
skew.250 <- sim_ss_over_lambda(lseq, 250, 1e4, "skew")

# ----- n = 500
sym.500 <- sim_ss_over_lambda(lseq, 500, 1e4, "n")
skew.500 <- sim_ss_over_lambda(lseq, 500, 1e4, "skew")

# ----- n = 1000
sym.1000 <- sim_ss_over_lambda(lseq, 1000, 1e4, "n")
skew.1000 <- sim_ss_over_lambda(lseq, 1000, 1e4, "skew")
