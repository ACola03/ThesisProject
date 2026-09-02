library(ggplot2)
library(checkPlotR)
library(patchwork)
library(MASS)
library(sp)
source("betaBootF.R")
source("betaSlugsF.R")

# default unless otherwise stated
set.seed(99)
n <- 2000; boots <- 1e4
conf.level <- 0.95

# ===== SlugPlots: Uniform

# ----- n = 50
n <- 50
shape1 <- 1; shape2 <- 1
paras.unif50 <- paraBoot(n, boots, shape1, shape2, conf.level)
wrap(paras.unif50, shape1, shape2, n, boots)

# ----- n = 100
n <- 100
shape1 <- 1; shape2 <- 1
paras.unif100 <- paraBoot(n, boots, shape1, shape2, conf.level)
wrap(paras.unif100, shape1, shape2, n, boots)

# ----- n = 250
n <- 250
shape1 <- 1; shape2 <- 1
paras.unif250 <- paraBoot(n, boots, shape1, shape2, conf.level)
wrap(paras.unif250, shape1, shape2, n, boots)

# ----- n = 500
n <- 500
shape1 <- 1; shape2 <- 1
paras.unif500 <- paraBoot(n, boots, shape1, shape2, conf.level)
wrap(paras.unif500, shape1, shape2, n, boots)

# ----- n = 1000
n <- 1000
shape1 <- 1; shape2 <- 1
paras.unif1000 <- paraBoot(n, boots, shape1, shape2, conf.level)
wrap(paras.unif1000, shape1, shape2, n, boots)

# ----- n = 2000
n <- 2000
shape1 <- 1; shape2 <- 1
paras.unif <- paraBoot(n, boots, shape1, shape2, conf.level)
wrap(paras.unif, shape1, shape2, n, boots)

# ===== SlugPlots: Non-Uniform

# ----- 1. Beta(0.9,0.9)

shape1 <- 0.9; shape2 <- 0.9
paras.bad1 <- paraBoot(n, boots, shape1, shape2, conf.level)
wrap(paras.bad1, shape1, shape2, n, boots)

# ----- 2. Beta(1.1,1.1)

shape1 <- 1.1; shape2 <- 1.1
paras.bad2 <- paraBoot(n, boots, shape1, shape2, conf.level)
wrap(paras.bad2, shape1, shape2, n, boots)

# ----- 3. Beta(1.05,0.95)

shape1 <- 1.05; shape2 <- 0.95
paras.bad3 <- paraBoot(n, boots, shape1, shape2, conf.level)
wrap(paras.bad3, shape1, shape2, n, boots)

# ----- 4. Beta(0.85,0.85)

shape1 <- 0.85; shape2 <- 0.85
paras.bad4 <- paraBoot(n, boots, shape1, shape2, conf.level)
wrap(paras.bad4, shape1, shape2, n, boots)

# ----- 5. Beta(1.15,1.15)

shape1 <- 1.15; shape2 <- 1.15
paras.bad5 <- paraBoot(n, boots, shape1, shape2, conf.level)
wrap(paras.bad5, shape1, shape2, n, boots)

# ----- 6. Beta(1.25,1)

shape1 <- 1.1; shape2 <- 1
paras.bad6 <- paraBoot(n, boots, shape1, shape2, conf.level)
wrap(paras.bad6, shape1, shape2, n, boots)

# ----- new

# lower miss range
c(250 - 1.96 * sqrt(10000 * 0.025 * 0.975), 
  250 + 1.96 * sqrt(10000 * 0.025 * 0.975))

# upper miss range
c(9750 - 1.96 * sqrt(10000 * 0.025 * 0.975), 
  9750 + 1.96 * sqrt(10000 * 0.025 * 0.975))

# total misses
c(500 - 1.96 * sqrt(10000 * 0.05 * 0.95), 
  500 + 1.96 * sqrt(10000 * 0.05 * 0.95))

2*0.03 - 0.03^2
0^0

# ----- new 2

lower.m <- quantile(paras.unif$mean, 0.025)
upper.m <- quantile(paras.unif$mean, 0.975)
lower.s <- quantile(paras.unif$shape, 0.025)
upper.s <- quantile(paras.unif$shape, 0.975)

m <- paras.unif$mean
s <- paras.unif$shape

l <- mean((m < lower.m) + (s < lower.s)) - mean((m < lower.m) & (s < lower.s))
u <- mean((m > upper.m) + (s > upper.s)) - mean((m > upper.m) & (s > upper.s))
l+u
mean((m < lower.m))

# ----- new 3

lower.m <- paras.unif$mean.lower
upper.m <- paras.unif$mean.upper
lower.s <- paras.unif$shape.lower
upper.s <- paras.unif$shape.upper

l <- mean((lower.m > 0.5) + (lower.s > 1)) - mean((lower.m > 0.5) & (lower.s > 1))
u <- mean((upper.m < 0.5) + (upper.s < 1)) - mean((upper.m < 0.5) & (upper.s < 1))
l+u
mean(lower.m > 0.5)

