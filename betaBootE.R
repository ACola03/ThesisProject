# ===== PARAMETRIC BOOTSTRAP EXECUTION

library(bbmle)
library(MASS)
source("betaBootF.R")

# save and load to bypass runtime (for n = 2000)
# saveRDS(kde.ab, file = "alphaBeta_kde.rds")
# kde_load <- readRDS("alphaBeta_kde.rds")

# ===== Confidence Regions Under Uniformity

set.seed(99)
n <- 2000
boots <- 1e4
conf.level <- 0.95
paras <- paraBoot(n, boots)

# ----- (alpha, beta)

alpha <- paras$shape1
beta <- paras$shape2
kde.ab <- kde2d(alpha, beta, n = 100)

dx <- diff(kde.ab$x[1:2])
dy <- diff(kde.ab$y[1:2])
sz <- sort(kde.ab$z)
c1 <- cumsum(sz) * dx * dy

level_cutoff <- approx(c1, sz, xout = 1 - conf.level)$y

plot(alpha, beta, pch = 19, col = "gray50", 
     main = paste0(sprintf("95%% 2D Confidence Region via kde2d \n(n = %d, N = %d)", n, boots)))
contour(kde.ab, levels = level_cutoff, labels = paste0(conf.level * 100, "%"), 
        col = "red", lwd = 2, add = TRUE)

# ----- (mean, shape)

mean <- paras$mean
shape <- paras$shape
kde.ms <- kde2d(mean, shape, n = 100)

dx <- diff(kde.ms$x[1:2])
dy <- diff(kde.ms$y[1:2])
sz <- sort(kde.ms$z)
c1 <- cumsum(sz) * dx * dy

level_cutoff <- approx(c1, sz, xout = 1 - conf.level)$y

plot(mean, shape, pch = 19, col = "gray50", 
     main = paste0(sprintf("95%% 2D Confidence Region via kde2d \n(n = %d, N = %d)", n, boots)))
contour(kde.ms, levels = level_cutoff, labels = paste0(conf.level * 100, "%"), 
        col = "red", lwd = 2, add = TRUE)


# -----

dx <- diff(kdeMS_load$x[1:2])
dy <- diff(kdeMS_load$y[1:2])
sz <- sort(kdeMS_load$z)
c1 <- cumsum(sz) * dx * dy

level_cutoff <- approx(c1, sz, xout = 1 - conf.level)$y










