# ===== PARAMETRIC BOOTSTRAP EXECUTION

library(bbmle)
library(MASS)
library(sp)
source("betaBootF.R")

# ===== Working With Uniform Contour

# load 
kdeAB_load <- readRDS("alphaBeta_kde.rds")
kdeMS_load <- readRDS("meanShape_kde.rds")

# compute the contour cut-offs
cutoff <- function(kde.dat){
  dx <- diff(kde.dat$x[1:2])
  dy <- diff(kde.dat$y[1:2])
  sz <- sort(kde.dat$z)
  c1 <- cumsum(sz) * dx * dy
  level_cutoff <- approx(c1, sz, xout = 1 - conf.level)$y
  return(level_cutoff)
}

level_cutoff.ab <- cutoff(kdeAB_load) 
level_cutoff.ms <- cutoff(kdeMS_load) 

contour.ab <- contourLines(kdeAB_load$x, kdeAB_load$y, kdeAB_load$z, levels = level_cutoff.ab)
contour.ms <- contourLines(kdeMS_load$x, kdeMS_load$y, kdeMS_load$z, levels = level_cutoff.ms)

contourPolygon.ms <- do.call(rbind, lapply(contour.ms, function(contour) data.frame(x = contour$x, y = contour$y)))
contourPolygon.ab <- do.call(rbind, lapply(contour.ab, function(contour) data.frame(x = contour$x, y = contour$y)))


# ===== Confidence Regions Under Non-Uniformity (Bad Tests)

set.seed(99)
n <- 2000
boots <- 1e4
conf.level <- 0.95

badTests <- function(n, boots, shape1, shape2, type, paraBad){

  if (type == "ms"){
    x <- paraBad$mean
    y <- paraBad$shape
    kde.load <- readRDS("meanShape_kde.rds")
  }
  else{
    x <- paraBad$shape1
    y <- paraBad$shape2
    kde.load <- readRDS("alphaBeta_kde.rds")
  }
  
  # chosen Beta
  kde <- kde2d(x, y, n = 100)
  dx <- diff(kde$x[1:2])
  dy <- diff(kde$y[1:2])
  sz <- sort(kde$z)
  c1 <- cumsum(sz) * dx * dy
  level_cutoff <- approx(c1, sz, xout = 1 - conf.level)$y
  
  # uniform
  level_cutoff.load <- cutoff(kde.load) 
  contour.load <- contourLines(kde.load$x, kde.load$y, kde.load$z, levels = level_cutoff.load)
  contourPolygon.load <- do.call(rbind, lapply(contour.load, function(contour) data.frame(x = contour$x, y = contour$y)))
 
  # plot
  title <- paste0(sprintf("Simulated Parameters of Beta(%.2f, %.2f) \n(n = %d, N = %d)", shape1, shape2, n, boots))
  plot(x, y, pch = 19, col = "gray50", main = title) # pert points
  contour(kde, levels = level_cutoff, labels = paste0("Perturbated " , conf.level * 100, "%"), 
          col = "red", lwd = 2, labcex = 1, add = TRUE) # perturbate contour
  contour(kde.load, levels = level_cutoff.load, labels = paste0("Uniform ", conf.level * 100, "%"),
          col = "green", lwd = 2, labcex = 1, add = TRUE) # unif contour
  
  # overlapping
  overlap <- mean(point.in.polygon(x,y,contourPolygon.load$x, contourPolygon.load$y) != 0)
  sprintf("Overlap within the Uniform 95%% 2D Region: %.4f", overlap)
}

badTests(n, boots, shape1, shape2, "ms", paras.bad1)
badTests(n, boots, shape1, shape2, "ab", paras.bad1)

# ----- 1. Beta(0.9,0.9)

shape1 <- 0.9; shape2 <- 0.9
paras.bad1 <- paraBoot(n, boots, shape1, shape2)
badTests(n, boots, shape1, shape2, "ms", paras.bad1)
badTests(n, boots, shape1, shape2, "ab", paras.bad1)

# ----- 2. Beta(1.1,1.1)

shape1 <- 1.1; shape2 <- 1.1
paras.bad2 <- paraBoot(n, boots, shape1, shape2)
badTests(n, boots, shape1, shape2, "ms", paras.bad2)
badTests(n, boots, shape1, shape2, "ab", paras.bad2)

# ----- 3. Beta(1.05,0.95)

shape1 <- 1.05; shape2 <- 0.95
paras.bad3 <- paraBoot(n, boots, shape1, shape2)
badTests(n, boots, shape1, shape2, "ms", paras.bad3)
badTests(n, boots, shape1, shape2, "ab", paras.bad3)

# ----- 4. Beta(0.85,0.85)

shape1 <- 0.85; shape2 <- 0.85
paras.bad4 <- paraBoot(n, boots, shape1, shape2)
badTests(n, boots, shape1, shape2, "ms", paras.bad4)
badTests(n, boots, shape1, shape2, "ab", paras.bad4)

# ----- 5. Beta(1.15,1.15)

shape1 <- 1.15; shape2 <- 1.15
paras.bad5 <- paraBoot(n, boots, shape1, shape2)
badTests(n, boots, shape1, shape2, "ms", paras.bad5)
badTests(n, boots, shape1, shape2, "ab", paras.bad5)
