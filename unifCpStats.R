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

# ----- Generate checkPlot Statistic Distributions

empDist <- function(numSims, numReps){
  
  var2s <- c()
  var3s <- c()
  
  for (x in 1:numReps){
    unif.dat <- data.frame(lambda = 1, p=runif(numSims))
    vars <- checkplotStats(dat = unif.dat, binwidth = 0.01, 0)
    var2s <- c(var2s, vars$var2[[1]])
    var3s <- c(var3s, vars$var3[[1]])
  }
  
  return(list("var2s" = var2s, "var3s" = var3s))
}

vars <- empDist(1e4, 1e4)

# ----- Spacing Statistics: Simulated

hist(vars$var2s)
mean(vars$var2s)
sqrt(var(vars$var2s))

# ----- Location Statistics: Simulated

hist(vars$var3s)
mean(vars$var3s)
sqrt(var(vars$var3s))

# ----- Spacing Statistics: Theoretical

(4*(n+1)^4) / ((n+4)*(n+3)*(n+2)^2*(n))

# ----- Location Statistics: Theoretical

i <- 1:n
mean.theory <- (sum(i^2) + sum(i))/((n+1)*(n+2)) - sum(i^2)/((n+1)^2)

n <- 1e4
c <- combn(1:n, 2, FUN = locationVar.aux2, simplify = FALSE)
t1 <- 2*sum(unlist(c))

c2 <- combn(1:n, 1, FUN = locationVar.aux3, simplify = FALSE)
t2 <- sum(unlist(c2))

(t1 + t2) - mean.theory^2

# ----- Scaled Location Statistics: Theoretical

i <- 1:n
vi <- i*(n-i+1)/((n+1)^2*(n+2))
mean.theory.scaled <-
  sum(i * (i + 1) / vi) / ((n + 1) * (n + 2)) -
  sum(i^2 / vi) / ((n + 1)^2)

n <- 1e4
c.scaled <- combn(1:n, 2, FUN = locationVar.aux2.scaled, simplify = FALSE)
t1.scaled <- 2*sum(unlist(c.scaled))

c2.scaled <- combn(1:n, 1, FUN = locationVar.aux3.scaled, simplify = FALSE)
t2.scaled <- sum(unlist(c2.scaled))

(t1.scaled + t2.scaled) - mean.theory.scaled^2

# -----

locationVar.aux2 <- function(combination){
  i <- combination[1]
  j <- combination[2]
  
  beta.m <- list( 
    "4i" = (i*(i+1)*(i+2)*(i+3)) / ((n+1)*(n+2)*(n+3)*(n+4)),
    "3i" = (i*(i+1)*(i+2)) / ((n+1)*(n+2)*(n+3)),
    "2i" = (i*(i+1)) / ((n+1)*(n+2)),
    "1i" = (i) / ((n+1)),
    
    "4j" = (j*(j+1)*(j+2)*(j+3)) / ((n+1)*(n+2)*(n+3)*(n+4)),
    "3j" = (j*(j+1)*(j+2)) / ((n+1)*(n+2)*(n+3)),
    "2j" = (j*(j+1)) / ((n+1)*(n+2)),
    "1j" = (j) / ((n+1))
  )
  
  diri.m <- list(
    "4" = 24/((n+4)*(n+3)*(n+2)*(n+1)),
    "3-1" = 6/((n+4)*(n+3)*(n+2)*(n+1)),
    "2-2" = 4/((n+4)*(n+3)*(n+2)*(n+1)),
    "2-1-1" = 2/((n+4)*(n+3)*(n+2)*(n+1)),
    "1-1-1-1" = 1/((n+4)*(n+3)*(n+2)*(n+1)),
    
    "3" = 6/((n+3)*(n+2)*(n+1)),
    "2-1" = 2/((n+3)*(n+2)*(n+1)),
    "1-1-1" = 1/((n+3)*(n+2)*(n+1)),
    
    "2" = 2/((n+2)*(n+1)),
    "1-1" = 1/((n+2)*(n+1))
  )
  
  ci <- i/(n+1)
  cj <- j/(n+1)
  
  # COMPONENT 1:
  C1c1 <- i*diri.m$`4` + (i*j - i)*diri.m$`2-2`
  C1c2 <- 2*i*(j-1)*diri.m$`3-1`+ (i*j*(j-1) - 2*i*(j-1))*diri.m$`2-1-1`
  C1c3 <- 2*i*(i-1)*diri.m$`3-1`+ (i*j*(i-1) - 2*i*(i-1))*diri.m$`2-1-1`
  C1c4 <- 2*i*(i-1)*diri.m$`2-2` + i*(i-1)*(j-2)*(j-3)*diri.m$`1-1-1-1` +
    (i*(i-1)*j*(j-1) - 2*i*(i-1) - i*(i-1)*(j-2)*(j-3))*diri.m$`2-1-1`
  C1 <- C1c1 + C1c2 + C1c3 + C1c4
  
  # COMPONENT 2:
  C2c1 <- i*diri.m$`3` + (i*j - i)*diri.m$`2-1`
  C2c2 <- 2*i*(i-1)*diri.m$`2-1` + (i*j*(i-1) - 2*i*(i-1))*diri.m$`1-1-1`
  C2 <- -2*cj*(C2c1 + C2c2)
  
  # COMPONENT 3:
  C3c1 <- beta.m$`2i`
  C3 <- cj^2 * C3c1
  
  # COMPONENT 4:
  C4c1 <- i*diri.m$`3` + (i*j - i)*diri.m$`2-1`
  C4c2 <- 2*i*(j-1)*diri.m$`2-1` + (i*j*(j-1) - 2*i*(j-1))*diri.m$`1-1-1`
  C4 <- -2*ci*(C4c1 + C4c2)
  
  # COMPONENT 5:
  C5c1 <- i*diri.m$`2`
  C5c2 <- (i*j - i)*diri.m$`1-1`
  C5 <- 4*ci*cj*(C5c1 + C5c2)
  
  # COMPONENT 6:
  C6 <- -2*ci*cj^2*beta.m$`1i`
  
  # COMPONENT 7:
  C7 <- ci^2*beta.m$`2j`
  
  # COMPONENT 8:
  C8 <- -2*ci^2*cj*beta.m$`1j`
  
  # COMPONENT 9:
  C9 <- ci^2*cj^2
  
  return(C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 + C9)
}

# =========

locationVar.aux3 <- function(combination){
  i <- combination[1]
  
  beta.m <- list( 
    "4" = (i*(i+1)*(i+2)*(i+3)) / ((n+1)*(n+2)*(n+3)*(n+4)),
    "3" = (i*(i+1)*(i+2)) / ((n+1)*(n+2)*(n+3)),
    "2" = (i*(i+1)) / ((n+1)*(n+2)),
    "1" = (i) / ((n+1))
  )
  
  diri.m <- list(
    "4" = 24/((n+4)*(n+3)*(n+2)*(n+1)),
    "3-1" = 6/((n+4)*(n+3)*(n+2)*(n+1)),
    "2-2" = 4/((n+4)*(n+3)*(n+2)*(n+1)),
    "2-1-1" = 2/((n+4)*(n+3)*(n+2)*(n+1)),
    "1-1-1-1" = 1/((n+4)*(n+3)*(n+2)*(n+1)),
    
    "3" = 6/((n+3)*(n+2)*(n+1)),
    "2-1" = 2/((n+3)*(n+2)*(n+1)),
    "1-1-1" = 1/((n+3)*(n+2)*(n+1)),
    
    "2" = 2/((n+2)*(n+1)),
    "1-1" = 1/((n+2)*(n+1))
  )
  
  ci <- i/(n+1)
  
  C <- beta.m$`4`- 4*ci*beta.m$`3` + 6*ci^2*beta.m$`2` - 4*ci^3*beta.m$`1` + ci^4
  
  return(C)
}

# =====

locationVar.aux2.scaled <- function(combination){
  i <- combination[1]
  j <- combination[2]
  
  beta.m <- list( 
    "4i" = (i*(i+1)*(i+2)*(i+3)) / ((n+1)*(n+2)*(n+3)*(n+4)),
    "3i" = (i*(i+1)*(i+2)) / ((n+1)*(n+2)*(n+3)),
    "2i" = (i*(i+1)) / ((n+1)*(n+2)),
    "1i" = (i) / ((n+1)),
    
    "4j" = (j*(j+1)*(j+2)*(j+3)) / ((n+1)*(n+2)*(n+3)*(n+4)),
    "3j" = (j*(j+1)*(j+2)) / ((n+1)*(n+2)*(n+3)),
    "2j" = (j*(j+1)) / ((n+1)*(n+2)),
    "1j" = (j) / ((n+1))
  )
  
  diri.m <- list(
    "4" = 24/((n+4)*(n+3)*(n+2)*(n+1)),
    "3-1" = 6/((n+4)*(n+3)*(n+2)*(n+1)),
    "2-2" = 4/((n+4)*(n+3)*(n+2)*(n+1)),
    "2-1-1" = 2/((n+4)*(n+3)*(n+2)*(n+1)),
    "1-1-1-1" = 1/((n+4)*(n+3)*(n+2)*(n+1)),
    
    "3" = 6/((n+3)*(n+2)*(n+1)),
    "2-1" = 2/((n+3)*(n+2)*(n+1)),
    "1-1-1" = 1/((n+3)*(n+2)*(n+1)),
    
    "2" = 2/((n+2)*(n+1)),
    "1-1" = 1/((n+2)*(n+1))
  )
  
  ci <- i/(n+1)
  cj <- j/(n+1)
  
  # COMPONENT 1:
  C1c1 <- i*diri.m$`4` + (i*j - i)*diri.m$`2-2`
  C1c2 <- 2*i*(j-1)*diri.m$`3-1`+ (i*j*(j-1) - 2*i*(j-1))*diri.m$`2-1-1`
  C1c3 <- 2*i*(i-1)*diri.m$`3-1`+ (i*j*(i-1) - 2*i*(i-1))*diri.m$`2-1-1`
  C1c4 <- 2*i*(i-1)*diri.m$`2-2` + i*(i-1)*(j-2)*(j-3)*diri.m$`1-1-1-1` +
    (i*(i-1)*j*(j-1) - 2*i*(i-1) - i*(i-1)*(j-2)*(j-3))*diri.m$`2-1-1`
  C1 <- C1c1 + C1c2 + C1c3 + C1c4
  
  # COMPONENT 2:
  C2c1 <- i*diri.m$`3` + (i*j - i)*diri.m$`2-1`
  C2c2 <- 2*i*(i-1)*diri.m$`2-1` + (i*j*(i-1) - 2*i*(i-1))*diri.m$`1-1-1`
  C2 <- -2*cj*(C2c1 + C2c2)
  
  # COMPONENT 3:
  C3c1 <- beta.m$`2i`
  C3 <- cj^2 * C3c1
  
  # COMPONENT 4:
  C4c1 <- i*diri.m$`3` + (i*j - i)*diri.m$`2-1`
  C4c2 <- 2*i*(j-1)*diri.m$`2-1` + (i*j*(j-1) - 2*i*(j-1))*diri.m$`1-1-1`
  C4 <- -2*ci*(C4c1 + C4c2)
  
  # COMPONENT 5:
  C5c1 <- i*diri.m$`2`
  C5c2 <- (i*j - i)*diri.m$`1-1`
  C5 <- 4*ci*cj*(C5c1 + C5c2)
  
  # COMPONENT 6:
  C6 <- -2*ci*cj^2*beta.m$`1i`
  
  # COMPONENT 7:
  C7 <- ci^2*beta.m$`2j`
  
  # COMPONENT 8:
  C8 <- -2*ci^2*cj*beta.m$`1j`
  
  # COMPONENT 9:
  C9 <- ci^2*cj^2
  
  sum <- C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 + C9
  scaled.sum <- sum/(beta.m$`1i` * beta.m$`1j`)
  
  return(scaled.sum)
}

locationVar.aux3.scaled <- function(combination){
  i <- combination[1]
  
  beta.m <- list( 
    "4" = (i*(i+1)*(i+2)*(i+3)) / ((n+1)*(n+2)*(n+3)*(n+4)),
    "3" = (i*(i+1)*(i+2)) / ((n+1)*(n+2)*(n+3)),
    "2" = (i*(i+1)) / ((n+1)*(n+2)),
    "1" = (i) / ((n+1))
  )
  
  diri.m <- list(
    "4" = 24/((n+4)*(n+3)*(n+2)*(n+1)),
    "3-1" = 6/((n+4)*(n+3)*(n+2)*(n+1)),
    "2-2" = 4/((n+4)*(n+3)*(n+2)*(n+1)),
    "2-1-1" = 2/((n+4)*(n+3)*(n+2)*(n+1)),
    "1-1-1-1" = 1/((n+4)*(n+3)*(n+2)*(n+1)),
    
    "3" = 6/((n+3)*(n+2)*(n+1)),
    "2-1" = 2/((n+3)*(n+2)*(n+1)),
    "1-1-1" = 1/((n+3)*(n+2)*(n+1)),
    
    "2" = 2/((n+2)*(n+1)),
    "1-1" = 1/((n+2)*(n+1))
  )
  
  ci <- i/(n+1)
  
  sum <- beta.m$`4`- 4*ci*beta.m$`3` + 6*ci^2*beta.m$`2` - 4*ci^3*beta.m$`1` + ci^4
  scaled.sum <- sum / (beta.m$`2`)
  
  return(scaled.sum)
}


