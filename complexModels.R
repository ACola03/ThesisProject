library(checkPlotR)
library(dplyr)
library(purrr)
library(reshape2)
library(shellpipes)
library(patchwork)
source("trianglePlot_Poisson.R")

generatePoissonGLM <- function(n,
                               beta = c(1, 0.5),
                               x = NULL,
                               x.dist = function(n) runif(n, -1, 1))
{
  # Generate covariate if not supplied
  if (is.null(x))
    x <- x.dist(n)
  
  # Linear predictor
  eta <- beta[1] + beta[2] * x
  
  # Mean
  mu <- exp(eta)
  
  # Response
  y <- rpois(n, mu)
  
  data.frame(
    y = y,
    x = x,
    mu = mu,
    eta = eta
  )
}

# -----

set.seed(1)

dat <- generatePoissonGLM(
  n = 30,
  beta = c(log(3), 0.8)
)

head(dat)

# -----

null <- glm(y ~ 1,
            family = poisson,
            data = dat)

alt <- glm(y ~ x,
           family = poisson,
           data = dat)

an <- anova(null, alt, test = "LRT")
an$`Pr(>Chi)`[2]

# -----

generateP <- function(n, numSims, beta1, beta2){
  
  p.vals <- c()
    
  for (i in 1:numSims){

    dat <- generatePoissonGLM(
      n,
      beta = c(beta1, beta2)
    )
    
    null <- glm(y ~ 1,
                family = poisson,
                data = dat)
    
    alt <- glm(y ~ x,
               family = poisson,
               data = dat)
    
    an <- anova(null, alt, test = "LRT")
    p <- an$`Pr(>Chi)`[2]
  
    p.vals <- c(p.vals, p)
  }
  return(p.vals)
}

p.vals <- generateP(
  n = 5,
  numSims = 1e4,
  beta1 = log(3),
  beta2 = 0
)

hist(p.vals)


length(unique(p.vals))


