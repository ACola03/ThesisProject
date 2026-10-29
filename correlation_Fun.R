cor.test2 <- function(x, y, r.null, p.out){
  
  r <- cor(x, y)
  n <- length(x)
  
  if (r.null == 0){
    stat <- r*sqrt((n-2)/(1-r^2)) # ~ t_{n-2}
    pval.two <- 2*min(pt(stat, n-2, lower.tail = TRUE), pt(stat, n-2, lower.tail = FALSE))
  } else {
    stat <- (atanh(r) - atanh(r.null))/(1/sqrt(n-3)) # approx. ~ N(0,1)
    pval.two <- 2*min(pnorm(stat, lower.tail = TRUE), pnorm(stat, lower.tail = FALSE))
  }
  
  if (p.out == "two") return(pval.two)
  
  if (stat < 0){
    pval <- pval.two/2
  } else{
    pval <- 1-pval.two/2
  }
  
  return(pval)
}


simCorr <- function(n, rho = 0.6, p.out = c("two", "signed")) {
  
  p.out <- match.arg(p.out)
  
  if (n < 4)
    stop("n must be at least 4.")
  
  if (!is.finite(rho) || rho == 0 || abs(rho) >= 1)
    stop("rho must be nonzero and strictly between -1 and 1.")
  
  # 1. Independent, jointly normal variables
  dat.null <- mvrnorm(
    n, mu = c(0, 0), Sigma = diag(2)
  )
  
  # 2. Independent variables, one squared
  dat.null.indep.sqr <- dat.null
  dat.null.indep.sqr[, 2] <- dat.null.indep.sqr[, 2]^2
  
  # 3. Nonlinear dependence: X versus X^2
  x <- rnorm(n)
  dat.null.nonlin <- cbind(x, x^2)
  
  # 4. Correlated, jointly normal variables
  Sigma.rho <- matrix(
    c(1, rho, rho, 1),
    nrow = 2,
    byrow = TRUE
  )
  
  dat.corr <- mvrnorm(
    n, mu = c(0, 0), Sigma = Sigma.rho
  )
  
  # 5. Correlated normal variables, then square Y
  dat.corr.sqr <- dat.corr
  dat.corr.sqr[, 2] <- dat.corr.sqr[, 2]^2
  
  # Calculate p-values under the appropriate null hypotheses
  pvals <- c(
    null = cor.test2(
      dat.null[, 1], dat.null[, 2],
      r.null = 0, p.out = p.out
    ),
    
    null.indep.sqr = cor.test2(
      dat.null.indep.sqr[, 1], dat.null.indep.sqr[, 2],
      r.null = 0, p.out = p.out
    ),
    
    null.nonlin = cor.test2(
      dat.null.nonlin[, 1], dat.null.nonlin[, 2],
      r.null = 0, p.out = p.out
    ),
    
    corr = cor.test2(
      dat.corr[, 1], dat.corr[, 2],
      r.null = rho, p.out = p.out
    ),
    
    corr.sqr = cor.test2(
      dat.corr.sqr[, 1], dat.corr.sqr[, 2],
      r.null = 0, p.out = p.out
    )
  )
  
  return(pvals)
}