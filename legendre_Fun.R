# =======
# Evaluating t-Test validity when assumptions aren't satisfied
# using the Neyman Smooth Test (Legendre Polynomials)
# =======

# ===== Perform Neyman Smoothness Test

tTest.legendre.chunk <- function(
    n.neyman, n.legendre, n.t,
    sim.dist = c("norm", "t", "lognorm", "gamma"),
    dist.args = list(),
    chunk.neyman = 10, # compute 10 Neyman stats at a time
    t.p.out = c("two", "left"), ney.p.out = c("two", "left")
){
  
  # Select distribution and default parameters
  sim.dist <- match.arg(sim.dist); t.p.out <- match.arg(t.p.out); ney.p.out <- match.arg(ney.p.out)
  
  default.args <- switch(
    sim.dist,
    norm    = list(sd = 1),
    t       = list(df = 3),
    lognorm = list(meanlog = 0, sdlog = 1),
    gamma   = list(shape = 2, rate = 1)
  )
  
  dist.args <- utils::modifyList(default.args, dist.args)
  
  # Select random generation function
  rfun <- switch(
    sim.dist,
    norm    = rnorm,
    t       = rt,
    lognorm = rlnorm,
    gamma   = rgamma
  )
  
  # Determine theoretical population mean
  pop.mean <- switch(
    sim.dist,
    norm    = 0,
    t       = 0,
    lognorm = exp(dist.args$meanlog + dist.args$sdlog^2 / 2),
    gamma   = dist.args$shape / dist.args$rate
  )
  
  # Check basic requirements
  if (n.t < 2 || n.legendre < 1 || n.neyman < 1 || chunk.neyman < 2)
    stop("Sample sizes and chunk.neyman must be positive, with n.t >= 2 and chunk.neyman >= 2.")
  
  if (n.neyman %% chunk.neyman == 1)
    stop("`colMeans()` requies multiple dimensions; check division of n.neyman by chunk.neyman ")
  
  if (sim.dist == "t" && dist.args$df <= 1)
    stop("Student-t df must exceed 1 for a finite population mean.")
  
  # Initialize output
  neyman.out <- matrix(
    NA_real_,
    nrow = n.neyman,
    ncol = 7,
    dimnames = list(
      NULL,
      c(paste0("theta", 1:6), "joint")
    )
  )
  
  starts <- seq(1, n.neyman, chunk.neyman)
  for (start in starts){
    
    end <- min(start + chunk.neyman - 1, n.neyman)
    k <- end - start + 1
    tot.pvals <- k * n.legendre
    
    # Generate and center observations
    t.data <- matrix(
      do.call(rfun, c(list(n = n.t * tot.pvals), dist.args)) - pop.mean,
      nrow = n.t,
      ncol = tot.pvals
    )
    
    # Calculate one-sample t-statistics
    xbar <- colMeans(t.data)
    
    s2 <- (colSums(t.data^2) - n.t * xbar^2) / (n.t - 1)
    
    test.statistic.temp <- xbar / sqrt(s2 / n.t)
    
    # Calculate left-tailed (or two-tailed) p-values
    pvals <- switch(
      t.p.out,
      two  = 2 * pt(-abs(test.statistic.temp), df = n.t - 1),
      left = pt(test.statistic.temp, df = n.t - 1, lower.tail = TRUE)
    )
    
    pval.matrix <- matrix(
      pvals,
      nrow = n.legendre,
      ncol = k
    )
    
    # Calculate PianoPlot statistics
    theta.mat <- legendre.coef(pval.matrix)
    neyman.out[start:end,] <- legendre.pval(n.legendre, theta.mat, ney.p.out)
  }
  
  return(neyman.out)
}

# -----

tTest.legendre.chunk.nlist <- function(
    n.neyman, n.legendre, n.t.list,
    sim.dist = c("norm", "t", "lognorm", "gamma"),
    dist.args = list(),
    chunk.neyman = 10, # compute 10 Neyman stats at a time
    t.p.out = c("two", "left"), ney.p.out = c("two", "left")
){
  
  # Select distribution and default parameters
  sim.dist <- match.arg(sim.dist); t.p.out <- match.arg(t.p.out); ney.p.out <- match.arg(ney.p.out)
  
  default.args <- switch(
    sim.dist,
    norm    = list(sd = 1),
    t       = list(df = 3),
    lognorm = list(meanlog = 0, sdlog = 1),
    gamma   = list(shape = 2, rate = 1)
  )
  
  dist.args <- utils::modifyList(default.args, dist.args)
  
  # Select random generation function
  rfun <- switch(
    sim.dist,
    norm    = rnorm,
    t       = rt,
    lognorm = rlnorm,
    gamma   = rgamma
  )
  
  # Determine theoretical population mean
  pop.mean <- switch(
    sim.dist,
    norm    = 0,
    t       = 0,
    lognorm = exp(dist.args$meanlog + dist.args$sdlog^2 / 2),
    gamma   = dist.args$shape / dist.args$rate
  )
  
  # Check basic requirements
  if (any(n.t.list < 2) || n.legendre < 1 || n.neyman < 1 || chunk.neyman < 2)
    stop("Sample sizes and chunk.neyman must be positive, with n.t >= 2 and chunk.neyman >= 2.")
  
  if (n.neyman %% chunk.neyman == 1)
    stop("`colMeans()` requies multiple dimensions; check division of n.neyman by chunk.neyman ")
  
  if (sim.dist == "t" && dist.args$df <= 1)
    stop("Student-t df must exceed 1 for a finite population mean.")
  
  neyman.out.all <- map(n.t.list, .f = function(n.t){
    
    # Initialize output
    neyman.out <- matrix(
      NA_real_,
      nrow = n.neyman,
      ncol = 7,
      dimnames = list(
        NULL,
        c(paste0("theta", 1:6), "joint")
      )
    )
    
    starts <- seq(1, n.neyman, chunk.neyman)
    for (start in starts){
      
      end <- min(start + chunk.neyman - 1, n.neyman)
      k <- end - start + 1
      tot.pvals <- k * n.legendre
      
      # Generate and center observations
      t.data <- matrix(
        do.call(rfun, c(list(n = n.t * tot.pvals), dist.args)) - pop.mean,
        nrow = n.t,
        ncol = tot.pvals
      )
      
      # Calculate one-sample t-statistics
      xbar <- colMeans(t.data)
      s2 <- (colSums(t.data^2) - n.t * xbar^2) / (n.t - 1)
      test.statistic.temp <- xbar / sqrt(s2 / n.t)
      
      # Calculate left-tailed (or two-tailed) p-values
      pvals <- switch(
        t.p.out,
        two  = 2 * pt(-abs(test.statistic.temp), df = n.t - 1),
        left = pt(test.statistic.temp, df = n.t - 1, lower.tail = TRUE)
      )
      
      pval.matrix <- matrix(
        pvals,
        nrow = n.legendre,
        ncol = k
      )
      
      # Calculate PianoPlot statistics
      theta.mat <- legendre.coef(pval.matrix)
      neyman.out[start:end,] <- legendre.pval(n.legendre, theta.mat, ney.p.out)
    }
    return(neyman.out)
  }
  )
  return(neyman.out.all)
}

legendre.coef <- function(p.vals){
  
  theta.mat <- cbind(
    theta1 = colMeans(sqrt(3)*(2*p.vals - 1)),
    theta2 = colMeans(sqrt(5)*(6*p.vals^2 - 6*p.vals + 1)),
    theta3 = colMeans(sqrt(7)*(20*p.vals^3 - 30*p.vals^2 + 12*p.vals - 1)),
    theta4 = colMeans(3*(70*p.vals^4 - 140*p.vals^3 + 90*p.vals^2 - 20*p.vals + 1)),
    theta5 = colMeans(sqrt(11)*(252*p.vals^5 - 630*p.vals^4 +
                                  560*p.vals^3 - 210*p.vals^2 + 30*p.vals - 1)),
    theta6 = colMeans(sqrt(13)*(924*p.vals^6 - 2772*p.vals^5 +
                                  3150*p.vals^4 - 1680*p.vals^3 +
                                  420*p.vals^2 - 42*p.vals + 1))
  )
  return(theta.mat)
}

legendre.pval <- function(n.legendre, theta.mat, ney.p.out){
  
  z.mat <- sqrt(n.legendre) * theta.mat
  
  pval.marginal <- if (ney.p.out == "two") {
    2 * pnorm(-abs(z.mat))
    
  } else {
    pval.two <- 2 * pnorm(-abs(z.mat)) 
    ifelse(theta.mat >= 0, (pval.two / 2), 1 - (pval.two / 2))
  }
  
  statistic.joint <- n.legendre * rowSums(theta.mat^2)
  pval.joint <- pchisq(statistic.joint, df = ncol(theta.mat),lower.tail = FALSE)
  
  return(cbind(
    pval.marginal,
    joint = pval.joint
    ))
}

# ===== Get p-values of the t-test itself

tTest.pvalRaw <- function(
    n.legendre, n.t,
    sim.dist = c("norm", "t", "lognorm", "gamma"),
    dist.args = list(), t.p.out = c("two", "left")
){
  
  # Select distribution and default parameters
  sim.dist <- match.arg(sim.dist); t.p.out <- match.arg(t.p.out)
  
  default.args <- switch(
    sim.dist,
    norm    = list(sd = 1),
    t       = list(df = 3),
    lognorm = list(meanlog = 0, sdlog = 1),
    gamma   = list(shape = 2, rate = 1)
  )
  
  dist.args <- utils::modifyList(default.args, dist.args)
  
  # Select random generation function
  rfun <- switch(
    sim.dist,
    norm    = rnorm,
    t       = rt,
    lognorm = rlnorm,
    gamma   = rgamma
  )
  
  # Determine theoretical population mean
  pop.mean <- switch(
    sim.dist,
    norm    = 0,
    t       = 0,
    lognorm = exp(dist.args$meanlog + dist.args$sdlog^2 / 2),
    gamma   = dist.args$shape / dist.args$rate
  )
  
  # Check basic requirements
  if (n.t < 2 || n.legendre < 1)
    stop("Sample sizes and chunk.neyman must be positive, with n.t >= 2 and chunk.neyman >= 2.")
  
  if (sim.dist == "t" && dist.args$df <= 1)
    stop("Student-t df must exceed 1 for a finite population mean.")
  
  t.data <- matrix(
    do.call(rfun, c(list(n = n.t * n.legendre), dist.args)) - pop.mean,
    nrow = n.t,
    ncol = n.legendre
  )
  
  # Calculate one-sample t-statistics
  xbar <- colMeans(t.data)
  s2 <- (colSums(t.data^2) - n.t * xbar^2) / (n.t - 1)
  test.statistic.temp <- xbar / sqrt(s2 / n.t)
  
  # Calculate left-tailed p-values
  pvals <- switch(
    t.p.out,
    two  = 2 * pt(-abs(test.statistic.temp), df = n.t - 1),
    left = pt(test.statistic.temp, df = n.t - 1, lower.tail = TRUE)
  )
  
  return(pvals)
}

# =====

tTest.pvalRaw.nlist <- function(
    n.legendre, n.t.list,
    sim.dist = c("norm", "t", "lognorm", "gamma"),
    dist.args = list(), t.p.out = c("two", "left")
){
  
  # Select distribution and default parameters
  sim.dist <- match.arg(sim.dist); t.p.out <- match.arg(t.p.out)
  
  default.args <- switch(
    sim.dist,
    norm    = list(sd = 1),
    t       = list(df = 3),
    lognorm = list(meanlog = 0, sdlog = 1),
    gamma   = list(shape = 2, rate = 1)
  )
  
  dist.args <- utils::modifyList(default.args, dist.args)
  
  # Select random generation function
  rfun <- switch(
    sim.dist,
    norm    = rnorm,
    t       = rt,
    lognorm = rlnorm,
    gamma   = rgamma
  )
  
  # Determine theoretical population mean
  pop.mean <- switch(
    sim.dist,
    norm    = 0,
    t       = 0,
    lognorm = exp(dist.args$meanlog + dist.args$sdlog^2 / 2),
    gamma   = dist.args$shape / dist.args$rate
  )
  
  # Check basic requirements
  if (any(n.t.list < 2) || n.legendre < 1)
    stop("Sample sizes and chunk.neyman must be positive, with n.t >= 2 and chunk.neyman >= 2.")
  
  if (sim.dist == "t" && dist.args$df <= 1)
    stop("Student-t df must exceed 1 for a finite population mean.")
  
  t.data.all <- map(n.t.list, .f = function(n.t){

    t.data <- matrix(
      do.call(rfun, c(list(n = n.t * n.legendre), dist.args)) - pop.mean,
      nrow = n.t,
      ncol = n.legendre
    )
  
    # Calculate one-sample t-statistics
    xbar <- colMeans(t.data)
    s2 <- (colSums(t.data^2) - n.t * xbar^2) / (n.t - 1)
    test.statistic.temp <- xbar / sqrt(s2 / n.t)
    
    # Calculate left-tailed p-values
    pvals <- switch(
      t.p.out,
      two  = 2 * pt(-abs(test.statistic.temp), df = n.t - 1),
      left = pt(test.statistic.temp, df = n.t - 1, lower.tail = TRUE)
    )
    
    return(pvals)
    }
  )
  
  return(t.data.all)
}

# ===== Plot the t-test p-values and the joint neyman p-values

pianoPanel <- function(pvalRaw, pvalNeyman, n.t.list, df){
  n <- length(n.t.list)
  par(mfrow = c(2, n))
  for (i in 1:n){
    title <- sprintf("Raw t-test Piano: t(df = %s, n = %s)", df, n.t.list[i])
    hist(pvalRaw[[i]], breaks = seq(0,1,length.out = 41), main = title, xlab = "p")
  }
  for (i in 1:n){
    title <- sprintf("Neyman Smoothness Test Piano")
    hist(pvalNeyman[[i]][,"joint"], breaks = seq(0,1,length.out = 41), main = title, xlab = "p")
  }
}

# ===== Basic warpper: only coded for t-distribution and df ... not logN or Gamma

legendre.wrap <- function(
    n.neyman, n.legendre, n.t.list,
    sim.dist = c("norm", "t", "lognorm", "gamma"),
    dist.args = list(),
    chunk.neyman = 10, # compute 10 Neyman stats at a time
    t.p.out = c("two", "left"), ney.p.out = c("two", "left"))
  {
  
  pvalNeyman <- tTest.legendre.chunk.nlist(n.neyman, n.legendre, n.t.list, sim.dist, 
                                          dist.args, chunk.neyman, 
                                          t.p.out, ney.p.out)
  pvalRaw <- tTest.pvalRaw.nlist(n.legendre, n.t.list, sim.dist, dist.args, t.p.out) 
  
  pianoPanel(pvalRaw, pvalNeyman, n.t.list, dist.args$df)
}

chiLimit <- function(stat, df, level = 0.95) {
	if (pchisq(stat, df) < 1-level) return(0)
	uniroot(
		function(l) pchisq(stat, df, ncp = l) - (1-level),
		c(0, stat),
		extendInt = "downX"
	)$root
}
