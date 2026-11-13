# =======
# Evaluating t-Test validity when assumptions aren't satisfied
# using the Neyman Smooth Test (Legendre Polynomials)
# =======

library(tidyr)
library(ggplot2)
library(dplyr)
library(purrr)
library(checkPlotR)

# ===== Perform Neyman Smoothness Test

tTest.legendre.chunk.nlist <- function(
    n.neyman, n.legendre, n.t.list,
    sim.dist = c("norm", "t", "lognorm", "gamma"),
    dist.args = list(), 
    which.basis = c(1,2,3,4,5,6),
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
  
  # Set names in n.t.list for nice map() outputs
  names(n.t.list) <- unique(n.t.list)
  
  neyman <- map(n.t.list, .f = function(n.t){
    
    # Initialize output
    neyman.stat <- matrix(
      NA_real_,
      nrow = n.neyman,
      ncol = length(which.basis) + 1,
      dimnames = list(
        NULL,
        c(paste0("theta", which.basis), "joint")
      )
    )
    neyman.p <- neyman.stat
    
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
      theta.mat <- legendre.coef(pval.matrix, which.basis)
      stats.pvals <- legendre.pval(n.legendre, theta.mat, ney.p.out)
      neyman.p[start:end,] <- stats.pvals$p
      neyman.stat[start:end,] <- stats.pvals$stat
    }
    return(list("stat" = neyman.stat, "p" = neyman.p))
  }
  )
  return(neyman)
}

legendre.coef <- function(p.vals, which.basis = c(1,2,3,4,5,6)){
  
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
  return(theta.mat[,which.basis])
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
  
  p <- cbind(pval.marginal, joint = pval.joint)
  stat <- cbind(z.mat, joint = statistic.joint)
  
  return(list("stat" = stat, "p" = p))
}

# ===== Get p-values of the t-test itself

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

diagnosticPanel <- function(pvalRaw, neyman.obs, neyman.null, n.t.list, sim.dist, dist.args, which.basis, level = 0.95, quantile = 0.95){ 
  n <- length(n.t.list)
  
  titlesRaw <- switch (
    sim.dist,
    norm = sprintf("Raw t-test: N(sd = %s, n = %s)", dist.args$sd, n.t.list),
    t = sprintf("Raw t-test: t(df = %s, n = %s)", dist.args$df, n.t.list),
    lognorm = sprintf("Raw t-test: logN(meanLog = %s, sdLog = %s, n = %s)", dist.args$meanlog, dist.args$sdlog, n.t.list),
    gamma = sprintf("Raw t-test: gamma(shape = %s, rate = %s, n = %s)", dist.args$shape, dist.args$rate, n.t.list)
  )
  
  ncp.null <- ncps.tibble(neyman.null, level)
  ncp.obs <- ncps.tibble(neyman.obs, level)
  
  plots <- list()
  
  for (i in 1:n){
    title <- titlesRaw[i]
    plots[[i]] <- pianoPlot(pvalRaw[[i]], titlesRaw[i])
  }
  for (i in 1:n){
    title <- sprintf("Neyman Smooth Test (b = %s)", toString(which.basis))
    plots[[n+i]] <- pianoPlot(neyman.obs[[i]][["p"]][,"joint"], title)
  }
  for (i in 1:n){
    ncp.threshold <- quantile(ncp.null[[i]]$upper, probs = quantile)
    plots[[2*n+i]] <- volcanoPlot(ncp.obs[[i]], target = ncp.threshold, orderFun = volcano)
  }
  gridExtra::grid.arrange(grobs = plots, nrow = 3, ncol = n)
}


# ===== Wrapper:

legendre.wrap <- function(
    n.neyman, n.legendre, n.t.list,
    sim.dist = c("norm", "t", "lognorm", "gamma"),
    dist.args = list(),
    which.basis = c(1,2,3,4,5,6),
    chunk.neyman = 10, # compute 10 Neyman stats at a time
    t.p.out = c("two", "left"), ney.p.out = c("two", "left"),
    level = 0.95, quantile = 0.95)
{
  
  neyman.null <- tTest.legendre.chunk.nlist(n.neyman, n.legendre, n.t.list, "norm", 
                                            dist.args = list(sd = 1), which.basis, 
                                            chunk.neyman, t.p.out, ney.p.out)
  
  neyman.obs <- tTest.legendre.chunk.nlist(n.neyman, n.legendre, n.t.list, sim.dist, 
                                           dist.args, which.basis, chunk.neyman, 
                                           t.p.out, ney.p.out)
  
  pvalRaw <- tTest.pvalRaw.nlist(n.legendre, n.t.list, sim.dist, dist.args, t.p.out) 
  
  diagnosticPanel(pvalRaw, neyman.obs, neyman.null, n.t.list, sim.dist, dist.args, which.basis, level, quantile)
}



# ===== Non-Centrality Parameter

chiLimit <- function(stat, df, level = 0.95) {
	if (pchisq(stat, df) < 1-level) return(0)
	uniroot(
		function(l) pchisq(stat, df, ncp = l) - (1-level),
		c(0, stat),
		extendInt = "downX"
	)$root
}


ncps.tibble <- function(neyman, level = 0.95){
  return(
    map(neyman, .f = function(x){
      stats <- x$stat
      df <- ncol(stats) - 1 # because `joint` included
      smooth.stat <- stats[,"joint"]
      tibble("stat" = smooth.stat, "upper" = sapply(smooth.stat, chiLimit, df, level))
    }
    )
  )
}


volcanoPlot <- function(tf
                        , target
                        , orderFun = slug
                        , conf = 0.95
                        , opacity = 0.2
                        , size = 0.1
                        , title = "Volcano plot"
                        , targNum = 1e3){
  
  thinner <- max(floor(length(tf$stat)/targNum), 1)
  thinned <- tf[seq(thinner, length(tf$stat), thinner),]
  
  return(ggplot(
    orderFun(thinned)
    , aes(x = quantile, y = upper, ymin = 0, ymax = upper)
  )
  + geom_pointrange(alpha = opacity, size = size
                    , aes(color = ifelse(upper > target, "red", "grey")))
  + geom_hline(yintercept = target, color = "blue")
  + geom_vline(
    xintercept = c((1 - conf), conf)
    , lty = 2, col = "red"
  )
  + xlab("index")
  + ylab("upper")
  + ggtitle(title)
  + scale_x_continuous(expand = c(0,0)
                       , breaks = c((1 - conf), 0.25, 0.5, 0.75, conf)
                       , labels = function(breaks){signif(breaks, 3)})
  + scale_color_manual(values = c("grey" = "grey", "red" = "red"))
  + guides(color = "none")
  + theme_classic()
  )
}


volcano <- function(tf){
  numEst <- nrow(tf)
  return(tf 
         %>% arrange(stat) 
         %>% mutate(quantile = ((1:numEst)-1/2)/numEst))
}


pianoPlot <- function(x, title) {
  ggplot(data.frame(p = x), aes(x = p)) +
    geom_histogram(
      breaks = seq(0, 1, length.out = 41),
      color = "black",
      fill = "grey"
    ) +
    labs(
      title = title,
      x = "p",
      y = "Frequency"
    ) +
    theme_classic()
}



















