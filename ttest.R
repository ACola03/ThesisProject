tTest.statSim.chunk <- function(
    n.pianos, n.pvals, n.t,
    sim.dist = c("norm", "t", "lognorm", "gamma"),
    dist.args = list(),
    chunk.pianos = 10
){
  
  # Select distribution and default parameters
  sim.dist <- match.arg(sim.dist)
  
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
  if (n.t < 2 || n.pianos < 1 || n.pvals < 1 || chunk.pianos < 1)
    stop("Sample sizes and chunk.pianos must be positive, with n.t >= 2.")
  
  if (sim.dist == "t" && dist.args$df <= 1)
    stop("Student-t df must exceed 1 for a finite population mean.")
  
  # Initialize output
  AD.out <- numeric(n.pianos)
  ZAD.out <- numeric(n.pianos)
  CVM.out <- numeric(n.pianos)
  ZCVM.out <- numeric(n.pianos)
  
  starts <- seq(1, n.pianos, chunk.pianos)
  
  for (start in starts){
    
    end <- min(start + chunk.pianos - 1, n.pianos)
    k <- end - start + 1
    tot.pvals <- k * n.pvals
    
    # Generate and centre observations
    t.data <- matrix(
      do.call(rfun, c(list(n = n.t * tot.pvals), dist.args)) - pop.mean,
      nrow = n.t,
      ncol = tot.pvals
    )
    
    # Calculate one-sample t-statistics
    xbar <- colMeans(t.data)
    
    s2 <- (colSums(t.data^2) - n.t * xbar^2) / (n.t - 1)
    
    test.statistic.temp <- xbar / sqrt(s2 / n.t)
    
    # Calculate left-tailed p-values
    pvals <- pt(test.statistic.temp, df = n.t - 1)
    
    pval.matrix <- matrix(
      pvals,
      nrow = n.pvals,
      ncol = k
    )
    
    # Calculate PianoPlot statistics
    AD.out[start:end] <- apply(
      pval.matrix, 2, zhang.AD, use.zhang = FALSE
    )
    
    ZAD.out[start:end] <- apply(
      pval.matrix, 2, zhang.AD, use.zhang = TRUE
    )
    
    CVM.out[start:end] <- apply(
      pval.matrix, 2, zhang.CVM, use.zhang = FALSE
    )
    
    ZCVM.out[start:end] <- apply(
      pval.matrix, 2, zhang.CVM, use.zhang = TRUE
    )
  }
  
  return(list(
    AD = AD.out,
    ZAD = ZAD.out,
    CVM = CVM.out,
    ZCVM = ZCVM.out
  ))
}

zhang.dist <- function(n.pianos, n.pvals, n.ts, df.ts,
                       lognorm.args = list(
                         list(meanlog = 0, sdlog = 0.25),
                         list(meanlog = 0, sdlog = 0.5),
                         list(meanlog = 0, sdlog = 1)
                       ),
                       gamma.args = list(
                         list(shape = 0.5, rate = 1),
                         list(shape = 1, rate = 1),
                         list(shape = 2, rate = 1),
                         list(shape = 5, rate = 1)
                       ),
                       chunk.pianos = 25) {
  
  ncols <- length(n.ts) * (
    1 + length(df.ts) +
      length(lognorm.args) + length(gamma.args)
  )
  
  AD.list <- vector("list", ncols)
  ZAD.list <- vector("list", ncols)
  CVM.list <- vector("list", ncols)
  ZCVM.list <- vector("list", ncols)
  
  config.names <- character(ncols)
  
  k <- 1
  
  for (n.t in n.ts) {
    
    # Normal baseline
    allStats.null <- tTest.statSim.chunk(
      n.pianos = n.pianos,
      n.pvals = n.pvals,
      n.t = n.t,
      sim.dist = "norm",
      chunk.pianos = chunk.pianos
    )
    
    config <- paste0("Normal(n = ", n.t, ")")
    
    AD.list[[k]] <- allStats.null$AD
    ZAD.list[[k]] <- allStats.null$ZAD
    CVM.list[[k]] <- allStats.null$CVM
    ZCVM.list[[k]] <- allStats.null$ZCVM
    
    config.names[k] <- config
    k <- k + 1
    
    # Student-t alternatives
    for (df.t in df.ts) {
      
      allStats.alt <- tTest.statSim.chunk(
        n.pianos = n.pianos,
        n.pvals = n.pvals,
        n.t = n.t,
        sim.dist = "t",
        dist.args = list(df = df.t),
        chunk.pianos = chunk.pianos
      )
      
      config <- paste0("t(n = ", n.t, ", df = ", df.t, ")")
      
      AD.list[[k]] <- allStats.alt$AD
      ZAD.list[[k]] <- allStats.alt$ZAD
      CVM.list[[k]] <- allStats.alt$CVM
      ZCVM.list[[k]] <- allStats.alt$ZCVM
      
      config.names[k] <- config
      k <- k + 1
    }
    
    # Lognormal alternatives
    for (args in lognorm.args) {
      
      allStats.alt <- tTest.statSim.chunk(
        n.pianos = n.pianos,
        n.pvals = n.pvals,
        n.t = n.t,
        sim.dist = "lognorm",
        dist.args = args,
        chunk.pianos = chunk.pianos
      )
      
      config <- paste0(
        "lognorm(n = ", n.t, ", ",
        paste(names(args), unlist(args),
              sep = " = ", collapse = ", "), ")"
      )
      
      AD.list[[k]] <- allStats.alt$AD
      ZAD.list[[k]] <- allStats.alt$ZAD
      CVM.list[[k]] <- allStats.alt$CVM
      ZCVM.list[[k]] <- allStats.alt$ZCVM
      
      config.names[k] <- config
      k <- k + 1
    }
    
    # Gamma alternatives
    for (args in gamma.args) {
      
      allStats.alt <- tTest.statSim.chunk(
        n.pianos = n.pianos,
        n.pvals = n.pvals,
        n.t = n.t,
        sim.dist = "gamma",
        dist.args = args,
        chunk.pianos = chunk.pianos
      )
      
      config <- paste0(
        "gamma(n = ", n.t, ", ",
        paste(names(args), unlist(args),
              sep = " = ", collapse = ", "), ")"
      )
      
      AD.list[[k]] <- allStats.alt$AD
      ZAD.list[[k]] <- allStats.alt$ZAD
      CVM.list[[k]] <- allStats.alt$CVM
      ZCVM.list[[k]] <- allStats.alt$ZCVM
      
      config.names[k] <- config
      k <- k + 1
    }
  }
  
  names(AD.list) <- config.names
  names(ZAD.list) <- config.names
  names(CVM.list) <- config.names
  names(ZCVM.list) <- config.names
  
  return(list(
    AD = AD.list,
    ZAD = ZAD.list,
    CVM = CVM.list,
    ZCVM = ZCVM.list
  ))
}

zhang.computePower <- function(stat.data, n, distribution, alpha){
  
  # Extract Normal baseline
  null.name <- paste0("Normal(n = ", n, ")")
  nullDist <- stat.data[[null.name]]
  
  if (is.null(nullDist))
    stop("Normal baseline not found.")
  
  # Extract alternatives for selected distribution and n
  pattern <- paste0("^", distribution, "\\(n = ", n, ",")
  
  alt.names <- grep(
    pattern,
    names(stat.data),
    value = TRUE
  )
  
  if (length(alt.names) == 0)
    stop("No matching distribution configurations found.")
  
  # Calculate critical value
  criticalValue <- quantile(nullDist, probs = 1 - alpha)
  
  # Calculate power for each configuration
  powerVec <- sapply(alt.names, function(config){
    mean(stat.data[[config]] > criticalValue)
  })
  
  # Extract configuration labels (parameters)
  config.labels <- sub(
    paste0("^", distribution, "\\(n = ", n, ",\\s*"),
    "",
    alt.names
  )
  
  config.labels <- sub("\\)$", "", config.labels)
  
  return(data.frame(
    n = n,
    config = config.labels,
    power = as.numeric(powerVec)
  ))
}

zhang.powerAnalysis <- function(data, statistic, distribution, alpha){
  
  statistic <- match.arg(
    statistic,
    choices = c("AD", "ZAD", "CVM", "ZCVM")
  )
  
  distribution <- match.arg(
    distribution,
    choices = c("t", "lognorm", "gamma")
  )
  
  stat.data <- data[[statistic]]
  
  # Identify available sample sizes from Normal baselines
  null.names <- grep(
    "^Normal\\(n = ",
    names(stat.data),
    value = TRUE
  )
  
  ns <- as.numeric(
    sub("^Normal\\(n = ([0-9]+)\\)$", "\\1", null.names)
  )
  
  # Calculate power across sample sizes
  power.list <- lapply(ns, function(n){
    
    zhang.computePower(
      stat.data = stat.data,
      n = n,
      distribution = distribution,
      alpha = alpha
    )
    
  })
  
  power.data <- do.call(rbind, power.list)
  
  # Treat parameter configurations as categorical
  power.data$config <- factor(
    power.data$config,
    levels = unique(power.data$config)
  )
  
  return(power.data)
}

powerPlot <- function(data, statistic, distribution, alpha){
  
  plot.data <- zhang.powerAnalysis(
    data = data,
    statistic = statistic,
    distribution = distribution,
    alpha = alpha
  )
  
  p_power <- ggplot(
    plot.data,
    aes(
      x = n,
      y = power,
      colour = config,
      group = config
    )
  ) +
    geom_line(linewidth = 1) +
    geom_point(size = 2) +
    scale_x_continuous(
      breaks = sort(unique(plot.data$n))
    ) +
    coord_cartesian(ylim = c(0, 1)) +
    theme_classic() +
    labs(
      title = paste0(statistic, " Power: ", distribution),
      x = "Sample Size (N)",
      y = bquote("Statistical Power (" * alpha == .(alpha) * ")"),
      colour = "Configuration"
    )
  
  print(p_power)
  return(invisible(p_power))
}

zhang.KS <- function(data, use.zhang = TRUE){
  n <- length(data)
  d <- sort(data)
  i <- 1:n
  
  if (use.zhang)
    zks <- max((i-0.5)*log((i-0.5)/(n*d)) + (n-i+0.5)*log((n-i+0.5)/(n*(1-d))))
  else 
    zks <- max(max(abs((1:n)/n) - d), max(abs(d - (1:n - 1)/n)))
  return(zks)
}

zhang.AD <- function(data, use.zhang = TRUE){
  n <- length(data)
  d <- sort(data)
  i <- 1:n
  
  if (use.zhang)
    ad <- -sum((log(d)/(n-i+0.5)) + (log(1-d)/(i-0.5)))
  else 
    ad <- (-2/n)*sum((i-0.5)*log(d) + (n-i+0.5)*log(1-d)) - n
  return(ad)
}

zhang.CVM <- function(data, use.zhang = TRUE){
  n <- length(data)
  d <- sort(data)
  i <- 1:n
  
  if (use.zhang)
    cvm <- sum((log((1/d - 1) / ((n-0.5)/(i-0.75) - 1)))^2)
  else 
    cvm <- sum((d - (i-0.5)/n)^2) + 1/(12*n)
  return(cvm)
}

results <- zhang.dist(
  n.pianos = 1000,
  n.pvals = 2000,
  n.ts = c(10,100),
  df.ts = c(2, 3, 6, 10, 30),
  
  lognorm.args = list(
    list(meanlog = 0, sdlog = 0.25),
    list(meanlog = 0, sdlog = 0.5),
    list(meanlog = 0, sdlog = 1)
  ),
  
  gamma.args = list(
    list(shape = 0.5, rate = 1),
    list(shape = 1, rate = 1),
    list(shape = 2, rate = 1)
  ),
  
  chunk.pianos = 25
)

# Student-t
powerPlot(results, "ZAD", "t", alpha = 0.05)

# Lognormal
powerPlot(results, "ZAD", "lognorm", alpha = 0.05)

# Gamma
powerPlot(results, "ZAD", "gamma", alpha = 0.05)