# T-Test Analysis Functions

# T-DISTRIBUTION SIMULATION: an entire piano of p-values
tTest.statSim <- function(n.t, n.unif, df = 3, use.normal = FALSE){
  
  if (!use.normal){
    t.data <- matrix(rt(n.t * n.unif, df = df), ncol = n.unif)
  } else if (use.normal){
    t.data <- matrix(rnorm(n.t * n.unif), ncol = n.unif)
  }
  
  test.statistic <- apply(t.data, MARGIN = 2, 
                          FUN = function(sim.dat){ 
                            (mean(sim.dat) - 0)/(sd(sim.dat)/sqrt(n.t))
                          })
  pvals <- pt(test.statistic, df = n.t - 1)
}

tTest.statSim.chunk <- function(n.pianos, n.pvals, n.t, df.t,  
                       chunk.pianos = 10, use.normal = FALSE){
  
  AD.out <- numeric(n.pianos); ZAD.out <- numeric(n.pianos)
  CVM.out <- numeric(n.pianos); ZCVM.out <- numeric(n.pianos)
  starts <- seq(1, n.pianos, chunk.pianos)
  
  for (start in starts){
    
    end <- min(start + chunk.pianos - 1, n.pianos)
    k <- (end - start + 1)
    tot.pvals <- k * n.pvals
    
    if (use.normal){
      t.data <- matrix(
        rnorm(n.t * tot.pvals),
        ncol = tot.pvals, nrow = n.t
      ) 
    } else if (!use.normal){
      t.data <- matrix(
        rt(n.t * tot.pvals, df = df.t),
        ncol = tot.pvals, nrow = n.t
      ) 
    }
    
    xbar <- colMeans(t.data)
    var <- (colSums(t.data^2) - n.t*xbar^2)/(n.t - 1)
    test.statistic.temp <- (xbar - 0)/(sqrt(var / n.t))
    
    pvals <- pt(test.statistic.temp, df = n.t - 1)
    pval.matrix <- matrix(pvals, nrow = n.pvals, ncol = k)
    
    AD.out[start:end] <- apply(pval.matrix, 2,
                               `zhang.AD`, use.zhang = FALSE)
    ZAD.out[start:end] <- apply(pval.matrix, 2,
                               `zhang.AD`, use.zhang = TRUE)
    CVM.out[start:end] <- apply(pval.matrix, 2,
                               `zhang.CVM`, use.zhang = FALSE)
    ZCVM.out[start:end] <- apply(pval.matrix, 2,
                                `zhang.CVM`, use.zhang = TRUE)
    
  }
  
  return(list("AD" = AD.out, "ZAD" = ZAD.out,
              "CVM" = CVM.out, "ZCVM" = ZCVM.out))
}

# ========== Beta-Fit Functions

nll <- function(mean, shape, x) {
  a <- 2*mean / shape
  b <- 2*(1 - mean) / shape
  -sum(dbeta(x, shape1 = a, shape2 = b, log = TRUE))
}

mle.est <- function(dat){
  mle2(
    nll,
    start = list(mean = 0.5, shape = 1),
    data = list(x = dat),
    method = "L-BFGS-B",
    lower = c(mean = 1e-4, shape = 1e-4),
    upper = c(mean = 1 - 1e-4, shape = Inf),
    trace = FALSE
  )
}

profileMean <- function(dat){
  mle2(
    nll,
    start = list(shape = 1),
    fixed = list(mean = 0.5),
    data = list(x = dat),
    method = "L-BFGS-B",
    lower = c(shape = 1e-4),
    upper = c(shape = Inf)
  )
}

profileShape <- function(dat){
  mle2(
    nll,
    start = list(mean = 0.5),
    fixed = list(shape = 1),
    data = list(x = dat),
    method = "L-BFGS-B",
    lower = c(mean = 1e-4),
    upper = c(mean = 1 - 1e-4)
  )
}

paraBoot.full <- function(n.unif, boots, sim.dist = "beta", shape1 = 1, shape2 = 1, n.t = 30, df.t = 3, conf = 0.95, meta.calc = FALSE){
  alphas <- numeric(boots); betas <- numeric(boots)
  means <- numeric(boots); shapes <- numeric(boots)
  means.se <- numeric(boots); shapes.se <- numeric(boots)
  
  means.lower.wald <- numeric(boots); means.upper.wald <- numeric(boots)
  means.lower.prof <- numeric(boots); means.upper.prof <- numeric(boots)
  means.pval.wald <- numeric(boots); means.pval.prof <- numeric(boots); 
  
  shapes.lower.wald <- numeric(boots); shapes.upper.wald <- numeric(boots)
  shapes.lower.prof <- numeric(boots); shapes.upper.prof <- numeric(boots)
  shapes.pval.wald <- numeric(boots); shapes.pval.prof <- numeric(boots); 
  
  success_count <- 0; success_flag <- logical(boots)
  z <- qnorm(1-(1-conf)/2)
  
  for (boot in 1:boots){
    
    if (sim.dist == "beta"){
      dat <- rbeta(n.unif, shape1, shape2)
    } else if (sim.dist == "t"){
      dat <- tTest.statSim(n.t, n.unif, df.t)
    }
    
    # for convergence failures
    fit_result <- tryCatch({
      suppressWarnings({
        fit <- mle.est(dat)
        # check if it converged (0 means success)
        if (fit@details$convergence != 0) stop("Non-zero convergence code")
        fit
      })
    }, error = function(e) {
      NULL # NULL if fails to converge
    })
    
    # if converged, store the values
    if (!is.null(fit_result)) {
      success_count <- success_count + 1; success_flag[boot] <- TRUE # increment successful fits
      coefs <- coef(fit_result); vcovs <- vcov(fit_result) # extract coefficients and variances
      pfint <- profile(fit_result); pconf <- confint(pfint) # profile intervals
      mean <- unname(coefs["mean"]); shape <- unname(coefs["shape"]) # assign coefficients
      mean.se <- sqrt(vcovs[1,1]); shape.se <- sqrt(vcovs[2,2]) # Wald Std.Err from inverting hessian
      
      # alpha, beta 
      alphas[boot] <- 2 * mean / shape
      betas[boot] <- 2 * (1 - mean) / shape
      
      # mean: estimate, stdErr, intervals (wald, profile)
      means[boot] <- mean; means.se[boot] <- mean.se 
      means.lower.wald[boot] <- mean-z*mean.se; means.upper.wald[boot] <- mean+z*mean.se
      means.lower.prof[boot] <- pconf["mean",1]; means.upper.prof[boot] <- pconf["mean",2]
      
      # shape: estimate, stdErr, intervals (wald, profile)
      shapes[boot] <- shape; shapes.se[boot] <- shape.se 
      shapes.lower.wald[boot] <- shape-z*shape.se; shapes.upper.wald[boot] <- shape+z*shape.se
      shapes.lower.prof[boot] <- pconf["shape",1]; shapes.upper.prof[boot] <- pconf["shape",2]
      
      # meta piano calculations
      if (meta.calc){
        # profile p-value
        restricted.fitMean <- profileMean(dat); restricted.fitShape <- profileShape(dat) 
        l.global <- logLik(fit_result); l.mean <- logLik(restricted.fitMean); l.shape <- logLik(restricted.fitShape) # likelihoods
        chi.mean <- 2*(l.global - l.mean); chi.shape <- 2*(l.global - l.shape) # likelihood ratio statistics
        
        mean.pval.prof.temp <- pchisq(chi.mean, df = 1, lower.tail = FALSE)/2;
        means.pval.prof[boot] <- ifelse(mean > 0.5, 1-mean.pval.prof.temp, mean.pval.prof.temp) # profile p-value: check p/2 or 1-p/2
        
        shape.pval.prof.temp <- pchisq(chi.shape, df = 1, lower.tail = FALSE)/2 
        shapes.pval.prof[boot] <- ifelse(shape > 1, 1-shape.pval.prof.temp, shape.pval.prof.temp) # profile p-value: check p/2 or 1-p/2
        
        # wald p-value
        norm.mean <- (mean - 0.5)/mean.se; norm.shape <- (shape - 1)/shape.se # convert to StdNorm   
        means.pval.wald[boot] <- pnorm(norm.mean); shapes.pval.wald[boot] <- pnorm(norm.shape) # wald p-value
      }
    } 
  }
  
  cat(sprintf("Successful fits: %d out of %d (%.1f%%)\n", 
              success_count, boots, (success_count/boots)*100))
  
  # filter out those which didn't converge
  success_indices <- which(success_flag)
  results <- list(
    "mean" = means[success_indices], 
    "mean.se" =  means.se[success_indices], 
    "mean.lower.wald" =  means.lower.wald[success_indices], 
    "mean.upper.wald" =  means.upper.wald[success_indices],
    "mean.lower.prof" =  means.lower.prof[success_indices],
    "mean.upper.prof" =  means.upper.prof[success_indices], 
    "mean.pval.wald" = means.pval.wald[success_indices],
    "mean.pval.prof" = means.pval.prof[success_indices],
    "shape" =  shapes[success_indices],
    "shape.se" =  shapes.se[success_indices], 
    "shape.lower.wald" =  shapes.lower.wald[success_indices], 
    "shape.upper.wald" =  shapes.upper.wald[success_indices], 
    "shape.lower.prof" =  shapes.lower.prof[success_indices], 
    "shape.upper.prof" =  shapes.upper.prof[success_indices], 
    "shape.pval.wald" = shapes.pval.wald[success_indices],
    "shape.pval.prof" = shapes.pval.prof[success_indices],
    "shape1" =  alphas[success_indices], 
    "shape2" =  betas[success_indices],
    "success_count" = success_count
  )
  
  return(results)
}

# ========== SlugPlot & Coverage

rangePlot <- function(tf
                      , target = "sampleMean"
                      , orderFun = slug
                      , conf = 0.95
                      , subplot.type = "full"
                      , opacity = 0.2
                      , size = 0.1
                      , title = "Range plot"
                      , targNum = 1e3){
  if(target == "sampleMean"){
    target <- mean(tf$est)
  }
  
  thinner <- max(floor(length(tf$est)/targNum), 1) # is this reducing the # of CIs plotted
  thinned <- tf[seq(thinner, length(tf$est), thinner),]
  plot.data <- orderFun(thinned)
  
  if (subplot.type == "left"){
    plot.data <- plot.data |> filter(quantile <= 0.05)
    breaks <- seq(0, 0.05, 0.01)
    xintercept <- c((1 - conf)/2)
    x_lab <- NULL; y_lab <- NULL; plot_title <- NULL
  } else if (subplot.type == "right"){
    plot.data <- plot.data |> filter(quantile >= 0.95)
    breaks <- seq(0.95, 1, 0.01)
    xintercept <- c(1 - (1 - conf)/2)
    x_lab <- NULL; y_lab <- NULL; plot_title <- NULL
  } else if (subplot.type == "full"){
    breaks <- c((1 - conf)/2, 0.25, 0.5, 0.75, conf + (1 - conf)/2)
    xintercept <- c((1 - conf)/2, 1 - (1 - conf)/2)
    x_lab <- "Quantile"; y_lab <- "Estimate"; plot_title <- title
  }
  
  slugPlot <- ggplot(plot.data, 
                     aes(x = quantile, y = est, ymin = lower, ymax = upper)) + 
    geom_pointrange(alpha = opacity, size = size, 
                    aes(color = ifelse(lower > target | upper < target, "red", "grey"))) +
    geom_hline(yintercept = target, color = "blue") +
    geom_vline(xintercept = xintercept,
               lty = 2, col = "red") +
    scale_x_continuous(expand = c(0,0.0005)
                       , breaks = breaks
                       , labels = function(breaks){signif(breaks, 3)}) + 
    scale_color_manual(values = c("grey" = "grey", "red" = "red")) +
    ggtitle(plot_title) +
    xlab(x_lab) +
    ylab(y_lab) +
    guides(color = "none") +
    theme_classic()
  
  return(slugPlot)
}

showSlugs <- function(param.dat, shape1, shape2, n, boots, alpha = 0.05, use.profile = TRUE, show.subplots = TRUE){
  
  # extract mean and shape data (wald or profile)
  means <- data.frame("est" = param.dat$mean, 
                      "lower" = if_else(rep(use.profile, param.dat$success_count), param.dat$mean.lower.prof, param.dat$mean.lower.wald),  
                      "upper" = if_else(rep(use.profile, param.dat$success_count), param.dat$mean.upper.prof, param.dat$mean.upper.wald))
  shapes <- data.frame("est" = param.dat$shape, 
                       "lower" = if_else(rep(use.profile, param.dat$success_count), param.dat$shape.lower.prof, param.dat$shape.lower.wald),
                       "upper" = if_else(rep(use.profile, param.dat$success_count), param.dat$shape.upper.prof, param.dat$shape.upper.wald))
  
  # slug plots of mean and shape
  interval.str <- ifelse(use.profile, "Profile Interval", "Wald Interval")
  
  meanPlot <- rangePlot(means, target = 0.5, orderFun = slug, title = paste0("SlugPlot of Mean: ", interval.str), targNum = boots/2)
  shapePlot <- rangePlot(shapes, target = 1, orderFun = slug, title = paste0("SlugPlot of Shape: ", interval.str), targNum = boots/2)
  
  if (show.subplots){
    meanPlot.left <- rangePlot(means, target = 0.5, orderFun = slug, title = "", targNum = boots/2, subplot.type = "left")
    meanPlot.right <- rangePlot(means, target = 0.5, orderFun = slug, title = "", targNum = boots/2, subplot.type = "right")
    meanPlot <- meanPlot / (meanPlot.left + meanPlot.right) + plot_layout(heights = c(4, 2))
    
    shapePlot.left <- rangePlot(shapes, target = 1, orderFun = slug, title = "", targNum = boots/2, subplot.type = "left")
    shapePlot.right <- rangePlot(shapes, target = 1, orderFun = slug, title = "", targNum = boots/2, subplot.type = "right")
    shapePlot <- shapePlot / (shapePlot.left + shapePlot.right) + plot_layout(heights = c(4, 2))
  }
  
  return(c(meanPlot, shapePlot))
}

coveragePlot <- function(coverage.list,
                         title = "Interval Coverage Analysis: Wald vs Profile", 
                         ymin = 0, ymax = 0.05, use.interval = TRUE){
  
  coverage.dat <- imap(coverage.list, .f = function(x, name){
    
    # Confidence Intervals: 
    mean.prof.interval.low <- mean(x$mean.upper.prof <= 0.5); mean.prof.interval.high <- mean(x$mean.lower.prof >= 0.5)
    mean.wald.interval.low <- mean(x$mean.upper.wald <= 0.5); mean.wald.interval.high <- mean(x$mean.lower.wald >= 0.5)
    shape.prof.interval.low <- mean(x$shape.upper.prof <= 1); shape.prof.interval.high <- mean(x$shape.lower.prof >= 1)
    shape.wald.interval.low <- mean(x$shape.upper.wald <= 1); shape.wald.interval.high <- mean(x$shape.lower.wald >= 1)
    
    # P-values:
    mean.prof.pval.low <- mean(x$mean.pval.prof <= 0.025); mean.prof.pval.high <- mean(x$mean.pval.prof >= 0.975)
    mean.wald.pval.low <- mean(x$mean.pval.wald <= 0.025); mean.wald.pval.high <- mean(x$mean.pval.wald >= 0.975)
    shape.prof.pval.low <- mean(x$shape.pval.prof <= 0.025); shape.prof.pval.high <- mean(x$shape.pval.prof >= 0.975)
    shape.wald.pval.low <- mean(x$shape.pval.wald <= 0.025); shape.wald.pval.high <- mean(x$shape.pval.wald >= 0.975)
    
    # Dataframe:
    mean.df <- data.frame(n = as.integer(name), "parameter" = "mean", 
                          "prof.interval.low" = mean.prof.interval.low, "prof.interval.high" = mean.prof.interval.high,
                          "wald.interval.low" = mean.wald.interval.low, "wald.interval.high" = mean.wald.interval.high,
                          "prof.pval.low" = mean.prof.pval.low, "prof.pval.high" = mean.prof.pval.high,
                          "wald.pval.low" = mean.wald.pval.low, "wald.pval.high" = mean.wald.pval.high)
    
    shape.df <- data.frame(n = as.integer(name), "parameter" = "shape", 
                           "prof.interval.low" = shape.prof.interval.low, "prof.interval.high" = shape.prof.interval.high,
                           "wald.interval.low" = shape.wald.interval.low, "wald.interval.high" = shape.wald.interval.high,
                           "prof.pval.low" = shape.prof.pval.low, "prof.pval.high" = shape.prof.pval.high,
                           "wald.pval.low" = shape.wald.pval.low, "wald.pval.high" = shape.wald.pval.high)
    
    rbind(mean.df, shape.df)
  }
  )
  
  coverage.dat <- tibble(do.call(rbind, coverage.dat))
  row.names(coverage.dat) <- NULL
  
  if (use.interval){ 
    coverage.dat_long <- coverage.dat %>%
      pivot_longer(
        cols = c(wald.interval.low, wald.interval.high, prof.interval.low, prof.interval.high),
        names_to = c("method", "metric", "direction"),
        names_sep = "\\.",
        values_to = "rate"
      )
  } else {
    coverage.dat_long <- coverage.dat %>%
      pivot_longer(
        cols = c(wald.pval.low, wald.pval.high, prof.pval.low, prof.pval.high),
        names_to = c("method", "metric", "direction"),
        names_sep = "\\.",
        values_to = "rate"
      )
  }
  
  ggplot(data = coverage.dat_long, aes(x = n, y = rate, color = direction, linetype = method)) +
    geom_point() +
    geom_line() +
    geom_hline(yintercept = c(0.02, 0.03), linetype = "dashed") +
    ylim(c(ymin, ymax)) +
    facet_wrap(~parameter, nrow = 2) +
    theme_minimal() +
    labs(
      y = "Coverage / Tail Rate",
      color = "Direction",
      linetype = "Method"
    ) + 
    ggtitle(title)
}

# ========== PianoPlot Statistic

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


zhang.dist <- function(n.pianos, n.pvals, n.ts, df.ts,
                       chunk.pianos = 25) {
  
  ncols <- length(n.ts) * (length(df.ts) + 1)
  
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
      df.t = NULL,
      chunk.pianos = chunk.pianos,
      use.normal = TRUE
    )
    
    config <- paste0("Normal(n = ", n.t, ")")
    
    AD.list[[k]] <- allStats.null$AD
    ZAD.list[[k]] <- allStats.null$ZAD
    CVM.list[[k]] <- allStats.null$CVM
    ZCVM.list[[k]] <- allStats.null$ZCVM
    
    config.names[k] <- config
    
    k <- k + 1
    
    # t-distribution alternatives
    for (df.t in df.ts) {
      
      allStats.alt <- tTest.statSim.chunk(
        n.pianos = n.pianos,
        n.pvals = n.pvals,
        n.t = n.t,
        df.t = df.t,
        chunk.pianos = chunk.pianos,
        use.normal = FALSE
      )
      
      config <- paste0("t(n = ", n.t, ", df = ", df.t, ")")
      
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
