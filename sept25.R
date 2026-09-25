# =====

# IN THIS FILE:
# 1. 4-Panel Histograms for Wald vs Profile P-Values (meta-analysis)
# 2. 2-Panel Line plots for coverage across n
# 3. 2-Panel Line plots for p-values across n
# 4. SlugPlot accompanied by a 2-panel focused view of each tail region 

library(ggplot2)
library(dplyr)
library(checkPlotR)
library(patchwork)
library(purrr)
library(tidyr)
library(bbmle)
library(MASS)

# =====

# MODIFIED SLUGPLOT:

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
  meanPlot <- rangePlot(means, target = 0.5, orderFun = slug, title = "SlugPlot of Mean", targNum = boots/2)
  shapePlot <- rangePlot(shapes, target = 1, orderFun = slug, title = "SlugPlot of Shape", targNum = boots/2)
  
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

# =====

# TABLES INTO PLOTS

coverage.250 <- readRDS("rdsFiles/coverage250_full.rds")
coverage.500 <- readRDS("rdsFiles/coverage500_full.rds")
coverage.750 <- readRDS("rdsFiles/coverage750_full.rds")
coverage.1000 <- readRDS("rdsFiles/coverage1000_full.rds")
coverage.1250 <- readRDS("rdsFiles/coverage1250_full.rds")
coverage.1500 <- readRDS("rdsFiles/coverage1500_full.rds")
coverage.1750 <- readRDS("rdsFiles/coverage1750_full.rds")
coverage.2000 <- readRDS("rdsFiles/coverage2000_full.rds")

coverage.list <- list("250" = coverage.250, "500" = coverage.500, "750" = coverage.750,
                      "1000" = coverage.1000, "1250" = coverage.1250, "1500" = coverage.1500, 
                      "1750" = coverage.1750, "2000" = coverage.2000)


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


# They are the same!
coveragePlot(coverage.list, ymin = 0.01, ymax = 0.04, use.interval = TRUE)
coveragePlot(coverage.list, ymin = 0.01, ymax = 0.04, use.interval = FALSE)

# =====

# META PIANOS

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

signif(coef(mle.est(coverage.2000$mean.pval.prof))[["shape"]], 4)

metaPianos <- function(coverage.list, n, show.flag = TRUE) {
  
  n.char <- as.character(n)
  plot.dat <- coverage.list[[n.char]]
  
  plot.df <- data.frame(
    mean.wald  = plot.dat$mean.pval.wald,
    shape.wald = plot.dat$shape.pval.wald,
    mean.prof  = plot.dat$mean.pval.prof,
    shape.prof = plot.dat$shape.pval.prof
  ) |> 
    pivot_longer(
      cols = everything(),
      names_to = c("parameter", "method"),
      names_sep = "\\."
    ) |>
    mutate(
      # Combine method and parameter to create 4 unique panel groups
      panel = paste0(method, "_", parameter),
      # Order them so Wald row is on top, Mean/Shape go left-to-right
      panel = factor(panel, levels = c("wald_mean", "wald_shape", "prof_mean", "prof_shape"))
    )
  
  # Compute bounds for mean/shape ... maybe replace with directly uniform (ADD)
  
  # Define completely custom titles for each of the 4 panels
  wald_mean.fit <- coef(mle.est(plot.dat$mean.pval.wald))
  wald_shape.fit <- coef(mle.est(plot.dat$shape.pval.wald))
  prof_mean.fit <- coef(mle.est(plot.dat$mean.pval.prof))
  prof_shape.fit <- coef(mle.est(plot.dat$shape.pval.prof))
  
  panel_titles <- c(
    "wald_mean"  = paste0("Wald Mean: Beta(", signif(wald_mean.fit[["mean"]],4), ", ", signif(wald_mean.fit[["shape"]],4), ")" ),
    "wald_shape" = paste0("Wald Shape: Beta(", signif(wald_shape.fit[["mean"]],4), ", ", signif(wald_shape.fit[["shape"]],4), ")" ),
    "prof_mean"  = paste0("Profile Mean: Beta(", signif(prof_mean.fit[["mean"]],4), ", ", signif(prof_mean.fit[["shape"]],4), ")" ),
    "prof_shape" = paste0("Profile Shape: Beta(", signif(prof_shape.fit[["mean"]],4), ", ", signif(prof_shape.fit[["shape"]],4), ")" )
    )
  
  metaPiano <- ggplot(plot.df, aes(x = value)) +
    geom_histogram(breaks = seq(0,1,length.out=41)) + 
    facet_wrap(
      ~ panel, 
      nrow = 2,
      labeller = as_labeller(panel_titles)
    ) +
    labs(
      title = paste0("P-Value Calibration / Piano Plots (n = ", n, ")"),
      x = "Left-Tail CDF P-Value",
      y = "Frequency"
    ) +
    xlim(0, 1) +
    theme_minimal() +
    theme(
      strip.text = element_text(face = "bold", size = 9),
      plot.title = element_text(face = "bold", hjust = 0.5)
    )
  
  print(metaPiano)
}

metaPianos(coverage.list, 2000)
