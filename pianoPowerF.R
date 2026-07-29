# ----- required

#library(tidyr)
#library(ggplot2)
#library(dplyr)
#library(patchwork)

# ----- Compute Kolmogorov-Smirnov

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

# ----- Compute Anderson-Darling

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

# ----- Compute Cramer-von Mises

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

# ----- Simulate an empirical distribution for 'stat', given shape parameters


shapes1 <- c(1, 0.90)
shapes2 <- c(1, 0.90)

zhang.dist <- function(N, n, shapes1, shapes2, stat, use.zhang){
  
  n.pairs <- length(shapes1)
  allStat_list <- vector("list", n.pairs) 
  
  statFunc <- switch(stat,
                     "KS" = `zhang.KS`,
                     "AD" = `zhang.AD`,
                     "CVM" = `zhang.CVM`)
  
  for (pair in seq_len(n.pairs)){
    shape1 <- shapes1[pair]
    shape2 <- shapes2[pair]
    
    beta_samples_matrix <- matrix(rbeta(n * N, shape1, shape2), nrow = n, ncol = N)
    allStat_list[[pair]] <- apply(beta_samples_matrix, 2, statFunc, use.zhang = use.zhang)
  }
  
  allStat.df <- do.call(cbind, allStat_list)
  colnames(allStat.df) <- paste0("Beta(", shapes1, ",", shapes2, ")")
  
  return(allStat.df)
}

# ----- Compute power for certain alpha

zhang.computePower <- function(data, alpha){
  nullDist <- data[,1]
  criticalValue <- quantile(nullDist, probs = alpha)
  powerVec <- numeric()
  
  for (col in 2:ncol(data)){
    powerVec[[col-1]] <- mean(data[,col] > criticalValue)
  }
  
  return(powerVec)
}

# ----- Power Analysis Wrapper: compute power for all shapes and alphas, for every statistic

zhang.powerAnalysis <- function(N, ns, shapes1, shapes2, alphas, add.plot){
  
  n.alphas <- length(alphas)
  n.betas <- length(shapes1) - 1
  
  # setup power list
  statPower_alpha <- lapply(vector("list", n.alphas), function(x) {
    vector("list", n.betas)
  })
  
  names(statPower_alpha) <- paste0(alphas)
  
  statPower_alpha <- lapply(statPower_alpha, function(sublist) {
    setNames(sublist, paste0("Beta", "(", shapes1[-1], ",", shapes2[-1], ")"))
  })
  
  # populate power list
  for (n in ns){
    
    ks.dist <- zhang.dist(N, n, shapes1, shapes2, "KS", FALSE)
    zks.dist <- zhang.dist(N, n, shapes1, shapes2, "KS", TRUE)
    ad.dist <- zhang.dist(N, n, shapes1, shapes2, "AD", FALSE)
    zad.dist <- zhang.dist(N, n, shapes1, shapes2, "AD", TRUE)
    cvm.dist <- zhang.dist(N, n, shapes1, shapes2, "CVM", FALSE)
    zcvm.dist <- zhang.dist(N, n, shapes1, shapes2, "CVM", TRUE)
    
    for (alpha in alphas){
      alpha.key <- as.character(alpha)
      ks.pow <- zhang.computePower(ks.dist, 1-alpha)
      zks.pow <- zhang.computePower(zks.dist, 1-alpha)
      ad.pow <- zhang.computePower(ad.dist, 1-alpha)
      zad.pow <- zhang.computePower(zad.dist, 1-alpha)
      cvm.pow <- zhang.computePower(cvm.dist, 1-alpha)
      zcvm.pow <- zhang.computePower(zcvm.dist, 1-alpha)
      
      for (shapePair in 1:n.betas){
        betaDist <- paste0("Beta", "(", shapes1[shapePair+1], ",", shapes2[shapePair+1], ")")
        
        row <- c("betaDist" = betaDist, "n" = n, 
                 "ks" = ks.pow[shapePair], "zks" = zks.pow[shapePair],
                 "ad" = ad.pow[shapePair], "zad" = zad.pow[shapePair],
                 "cvm" = cvm.pow[shapePair], "zcvm" = zcvm.pow[shapePair])
        
        statPower_alpha[[alpha.key]][[shapePair]] <- rbind(
          statPower_alpha[[alpha.key]][[shapePair]], 
          row, deparse.level = 0
        )
        
      }
    }
  }
  return(statPower_alpha)
}

# ----- Power Analysis Plot: shows power across increasing n for every shape parameter pair

powerPlot <- function(data, n.lim, alpha, show.piano = FALSE){
  alpha.key <- as.character(alpha)
  alpha.data <- data[[alpha.key]]
  
  tibble.data <- as_tibble(do.call(rbind, alpha.data))
  
  plot.data <- pivot_longer(tibble.data, cols = -c("n", "betaDist"), names_to = "stat", values_to = "power") |>
    mutate(stat = factor(stat, levels = c("ad", "cvm", "ks", "zad", "zcvm", "zks")), betaDist = factor(betaDist, levels = names(alpha.data)),
           n = as.numeric(n), power = as.numeric(power)) |>
    filter(n <= n.lim)
  
  p_power <- ggplot(data = plot.data) +
    geom_line(aes(x = n, y = power, group = stat, linetype = stat, colour = stat), linewidth = 1) +
    geom_point(aes(x = n, y = power, colour = stat)) + 
    scale_linetype_manual(values = c("solid", "solid", "solid", "dashed", "dashed", "dashed")) +
    scale_color_viridis_d(option = "turbo") +
    facet_wrap(~betaDist, nrow = 1) + 
    theme_minimal() +
    theme(strip.text = element_text(size = 12)) + 
    labs(y = bquote("Statistical Power (" * alpha == .(alpha) * ")"), x = "Sample Size (n)")
  
  if (!show.piano) {
    print(p_power)
    return(invisible(p_power))
  }
  
  unique_betas <- names(alpha.data)
  
  density_plots <- lapply(unique_betas, function(b_str) {
    params <- as.numeric(unlist(regmatches(b_str, gregexpr("[0-9.]+", b_str))))
    s1 <- params[1]
    s2 <- params[2]
    
    beta_vals <- rbeta(1e4, s1, s2)
    beta_df <- data.frame(p = beta_vals)
    
    ggplot(beta_df, aes(x = beta_vals)) +
      geom_histogram(breaks = seq(0, 1, length.out=40+1)) + 
      theme_minimal() +
      theme(
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title = element_blank(),
        plot.title = element_blank()
      ) +
      labs(title = b_str)
  })
  
  # Combine density plots horizontally using patchwork
  p_densities <- wrap_plots(density_plots, nrow = 1)
  
  # Stack power plot on top of density plots
  final_plot <- p_power / p_densities + plot_layout(heights = c(3, 1))
  
  print(final_plot)
  return(invisible(final_plot))
}

# ----- Plot Empirical Statistic Distribution: as a violin plot

statDistPlot <- function(data, stat, N, n){
  
  # visualize in boxplots
  tib <- as_tibble(data)
  plotData <- pivot_longer(tib, cols = everything(), names_to = "Beta_Type", values_to = "Stat") 
  axisText <- paste0(stat, " Statistic")
  titleText <- paste0(stat, " Statistic Comparisons: (N,n) = (", N, ",", n, ")")
  
  ggplot(plotData, aes(x = reorder(Beta_Type, Stat, mean), y = Stat, fill = Beta_Type)) +
    geom_violin(draw_quantiles = c(0.25, 0.5, 0.75)) +
    labs(x = "Beta Distribution", y = axisText, 
         title = titleText) +
    theme_minimal() +
    theme(legend.position = "none")
}

