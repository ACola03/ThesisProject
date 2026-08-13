# ===== sim_ss 

# Description: 
# - Compute mean and sd of CvM statistics using the expectation and midpoint

# Arguments:
# 1. lambda: controls the balance between the two shape parameters
# 2. samp: number of beta samples to simulate
# 3. reps: how many simulations are conducted to compute the avg & sd
# 4. shape: controls whether we examine u/n or a shifted shape

sim_ss <- function(lambda, samp, reps, shape){
  # some setup
  i <- 1:samp
  c2 <- (3-lambda)/2
  
  # two targets
  midpointTarget <- (2*i-1)/(2*samp)
  expectTarget <- i/(samp+1)
  
  # compute a CvM or AD statistic, using a specified target
  ss <- function(d, t, ad.weight){
    sq.diffs <- (d-t)^2
    weights <- 1/(d*(1-d))
    stat <- ifelse(ad.weight, sum(weights*sq.diffs), sum(sq.diffs))
    return(stat)
  }
  
  # store statistics
  cvm.midpoint <- numeric(reps); cvm.expect <- numeric(reps)
  ad.midpoint <- numeric(reps); ad.expect <- numeric(reps)
  
  # controls the beta shape (u/n or shift)
  c1 <- ifelse(shape %in% c("u","n"), c2, 1)
  
  # compute the collection of statistics
  for(iter in 1:reps){
    dist <- sort(rbeta(samp, c1, c2))
    cvm.midpoint[iter] <- ss(dist, midpointTarget, FALSE)
    cvm.expect[iter] <- ss(dist, expectTarget, FALSE)
    ad.midpoint[iter] <- ss(dist, midpointTarget, TRUE)
    ad.expect[iter] <- ss(dist, expectTarget, TRUE)
  }
  
  
  # returned object of stats across lambda
  out <- c(lambda, 
           mean_cvm.midpoint = mean(cvm.midpoint), 
           sig_cvm.midpoint=sd(cvm.midpoint),
           mean_cvm.expect = mean(cvm.expect), 
           sig_cvm.expect=sd(cvm.expect),
           mean_ad.midpoint = mean(ad.midpoint), 
           sig_ad.midpoint=sd(ad.midpoint),
           mean_ad.expect = mean(ad.expect), 
           sig_ad.expect=sd(ad.expect)
           )
  
  return(out)
}

# ===== sim_ss_over_lambda

# Description:
# - Compute the statistics described above across a range of shape parameters 
#   controlled by lambda.

sim_ss_over_lambda <- function(lambda, samp, reps, shape){
  results <- lapply(lambda, function(lam){
    sim_ss(lam, samp, reps, shape)
    })
  
  df <- as.data.frame(do.call(rbind, results))
  names(df)[[1]] <- "lambda"
  
  title.shape <- ifelse(shape %in% c("u","n"), "Symmetric", "Skewed")
  
  # plot cvm: midpoint vs expected; showing where the minimum occurs
  min_cvm.midpoint <-df$lambda[which.min(df$mean_cvm.midpoint)]
  min_cvm.expect <-df$lambda[which.min(df$mean_cvm.expect)]
  cvm.plot <- ggplot(df) + 
    aes(lambda) +
    geom_line(aes(y=mean_cvm.midpoint, color = "Midpoint")) + 
    geom_point(aes(y=mean_cvm.midpoint, color = "Midpoint")) +
    geom_vline(xintercept = min_cvm.midpoint, linetype = "dashed") +
    geom_line(aes(y=mean_cvm.expect, color = "Expectation")) +
    geom_point(aes(y=mean_cvm.expect, color = "Expectation")) +
    geom_vline(xintercept = min_cvm.expect, color = "blue", linetype = "dashed") +
    scale_color_manual(
      name = "Target Type", 
      values = c("Midpoint" = "black", "Expectation" = "blue")
    ) +
    labs(title = paste0("Identifying the Uniform Distribution: CvM Statistic", "\n",
                        "Beta Perturbation Shape: ", title.shape, " (n=",samp, ")"), 
         y = "CvM Statistic")
  print(cvm.plot)
  
  # plot ad: midpoint vs expected; showing where the minimum occurs
  min_ad.midpoint <-df$lambda[which.min(df$mean_ad.midpoint)]
  min_ad.expect <-df$lambda[which.min(df$mean_ad.expect)]
  ad.plot <- ggplot(df) + 
    aes(lambda) +
    scale_y_log10() +
    geom_line(aes(y=mean_ad.midpoint, color = "Midpoint")) +
    geom_point(aes(y=mean_ad.midpoint, color = "Midpoint")) +
    geom_vline(xintercept = min_ad.midpoint, linetype = "dashed") +
    geom_line(aes(y=mean_ad.expect, color = "Expectation")) +
    geom_point(aes(y=mean_ad.expect, color = "Expectation")) +
    geom_vline(xintercept = min_ad.expect, color = "blue", linetype = "dashed") +
    scale_color_manual(
      name = "Target Type", 
      values = c("Midpoint" = "black", "Expectation" = "blue")
    ) +
    labs(title = paste0("Identifying the Uniform Distribution: AD Statistic", "\n",
                        "Beta Perturbation Shape: ", title.shape, " (n=",samp, ")"), 
         y = "AD Statistic")
  print(ad.plot)
  
  return(df)
}
