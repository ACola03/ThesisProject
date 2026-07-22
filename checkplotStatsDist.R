library(dplyr)
library(purrr)
library(reshape2)
library(tidyr)

# functions used by powerAnalysis.R to compute checkplot statistics

# ==========
# Bin Variance Statistic

binvar_dist <- function(N, n, shapes1, shapes2, binwidth = 0.025){
  
  n.pairs <- length(shapes1)
  allStat_list <- vector("list", n.pairs)
  
  for (pair in seq_len(n.pairs)){
    shape1 <- shapes1[pair]
    shape2 <- shapes2[pair]
    
    beta_samples_matrix <- matrix(rbeta(n * N, shape1, shape2), nrow = n, ncol = N)
    allStat_list[[pair]] <- apply(beta_samples_matrix, 2, function(x) floor(x/binwidth) |> table() |> var())
  }
  
  allStat.df <- do.call(cbind, allStat_list)
  colnames(allStat.df) <- paste0("Beta(", shapes1, ",", shapes2, ")")
  
  return(allStat.df)
}

# ==========
# Spacing Statistic

space_dist <- function(N, n, shapes1, shapes2){
  
  n.pairs <- length(shapes1)
  allStat_list <- vector("list", n.pairs) # Pre-allocate list for memory efficiency
  
  for (pair in seq_len(n.pairs)){
    shape1 <- shapes1[pair]
    shape2 <- shapes2[pair]
    
    beta_samples_matrix <- rbind(matrix(rbeta(n * N, shape1, shape2), nrow = n, ncol = N), rep(0,N), rep(1,N))
    allStat_list[[pair]] <- apply(beta_samples_matrix, 2, function(x) var((n+1)*diff(sort(x))))
  }
  
  allStat.df <- do.call(cbind, allStat_list)
  colnames(allStat.df) <- paste0("Beta(", shapes1, ",", shapes2, ")")
  
  return(allStat.df)
}

# ==========
# Location Statistic

location_dist_aux <- function(data_matrix, beta.exp, beta.var, reweight, mean){
  pvec_matrix <- apply(data_matrix, 2, sort)
  summary.func <- ifelse(mean, `colMeans`, `colSums`)
  
  if (reweight){
    ssq.err <- summary.func(((pvec_matrix - beta.exp)^2) / beta.var)
  } else{
    ssq.err <- summary.func((pvec_matrix - beta.exp)^2)
  }
  
  return(ssq.err)
}

location_dist <- function(N, n, shapes1, shapes2, reweight, mean){
  beta.exp <- (1:n)/(n+1)
  beta.var <- ((1:n)*(n - (1:n) + 1)) / ((n+1)^2*(n+2))
  
  n.pairs <- length(shapes1)
  allStat_list <- vector("list", n.pairs) # Pre-allocate list for memory efficiency
  
  for (pair in seq_len(n.pairs)){
    shape1 <- shapes1[pair]
    shape2 <- shapes2[pair]
    
    beta_samples_matrix <- matrix(rbeta(n * N, shape1, shape2), nrow = n, ncol = N)
    
    # calculate statistics for all N simulations simultaneously
    allStat_list[[pair]] <- location_dist_aux(beta_samples_matrix, beta.exp, beta.var, reweight, mean)
  }
  
  allStat.df <- do.call(cbind, allStat_list)
  colnames(allStat.df) <- paste0("Beta(", shapes1, ",", shapes2, ")")
  
  return(allStat.df)
}

# ==========
# DTS Statistic

dts_dist <- function(N, n, shapes1, shapes2){
  
  n.pairs <- length(shapes1)
  allStat_list <- vector("list", n.pairs) 
  
  for (pair in seq_len(n.pairs)){
    shape1 <- shapes1[pair]
    shape2 <- shapes2[pair]
    
    beta_samples_matrix <- matrix(rbeta(n * N, shape1, shape2), nrow = n, ncol = N)
    unif_samples_matrix <- matrix(runif(n * N), nrow = n, ncol = N)
    
    allStat_list[[pair]] <- sapply(1:N, function(i) dts_stat(beta_samples_matrix[, i], unif_samples_matrix[, i]))
  }
  
  allStat.df <- do.call(cbind, allStat_list)
  colnames(allStat.df) <- paste0("Beta(", shapes1, ",", shapes2, ")")
  
  return(allStat.df)
}