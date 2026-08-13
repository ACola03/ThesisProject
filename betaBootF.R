# ------ PARAMETRIC BOOTSTRAP FUCNTIONS

library(bbmle)
library(MASS)

null.llk <- function(mean, shape, x) {
  a <- 2*mean / shape
  b <- 2*(1 - mean) / shape
  -sum(dbeta(x, shape1 = a, shape2 = b, log = TRUE))
}

mle.est <- function(dat){
  mle2(
    null.llk,
    start = list(mean = 0.5, shape = 1),
    data = list(x = dat),
    method = "L-BFGS-B",
    lower = c(mean = 1e-4, shape = 1e-4),
    upper = c(mean = 1 - 1e-4, shape = Inf),
    trace = FALSE
  )
}

paraBoot <- function(n.unif, boots, shape1 = 1, shape2 = 1){
  alphas <- numeric(boots); betas <- numeric(boots)
  means <- numeric(boots); shapes <- numeric(boots)
  success_count <- 0
  
  for (boot in 1:boots){
    dat <- rbeta(n.unif, shape1, shape2)
    
    # for convergence failures
    fit_result <- tryCatch({
      suppressWarnings({
        fit <- mle.est(dat)
        # check if it converged (0 means success)
        if (fit@details$convergence != 0) stop("Non-zero convergence code")
        coefs <- coef(fit)
        c(mean = unname(coefs["mean"]), shape = unname(coefs["shape"]))
      })
    }, error = function(e) {
      NULL # NULL if fails to converge
    })
    
    # if converged, store the values
    if (!is.null(fit_result)) {
      success_count <- success_count + 1
      mean <- fit_result["mean"]
      shape <- fit_result["shape"]
      
      means[boot] <- mean
      shapes[boot] <- shape
      alphas[boot] <- 2 * mean / shape
      betas[boot] <- 2 * (1 - mean) / shape
    } else {
      # temporarily store the non-converged
      means[boot] <- NA
      shapes[boot] <- NA
      alphas[boot] <- NA
      betas[boot] <- NA
    }
  }
  
  cat(sprintf("Successful fits: %d out of %d (%.1f%%)\n", 
              success_count, boots, (success_count/boots)*100))
  
  # filter out those which didn't converge
  results <- list(
    "mean" = na.omit(means), 
    "shape" = na.omit(shapes),
    "shape1" = na.omit(alphas), 
    "shape2" = na.omit(betas)
  )
  
  return(results)
}

# ----- Bad Tests





