library(shellpipes)


# =======================

# DATA GENERATION

generatePoisson <- function(lambdas, dat = NULL,  numSims = 1e4, 
                            numReps = 1, testv = "poisson.test", seed=1){
  set.seed(seed)
  ## Make these as parallel as possible, and explain why they even have to exist with a brief comment
  if (testv == "poisson.test"){
    poisson.data <- 
      purrr::map(lambdas, function(lambda){
        if (is.null(dat)){
          dat <- rpois(numSims, lambda)
        }
        df <- data.frame(multPois(dat, lambda, testv), lambda)
      }) |> list_rbind()
    
  } else if (testv == "wald.intercept"){
    poisson.data <- 
      purrr::map(lambdas, function(lambda){
        dat <- purrr::map(c(1:numSims), function(dummy){
          rpois(numReps, lambda)}) %>% rbind()
        
        df <- data.frame(multPois(dat, lambda, testv), lambda)
      }) %>% list_rbind()
  }
  return(poisson.data)
}

# -----

multPois <- function(dat, lambda0, testv){
  if (testv == "poisson.test"){
    df <- 
      purrr::map(dat, function(d){
        bt <- stats::poisson.test(x = d, T = 1, r = lambda0, alternative = "less")
        gt <- stats::poisson.test(x = d, T = 1, r = lambda0, alternative = "greater")
        gt.adj <- ppois(q = d, lambda = lambda0, lower.tail = FALSE)
        bt.adj <- ppois(q = d, lambda = lambda0, lower.tail = TRUE) - dpois(d, lambda0) 
        cp <- bt$p.value
        rp <- bt$p.value + runif(1) * (1-gt$p.value - cp)
        ci <- stats::poisson.test(d, T = 1, r = lambda0, alternative = "two.sided")
        return(data.frame(est = d,
                          left.exact = cp,
                          right.exact = gt$p.value,
                          left.adj = bt.adj,
                          right.adj = gt.adj,
                          p = rp,
                          lower = ci$conf.int[1],
                          upper = ci$conf.int[2],
                          test = testv)
        )
      }) %>%
      list_rbind()
    
  } else if (testv == "wald.intercept"){
    df <-
      purrr::map(dat, function(d){
        g <- glm(d ~ 1 + offset(log(0*d+lambda0)), family = poisson())
        est <- coef(g)[[1]] # N(0, s2) ... mean intercept 0 
        se <- sqrt(vcov(g))[[1,1]]
        w <- est/se # N(0,1) ... since divide by se -> s2
        bt <- pnorm(w, mean = 0, sd = 1, lower.tail = TRUE)
        gt <- pnorm(w, mean = 0, sd = 1, lower.tail = FALSE)
        ci.lower <- est - 1.96*se # double-check
        ci.upper <- est + 1.96*se # double-check
        # rp <- UNKNOWN IF POSSIBLE
        return(data.frame(pois.mean = mean(d),
                          int.est = est,
                          int.se = se,
                          z.value = w,
                          left = bt, # not really needed
                          right = gt, # not really needed
                          p = bt, 
                          lower = ci.lower,
                          upper = ci.upper,
                          testv = testv)
        )
      }) %>% 
      list_rbind()
  }
}


# =======================

# TRIANGLE PLOT PREPARATION


triangleData.poisson <- function(dat, plot = "left", point.mass = TRUE, 
                                 fuzz.x = FALSE, use.left = TRUE, testv = "poisson.test"){
  
  # for 1 lambda
  if (fuzz.x | testv != "poisson.test"){ # any continuous approximation
    thresholds <- sort(unique(round(c(0, 0.5, 1, dat$p),9)))
    
  } else if (use.left){ # left tail data
    if (plot == "left"){
      thresholds <- sort(unique(c(0, 0.5, 1, dat$left.exact)))
      dat$p <- dat$left.exact
    } else{
      thresholds <- sort(unique(c(0, 0.5, 1, dat$right.adj)))
      dat$p <- dat$right.adj
    }
    
  } else if (!use.left){ # right tail data
    if (plot == "left"){
      thresholds <- sort(unique(c(0, 0.5, 1, dat$right.exact)))
      dat$p <- dat$right.exact
    } else {
      thresholds <- sort(unique(c(0, 0.5, 1, dat$left.adj)))
      dat$p <- dat$left.adj
    }
  }
  
  if (plot == "left"){
    left.data <- triangleData.poisson.aux(dat, thresholds, point.mass, plot = "left")
    return(left.data)
  } else{
    right.data <- triangleData.poisson.aux(dat, thresholds, point.mass, plot = "right")
    return(right.data)
  }
}

# -----

triangleData.poisson.aux <- function(dat, thresholds, point.mass = TRUE, plot = "left"){
  
  if (point.mass){
    cmp <- if (plot == "left") `<=` else `>=`
    
    data.aux <-
      dat %>% split(dat$lambda) %>%
      purrr::map(function(dat){
        purrr::map(thresholds, function(f){
          sum <- sum(cmp(dat$p, f)) # HERE
          setNames(data.frame(sum), as.character(f))
        }) %>% list_cbind()
      }) %>% list_rbind()
  } else {
    cmp <- if (plot == "left") `<` else `>`
    
    data.aux <-
      dat %>% split(dat$lambda) %>%
      purrr::map(function(dat){
        purrr::map(thresholds, function(f){
          sum <- sum(cmp(dat$p, f)) # HERE
          setNames(data.frame(sum), as.character(f))
        }) %>% list_cbind()
      }) %>% list_rbind()
  }
  
  data.aux <- data.aux %>%
    mutate(lambda = unique(dat$lambda)) %>%
    select(lambda, everything())
  
  return(data.aux)
}

# =======================

# MAKE TRIANGLE PLOTS


trianglePlot.poisson <- function(dat.left, dat.right, testv, plot = "left", 
                                 fuzz.x = TRUE, point.mass = TRUE, add.points = FALSE){
  
  prep <- function(dat){
    dat %>% 
      melt(id.vars = "lambda", variable.name = "p", value.name = "count") %>%
      mutate(p = as.double(as.character(p)),
             cum.prop = count/max(count)
      )
  }
  
  if (plot == "left"){
    plot.data.left <- prep(dat.left)
    triangle.plot <- trianglePlot.poisson.aux(plot.data.left, NULL, 
                                              testv, plot = "left", fuzz.x, 
                                              point.mass, add.points)
    
  } else if (plot == "right"){
    plot.data.right <- prep(dat.right)
    triangle.plot <- trianglePlot.poisson.aux(NULL, plot.data.right, 
                                              testv, plot = "right", fuzz.x, 
                                              point.mass, add.points)
    
  } else if (plot == "full"){  
    plot.data.left <- prep(dat.left)
    plot.data.right <- prep(dat.right)
    
    # take the correct half of the data corresponding to triangle side
    plot.data.left <- plot.data.left %>% filter(p <= 0.5) 
    plot.data.right <- plot.data.right %>% filter(p >= 0.5)
    
    triangle.plot <- trianglePlot.poisson.aux(plot.data.left, plot.data.right, 
                                              testv, plot = "full", fuzz.x, 
                                              point.mass, add.points)
  }
  
  return(triangle.plot)
}

# -----

trianglePlot.poisson.aux <- function(dat.left, dat.right, testv, plot = "left", 
                                     fuzz.x = TRUE, point.mass = TRUE, add.points = FALSE){
  
  inclusive <- point.mass | testv != "poisson.test"
  step.left <- ifelse(inclusive | fuzz.x, "hv", "vh")
  step.right <- ifelse(inclusive| fuzz.x, "vh", "hv")
  
  if (plot == "left"){
    triangle.plot <- ggplot(data = dat.left) + 
      geom_step(aes(x = p, y = cum.prop), direction = step.left) +
      geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
      facet_grid(~lambda) +
      theme_classic() +
      theme(panel.spacing.x = unit(2, "lines"), 
            axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
      labs(x = "nominal p value",
           y = "cumulative prop. more extreme",
           title = "Triangle Plot: Poisson")
    if (add.points){
      triangle.plot <- triangle.plot + geom_point(aes(x = p, y = cum.prop))
    }
    
  } else if (plot == "right"){
    triangle.plot <- ggplot(data = dat.right) + 
      geom_step(aes(x = p, y = cum.prop), direction = step.right) +
      geom_abline(intercept = 1, slope = -1, linetype = "dashed") +
      facet_grid(~lambda) +
      theme_classic() +
      theme(panel.spacing.x = unit(2, "lines"), 
            axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
      labs(x = "nominal p value",
           y = "cumulative prop. more extreme",
           title = "Triangle Plot: Poisson") 
    
    if (add.points){
      triangle.plot <- triangle.plot + geom_point(aes(x = p, y = cum.prop))
    } 
    
  } else if (plot == "full"){
    # properly connects the left and right triangles ... facet on lambda
    joinData <- data.frame(lambda = rep(unique(dat.left$lambda), 2),
                           p = rep(0.5, 2*length(unique(dat.left$lambda))),
                           cum.prop = c(dat.left %>% filter(p == 0.5) %>% pull(cum.prop),
                                        dat.right %>% filter(p == 0.5) %>% pull(cum.prop)))
    
    # triangle reference lines ... facet on lambda
    referenceData.left <- data.frame(lambda = sort(rep(unique(dat.left$lambda), 2)),
                                     p = rep(c(0, 0.5), length(unique(dat.left$lambda))),
                                     cum.prop = rep(c(0, 0.5), length(unique(dat.left$lambda))))
    
    referenceData.right <- data.frame(lambda = sort(rep(unique(dat.left$lambda), 2)),
                                      p = rep(c(0.5, 1), length(unique(dat.left$lambda))), 
                                      cum.prop = rep(c(0.5, 0), length(unique(dat.left$lambda)))) 
    
    # plot the full triangle
    triangle.plot <- 
      ggplot() +
      xlim(c(0,1)) + 
      ylim(c(0,1)) +
      geom_step(data = dat.left, aes(x = p, y = cum.prop), direction = step.left) + 
      geom_step(data = dat.right, aes(x = p, y = cum.prop), direction = step.right) +
      geom_line(data = joinData, aes(x = p, y = cum.prop), col = "red") + 
      geom_line(data = referenceData.left, aes(x = p, y = cum.prop), col = "black", linetype = "dashed") +
      geom_line(data = referenceData.right, aes(x = p, y = cum.prop), col = "black", linetype = "dashed") +
      facet_grid(~lambda) + 
      theme_classic() +
      theme(panel.spacing.x = unit(2, "lines"),
            plot.title = element_text(size = 15))+
      labs(x="nominal p-value", 
           y = "cum. prop. more extreme",
           title = sprintf("Triangle Plots: Poisson"))
    
    if (add.points){
      triangle.plot <- triangle.plot + 
        geom_point(data = dat.left, aes(x = p, y = cum.prop)) + 
        geom_point(data = dat.right, aes(x = p, y = cum.prop))
    }
    
  }
  
  return(triangle.plot)
}

# ==========

# WRAPPER

trianglePlot <- function(lambdas, numSims = 1e4, numReps = 1, # defaulting to a discrete fuzzed approach
                         family = "poisson", testv = "poisson.test",
                         plot = "full", dat = NULL, use.left = TRUE, fuzz.x = TRUE,
                         point.mass = TRUE, add.points = FALSE, add.checkplot = FALSE){
  
  if (family == "poisson"){ # intention of adding more families
    # generate data
    dat <- generatePoisson(lambdas, dat, numSims, numReps, testv)
    
    if (add.checkplot){
      checkplot <- checkPlot(dat, facets = length(unique(dat$lambda))) + facet_grid(~lambda)
      print(checkplot)
    }
    
    # prepare plotting data accordingly & plot it
    if (plot == "left" | plot == "full"){
      triangle.data.left <- triangleData.poisson(dat, plot = "left", point.mass, 
                                                 fuzz.x, use.left, testv)
      
      left.plot <- trianglePlot.poisson(triangle.data.left, NULL, testv, 
                                        plot = "left", fuzz.x, point.mass,
                                        add.points)  
      print(left.plot)
    }
    
    if (plot == "right" | plot == "full"){
      triangle.data.right <- triangleData.poisson(dat, plot = "right", point.mass, 
                                                  fuzz.x, use.left, testv)
      
      right.plot <- trianglePlot.poisson(NULL, triangle.data.right, testv, 
                                         plot = "right", fuzz.x, point.mass, 
                                         add.points)
      print(right.plot)
    }
    
    if (plot == "full"){
      full.plot <- trianglePlot.poisson(triangle.data.left, triangle.data.right, 
                                        testv, plot, fuzz.x, point.mass, add.points)
      print(full.plot)
    }
  }
}

saveEnvironment()
