library(shellpipes)
rpcall("trianglePlot_Poisson.Rout trianglePlot_Poisson.pipestar trianglePlot_Poisson.R")

# =======================

# DATA GENERATION

# Changes:
# i) Made the cases as parallel as possible with brief descriptions
# ii) Using the 'proper' pipes (|>)

# Description: For each lambda, create a dataframe of numSims poisson responses
# i) poisson.test is the simplest case where counts don't require the consideration below
# ii) since GLMs improve with a larger sample size per experiment, we enable this with numReps
# iii) users can provide their own data for specific testing, given they generate it correctly
#      by providing their own data, they can evaluate it against candidate lambdas, which is useful
# iv) each map() call returns a dataframe per lambda, which are then bound by row

generatePoisson <- function(lambdas, dat = NULL,  numSims = 1e4, 
                            numReps = 1, testv = "poisson.test"){
  
  poisson.data <- purrr::map(lambdas, function(lambda){
    if (testv == "poisson.test" & is.null(dat)){
      dat <- rpois(numSims, lambda)
      
    } else if (testv == "wald.intercept" & is.null(dat)){
      dat <- purrr::map(c(1:numSims), function(dummy){
        rpois(numReps, lambda)
      }) |> rbind()
    }
    
    df <- data.frame(multPois(dat, lambda, testv), lambda)}) |> list_rbind()
  
  return(poisson.data)
}

# -----

# Changes:
# i) used ppois instead of poisson.test
# ii) added proper two-tailed p-values
# iii) an attempt at fuzzing wald.intercept p-values 
# iv) deciding to keep the confidence interval for slugPlots
# v) improved variable names (as requested) -> left.strict

# Description: calculate necessary quantities for p-value diagnostics
# i) the cases are separate due to different distributions (pois vs norm)
# ii) the wald.intercept case can take multiple poisson counts in a single experiment

multPois <- function(dat, lambda0, testv){
  if (testv == "poisson.test"){
    df <- 
      purrr::map(dat, function(d){
        bt <- ppois(q = d + 0.5, lambda = lambda0, lower.tail = TRUE)
        gt <- 1 - ppois(q = d - 0.5, lambda = lambda0, lower.tail = TRUE)
        bt.strict <- ppois(q = d, lambda = lambda0, lower.tail = TRUE) - dpois(d, lambda0)
        gt.strict <- ppois(q = d, lambda = lambda0, lower.tail = FALSE)
        two.tail <- min(2*min(bt,gt), 1) # instead of poisson.test exact two-tail
        rp <- bt + runif(1) * (1-gt - bt)
        r.two.tail <- min(2*min(rp,1-rp), 1)
        ci <- stats::poisson.test(d, T = 1, r = lambda0, alternative = "two.sided")
        return(data.frame(est = d,
                          left.exact = bt,
                          right.exact = gt,
                          two.exact = two.tail,
                          left.strict = bt.strict,
                          right.strict = gt.strict,
                          p = rp,
                          two.p = r.two.tail,
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
        two.tail <- min(2*min(bt,gt), 1) 
        ci.lower <- est - 1.96*se
        ci.upper <- est + 1.96*se 
        # rp <- I THINK THIS IS POSSIBLE ... COMING NEXT
        return(data.frame(pois.mean = mean(d),
                          int.est = est, # int meaning intercept
                          int.se = se,
                          z.value = w,
                          left = bt, # not really needed
                          right = gt, # not really needed
                          p = bt, # this currently isn't fuzzed but is treated as such
                          two.p = two.tail,
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

# Changes / Comments:
# i) use.left is to determine which is the primary tail (left or right)
#    this will determine which tail is measured in the left (primary) plot
#    however, I don't really see much use for this, which we can discuss

# ii) point.mass controls whether we calculate >= or >
#     I recommend setting this to TRUE since it is the most logical choice
#     however, the plots will look identical regardless of the input since 
#     the step drawing choice (hv or vh) is chosen appropriately 

# iii) upon construction, any 'p' or 'two.p' column represents fuzzed or continuous p-values
#     however, upon writing this, I realize that the wald.intercept 'p' represents a 
#     continuous p-value but not a fuzzed p-value so if I fuzz it, I need to change names
#     ... actually, I'll have 'p' represent the fuzz and 'left.exact' represent the non-fuzz
#     or I'll use 'rp' instead of 'p' to be super clear that it is randomized .. keep convention

# iv) I like using 'else if' instead of 'else' so it's clear what the other condition is

# v) technically, for the wald case, they aren't fuzzed and can seem discrete, but until 
#    they are fuzzed we assume they have such behaviour so the fuzz condition may not 
#    align with the choice of using continuous approximation so logicals may need to be changed 
#    (in reference to fuzz.x | testv != ... )

# Description: prepare data for triangle plot creation
# i) the acceptable 'plot' input is 'one' or 'two' to denote which tail types we want to plot
#    'one' will plot one-sided triangles, while 'two' will plot the combined triangles
#    the choice determines how the data will be pre-processed

# ii) eventually, the 'p' column is modified as desired so I can avoid considering
#     multiple cases in the auxiliary function

triangleData.poisson <- function(dat, plot = "one", point.mass = TRUE, 
                                 fuzz.x = TRUE, use.left = TRUE, testv = "poisson.test"){
  
  if (fuzz.x | testv != "poisson.test"){ # any continuous approximation
    
    p.col <- if (plot == "one") dat$p else dat$two.p
    thresholds <- sort(unique(round(c(0, 0.5, 1, p.col), 6)))
    dat$p <- p.col
    
    left.data <- triangleData.poisson.aux(dat, thresholds, point.mass, plot="left")
    right.data <- triangleData.poisson.aux(dat, thresholds, point.mass, plot="right")
    
  } else if (plot == "one"){
    thresholds <- sort(unique(c(0, 0.5, 1, dat$left.exact)))
    dat$p <- dat$left.exact
    left.data <- triangleData.poisson.aux(dat, thresholds, point.mass, plot = "left")
    
    thresholds <- sort(unique(c(0, 0.5, 1, dat$right.strict)))
    dat$p <- dat$right.strict
    right.data <- triangleData.poisson.aux(dat, thresholds, point.mass, plot = "right")
    
  } else if (plot == "two"){
    thresholds <- sort(unique(c(0, 0.5, 1, dat$two.exact)))
    dat$p <- dat$two.exact
    
    left.data <- triangleData.poisson.aux(dat, thresholds, point.mass, plot = "left")
    right.data <- triangleData.poisson.aux(dat, thresholds, point.mass, plot = "right")
  }
  
  return(list("left.data" = left.data, "right.data" = right.data))
  
}

# -----

triangleData.poisson.aux <- function(dat, thresholds, point.mass = TRUE, plot = "left"){
  
  if (point.mass){
    cmp <- if (plot == "left") `<=` else `>=`
  } else {
    cmp <- if (plot == "left") `<` else `>`
  }
  
  data.aux <-
    dat %>% split(dat$lambda) %>%
    purrr::map(function(dat){
      purrr::map(thresholds, function(f){
        sum <- sum(cmp(dat$p, f))
        setNames(data.frame(sum), as.character(f))
      }) %>% list_cbind()
    }) %>% list_rbind()
  
  data.aux <- data.aux %>%
    mutate(lambda = unique(dat$lambda)) %>%
    select(lambda, everything())
  
  return(data.aux)
}

# =======================

# MAKE TRIANGLE PLOTS

# Comments:
# i) Cannot construct both 'one' and 'two' plots at the same time since 
#    'dat' is specific to 'one' or 'two'. Will be handled by wrapper.

trianglePlot.poisson <- function(dat, testv, plot = "one", 
                                 fuzz.x = TRUE, point.mass = TRUE, add.points = FALSE){
  
  dat.left <- dat$left.data
  dat.right <- dat$right.data
  
  prep <- function(dat){
    dat %>% 
      melt(id.vars = "lambda", variable.name = "p", value.name = "count") %>%
      mutate(p = as.double(as.character(p)),
             cum.prop = count/max(count)
      )
  }
  
  if (plot == "one"){
    plot.data.left <- prep(dat.left)
    triangle.plot.left <- trianglePlot.poisson.aux(plot.data.left, NULL, 
                                                   testv, plot = "left", fuzz.x, 
                                                   point.mass, add.points)
    
    plot.data.right <- prep(dat.right)
    triangle.plot.right <- trianglePlot.poisson.aux(NULL, plot.data.right, 
                                                    testv, plot = "right", fuzz.x, 
                                                    point.mass, add.points)
    
    triangle.plot.two <- NULL
    
  } else if (plot == "two"){  
    plot.data.left <- prep(dat.left)
    plot.data.right <- prep(dat.right)
    
    # take the correct half of the data corresponding to triangle side
    plot.data.left <- plot.data.left %>% filter(p <= 0.5) 
    plot.data.right <- plot.data.right %>% filter(p >= 0.5)
    
    triangle.plot.two <- trianglePlot.poisson.aux(plot.data.left, plot.data.right, 
                                                  testv, plot = "two", fuzz.x, 
                                                  point.mass, add.points)
    
    triangle.plot.left <- NULL
    triangle.plot.right <- NULL
  }
  
  return(list("left.plot" = triangle.plot.left, 
              "right.plot" = triangle.plot.right,
              "two.plot" = triangle.plot.two))
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
    
  } else if (plot == "two"){
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
      geom_line(data = joinData, aes(x = p, y = cum.prop), col = "red", linetype = "dashed") + 
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
                         plot = "both", dat = NULL, use.left = TRUE, fuzz.x = TRUE,
                         point.mass = TRUE, add.points = FALSE, add.checkplot = FALSE){
  
  if (family == "poisson"){ # intention of adding more families
    # generate data
    dat <- generatePoisson(lambdas, dat, numSims, numReps, testv)
    
    if (add.checkplot){
      checkplot <- checkPlot(dat, facets = length(unique(dat$lambda))) + facet_grid(~lambda)
      print(checkplot)
      #hist(dat$two.exact)
      hist(dat$two.p)
    }
    
    # prepare plotting data & plot it
    if (plot == "one" | plot == "both"){
      triangle.data <- triangleData.poisson(dat, plot = "one", point.mass, 
                                            fuzz.x, use.left, testv)
      
      triangle.plots <- trianglePlot.poisson(triangle.data, testv, 
                                             plot = "one", fuzz.x, point.mass,
                                             add.points)  
      
      left.plot <- triangle.plots$left.plot
      right.plot <- triangle.plots$right.plot
      
      print(left.plot)
      print(right.plot)
      
    }
    
    if (plot == "two" | plot == "both"){
      triangle.data <- triangleData.poisson(dat, plot = "two", point.mass, 
                                            fuzz.x, use.left, testv)
      
      triangle.plots <- trianglePlot.poisson(triangle.data, testv, 
                                             plot = "two", fuzz.x, point.mass,
                                             add.points)  
      
      two.plot <- triangle.plots$two.plot
      print(two.plot)
      
    }
  }
}

saveEnvironment()
