library(shellpipes)
rpcall("trianglePlot_Poisson.Rout trianglePlot_Poisson.pipestar trianglePlot_Poisson.R")

# =======================

# DATA GENERATION

# Description: For each lambda, create a dataframe of numSims poisson responses
# i) poisson.test is the simplest case where counts don't require the consideration below
# ii) since GLMs improve with a larger sample size per experiment, we enable this with numReps
# iii) users can provide their own data for specific testing, given they generate it correctly
#      by providing their own data, they can evaluate it against candidate lambdas, which is useful
# iv) each map() call returns a dataframe per lambda, which are then bound by row

generatePoisson <- function(lambdas, dat = NULL,  numSims = 1e4, 
                            numReps = 1, testv = "poisson.test"){
  dat.null <- is.null(dat)
  
  poisson.data <- purrr::map(lambdas, function(lambda){
    if (testv == "poisson.test" & dat.null){
      dat <- rpois(numSims, lambda)
      
    } else if (testv == "wald.intercept" & dat.null){
      dat <- purrr::map(c(1:numSims), function(dummy){
        rpois(numReps, lambda)
      }) |> rbind()
    }
    
    df <- data.frame(multPois(dat, lambda, testv), lambda)}) |> list_rbind()
  
  return(poisson.data)
}

# -----

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
        rp <- bt + runif(1) * (1-gt - bt)
        ci <- stats::poisson.test(d, T = 1, r = lambda0, alternative = "two.sided")
        return(data.frame(est = d,
                          left.exact = bt,
                          right.exact = gt,
                          left.strict = bt.strict,
                          right.strict = gt.strict,
                          p = round(rp,6),
                          lower = ci$conf.int[1],
                          upper = ci$conf.int[2],
                          test = testv)
        )
      }) |>
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
        ci.lower <- est - 1.96*se
        ci.upper <- est + 1.96*se 
        # rp <- I THINK THIS IS POSSIBLE ... COMING NEXT
        return(data.frame(pois.mean = mean(d),
                          int.est = est, # int meaning intercept
                          int.se = se,
                          z.value = w,
                          left.exact = bt, 
                          right.exact = gt, 
                          p = round(bt,6), # this currently isn't fuzzed but is treated as such
                          lower = ci.lower,
                          upper = ci.upper,
                          testv = testv)
        )
      }) |> 
      list_rbind()
  }
}

# =======================

# TRIANGLE PLOT PREPARATION

# Description: prepare data for triangle plot creation
# i) the acceptable 'plot' input is 'one' or 'two' to denote which tail types we want to plot
#    'one' will plot one-sided triangles, while 'two' will plot the combined triangles
#    the choice determines how the data will be pre-processed

# ii) eventually, the 'p' column is modified as desired so I can avoid considering
#     multiple cases in the auxiliary function

triangleData.poisson <- function(dat, point.mass = TRUE, fuzz.x = FALSE, 
                                 testv = "poisson.test"){
  
  if (fuzz.x | testv != "poisson.test"){ # any continuous approximation
    thresholds <- sort(unique(c(0, 0.5, 1, dat$p)))
    left.data <- triangleData.poisson.aux(dat, thresholds, point.mass, plot = "left")
    
    dat$p <- 1-dat$p
    thresholds <- sort(unique(c(0, 0.5, 1, dat$p)))
    right.data <- triangleData.poisson.aux(dat, thresholds, point.mass, plot = "right")
    
    # maybe i can simplify this by taking thresholds from dat$p and always changing $p to what i want
    
  } else {
    thresholds <- sort(unique(c(0, 0.5, 1, dat$left.exact)))
    dat$p <- dat$left.exact
    left.data <- triangleData.poisson.aux(dat, thresholds, point.mass, plot = "left")
    
    thresholds <- sort(unique(c(0, 0.5, 1, dat$right.exact)))
    dat$p <- dat$right.exact
    right.data <- triangleData.poisson.aux(dat, thresholds, point.mass, plot = "right")
  }
  
  return(list("left.data" = left.data, "right.data" = right.data))
}

# -----

triangleData.poisson.aux <- function(dat, thresholds, point.mass = TRUE, plot = "left"){
  
  if (point.mass){
    cmp <- `<=` 
  } else {
    cmp <- `<` 
  }
  
  data.aux <-
    dat |> split(dat$lambda) |>
    purrr::map(function(dat){
      purrr::map(thresholds, function(f){
        sum <- sum(cmp(dat$p, f))
        setNames(data.frame(sum), as.character(f))
      }) |> list_cbind()
    }) |> list_rbind()
  
  data.aux <- data.aux |>
    mutate(lambda = unique(dat$lambda)) |>
    select(lambda, everything())
  
  return(data.aux)
}

# =======================

# MAKE TRIANGLE PLOTS

trianglePlot.poisson <- function(dat, testv, plot = "one", 
                                 fuzz.x = TRUE, point.mass = TRUE, add.points = FALSE){
  
  prep <- function(dat, plot = "left"){
    dat |> 
      melt(id.vars = "lambda", variable.name = "p", value.name = "count") |>
      mutate(p = as.double(as.character(p)),
             p = if (plot == "right") 1 - p else p,
             cum.prop = count/max(count)
      )
  }
  
  dat.left <- dat$left.data
  dat.right <- dat$right.data
  
  plot.data.left <- prep(dat.left, plot = "left")
  plot.data.right <- prep(dat.right, plot = "right")
  
  triangle.plot.left <- NULL
  triangle.plot.right <- NULL
  triangle.plot.two <- NULL
  
  if (plot %in% c("one", "both")) {
    triangle.plot.left <- trianglePlot.poisson.aux(
      plot.data.left, NULL, testv, plot = "left", 
      fuzz.x, point.mass, add.points)
    
    triangle.plot.right <- trianglePlot.poisson.aux(
      NULL, plot.data.right, testv, plot = "right",
      fuzz.x, point.mass, add.points)
  }
  
  if (plot %in% c("two", "both")) {
    plot.data.left <- plot.data.left |> filter(p <= 0.5) 
    plot.data.right <- plot.data.right |> filter(p >= 0.5)
    
    triangle.plot.two <- trianglePlot.poisson.aux(
      plot.data.left, plot.data.right, testv,
      plot = "two", fuzz.x, point.mass, add.points)
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
  
  base.theme <- theme_classic() +
    theme(panel.spacing.x = unit(2, "lines"), 
          axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
  
  base.labs <- labs(x = "nominal p value",
                    y = "cumulative prop. more extreme",
                    title = "Triangle Plot: Poisson")
  
  if (plot %in% c("left", "right")){
    
    dat <- if (plot == "left") dat.left else dat.right
    step <- if (plot == "left") step.left else step.right
    ref.line <- 
      if (plot == "left") geom_abline(intercept = 0, slope = 1, linetype = 2) 
    else geom_abline(intercept = 1, slope = -1, linetype = 2) 
    
    triangle.plot <- ggplot(data = dat) +
      geom_step(aes(x = p, y = cum.prop), direction = step) +
      ref.line + 
      facet_grid(~lambda) +
      base.theme +
      base.labs
    
    if (add.points)
      triangle.plot <- triangle.plot + geom_point(aes(x = p, y = cum.prop))
    
  } else if (plot == "two"){
    # properly connects the left and right triangles ... facet on lambda
    joinData <- data.frame(lambda = rep(unique(dat.left$lambda), 2),
                           p = rep(0.5, 2*length(unique(dat.left$lambda))),
                           cum.prop = c(dat.left |> filter(p == 0.5) |> pull(cum.prop),
                                        dat.right |> filter(p == 0.5) |> pull(cum.prop)))
    
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
      base.theme +
      base.labs
    
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

# Comments:
# i) The 'family == "poisson"' condition exists for expsanion into more families
# ii) The default parameters behave for a non-fuzzed approach

trianglePlot <- function(lambdas, numSims = 1e4, numReps = 1, 
                         family = "poisson", testv = "poisson.test",
                         plot = "both", dat = NULL, fuzz.x = FALSE,
                         point.mass = TRUE, add.points = FALSE, add.checkplot = FALSE){
  
  if (family == "poisson"){ 
    dat <- generatePoisson(lambdas, dat, numSims, numReps, testv)
    
    if (add.checkplot){
      dat.temp <- dat
      
      if (!fuzz.x) 
        dat.temp$p <- dat.temp$left.exact
      
      checkplot <- checkPlot(dat.temp, facets = length(unique(dat$lambda))) + facet_grid(~lambda)
      print(checkplot) 
    }
    
    triangle.data <- triangleData.poisson(dat, point.mass, fuzz.x, testv)
    triangle.plots <- trianglePlot.poisson(triangle.data, testv, plot, 
                                           fuzz.x, point.mass, add.points)  
    
    left.plot <- triangle.plots$left.plot
    right.plot <- triangle.plots$right.plot
    two.plot <- triangle.plots$two.plot
    
    if (!is.null(left.plot)) print(left.plot)
    if (!is.null(right.plot)) print(right.plot)
    if (!is.null(two.plot)) print(two.plot)
  }
}

saveEnvironment()
