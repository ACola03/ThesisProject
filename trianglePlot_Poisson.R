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

generatePoisson <- function(lambdas, dat = NULL,  
                            numSims = 1e4, numReps = 1, testv = "poisson.test"){
  
  dat.null <- is.null(dat)
  
  poisson.data <- purrr::map(lambdas, function(lambda){
    if (testv == "poisson.test" & dat.null){
      dat <- rpois(numSims, lambda)
      
    } else if (testv %in% c("wald.intercept", "lrt") & dat.null){
      dat <- purrr::map(c(1:numSims), function(dummy){
        rpois(numReps, lambda)
      }) |> rbind()
    }
    
    df <- data.frame(testv, lambda, multPois(dat, lambda, testv))}
  ) |> list_rbind()
  
  # for stability
  if (testv %in% c("wald.intercept", "lrt"))
    poisson.data <- poisson.data |> group_by(pois.mean, lambda) |> mutate(cp = round(max(cp),6)) %>% ungroup() # !
  
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
        return(data.frame(pois.mean = d,
                          left.exact = bt,
                          right.exact = gt,
                          cp = bt)
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
        return(data.frame(pois.mean = mean(d),
                          int.est = est, 
                          int.se = se,
                          z.value = w,
                          left.exact = bt, 
                          right.exact = gt, 
                          cp = bt)
        )
      }) |> 
      list_rbind()
    
  } else if (testv == "lrt"){
    df <-
      purrr::map(dat, function(d){
        null.model <- glm(d ~ -1 + offset(log(0*d+lambda0)), family = poisson())
        full.model <- glm(d ~ 1 + offset(log(0*d+lambda0)), family = poisson())
        est <- coef(full.model)[[1]] # N(0, s2) ... mean intercept 0 
        anova <- anova(null.model, full.model, test = "LRT")
        deviance <- anova$Deviance[2]
        p <- ifelse(deviance < 0, 1, anova$`Pr(>Chi)`[2]) # sometimes negative -> NA
        return(data.frame(pois.mean = mean(d),
                          int.est = est,
                          deviance = deviance,
                          left.exact = 1-p, 
                          right.exact = p, 
                          cp = p)
        )
      }) |> 
      list_rbind()
  } 
}

# -----

# FUZZING

fuzzPoisson <- function(dat, testv, fuzz.type = "supervised", use.fuzz = TRUE){
  
  if (fuzz.type == "supervised"){
    dat <- supervised.fuzz(dat, testv)
  } 
  else if (fuzz.type == "unsupervised"){
    dat <- unsupervised.fuzz(dat)
  }
  
  if (use.fuzz) dat$p <- dat$rp
  else dat$p <- dat$cp
  
  return(dat)
}

supervised.fuzz <- function(dat, testv){
  
  if (testv == "poisson.test"){
    supervised <- dat%>%
      mutate(rp.lower = 1 - right.exact,
             rp.upper = left.exact,
             rp = runif(n(), rp.lower, rp.upper))
  }
  else if (testv == "lrt"){
    lrt.intervals <- dat %>%
      distinct(lambda, pois.mean, cp) %>%
      arrange(lambda, pois.mean) %>%
      group_by(lambda) %>%
      mutate(rp.lower  = ifelse(pois.mean <= lambda, lag(cp, default = 0), cp),
             rp.upper = ifelse(pois.mean <= lambda, cp, lag(cp, default = 0))) %>%
      ungroup()
    
    supervised <- dat %>% 
      left_join(lrt.intervals, by  = c("lambda", "pois.mean", "cp")) %>%
      mutate(rp = runif(n(), rp.lower, rp.upper))
  }
  else if (testv == "wald.intercept"){ # AN IDEA I HAD
    wald.intervals <- dat %>%
      distinct(lambda, pois.mean, p) %>%
      arrange(lambda, pois.mean) %>%
      mutate(p = ifelse(pois.mean == 0, 0, p)) %>%
      group_by(lambda) %>%
      mutate(rp.lower  = ifelse(pois.mean <= lambda, p, lag(p)),
             rp.upper = ifelse(pois.mean <= lambda, lead(p, default = 1), p)) %>%
      ungroup() %>%
      select(-p)
    
    supervised <- dat %>% 
      left_join(wald.intervals, by  = c("lambda", "pois.mean")) %>%
      mutate(rp = runif(n(), rp.lower, rp.upper))
  }
  
  return(supervised)
}

unsupervised.fuzz <- function(dat){
  
  fuzz.intervals <- dat %>% 
    distinct(lambda, pois.mean, cp) %>% 
    arrange(lambda, cp) %>% # if you sort based on pois.mean it would be supervised, otherwise arrange on cp
    group_by(lambda) %>%
    mutate(
      rp.lower = lag(cp, default = 0),
      rp.upper = cp
    ) %>%
    ungroup()
  
  fuzz.df <- dat %>% 
    left_join(fuzz.intervals, by = c("lambda", "cp", "pois.mean")) %>% 
    mutate(rp = runif(n(), min = rp.lower, max = rp.upper))
  
  return(fuzz.df)
}

# =======================

# TRIANGLE PLOT PREPARATION

# Description: prepare data for triangle plot creation
# i) the acceptable 'plot' input is 'one' or 'two' to denote which tail types we want to plot
#    'one' will plot one-sided triangles, while 'two' will plot the combined triangles
#    the choice determines how the data will be pre-processed

triangleData.poisson <- function(dat, point.mass = TRUE, use.fuzz = FALSE, 
                                 testv = "poisson.test"){
  
  dat$p <- if (use.fuzz | testv != "poisson.test") dat$p else dat$left.exact
  thresholds <- sort(unique(c(0, 0.5, 1, dat$p)))
  left.data <- triangleData.poisson.aux(dat, thresholds, point.mass)
  
  dat$p <- if (use.fuzz | testv != "poisson.test") 1 - dat$p else dat$right.exact
  thresholds.right <- sort(unique(c(0, 0.5, 1, dat$p)))
  right.data <- triangleData.poisson.aux(dat, thresholds, point.mass)
  
  return(list("left.data" = left.data, "right.data" = right.data))
}

# -----

triangleData.poisson.aux <- function(dat, thresholds, point.mass = TRUE){
  
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
                                 use.fuzz = TRUE, point.mass = TRUE, add.points = FALSE){
  
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
      use.fuzz, point.mass, add.points)
    
    triangle.plot.right <- trianglePlot.poisson.aux(
      NULL, plot.data.right, testv, plot = "right",
      use.fuzz, point.mass, add.points)
  }
  
  if (plot %in% c("two", "both")) {
    plot.data.left <- plot.data.left |> filter(p <= 0.5) 
    plot.data.right <- plot.data.right |> filter(p >= 0.5)
    
    triangle.plot.two <- trianglePlot.poisson.aux(
      plot.data.left, plot.data.right, testv,
      plot = "two", use.fuzz, point.mass, add.points)
  }
  
  return(list("left.plot" = triangle.plot.left, 
              "right.plot" = triangle.plot.right,
              "two.plot" = triangle.plot.two))
}

# -----

trianglePlot.poisson.aux <- function(dat.left, dat.right, testv, plot = "left", 
                                     use.fuzz = TRUE, point.mass = TRUE, add.points = FALSE){
  
  inclusive <- point.mass | testv != "poisson.test"
  step.left <- ifelse(inclusive | use.fuzz, "hv", "vh")
  step.right <- ifelse(inclusive| use.fuzz, "vh", "hv")
  
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

# CheckPlots:
# A wrapper for the checkPlot() function that add titles that feature
# some meaningful statistics on the p-values of interest

# i) varStat = 0 is for testing and direct comparison of variance methods

checkplotWrapper <- function(dat, numSims, numReps, testv, binwidth, varStat = 0){
  
  checkPlot.title <- paste0(
    sprintf("CheckPlot: %s (N = %d, n = %d)",
            testv, numSims, numReps))
  
  checkplot.var <- checkplotStats(dat, binwidth, varStat)
  var1 <- checkplot.var$var1
  var2 <- checkplot.var$var2
  var3 <- checkplot.var$var3
  
  checkPlot.subtitle <- switch(
    as.character(varStat),
    "1" = paste0(sprintf("Lambda = %-3d | Bar = %-6.3f", sort(unique(dat$lambda)), var1), collapse = "\n"),
    "2" = paste0(sprintf("Lambda = %-3d | Space = %-3.3f", sort(unique(dat$lambda)), var2), collapse = "\n"),
    "3" = paste0(sprintf("Lambda = %-3d | CDF = %-3.3f", sort(unique(dat$lambda)), var3), collapse = "\n"),
    "0" = paste0(sprintf("Lambda = %-3d | Bar = %-6.3f | Space = %-3.3f | CDF = %-3.3f", sort(unique(dat$lambda)), var1, var2, var3), collapse = "\n"),
    NULL
    )  # what if i have (10, 1) is it sorted ... split and facet will return it as sorted, so var is in ascending lambda order
  
  checkplot <- checkPlot(dat, breaks = seq(0,1,binwidth), facets = length(unique(dat$lambda))) + 
    facet_grid(~lambda) + 
    theme_classic() + 
    ggtitle(checkPlot.title, subtitle = checkPlot.subtitle) +
    theme(
      plot.title = element_text(family = "mono"), 
      plot.subtitle = element_text(family = "mono"),
      panel.spacing.x = unit(1, "cm")) # required for aligned titles
  
  print(checkplot)
}


checkplotStats <- function(dat, binwidth, varStat = 3){
  
  var1 <- NA
  var2 <- NA
  var3 <- NA
  
  # stat 1 bins: like the checkPlot
  if (varStat == 1 | varStat == 0){
    counts <- dat %>%
      ungroup() %>%
      mutate(binID = floor(p/binwidth)) %>%
      count(lambda, binID) %>%
      tidyr::complete(lambda, binID, fill = list(n = 0))
    
    var1 <- counts %>%
      split(counts$lambda) %>%
      purrr::map(function(dat.split){
        var(dat.split$n)
      })
  }
  
  # Stat2 Spacing: relies on the spacings of adjacent p-values
  if (varStat == 2 | varStat == 0){
    var2 <- dat %>%
      split(dat$lambda) %>%
      purrr::map(function(dat.split){
        pvals <- sort(c(0,dat.split$p,1))
        gaps <- diff(pvals)
        var.gaps <- var(length(gaps)*gaps) # (n+1)*(p_{n+1} - p_{n})
      }) 
  }
  
  # Stat3 CDF: compares empirical p-value to expected order statistic
  if (varStat == 3 | varStat == 0){
    var3 <- dat %>%
      split(dat$lambda) %>%
      purrr::map(function(dat.split){
        n <- nrow(dat.split)
        pvec <- sort(dat.split$p)
        lvec <- (1:n)/(n+1)
        var <- ((1:n)*(n - (1:n) + 1)) / ((n+1)^2*(n+2))
        sign <- sign(sum(pvec - lvec))
        sq.err <- (pvec - lvec)^2
        location <- sign * sum(sq.err/var)
        return(location) 
      })
  }
  
  return(list("var1" = var1, "var2" = var2, "var3" = var3))
}

# ==========

# WRAPPER

# Comments:
# i) The 'family == "poisson"' condition exists for expansion into more families
# ii) The default parameters behave for a non-fuzzed approach

trianglePlot <- function(lambdas, numSims = 1e4, numReps = 1, 
                         family = "poisson", testv = "poisson.test",
                         fuzz.type = "supervised", use.fuzz = FALSE,
                         binwidth = 0.05, varStat = 3,
                         plot = "both", dat = NULL, 
                         point.mass = TRUE, add.points = FALSE, add.checkplot = FALSE){
  
  if (family == "poisson"){ 
    dat <- generatePoisson(lambdas, dat, numSims, numReps, testv)
    dat <- fuzzPoisson(dat, testv, fuzz.type, use.fuzz)
    
    if (add.checkplot){
      checkplotWrapper(dat, numSims, numReps, testv, binwidth, varStat)
    }
    
    
    if (plot != "none") { # for checkPlot priority
      triangle.data <- triangleData.poisson(dat, point.mass, use.fuzz, testv)
      triangle.plots <- trianglePlot.poisson(triangle.data, testv, plot, 
                                             use.fuzz, point.mass, add.points)  
      
      left.plot <- triangle.plots$left.plot
      right.plot <- triangle.plots$right.plot
      two.plot <- triangle.plots$two.plot
      
      if (!is.null(left.plot)) print(left.plot)
      if (!is.null(right.plot)) print(right.plot)
      if (!is.null(two.plot)) print(two.plot)
    }
  }
}


# ==========

# WALD FUZZ

# Notes:
# i) filter.zero controls whether the count is conditional on being greater than 0
# ii) currently, due to decimal precision p-values for identical counts can be very
#     slightly different, but identical under print(). This causes some issue when 
#     naming columns by the p-value, which is needed to determine the empirical 
#     extremeness. Because of this, I just round to 4 (could change to 6) decimal
#     places. However, a better way to do this is to just take the maximum p-value
#     for each count, since they should only be slightly different (unnoticable). 

wald.fuzz <- function(wald.data, filter.zero = FALSE){
  
  if (filter.zero) wald.data <- wald.data %>% filter(pois.mean > 0)
  
  wald.intervals <- wald.data %>%
    distinct(lambda, pois.mean, p) %>% 
    arrange(lambda, p) %>% # if you sort based on pois.mean it would be supervised, otherwise arrange on p
    group_by(lambda) %>%
    mutate(
      p.lower = lag(p, default = 0),
      p.upper = p
    ) %>%
    ungroup()
  
  
  fuzz.df <- wald.data %>% 
    left_join(wald.intervals, by = c("lambda", "p", "pois.mean")) %>% # but joining on pois.mean doesn't invoke supervision
    mutate(rp = runif(n(), min = p.lower, max = p.upper))
  
  return(fuzz.df)
}

saveEnvironment()
