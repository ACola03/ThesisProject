# ==== Small rangePlot Issue

# Test: A dominating target that is completely above or below all CIs
# Issue: All intervals should be colored red, but actually appear grey
# Reason: if the "grey" group never shows up in the geom_pointrange, then the 
#         first alphabetical group is "red" which will be assigned to "grey"
# Solution: manually specify the colors according to each group

# ===== Minor rangePlot Modification

rangePlot <- function(tf
                      , target = "sampleMean"
                      , orderFun = slug
                      , conf = 0.95
                      , opacity = 0.2
                      , size = 0.1
                      , title = "Range plot"
                      , targNum = 1e3){
  if(target == "sampleMean"){
    target <- mean(tf$est)
  }
  thinner <- max(floor(length(tf$est)/targNum), 1) # is this reducing the # of CIs plotted
  thinned <- tf[seq(thinner
                    , length(tf$est)
                    , thinner)
                ,
  ]
  return(ggplot(
    orderFun(thinned)
    , aes(x = quantile, y = est, ymin = lower, ymax = upper)
  )
  + geom_pointrange(alpha = opacity, size = size
                    , aes(color = ifelse(lower > target | upper < target
                                         , "red", "grey")))
  + geom_hline(yintercept = target, color = "blue")
  + geom_vline(
    xintercept = c((1 - conf)/2, 1 - (1 - conf)/2)
    , lty = 2, col = "red"
  )
  + xlab("index")
  + ylab("estimate")
  + ggtitle(title)
  + scale_x_continuous(expand = c(0,0)
                       , breaks = c((1 - conf)/2, 0.25, 0.5, 0.75
                                    , conf + (1 - conf)/2)
                       , labels = function(breaks){signif(breaks, 3)})
  + scale_color_manual(values = c("grey" = "grey", "red" = "red"))
  + guides(color = "none")
  + theme_classic()
  )
}

# ===== SlugPlot Related

# ----- tailCoverage()
# Description:
# - computes the proportion of CIs that miss the target value
# - used for understanding what 'n' is needed to obtain 2-3% in each tail

tailCoverage <- function(param.dat, type, n, boots){
  if (type == "mean"){
    target <- 0.5
    lower <- param.dat$mean.lower
    upper <- param.dat$mean.upper
  } else{
    target <- 1
    lower <- param.dat$shape.lower
    upper <- param.dat$shape.upper
  }
  
  u <- mean(lower > target)
  l <- mean(upper < target)
  c.message <- paste0("Tail Coverage of ", paste0(toupper(substring(type,1,1)), substring(type, 2, )), 
                      " Parameter: n = ", n, ", boots = ", boots)
  u.message <- paste0("Upper-Tail Red Zone: ", signif(u,3))
  l.message <- paste0("Lower-Tail Red Zone: ", signif(l,3))
  t.message <- paste0("Total Red Zone: ", signif(l + u,3))
  f.message <- paste(c.message, u.message, l.message, t.message, sep = "\n")
  message(f.message)
}

# ----- showSlugs()
# Description:
# - a slugPlot wrapper for rangePlot() and tailCoverage()

showSlugs <- function(param.dat, shape1, shape2, n, boots){
  # compute how many times the CIs miss
  tailCoverage(param.dat, "mean",  n, boots)
  tailCoverage(param.dat, "shape",  n, boots)
  
  # extract mean and shape data
  means <- data.frame("est" = param.dat$mean, "se" = param.dat$mean.se, 
                      "lower" = param.dat$mean.lower, "upper" = param.dat$mean.upper)
  shapes <- data.frame("est" = param.dat$shape, "se" = param.dat$shape.se, 
                       "lower" = param.dat$shape.lower, "upper" = param.dat$shape.upper)
  
  # slug plots of mean and shape
  meanPlot <- rangePlot(means, target = 0.5, orderFun = slug, title = "SlugPlot of Mean")
  shapePlot <- rangePlot(shapes, target = 1, orderFun = slug, title = "SlugPlot of Size")
  
  # combined plots
  combinedPlot <- (meanPlot + shapePlot)

  return(combinedPlot)
}

# ===== ContourPlot Related

# ----- cutoff()
# Description:
# - helper function to determine where the low-density region begins

cutoff <- function(kde.dat){
  dx <- diff(kde.dat$x[1:2])
  dy <- diff(kde.dat$y[1:2])
  sz <- sort(kde.dat$z)
  c1 <- cumsum(sz) * dx * dy
  level_cutoff <- approx(c1, sz, xout = 1 - conf.level)$y
  return(level_cutoff)
}

# ----- showContours()
# Description:
# - shows the 95% high-density region of the joint parameter distribution
# - useful in comparing the joint dist between a uniform and another beta

showContours <- function(param.dat, shape1, shape2, n, boots){
  # extract (mean,shape) data
  x <- param.dat$mean
  y <- param.dat$shape
  df_pts <- data.frame(x = x, y = y)
  
  kde.load <- readRDS("meanShape_kde.rds")
  
  # cutoff contours for chosen and loaded beta
  kde <- MASS::kde2d(x, y, n = 100)
  level_cutoff <- cutoff(kde)
  level_cutoff.load <- cutoff(kde.load) 
  
  # extract contour lines for the loaded (uniform) distribution
  cl.load <- contourLines(kde.load$x, kde.load$y, kde.load$z, levels = level_cutoff.load)
  df_contour.load <- do.call(rbind, lapply(seq_along(cl.load), function(i) {
    data.frame(x = cl.load[[i]]$x, y = cl.load[[i]]$y, group = i)
  }))
  
  unif_label <- paste0("Uniform ", conf.level * 100, "%")
  
  # build the contour plot with the joint distribution
  p <- ggplot() +
    geom_point(data = df_pts, aes(x = x, y = y), color = "gray65", shape = 1) +
    theme_classic() +
    labs(
      title = "Joint Parameter Distribution",
      x = "mean",
      y = "size",
      color = "Contour"
    )
  
  # add the perturbed contour only if it's non-uniform
  if (!all(c(shape1, shape2) == c(1, 1))){
    cl.pert <- contourLines(kde$x, kde$y, kde$z, levels = level_cutoff)
    df_contour.pert <- do.call(rbind, lapply(seq_along(cl.pert), function(i) {
      data.frame(x = cl.pert[[i]]$x, y = cl.pert[[i]]$y, group = i)
    }))
    
    pert_label <- paste0("Perturbated ", conf.level * 100, "%")
    
    p <- p + geom_path(data = df_contour.pert, aes(x = x, y = y, group = group, color = pert_label), 
                       linewidth = 1)
    
    # overlap calculation in the 95% uniform CI
    contourPolygon.load <- do.call(rbind, lapply(cl.load, function(contour) data.frame(x = contour$x, y = contour$y)))
    overlap <- mean(sp::point.in.polygon(x, y, contourPolygon.load$x, contourPolygon.load$y) != 0)
    
    message(sprintf("Overlap within the Uniform 95%% 2D Region: %.4f", overlap))
  }
  
  # now overlap the uniform so that the plot is focused on the perturbed region
  p <- p + geom_path(data = df_contour.load, 
                     aes(x = x, y = y, group = group, color = unif_label), # will the 'color' in labs work
                     linewidth = 1) 
  
  # green: uniform; red: perturbed
  p <- p + scale_color_manual(values = setNames(
    c("#1cd44d", "red"), 
    c(paste0("Uniform ", conf.level * 100, "%"), paste0("Perturbated ", conf.level * 100, "%"))
  ))
  
  return(p)
}

# ===== Wrapper

wrap <- function(param.dat, shape1, shape2, n, boots){
  
  slugPlot <- showSlugs(param.dat, shape1, shape2, n, boots)
  contourPlot <- showContours(param.dat, shape1, shape2, n, boots)
  pianoPlot <- ggplot(data = data.frame(x = rbeta(boots, shape1, shape2))) +
    geom_histogram(aes(x), breaks = seq(0,1,length.out = 41)) +
    labs(title = "PianoPlot")
  
  mainTitle <- paste0("Beta(", shape1, ",", shape2, ") Parameter Analysis\n(N,n) = (", boots, ",", n, ")")
  
  # combined plots: weird warning from interaction of ggplot and patchwork
  suppressWarnings({
    combinedPlot <- (slugPlot) / (contourPlot + pianoPlot) +
      plot_layout(guides = "collect") +
      patchwork::plot_annotation(
        title = mainTitle,
        theme = theme(plot.title = element_text(size = 14))
      )
    print(combinedPlot)
  })
  
}

# Notes:
# 1. Eventually fix title: bigger font
