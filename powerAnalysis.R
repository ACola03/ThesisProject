library(dplyr)
library(purrr)
library(reshape2)
library(shellpipes)
library(tidyr)
library(ggplot2)
library(twosamples)
library(knitr)
library(kableExtra)
source("checkplotStatsDist.R")

# Used to produce July22.md 
# Demonstrates the checkplot statistics under beta perturbations

# -----

plotAnalysis <- function(data, stat, N, n, showStats = FALSE){
  
  # add summary statistics to column names
  if (showStats){
    means <- numeric()
    sds <- numeric()
    ranges <- character()
    
    for (col in 1:ncol(data)){
      dat <- data[,col]
      means[[col]] <- mean(dat)
      sds[[col]] <- sd(dat)
      rng_vals    <- range(dat)
      ranges[col] <- paste0("[", round(rng_vals[1], 3), ", ", round(rng_vals[2], 3), "]")
    }
    colnames(data) <- paste0(colnames(data),"\n" , ranges, "\nMean:", round(means,3), "\nSd:", round(sds,3))
  }
  
  # visualize in boxplots
  tib <- as_tibble(data)
  plotData <- pivot_longer(tib, cols = everything(), names_to = "Beta_Type", values_to = "Stat") 
  axisText <- paste0(stat, " Statistic")
  titleText <- paste0(stat, " Statistic Comparisons: (N,n) = (", N, ",", n, ")")
  
  ggplot(plotData, aes(x = reorder(Beta_Type, Stat, mean), y = Stat, fill = Beta_Type)) +
    geom_boxplot() +
    labs(x = "Beta Distribution", y = axisText, 
         title = titleText) +
    theme_minimal() +
    theme(legend.position = "none")
}

plotAllAnalysis <- function(bin, space, location, dts, N, n, showStats = FALSE){
  # combine into one dataframe
  statType <- factor(c(rep("Bin",N), rep("Space",N), rep("Location",N), rep("DTS",N)), 
                     levels = c("Bin", "Space", "Location", "DTS"))
  allStats <- data.frame("statType" = statType, rbind(bin,space,location,dts), check.names = FALSE) |> 
    pivot_longer(cols = -c("statType"), names_to = "Beta_Type", values_to = "Value") |>
    mutate(Beta_Type = as.factor(Beta_Type)) 
  
  # visualize in boxplots
  titleText <- paste0("Checkplot Statistic Comparisons: (N,n) = (", N, ",", n, ")")
  ggplot(data = allStats) +
    geom_boxplot(aes(x = reorder(Beta_Type, Value, mean), y = Value, fill = Beta_Type)) + 
    facet_wrap(~statType, scales = "free_y", nrow = 4) +
    theme_bw() + 
    theme(legend.position = "none") + 
    labs(y = "Statistic Value", x = "Beta Distribution", title = titleText)
}

# Note: Requires Beta(1,1) to be the first column
powerAnalysis <- function(data, stat, N, n){
  nullMin <- min(data[,1])
  nullMax <- max(data[,1])
  powerVec <- numeric()
  
  for (col in 2:ncol(data)){
    powerVec[[col-1]] <- 1 - mean(between(data[,col], nullMin, nullMax))
  }
  
  out <- paste(paste0(sprintf("%s Power Analysis: (N,n) = (%d,%d)", stat, N, n)),
               paste0(sprintf("%s: %.4f", colnames(data)[2:ncol(data)], powerVec), collapse = "\n"),
               sep = "\n")
  cat(out)
  return(powerVec)
}

makeTablePower <- function(dat, stat, dat.n1, dat.n2, dat.n3){
  tableData <- data.frame(
    Distribution = colnames(dat)[2:ncol(dat)],
    powerLow = dat.n1, 
    powerMed = dat.n2,
    powerHigh = dat.n3)
  
  colnames(tableData) <- c("Distribution","Power","Power","Power")
  caption <- paste(stat, "Power Analysis Across Sample Sizes (N = 10000)")
  
  ltx <- kbl(tableData, 
             format = "latex", 
             booktabs = TRUE, 
             align = "lccc",
             caption = caption) |>
    kable_styling(
      latex_options = c("hold_position"), 
      font_size = 12 
    ) |>
    add_header_above(c(" " = 1, "n = 100" = 1, "n = 1000" = 1, "n = 10000" = 1), bold = TRUE)
  
  print(ltx)
}

makeTableStats <- function(bin, space, location, dts, n){
  tableData <- data.frame(
    Distribution = colnames(bin),
    meanBin = round(unname(apply(bin, 2, mean)),3), 
    sdBin = round(unname(apply(bin, 2, sd)),3), 
    meanSpace = round(unname(apply(space, 2, mean)),3),
    sdSpace = round(unname(apply(space, 2, sd)),3),
    meanLoc = round(unname(apply(location, 2, mean)),3),
    sdLoc = round(unname(apply(location, 2, sd)),3),
    meanDts = round(unname(apply(dts, 2, mean)),3),
    sdDts = round(unname(apply(dts, 2, sd)),3)
    )
  
  colnames(tableData) <- c("Distribution","Mean","Sd","Mean","Sd","Mean","Sd","Mean","Sd")
  caption <- paste("Summary Statistics (n = ", n, ")", sep = "")
  
  ltx <- kbl(tableData, 
             format = "latex", 
             booktabs = TRUE, 
             align = "lcccccccc",
             caption = caption) |>
    kable_styling(
      latex_options = c("hold_position"), 
      font_size = 12 
    ) |>
    add_header_above(c(" " = 1, "Bin" = 2, "Space" = 2, "Location" = 2, "DTS" = 2), bold = TRUE)
  
  print(ltx)
}


# -----

# configurations
set.seed(1)
N <- 1e4; n1 <- 1e2; n2 <- 1e3; n3 <- 1e4
shapes1 <- c(1.0, 1.25, 2.0,  0.75)
shapes2 <- c(1.0, 0.75, 2.0,  0.75)

# compute
binvarAnalysis1 <- binvar_dist(N, n1, shapes1, shapes2, 0.025)
binvarAnalysis2 <- binvar_dist(N, n2, shapes1, shapes2, 0.025)
binvarAnalysis3 <- binvar_dist(N, n3, shapes1, shapes2, 0.025)

locationAnalysis1 <- location_dist(N, n1, shapes1, shapes2, TRUE, TRUE)
locationAnalysis2 <- location_dist(N, n2, shapes1, shapes2, TRUE, TRUE)
locationAnalysis3 <- location_dist(N, n3, shapes1, shapes2, TRUE, TRUE)

spaceAnalysis1 <- space_dist(N, n1, shapes1, shapes2)
spaceAnalysis2 <- space_dist(N, n2, shapes1, shapes2)
spaceAnalysis3 <- space_dist(N, n3, shapes1, shapes2)

dtsAnalysis1 <- dts_dist(N, n1, shapes1, shapes2)
dtsAnalysis2 <- dts_dist(N, n2, shapes1, shapes2)
dtsAnalysis3 <- dts_dist(N, n3, shapes1, shapes2)

# plot
plotAllAnalysis(binvarAnalysis1, spaceAnalysis1, locationAnalysis1, dtsAnalysis1, N, n1)
plotAllAnalysis(binvarAnalysis2, spaceAnalysis2, locationAnalysis2, dtsAnalysis2, N, n2)
plotAllAnalysis(binvarAnalysis3, spaceAnalysis3, locationAnalysis3, dtsAnalysis3, N, n3)

# power
binvar1.p <- powerAnalysis(binvarAnalysis1, "Bin (0.025) Variance", N, n1)
binvar2.p <- powerAnalysis(binvarAnalysis2, "Bin (0.025) Variance", N, n2)
binvar3.p <- powerAnalysis(binvarAnalysis3, "Bin (0.025) Variance", N, n3)

loc1.p <- powerAnalysis(locationAnalysis1, "Location", N, n1)
loc2.p <- powerAnalysis(locationAnalysis2, "Location", N, n2)
loc3.p <- powerAnalysis(locationAnalysis3, "Location", N, n3)

sp1.p <- powerAnalysis(spaceAnalysis1, "Space", N, n1)
sp2.p <- powerAnalysis(spaceAnalysis2, "Space", N, n2)
sp3.p <- powerAnalysis(spaceAnalysis3, "Space", N, n3)

dts1.p <- powerAnalysis(dtsAnalysis1, "DTS", N, n1)
dts2.p <- powerAnalysis(dtsAnalysis2, "DTS", N, n2)
dts3.p <- powerAnalysis(dtsAnalysis3, "DTS", N, n3)

# tabulated (power)
makeTablePower(binvarAnalysis1, "Bin (0.025) Variance", binvar1.p, binvar2.p, binvar3.p)
makeTablePower(locationAnalysis1, "Location", loc1.p, loc2.p, loc3.p)
makeTablePower(spaceAnalysis1, "Space", sp1.p, sp2.p, sp3.p)
makeTablePower(dtsAnalysis1, "DTS", dts1.p, dts2.p, dts3.p)

# tabulated (summary stats)
makeTableStats(binvarAnalysis1, spaceAnalysis1, locationAnalysis1, dtsAnalysis1, n1)
makeTableStats(binvarAnalysis2, spaceAnalysis2, locationAnalysis2, dtsAnalysis2, n2)
makeTableStats(binvarAnalysis3, spaceAnalysis3, locationAnalysis3, dtsAnalysis3, n3)

# -----

N <- 10000

scenario_labels <- c(
  "Beta(1, 1)", 
  "Beta(1.25, 0.75)", 
  "Beta(2, 2)", 
  "Beta(0.75, 0.75)"
)

beta_data <- data.frame(
  "Beta(1, 1)"         = rbeta(N, 1, 1),
  "Beta(1.25, 0.75)" = rbeta(N, 1.25, 0.75),
  "Beta(2, 2)"   = rbeta(N, 2, 2),
  "Beta(0.75, 0.75)" = rbeta(N, 0.75, 0.75),
  check.names = FALSE
)

plot_df <- beta_data %>%
  pivot_longer(cols = everything(), names_to = "Scenario", values_to = "Value") %>%
  mutate(Scenario = factor(Scenario, levels = scenario_labels))

ggplot(plot_df, aes(x = Value, fill = Scenario)) +
  geom_histogram(bins = 40, color = "white", alpha = 0.85, boundary = 0) +
  facet_wrap(~Scenario, scales = "free_y", ncol = 2) +
  theme_minimal(base_size = 13) +
  theme(
    legend.position = "none",
    strip.text = element_text(face = "bold", size = 11),
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold", hjust = 0.5)
  ) +
  labs(
    title = paste0("Beta Distribution Histograms (N = ", format(N, scientific=FALSE), ")"),
    x = "Random Variate Value",
    y = "Frequency"
  )

# -----

