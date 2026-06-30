library(checkPlotR)
library(dplyr)
library(purrr)
library(reshape2)
library(shellpipes)
library(patchwork)
source("trianglePlot_Poisson.R")

fuzzingChi <- function(sims, lambda){
  
  pois <- rpois(sims, lambda) 
  p.pois <- ppois(pois, lambda)
  pois.chi <- ((pois  - lambda)/sqrt(lambda))^2
  p.pois.chi <- pchisq(pois.chi, df = 1, lower.tail = FALSE)
  pois.df <- data.frame(lambda = lambda, x = pois, x.chi = pois.chi, p = p.pois, p.chi = p.pois.chi)
  
  ggplot(data = pois.df) + 
    geom_histogram(aes(x=p.pois.chi), binwidth = 0.01)
  
  lrt.intervals.og <- pois.df %>%
    distinct(x, p.chi) %>%
    arrange(x) %>%
    mutate(rp.lower  = ifelse(x <= lambda, lag(p.chi, default = 0), p.chi),
           rp.upper = ifelse(x <= lambda, p.chi, lag(p.chi, default = 0))) 
  
  lrt.intervals.new <- pois.df %>%
    distinct(x, p.chi) %>%
    arrange(x) %>%
    mutate(rp.lower  = ifelse(x <= lambda, lag(p.chi, default = 0), lead(p.chi, default = 0)),
           rp.upper = p.chi)
  
  supervised.og <- pois.df %>% 
    left_join(lrt.intervals.og, by  = c("x", "p.chi")) %>%
    mutate(rp = runif(n(), rp.lower, rp.upper), p = rp)
  
  supervised.new <- pois.df %>% 
    left_join(lrt.intervals.new, by  = c("x", "p.chi")) %>%
    mutate(rp = runif(n(), rp.lower, rp.upper), p = rp)
  
  stats.og <- checkplotStats(supervised.og, binwidth = 0.01, varStat = 3)
  stats.new <- checkplotStats(supervised.new, binwidth = 0.01, varStat = 3)
  
  new.plot <- ggplot(data = supervised.new) + 
    geom_histogram(aes(x=rp), binwidth = 0.01) +
    theme_classic() +
    ggtitle("Alternate Method", 
            subtitle = paste0(sprintf("Lambda = %-3d | CDF = %-3.3f", lambda, stats.new$var3), collapse = "\n")) +
    xlim(0, 1)
  
  og.plot <- ggplot(data = supervised.og) + 
    geom_histogram(aes(x=rp), binwidth = 0.01) +
    theme_classic() +
    ggtitle("Chosen Method", 
            subtitle = paste0(sprintf("Lambda = %-3d | CDF = %-3.3f", lambda, stats.og$var3), collapse = "\n")) +
    xlim(0, 1)
  
  new.plot + og.plot
}

set.seed(1)
fuzzingChi(1e4, 5)
