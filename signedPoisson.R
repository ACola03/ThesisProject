# Signed Poisson Data:

# i) Calculate the signed p-value for poisson data
signedP.pois <- function(x, lambda){
  pl <- poisson.test(x, r = lambda, alternative = "less")$p.value
  pg <- poisson.test(x, r = lambda, alternative = "greater")$p.value
  return(ifelse(pl < pg , pl, -pg))
}

# ii) Aggregate poisson p-values into a dataframe
aggP.pois <- function(lambda, min, max){  
  plist.sign <- numeric()
  plist.two <- numeric()
  counts <- max(min-1,0):(max+1)

  for (x in 1:length(counts)){
    c <- counts[x]
    plist.sign[[x]] <- (signedP.pois(c, lambda))
    plist.two[[x]] <- poisson.test(c, r = lambda, alternative = "two.sided")$p.value
  }
  
  df <- data.frame(x = counts
                   , d = dpois(counts, lambda)
                   , p.left = ppois(counts, lambda, lower.tail = TRUE)
                   , p.right = ppois(counts, lambda, lower.tail = FALSE) + dpois(counts, lambda)
                   , p.sign = plist.sign
                   , p.two = plist.two)
  return(df)
}

# iii) Fuzz the signed p-values (toward 0)
signedFuzz <- function(dat, agg){
  
  signedFuzz.intervals <- signedFuzzInt(agg)
  signed <- dat %>% 
    left_join(signedFuzz.intervals, by = c("x")) %>%
    mutate(p.fuzz = runif(n(), min = p.min, max = p.max)) 
  
  return(signed)
}


# iv) Compute the signed fuzz intervals (auxiliary to signedFuzz)
signedFuzzInt <- function(agg){
  
  # Positive Signed
  pos <- agg %>% 
    filter(p.sign > 0) %>%
    distinct(p.sign, .keep_all = TRUE) %>%
    arrange(p.sign) %>%
    mutate(p.max = p.sign,
           p.min = lag(p.sign, default = 0)) %>%
    select(x, d, p.left, p.right, p.sign, p.min, p.max)
  
  # Negative Signed
  neg <- agg %>%
    filter(p.sign < 0) %>%
    distinct(p.sign, .keep_all = TRUE) %>%
    arrange(p.sign) %>%
    mutate(p.min = p.sign,
           p.max = lead(p.sign, default = 0)) %>%
    select(x, d, p.left, p.right, p.sign, p.min, p.max)
  
  signedFuzz.intervals <- rbind(neg, pos)
  return(signedFuzz.intervals)
}

# v) Poisson Wrapper: simulate signed poisson data and fuzz p-values
poissonData <- function(N, lambda, recover = c("left", "right", "none")){
  
  # Simulate poisson data
  poisson.sim <- data.frame(x = rpois(N, lambda))

  # Compute p-data (w/ signs) for all possible events
  poisson.agg <- aggP.pois(lambda, min(poisson.sim), max(poisson.sim))

  # Compute all fuzz intervals and fuzz the signed p-values
  poisson.fuzz <- signedFuzz(poisson.sim, poisson.agg)
  
  # Recover tails if desired
  if (recover != "none")
    poisson.fuzz <- recoverTails(poisson.fuzz, recover)
  
  return(poisson.fuzz)
}

# vi) Recover left or right tail if needed
recoverTails <- function(dat, recover = c("left", "right", "none")){
  
  if (recover == "left"){
    dat <- dat %>% 
      mutate(p.fuzz = ifelse(dat$p.fuzz < 0, dat$p.fuzz + 1, dat$p.fuzz))
  }
  else if (recover == "right"){
    dat <- dat %>%
      mutate(p.fuzz = ifelse(dat$p.fuzz > 0, -1*(dat$p.fuzz-1) , -1*dat$p.fuzz))
  }
 
  return(dat) 
}

# -----

N <- 1e5; lambda <- 100
dat <- poissonData(N, lambda, "none") 

range(dat$p.fuzz)

min_val <- min(dat$p.fuzz)
max_val <- max(dat$p.fuzz)

ggplot(data = dat) +
  geom_histogram(aes(x = p.fuzz), binwidth = 0.025, boundary = min_val) +
  coord_cartesian(xlim = c(min_val, max_val))

# -----
