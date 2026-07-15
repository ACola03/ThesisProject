# Signed Binomial Data:

library(dplyr)
library(ggplot2)

# i) Calculate the signed p-value for binomial data
signedP.binom <- function(x, n, p){
  pl <- binom.test(x, n, p, alternative = "less")$p.value
  pg <- binom.test(x, n, p, alternative = "greater")$p.value
  return(ifelse(pl < pg , pl, -pg))
}

# ii) Aggregate binomial p-values into a dataframe
aggP.binom <- function(n, p){  
  plist.sign <- numeric()
  counts <- 0:n
  
  for (x in counts){
    plist.sign[[x+1]] <- (signedP.binom(x, n, p))
    ## plist.two[[x+1]] <- binom.test(x, n, p, alternative = "two.sided")$p.value
  }
  
  df <- tibble(x = counts
                   , d = dbinom(counts, n, p)
                   , p.left = pbinom(counts, n, p, lower.tail = TRUE)
                   , p.right = pbinom(counts, n, p, lower.tail = FALSE) + dbinom(counts, n, p)
                   , p.sign = plist.sign
  )
  return(df)
}

# iii) Fuzz the signed p-values (toward 0)
signedFuzz <- function(dat, agg){
  
  signedFuzz.intervals <- signedFuzzInt(agg)
  signed <- dat |> 
    left_join(signedFuzz.intervals, by = c("x")) |>
    mutate(sp.fuzz = runif(n(), min = p.min, max = p.max)
	 	, p.fuzz = if_else(sp.fuzz<0, 1+sp.fuzz, sp.fuzz)
	 ) 
  
  return(signed)
}


# iv) Compute the signed fuzz intervals (auxiliary to signedFuzz)
signedFuzzInt <- function(agg){
  
  # Positive Signed
  pos <- agg |> 
    filter(p.sign > 0) |>
    distinct(p.sign, .keep_all = TRUE) |>
    arrange(p.sign) |>
    mutate(p.max = p.sign,
           p.min = lag(p.sign, default = 0)) |>
    select(x, d, p.left, p.right, p.sign, p.min, p.max)
  
  # Negative Signed
  neg <- agg |>
    filter(p.sign < 0) |>
    distinct(p.sign, .keep_all = TRUE) |>
    arrange(p.sign) |>
    mutate(p.min = p.sign,
           p.max = lead(p.sign, default = 0)) |>
    select(x, d, p.left, p.right, p.sign, p.min, p.max)
  
  signedFuzz.intervals <- rbind(neg, pos)
  return(signedFuzz.intervals)
}

# v) Binomial Wrapper: simulate signed binomial data and fuzz p-values
binomialData <- function(N, n, p){
  
  # Simulate binomial data
  binomial.sim <- data.frame(x = rbinom(N, n, p))
  
  # Compute p-data (w/ signs) for all possible events
  binomial.agg <- aggP.binom(n, p)

  # Compute all fuzz intervals and fuzz the signed p-values
  binomial.fuzz <- signedFuzz(binomial.sim, binomial.agg)
  
  return(binomial.fuzz)
}

# -----

N <- 1e5; n <- 6; p <- 0.5
dat <- binomialData(N, n, p)

bins <- 200

print(ggplot(data = dat)
	+ geom_histogram(aes(x=p.fuzz), breaks = seq(0, 1, length.out=bins+1))
)

print(ggplot(data = dat)
	+ geom_histogram(aes(x=sp.fuzz), bins=bins)
)

range(dat$p.fuzz)

# -----
