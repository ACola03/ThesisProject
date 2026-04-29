library(purrr)
library(checkPlotR)

lam <- 12
N <- 10
reps <- 1e2

dat <- (map_dfr( c(rep=1:reps), function(dummy){rpois(N, lam)}))
print(dat)
d <- dat[[1]]

g <- glm(d~1+offset(log(0*d+lam)), family=poisson())
summary(g)

data.frame(t(rnorm(10)))

q()

tests <- map_dfr(1:nrow(dat), function(samp){
  p<-t.test(datNorm[samp, -1], mu=datNorm[samp,1], alternative="l")$p.value
  ci<-t.test(datNorm[samp, -1], mu=datNorm[samp,1])
  lower<-ci$conf.int[1]
  upper<-ci$conf.int[2]
  est<-ci$estimate
  return(data.frame(p, lower, upper, est, tm=datNorm[samp,1]))
})

