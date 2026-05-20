
set.seed(23)
lambda <- 2

d <- rpois(1, lambda)
d <- 0:8

data.frame(
	lower=ppois(q=d+0.5, lambda = lambda)
	, upper=1-ppois(q=d-0.5, lambda = lambda)
	, one=1
)

