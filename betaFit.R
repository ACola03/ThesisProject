library(bbmle)

set.seed(99)
dat <- runif(2000)

nll <- function(mean, shape, x) {
	a <- 2*mean / shape
	b <- 2*(1 - mean) / shape
	-sum(
		dbeta(x, shape1 = a, shape2 = b, log = TRUE)
	)
}

fit <- mle2(
	nll,
	start = list(mean = 0.5, shape = 1),
	data = list(x = dat),
	method = "L-BFGS-B",
	lower = c(mean = 1e-4, shape = 1e-4),
	upper = c(mean = 1 - 1e-4, shape = Inf)
)

summary(fit)

pfit <- profile(fit)
confint(pfit)
plot(pfit)
