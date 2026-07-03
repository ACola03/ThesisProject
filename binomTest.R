
binom.test(1, 4, alternative = "less")$p.value

## Calculate a signed, conservative P value (
binomP <- function(x, n, p){
	pl <- binom.test(x, n, p, alternative = "less")$p.value
	pg <- binom.test(x, n, p, alternative = "greater")$p.value
	return(ifelse(pl < pg , pl, -pg))
}

n <- 5; p <- 0.37
v <- 0:n

plist <- numeric()

for (x in v){
	plist[[x+1]] <- (binomP(x, n, p))
}

data.frame(x = v
	, P = plist
	, d = dbinom(v, n, p)
)
