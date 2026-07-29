library(shellpipes)

set.seed(729)

sim_ss <- function(lambda, samp, reps){
	i <- 1:samp
	c <- (3-lambda)/2
	cvmTarget <- (2*i-1)/(2*samp)
	expectTarget <- i/(samp+1)
	ss <- function(d, t){
		return(sum((d-t)^2))
	}
	cvm_ss <- numeric(reps)
	expect_ss <- numeric(reps)
	for(iter in 1:reps){
		dist <- sort(rbeta(samp, c, c))
		cvm_ss[iter] <- ss(dist, cvmTarget)
		expect_ss[iter] <- ss(dist, expectTarget)
	}
	return(c(lambda
		, cvm=mean(cvm_ss)
		, sig_cvm=sd(cvm_ss)
		, expect=mean(expect_ss)
		, sig_expect=sd(expect_ss)
	))
}

sim_ss_over_lambda <- function(lambda, samp, reps){
	results <- lapply(lambda, function(lam){
		sim_ss(lam, samp, reps)
	})
	df <- as.data.frame(do.call(rbind, results))
	names(df)[[1]] <- "lambda"

	return(df)
}

lseq <- exp(seq(-3, 1, length.out=21))

out <- sim_ss_over_lambda(lseq, 3, 1e4)
print(out)

rdsSave(out)


