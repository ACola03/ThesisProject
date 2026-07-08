
n <- 1000
reps <- 1000

ls <- numeric()
lvec <- (1:n)/(n+1)

set.seed(708)
for (r in 1:reps){
	pvec <- sort(runif(n))
	ls[[r]] <- sum((pvec-lvec)^2)
	## ls[[r]] <- sum(abs(pvec-lvec))
}

print(mean(ls))
print(sd(ls))
print(range(ls))

silent <- hist(ls)
