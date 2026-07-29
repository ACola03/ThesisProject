
set.seed(729)

n <- 3
i <- 1:n
lam <- 1
c <- (3-lam)/2
dist <- sort(rbeta(n, c, c))

cvmTarget <- (2*i-1)/(2*n)
expectTarget <- i/(n+1)
