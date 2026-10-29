corSignedP <- function(s1, s2, alt){
  test <- cor.test(s1, s2, alternative = alt)
  P <- test$p.value
  if (test$statistic <0) return(P)
  return(1-P)
}

sampSize <- 10

s1 <- rnorm(sampSize)
s2 <- rnorm(sampSize)
corSignedP(s1, s2)

s1 <- rnorm(sampSize)
s2 <- rnorm(sampSize)
corSignedP(s1, s2)

s1 <- rnorm(sampSize)
s2 <- rnorm(sampSize)
corSignedP(s1, s2)

s1 <- rnorm(sampSize)
s2 <- rnorm(sampSize)
cor.test(s1, s2, alternative = "less")
cor.test(s1, s2, alternative = "greater")
cor.test(s1, s2, alternative = "two.sided")