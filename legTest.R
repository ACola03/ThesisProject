library(shellpipes)
library(magrittr)

sourceFiles()

numPianos <- 100
numReps <- 2000
(numReps
	%>% rbeta(1, 1)
	%>% matrix(ncol=numPianos)
	%>% legendre.coef
	%>% extract(, 1:4)
	%>% raise_to_power(2)
	%>% rowSums()
	|> mean()
)
