library(shellpipes)
startGraphics()

library(ggplot2); theme_set(theme_bw())

dat <- rdsRead()

print(ggplot(dat)
	+ aes(lambda)
	+ scale_x_log10()
	+ geom_line(aes(y=cvm))
	+ geom_line(aes(y=expect), color="blue")
)
