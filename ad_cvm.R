# ===== ONE SAMPLE PLOTS: CvM & AD

# ----- Data

set.seed(30)
n <- 20
u <- sort(runif(n))
i <- 1:n

# Correct segment coordinates using the exact right-continuous ECDF value (i/n)
segments_df <- data.frame(
  x = u,
  y_cdf = u,          # theoretical y = x
  y_ecdf = i / n,      # exact value of right-continuous ECDF at u[i]
  ad_weight = 1/(u*(1-u))
)

# ----- CvM

plot(c(0, 1), c(0, 1), type = "n", xlab = "u", ylab = "F(u)", 
     main = paste0("Visual CvM Calculation\nECDF Deviations at Order Statistics"))

# Theoretical line
abline(0, 1, col = "gray45", lty = 2, lwd = 2)

# ECDF step function
lines(stepfun(u, (0:n)/n), col = "black", lwd = 2, pch = 1)

# AD vertical line segments
segments(x0 = segments_df$x, y0 = segments_df$y_cdf, 
         x1 = segments_df$x, y1 = segments_df$y_ecdf, 
         col = "grey45", lwd = 1.5)
points(segments_df$x, segments_df$y_ecdf, pch = 16, col = "grey45", cex = 1)


# ----- Anderson-Darling

# Base plot
plot(c(0, 1), c(0, 1), type = "n", xlab = "u", ylab = "F(u)", 
     main = paste0("Visual Anderson-Darling Calculation\nECDF Deviations at Order Statistics"))

# Theoretical line
abline(0, 1, col = "gray40", lty = 2, lwd = 2)

# ECDF step function
lines(stepfun(u, (0:n)/n), col = "black", lwd = 2)

# Scale weights
weights <- segments_df$ad_weight
log_w <- log(weights)
norm_w <- (log_w - min(log_w)) / (max(log_w) - min(log_w))

# Color map by AD weight
pal <- colorRampPalette(rev(c("gray10", "gray45","grey85")))
color_indices <- findInterval(norm_w, seq(0, 1, length.out = 100))
line_colors <- pal(100)[color_indices]

# AD vertical line segments
segments(x0 = segments_df$x, y0 = segments_df$y_cdf, 
         x1 = segments_df$x, y1 = segments_df$y_ecdf, 
         col = line_colors, lwd = 3)
points(segments_df$x, segments_df$y_ecdf, pch = 16, col = line_colors, cex = 1)

# ----- Highlight Differences
















