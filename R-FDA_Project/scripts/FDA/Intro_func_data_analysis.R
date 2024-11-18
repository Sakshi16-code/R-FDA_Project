# Load required libraries
library(fda)
library(tidyverse)
library(MASS) # for ginv()

# Set seed for reproducibility
set.seed(999)
n_obs <- 80
time_span <- 100
time <- sort(runif(n_obs, 0, time_span))
Wiener <- cumsum(rnorm(n_obs)) / sqrt(n_obs)
y_obs <- Wiener + rnorm(n_obs, 0, .05)

# Setting up the basis
times_basis = seq(0, time_span, 1)
knots = c(seq(0, time_span, 5))
n_knots = length(knots)
n_order = 4
n_basis = length(knots) + n_order - 2
basis = create.bspline.basis(c(min(times_basis), max(times_basis)), n_basis, n_order, knots)
n_basis

PHI = eval.basis(time, basis)
dim(PHI)

# Plot the basis functions and locations of the knots
jpeg('~/path/R-FDA_Project/outputs/basis_vs_time_plot.jpg')
matplot(time, PHI, type = 'l', lwd = 1, lty = 1, xlab = 'time', ylab = 'basis', cex.lab = 1, cex.axis = 1)
for (i in 1:n_knots) {
  abline(v = knots[i], lty = 2, lwd = 1)
}
dev.off()

# Estimating the Basis Coefficients
M = ginv(t(PHI) %*% PHI) %*% t(PHI)
c_hat = M %*% Wiener
y_hat = PHI %*% c_hat

# Augment data frame for plotting
df <- data.frame(time = seq_along(Wiener), Wiener = Wiener, y_obs = y_obs)
df$y_hat <- as.numeric(y_hat)

# Plotting least squares estimate
p2 <- df %>%
  ggplot() +
  geom_line(aes(x = time, y = Wiener), col = "grey") +
  geom_point(aes(x = time, y = y_obs)) +
  geom_line(aes(x = time, y = y_hat), col = "red") +
  ggtitle("Original curve and least squares estimate") +
  xlab("time") + ylab("f(time)")

# Save the plot
ggsave('~/path/R-FDA_Project/outputs/least_sq_est.jpg', plot = p2)

# Estimate the variance of noise
SSE <- t(y_hat - y_obs) %*% (y_hat - y_obs)
sigma2 <- SSE / (n_obs - n_basis)

# Estimate the variance of the fitted curve
H = PHI %*% M
varYhat = diag(H %*% H * matrix(sigma2, n_obs, n_obs))

# 95% confidence interval
y_hat025 = y_hat - 1.96 * sqrt(varYhat)
y_hat975 = y_hat + 1.96 * sqrt(varYhat)

# Augment data frame for plotting confidence intervals
df <- mutate(df, y_hat025 = y_hat025, y_hat975 = y_hat975)

# Plotting curve with error bars
p3 <- df %>%
  ggplot() +
  geom_line(aes(x = time, y = Wiener), col = "grey") +
  geom_point(aes(x = time, y = y_obs)) +
  geom_line(aes(x = time, y = y_hat), col = "red") +
  geom_line(aes(x = time, y_hat025), col = "green") +
  geom_line(aes(x = time, y_hat975), col = "green") +
  ggtitle("Estimated curve with error bars") +
  xlab("time") + ylab("f(time)")

# Save the plot
ggsave('~/path/R-FDA_Project/outputs/curve_with_error_bars.jpg', plot = p3)

# Smoothing using FDA package
Wiener_obj <- smooth.basis(argvals = time, y = y_obs, fdParobj = basis)

# Comparing smooth.basis() with our hand calculated curve
jpeg('~/path/outputs/comparison_plot.jpg')
plot(time, Wiener, type = "l", xlab = "time", ylab = "f(time)", main = "Comparison of fda package and naive smoothing estimates", col = "grey")
lines(time, y_hat, type = "l", col = "red")
lines(Wiener_obj, lwd = 1, col = "blue")
dev.off()
