# Ensure required packages from FDA, tidyverse, and admiral are installed and loaded
install_if_missing <- function(package) {
  if (!require(package, character.only = TRUE)) {
    install.packages(package, dependencies = TRUE)
    library(package, character.only = TRUE)
  }
}
packages <- c("fda", "tidyverse", "admiral")
lapply(packages, install_if_missing)

# Function to simulate data
fake_curves <- function(n_curves = 100, n_points = 80, max_time = 100){
  ID <- 1:n_curves
  x <- vector(mode = "list", length = n_curves)
  t <- vector(mode = "list", length = n_curves)
  
  for (i in 1:n_curves){
    t[i] <- list(sort(runif(n_points, 0, max_time)))
    x[i] <- list(cumsum(rnorm(n_points)) / sqrt(n_points))
  }
  df <- tibble(ID, t, x)
  names(df) <- c("ID", "Time", "Curve")
  return(df)
}

# Generate the data
set.seed(123)
n_curves <- 40
n_points <- 80
max_time <- 100
df <- fake_curves(n_curves = n_curves, n_points = n_points, max_time = max_time)

# Data preparation for ggplot
df_1 <- df %>% select(!c(ID, Curve)) %>% unnest_longer(Time)
df_2 <- df %>% select(!c(ID, Time)) %>% unnest_longer(Curve)
ID <- sort(rep(1:n_curves, n_points))
df_l <- cbind(ID, df_1, df_2)

# Plot the simulated data with ggplot
p <- ggplot(df_l, aes(x = Time, y = Curve, group = ID, col = as.factor(ID))) +
  geom_line() +
  labs(title = "Simulated Functional Data", x = "Time", y = "Curve Value") +
  theme_minimal()
print(p)

# Save the ggplot as PNG in the current working directory
ggsave(file.path(getwd(), "curve_plot.png"), plot = p, width = 8, height = 5)

# Create B-spline basis
knots <- seq(0, max_time, by = 5)
n_order <- 4
n_basis <- length(knots) + n_order - 2
basis <- create.bspline.basis(rangeval = c(0, max_time), nbasis = n_basis)

# Plot and save the basis plot in the current working directory
png(file.path(getwd(), "basis_plot.png"), width = 800, height = 600)
plot(basis, main = "B-Spline Basis Functions")
dev.off()

# Convert data into functional data object
argvals <- matrix(df_l$Time, nrow = n_points, ncol = n_curves)
y_mat <- matrix(df_l$Curve, nrow = n_points, ncol = n_curves)
W.obj <- Data2fd(argvals = argvals, y = y_mat, basisobj = basis, lambda = 0.5)

# Calculate mean and standard deviation for the functional data
W_mean <- mean.fd(W.obj)
W_sd <- std.fd(W.obj)

# Create objects for upper and lower standard deviations
SE_u <- fd(coef = W_mean$coefs + 1.96 * W_sd$coefs / sqrt(n_curves), basisobj = basis)
SE_l <- fd(coef = W_mean$coefs - 1.96 * W_sd$coefs / sqrt(n_curves), basisobj = basis)

# Plot the smoothed curves with confidence intervals and save in the current working directory
png(file.path(getwd(), "smoothed_curves.png"), width = 800, height = 600)
plot(W.obj, xlab = "Time", ylab = "Value", main = "Smoothed Functional Data", lty = 1)
lines(SE_u, col = "blue", lwd = 2, lty = 2)
lines(SE_l, col = "red", lwd = 2, lty = 2)
lines(W_mean, col = "darkgreen", lwd = 3)
legend("topright", legend = c("Mean", "Upper Bound", "Lower Bound"),
       col = c("darkgreen", "blue", "red"), lty = c(1, 2, 2), lwd = 2)
dev.off()

# Create the covariance surface plot using admiral and plotly functions (if necessary)
days <- seq(0, max_time, by = 2)
cov_W <- var.fd(W.obj)
var_mat <- eval.bifd(days, days, cov_W)

# Save covariance surface plot as PNG using admiral or base plotting functions (not interactive)
png(file.path(getwd(), "covariance_surface.png"), width = 800, height = 600)
image(days, days, var_mat, main = "Covariance Surface")
dev.off()

# Display the saved plots
