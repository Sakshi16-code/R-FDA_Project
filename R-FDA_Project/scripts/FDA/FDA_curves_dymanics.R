# Function to install and load packages if not already installed
install_if_missing <- function(package) {
  if (!require(package, character.only = TRUE)) {
    install.packages(package, dependencies = TRUE)
    library(package, character.only = TRUE)
  }
}

# List of required packages
packages <- c("fda", "tidyverse", "admiral", "fdaoutlier", "gganimate", "gridExtra")
lapply(packages, install_if_missing)

# Load the necessary libraries
library(fdaoutlier)
library(tidyverse)
library(fda)
library(gganimate)
library(gridExtra)

# Redirect console output to a text file
sink("output_log.txt")

# Generate simulated data using model 4
set.seed(99999)
n_curves <- 100
n_obs <- 50
mod4 <- simulation_model4(n = n_curves, p = n_obs, outlier_rate = .5, seed = 50, plot = FALSE)
print("Simulation data generated")

index <- 1:n_curves
index1 <- mod4$true_outliers

# Print the indices of outliers
cat("Outlier indices:", index1, "\n")

# Save the outliers indices to a CSV
write.csv(index1, "outliers_indices.csv", row.names = FALSE)

# Extract curves matrix
curves_mat <- mod4$data
print("Curves matrix extracted")

# Add treatment groups
treatment <- rep(2, n_obs)
curves <- data.frame(index, treatment, curves_mat)
curves <- curves %>% mutate(treatment = if_else((index %in% index1), 1, 2))
print("Treatment groups assigned")

# Reshape data for plotting
time <- 1:n_obs
curves_l <- pivot_longer(curves, cols = !c("index", "treatment"), names_to = "Xval") %>%
  mutate(time = rep(time, 100), .before = "treatment", treatment = as.factor(treatment)) %>%
  dplyr::select(-Xval)
print("Data reshaped for plotting")

# Plot original curves
p <- curves_l %>% ggplot(aes(time, value, color = treatment)) +
  geom_line(aes(group = index)) + 
  scale_color_manual(values = c("navy blue", "dark grey")) +
  ggtitle("Model 4 Curves")
print("Plotting Model 4 Curves")
print(p)
ggsave("model4_curves.png", plot = p, width = 8, height = 6)

# Create B-spline basis functions
knots <- seq(0, n_obs, 5)
n_knots <- length(knots)
n_order <- 4
n_basis <- length(knots) + n_order - 2
spline_basis <- create.bspline.basis(rangeval = c(0, n_obs), nbasis = n_basis, norder = n_order)
print("B-spline basis created")
png("bspline_basis.png")
plot(spline_basis)
dev.off()

# Split data into two treatment groups
df1 <- curves_mat[index1,] # treatment 1
index2 <- index[!(index %in% index1)]
df2 <- curves_mat[index2,] # treatment 2
print("Data split into two treatment groups")

# Represent curves in the function space using B-spline basis
df1_obj <- Data2fd(argvals = 1:n_obs, y = t(df1), basisobj = spline_basis, lambda = 0.5)
df2_obj <- Data2fd(argvals = 1:n_obs, y = t(df2), basisobj = spline_basis, lambda = 0.5)
print("Data converted to functional data objects")

# Evaluate derivatives
tfine <- seq(0, 50, by = 0.5)
pos1 <- as.vector(eval.fd(tfine, df1_obj))
pos2 <- as.vector(eval.fd(tfine, df2_obj))
vel1 <- as.vector(eval.fd(tfine, df1_obj, 1))
vel2 <- as.vector(eval.fd(tfine, df2_obj, 1))
acc1 <- as.vector(eval.fd(tfine, df1_obj, 2))
acc2 <- as.vector(eval.fd(tfine, df2_obj, 2))
print("Derivatives evaluated")

# Save derivatives data to CSV
derv1 <- data.frame(time = rep(tfine, 50), id = rep(1:50, each = 101), pos1, vel1, acc1)
derv2 <- data.frame(time = rep(tfine, 50), id = rep(51:100, each = 101), pos2, vel2, acc2)
write.csv(derv1, "derivatives_treatment1.csv", row.names = FALSE)
write.csv(derv2, "derivatives_treatment2.csv", row.names = FALSE)

# Plot velocity and acceleration
pv1 <- derv1 %>% ggplot(aes(time, vel1, col = id)) + geom_line() + ggtitle("Velocity Treatment 1")
pv2 <- derv2 %>% ggplot(aes(time, vel2, col = id)) + geom_line() + ggtitle("Velocity Treatment 2")
pa1 <- derv1 %>% ggplot(aes(time, acc1, col = id)) + geom_line() + ggtitle("Acceleration Treatment 1")
pa2 <- derv2 %>% ggplot(aes(time, acc2, col = id)) + geom_line() + ggtitle("Acceleration Treatment 2")
print("Velocity and acceleration plots")
grid.arrange(pv1, pa1, pv2, pa2, nrow = 2, ncol = 2)
ggsave("velocity_acceleration_plots.png", arrangeGrob(pv1, pa1, pv2, pa2), width = 12, height = 10)

# Perform two-sample t-test on derivatives
dfdt1 <- deriv.fd(df1_obj, 1)
dfdt2 <- deriv.fd(df2_obj, 1)
tres <- tperm.fd(dfdt1, dfdt2, plotres = FALSE)
print("Two-sample t-test performed")


# Plot t-test results
tres_dat <- tibble(time = tres$argvals, t_values = tres$Tvals, q_vals = tres$qvals.pts)
max_q <- tres$qval
ttest <- tres_dat %>% ggplot(aes(time, t_values, colour = "t_value")) + geom_line() +
  geom_line(aes(time, q_vals, colour = "q_vals")) +
  geom_line(aes(time, max_q, colour = "max_q"), linetype = "dashed") +
  labs(x = "time", y = "") +
  ggtitle("Statistics for Pointwise t-test")
print("Plotting t-test statistics")
print(ttest)
ggsave("t_test_statistics.png", plot = ttest, width = 8, height = 6)

# Create phase space animation
pos <- eval.fd(tfine, df2_obj[1])
vel <- eval.fd(tfine, df2_obj[1], 1)
acc <- eval.fd(tfine, df2_obj[1], 2)
phase_dat <- tibble(tfine, vel, acc)
phase_plot <- ggplot(phase_dat, aes(vel, acc)) + geom_point() + ggtitle("Trajectory in Phase Space")
anim <- phase_plot + transition_time(tfine) + shadow_mark()
print("Creating phase space animation")
anim_save("phase_animation.gif", anim)

# End the sink
sink()

cat("All outputs saved successfully.\n")
