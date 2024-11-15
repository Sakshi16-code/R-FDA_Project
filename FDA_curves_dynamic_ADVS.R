# Function to install and load packages if missing
install_if_missing <- function(package) {
  if (!require(package, character.only = TRUE)) {
    install.packages(package, dependencies = TRUE)
    library(package, character.only = TRUE)
  }
}

# List of packages required
packages <- c("tidyverse", "fda", "plotly", "admiral", "haven", "pharmaverseadam")

# Install each package in the list
lapply(packages, install_if_missing)

data("advs")
unique(advs$PARAM)
unique(advs$PARAMCD)
unique(advs$AVISIT)
# Filter for systolic blood pressure measurements
advs_filtered <- advs %>%
  filter(PARAMCD == "DIABP", 
         AVISIT %in% c("Baseline", "Week 2", "Week 4", "Week 6", 
                                                           "Week 8", "Week 12", "End of Treatment", 
                                                           "Week 16", "Week 20", "Week 24", "Week 26")) %>%
  select(USUBJID, AVISIT, AVAL, TRT01A) %>%
  mutate(AVISIT = factor(AVISIT, levels = c("Baseline", "Week 2", "Week 4", "Week 6", 
                                            "Week 8", "Week 12", "End of Treatment", 
                                            "Week 16", "Week 20", "Week 24", "Week 26")))


duplicates <- advs_filtered %>%
  group_by(USUBJID, TRT01A, AVISIT) %>%
  summarise(n = n(), .groups = 'drop') %>%
  filter(n > 1)
print(duplicates)

advs_filtered_summary <- advs_filtered %>% 
  group_by(USUBJID, TRT01A, AVISIT) %>% 
  summarise(AVAL = mean(AVAL, na.rm = TRUE), .groups = 'drop')
advs_filtered_summary

advs_wide <- advs_filtered %>%
  pivot_wider(names_from = AVISIT, values_from = AVAL, values_fill = NA) 
advs_wide <- advs_wide %>%  
  drop_na()

#Extracting a matrix 
n_curves <- nrow(advs_wide)
n_obs <- ncol(advs_wide)
n_obs
print(paste("Number of Curves:", n_curves)) # Should be the number of unique subjects
print(paste("Number of Observations per Curve:", n_obs))

curves_mat <- advs_wide %>% select(-c(USUBJID, TRT01A)) %>% as.matrix()
treatment <- as.numeric(factor(advs_wide$TRT01A))
# Create a data frame for plotting
curves <- data.frame(
  index = 1:n_curves,
  treatment = treatment,
  curves_mat
)
head(curves)
print(dim(curves_mat))
time <- 1:n_obs
head(time)
# Pivot the data for plotting
curves_long <- curves %>%
  pivot_longer(cols = -c(index, treatment), names_to = "Week", values_to = "value") %>%
  mutate(time = rep(time, times = n_curves), .before = "treatment") %>%
  dplyr::select(-Week)

# Plot the ADVS Curves
p_curves <- ggplot(curves_long, aes(x = time, y = value, color = as.factor(treatment))) +
  geom_line(aes(group = index)) +
  scale_color_manual(values = c("navy blue", "dark grey")) +
  ggtitle("ADVS Curves by Treatment") +
  labs(x = "Time (Weeks)", y = "Diastolic BP") +
  theme_minimal()
print(p_curves)

# B-Spline Basis Setup
knots <- seq(0, n_obs, by = 1)
n_order <- 4
n_basis <- length(knots) + n_order - 2
spline_basis <- create.bspline.basis(rangeval = c(0, n_obs), nbasis = n_basis, norder = n_order)

# Create Functional Data Objects for Each Treatment Group
treatment1_idx <- which(treatment == 1)
treatment2_idx <- which(treatment == 2)

df1 <- curves_mat[treatment1_idx, ]
df2 <- curves_mat[treatment2_idx, ]

# Functional Data Objects
df1_fd <- Data2fd(argvals = 1:n_obs, y = t(df1), basisobj = spline_basis, lambda = 0.5)
df2_fd <- Data2fd(argvals = 1:n_obs, y = t(df2), basisobj = spline_basis, lambda = 0.5)

# Evaluate Functional Data for Position, Velocity, and Acceleration
tfine <- seq(0, n_obs, by = 0.1)
pos1 <- eval.fd(tfine, df1_fd)
vel1 <- eval.fd(tfine, df1_fd, 1)
acc1 <- eval.fd(tfine, df1_fd, 2)

pos2 <- eval.fd(tfine, df2_fd)
vel2 <- eval.fd(tfine, df2_fd, 1)
acc2 <- eval.fd(tfine, df2_fd, 2)

# Data Preparation for Derivative Plots
time <- rep(tfine, each = n_obs)
id1 <- rep(1:nrow(df1), each = length(tfine))
id2 <- rep(1:nrow(df2), each = length(tfine))

deriv1 <- data.frame(time = rep(tfine, nrow(df1)), id1, pos1, vel1, acc1)
deriv2 <- data.frame(time = rep(tfine, nrow(df2)), id2, pos2, vel2, acc2)

# Plot Velocity and Acceleration for Treatment Groups
p_vel_acc_1 <- ggplot(deriv1, aes(x = time, y = vel1, col = factor(id1))) +
  geom_line() + ggtitle("Velocity for Treatment 1")

p_vel_acc_2 <- ggplot(deriv2, aes(x = time, y = vel2, col = factor(id2))) +
  geom_line() + ggtitle("Velocity for Treatment 2")

grid.arrange(p_vel_acc_1, p_vel_acc_2, ncol = 2)

# Phase Space Analysis
phase1 <- deriv1 %>% ggplot(aes(vel1, acc1, col = factor(id1))) + 
  geom_point() + ggtitle("Phase Space for Treatment 1")

phase2 <- deriv2 %>% ggplot(aes(vel2, acc2, col = factor(id2))) + 
  geom_point() + ggtitle("Phase Space for Treatment 2")

grid.arrange(phase1, phase2, ncol = 2)

# Animation for Phase Space Trajectory
phase_dat <- data.frame(tfine = tfine, vel = vel2, acc = acc2)
p_anim <- ggplot(phase_dat, aes(vel, acc)) +
  geom_point() + ggtitle("Trajectory in Phase Space")

animation <- p_anim + transition_time(tfine) + shadow_mark()
anim_save("advs_animation.gif", animation)
animation




# Create a separate dataset with unique subjects
#advs_long <- advs_filtered %>% 
#  group_by(USUBJID, AVISIT) %>%
#  summarise(n = n()) %>%
#  filter(n > 1)

#advs_long
# Convert data to a wide format to create a matrix for FDA
e



# Remove columns with missing values
advs_matrix <- advs_wide %>%
  select(-USUBJID) %>%
  na.omit() %>%
  as.matrix()
advs_matrix

advs
advs_filtered
advs_long
advs_wide
advs_matrix

# Define the time points (days) and subjects
time_points <- unique(advs_long$ADY)
n_basis <- 15 # Number of basis functions
basis <- create.bspline.basis(rangeval = c(min(time_points), max(time_points)), nbasis = n_basis)

# Convert the wide data into a functional data object
advs_fd <- Data2fd(argvals = time_points, y = t(advs_matrix), basisobj = basis)

# Plot the raw curves
plot(advs_fd, xlab = "Day", ylab = "Systolic Blood Pressure", main = "FDA of Systolic Blood Pressure")


# Conduct fPCA on the functional data object
fPCA_result <- pca.fd(advs_fd, nharm = 3)

# Plot the mean function and the first few principal components
plot(fPCA_result$harmonics, xlab = "Day", ylab = "Harmonics", main = "fPCA Components")

# Explained variance
explained_variance <- fPCA_result$varprop * 100
print(paste("Variance explained by each component:", explained_variance))

# Calculate the covariance surface
days <- seq(min(time_points), max(time_points), by = 1)
cov_fd <- var.fd(advs_fd)
cov_matrix <- eval.bifd(days, days, cov_fd)

# 3D surface plot using Plotly
plot_ly(x = days, y = days, z = ~cov_matrix) %>%
  add_surface(contours = list(
    z = list(show = TRUE, usecolormap = TRUE, highlightcolor = "#ff0000", project = list(z = TRUE))
  )) %>%
  layout(scene = list(camera = list(eye = list(x = 1.87, y = 0.88, z = -0.64))),
         title = "Variance-Covariance Surface of Systolic Blood Pressure")

# Save FDA plots
ggsave("systolic_bp_fda.png")

# Save the functional data object and fPCA results
saveRDS(advs_fd, file = "advs_fd.rds")
saveRDS(fPCA_result, file = "fPCA_result.rds")
