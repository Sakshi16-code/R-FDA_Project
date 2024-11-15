install_if_missing <- function(package) {
  if (!require(package, character.only = TRUE)) {
    install.packages(package, dependencies = TRUE)
    library(package, character.only = TRUE)
  }
}
packages <- c("fda", "tidyverse", "admiral", "fdaoutlier", "gganimate", "gridExtra")
lapply(packages, install_if_missing)

library(fdaoutlier)
library(tidyverse)
library(fda)
library(gganimate)
library(gridExtra)

set.seed(99999)
n_curves <- 100
n_obs <- 50
mod4 <- simulation_model4(n = n_curves, p = n_obs, outlier_rate = .5, seed = 50, plot = FALSE)
index <- 1:n_curves
index1 <- mod4$true_outliers
# curves_mat is an n_curves x n_obs matrix
curves_mat <- mod4$data
treatment <- rep(2,n_obs)
curves <- data.frame(index, treatment, curves_mat)
curves <- curves %>% mutate(treatment = if_else((index %in% index1),1,2))

time <- 1:n_obs

curves_l <- pivot_longer(curves, cols = !c("index", "treatment"), names_to = "Xval") %>%
  mutate(time = rep(time,100), .before = "treatment", treatment = as.factor(treatment)) %>%
  dplyr::select(-Xval)

p <- curves_l %>% ggplot(aes(time,value, color = treatment)) +
  geom_line(aes(group = index)) + 
  scale_color_manual(values=c("navy blue", "dark grey")) +
  ggtitle("Model 4 Curves")
p

knots    = c(seq(0,n_obs,5)) #Location of knots
n_knots   = length(knots) #Number of knots
n_order   = 4 # order of basis functions: for cubic b-splines: order = 3 + 1
n_basis   = length(knots) + n_order - 2;
spline_basis = create.bspline.basis(rangeval = c(0,n_obs), nbasis = n_basis, norder = n_order)
#plot(spline_basis)

# df1 is an (n_curves/2) x (n_obs) matrix 
df1 <- curves_mat[index1,] # data for treatment 1
index2 <- index[!(index %in% index1)]
df2 <- curves_mat[index2,]  # data for treatment 2
# Use the b-spline basis to create represent the curves as vectors in the function space
df1_obj <- Data2fd(argvals = 1:n_obs, y = t(df1), basisobj = spline_basis, lambda = 0.5)
df2_obj <- Data2fd(argvals = 1:n_obs, y = t(df2), basisobj = spline_basis, lambda = 0.5)

tfine <- seq(0,50,by=.5)
# Each matrix is 101 x 50 rows are different times, columns are curves
pos1 <- as.vector(eval.fd(tfine, df1_obj));     pos2 <- as.vector(eval.fd(tfine, df2_obj))
vel1 <- as.vector(eval.fd(tfine, df1_obj,1));   vel2 <- as.vector(eval.fd(tfine, df2_obj,1)) 
acc1 <- as.vector(eval.fd(tfine, df1_obj,2));   acc2 <- as.vector(eval.fd(tfine, df2_obj,2)) 

time <- rep(tfine,50)
id1 <- rep(1:50,each=101)
id2 <- rep(51:100,each=101)

derv1 <- data.frame(time, id1, pos1, vel1, acc1)
derv2 <- data.frame(time, id2, pos2, vel2, acc2)

pv1 <- derv1 %>% ggplot(aes(time,vel1,col=id1)) + geom_line(aes(group = id1)) + ggtitle("Velocity Treatment 1")
pv2 <- derv2 %>% ggplot(aes(time,vel2,col=id2)) + geom_line(aes(group = id2)) + ggtitle("Velocity Treatment 2")
pa1 <- derv1 %>% ggplot(aes(time,acc1,col=id1)) + geom_line(aes(group = id1)) +ggtitle("Acceleration Treatment 1")
pa2 <- derv2 %>% ggplot(aes(time,acc2,col=id2)) + geom_line(aes(group = id2)) + ggtitle("Acceleration Treatment 2")

grid.arrange(pv1, pa1, pv2, pa2, nrow = 2,ncol = 2, padding = unit(1, "line"))

dfdt1 <- deriv.fd(df1_obj,1)
dfdt2 <- deriv.fd(df2_obj,1)
tres <- tperm.fd(dfdt1,dfdt2,plotres=FALSE)

max_q <- tres$qval
tres_dat <-tibble(time = tres$argvals, t_values = tres$Tvals,
                  q_vals = tres$qvals.pts)

p <- tres_dat %>% ggplot(aes(time,t_values, colour = "t_value")) + geom_line() +
  geom_line(aes(time, q_vals, colour = "q_vals")) +
  geom_line(aes(time, max_q,colour = "max_q"), linetype= "dashed") + 
  labs(x = "time", y = "") +
  ggtitle("Statistics for Pointwise t-test")
p

phase1 <- derv1 %>% filter(id1 %in% 15:17)
phase2 <- derv2 %>% filter(id2 %in% 51:53)
pph1 <- phase1 %>% ggplot(aes(vel1,acc1,col=id1)) + geom_point() + ggtitle("Treatment 1 v vs. acc")
pph2 <- phase2 %>% ggplot(aes(vel2,acc2,col=id2)) + geom_point() + ggtitle("Treatment 2 v vs. acc")

grid.arrange(pph1, pph2, ncol = 2, padding = unit(1, "line"))

pos <- eval.fd(tfine, df2_obj[1])
vel <- eval.fd(tfine, df2_obj[1],1)
acc <- eval.fd(tfine, df2_obj[1],2)
phase_dat <- tibble(tfine, vel, acc)

p <- ggplot(phase_dat, aes(vel,acc)) +  geom_point() + ggtitle("Trajectory in Phase Space")

anim <- p + transition_time(tfine) + shadow_mark()
anim_save("anim.gif", anim)
anim
