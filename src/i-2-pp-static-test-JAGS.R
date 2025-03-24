# A JAGS model to estimate parameters of the pp static model.
# testing the model on simulated data

library(rjags)
library(MCMCvis)

nt <- 1 # number of time steps
# simulate the data
source("src/i-1-simulate-pp.R")

# Data
max.c <- 10 # maximum colonies (for data augmentation)
n.s <- tapply(sim.dat[,"time"], sim.dat[,"time"], length) # n surveys in each time period

data.list <- list(
  s.0 = sim.dat[, c("x", "y")], # survey site locations (matrix of X and Y coordinates)
  sur.lev.var = sim.dat[, "det.var"], # detection covariate (vector of length II)
  obs.i = sim.dat[, "obs"], # presence/absence observations
  u.0 = 300, # as data for now
  x.min = 0,
  y.min = 0, # possible spatial extent of the species across all time, bounding box
  x.max = 20000,
  y.max = 20000, # in metres
  II = nrow(sim.dat), # total number of surveys
  M = max.c # maximum number of colonies (data-augmentation approach)
)

# initials
init.list <- list(
  psi = 0.1,
  alpha.sig = log(100),
  beta.sig = log(1),
  Z = rep(1, max.c))

# parameters to monitor
params <- c("psi",
            "alpha.sig",
            "beta.sig",
            "loc")

# mcmc settings
nb <- 5000
ni <- 2000
nc <- 3

# the model
a <- jags.model(file = "src/model-files/pp-static-simtest-JAGS.txt", 
                data = data.list, 
                inits = init.list,
                n.chains = nc)

update(a, n.iter = nb) # burn in
b<-coda.samples(a, 
                variable.names = params, 
                n.iter = ni, 
                thin = 5)

gelman.diag(b, multivariate = FALSE)


summary(b)

## make some plots
### Raw data
pdf(file = "out/sim-static-map.pdf")
  plot(y~x, col = (obs*time+1), data = sim.dat)
  points(c.0.y~c.0.x, col = (gen+1), data = c.mat, pch = 12, cex =2)
  legend('topleft', pch = 21, col = 1:(nt+1), legend = c("absent", 1:nt), title = "Time step")
dev.off()

### Non-loc parameters
temp <- MCMCchains(b, params = c("psi",
                                 "alpha.sig",
                                 "beta.sig"))
pdf(file = "out/sim-static-pairs-parameters.pdf")
  pairs(temp)
dev.off()

### density map
temp <- MCMCchains(b, params = c("loc"))
x <- as.vector(temp[,1:(ncol(temp)/2)])
y <- as.vector(temp[,(ncol(temp)/2+1):ncol(temp)])
point.data <- data.frame(x = x, y = y) |> subset(x!=0)
# Estimate 2D kernel density
density_est <- MASS::kde2d(point.data$x, point.data$y, n = 100, lims = c(x.min, x.max, y.min, y.max))  # 100x100 grid
# Convert to a data frame for ggplot
density_df <- data.frame(
  x = rep(density_est$x, each = length(density_est$y)),
  y = rep(density_est$y, times = length(density_est$x)),
  density = as.vector(density_est$z)
)

# Plot the density with ggplot2
sim.dat <- as.data.frame(sim.dat)
library(ggplot2)
ggplot(density_df, aes(x=y, y=x, fill = density)) +
  geom_raster(interpolate = TRUE) +
  geom_contour(aes(z = density), color = "black", alpha = 0.5) +
  scale_fill_viridis_c(aes(fill = density), alpha = 0.4) +
  geom_point(data = c.mat, aes(x = c.0.x, y = c.0.y), shape = 12, colour = "red", size = 6) +
  geom_point(data = sim.dat, aes(x = x, y = y), colour = (sim.dat$obs*sim.dat$time+1)) +
  theme_minimal() +
  labs(title = "2D Kernel Density Estimation",
       x = "X Coordinate",
       y = "Y Coordinate",
       fill = "Density")
ggsave("out/sim-static-density.map.pdf")
