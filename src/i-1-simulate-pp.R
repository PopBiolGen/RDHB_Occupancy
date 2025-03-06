# Script to simulate data for proposed point-process model

source("src/i-0-simulate-pp_functions.R")

# some parameters
density <- 1e-7 # colonies per m^2
s.rad <- 3 # spot search effective radius
alpha.det <- 0 # intercept for log-odds of detection
beta.det <- 2 # slope of some linear effect on detection
#alpha.extent <- 50 # intercept for change in extent per time period
#beta.extent <- -10 # effect of control on extent growth per time period
#sd.extent <- 30 # standard deviation of random effect on change in extent per time period
sigma.u <- 500 # parameter affecting spread of the density kernel
u.0 <- 300 # expected number of encounters at 0 distance from a colony
alpha.lambda <- log(1.5) # intercept of growth rate function
beta.lambda <- 0.1 # slope of the growth rate function
sigma.d <- 900

# some setup variables
nt <- 3 # number of time steps
x.min <- y.min <- 0 # possible spatial extent of the species across all time, bounding box
x.max <- y.max <- 20000 # in metres

g0.x.min <- g0.y.min <- 7500 # bounding box for possible location of ground zero
g0.x.max <- g0.y.max <- 12500

survey.density <- 30*density # density of survey points

#mask.raster <- matrix(c(1,0,1,1), nrow = 2) # make a "raster" denoting available habitat (a triangle in this case)
#raster.scale <- x.max/nrow(mask.raster)

# place ground zero
g0.x <- 0.5*x.max # location of centre of invasion (putative origin)
g0.y <- 0.5*y.max

# initialisation for dynamic state variables
r0 <- 0.2*max(c(x.max, y.max)) # radius of invasion extent at time = 0, in m
c.n <- density*x.max*y.max # number of "colonies"
c.0.x <- runif(n = c.n, min = x.min, max = x.max) # place initial colonies at time 0
c.0.y <- runif(n = c.n, min = y.min, max = y.max)
c.0 <- cbind(c.0.x, c.0.y) # matrix of colony locations
z0 <- pairwise_distances(c.0, cbind(g0.x, g0.y)) < r0 # distances from g0 < r0 (i.e. which colonies are real.)
#z0 <- z0 * raster_to_points(c.0, mask.raster, rast.scale = x.max/2)
c.0 <- c.0[z0==1,]
s.n <- survey.density*x.max*y.max # number of surveys

# Time step 1
sim.dat <- generate_data(generate_surveys(), c.0)

# Time step > 1
if (nt > 1){
  c.t <- c.0
  print(c.t)
  lambda.x <- rnorm(nt) # some covariate of lambda
  lambda.t <- exp(alpha.lambda + beta.lambda*lambda.x) #lambda as a function of time step
  for (tt in 2:nt){
     c.t <- pp_growth(c.t, lambda.t[tt], sigma.d = sigma.d)
     print(c.t)
     s.t <- generate_surveys()
     sim.dat <- rbind(sim.dat, generate_data(s.t, c.t, t = tt))
  }  
}

plot(y~x, col = (obs*time+1), data = sim.dat)
legend('topleft', pch = 21, col = 1:(nt+1), legend = c("absent", 1:nt), title = "Time step")
