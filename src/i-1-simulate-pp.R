# Script to simulate data for proposed point-process model

source("src/i-0-simulate-pp_functions.R")

# some parameters
density <- 1e-7 # colonies per m^2
s.rad <- 3 # spot search effective radius
alpha.sig <- log(500) # log of mean sigma for the density kernel
beta.sig <- log(2)/2 # slope of some linear effect on sigma
u.0 <- 300 # expected number of encounters at 0 distance from a colony
lambda <- 1.5 # constant growth rate
sigma.d <- 900

# some setup variables
x.min <- y.min <- 0 # possible spatial extent of the species across all time, bounding box
x.max <- y.max <- 20000 # in metres

g0.x.min <- g0.y.min <- 7500 # bounding box for possible location of ground zero
g0.x.max <- g0.y.max <- 12500

survey.density <- 15*density # density of survey points

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
c.0 <- c.0[z0==1, , drop = FALSE]
s.n <- survey.density*x.max*y.max # number of surveys

# Time step 1
sim.dat <- generate_data(generate_surveys(), c.0, alpha.sig, beta.sig)
c.mat <- data.frame(cbind(c.0, gen = 1))

# Time step > 1
if (nt > 1){
  c.t <- c.0
  for (tt in 2:nt){
     c.t <- pp_growth(c.t, lambda, sigma.d = sigma.d)
     s.t <- generate_surveys()
     sim.dat <- rbind(sim.dat, generate_data(s.t, c.t, alpha.sig, beta.sig, t = tt))
     c.mat <- rbind(c.mat, cbind(c.t, gen = tt))
  }  
}

print(c.mat)
plot(y~x, col = (obs*time+1), data = sim.dat)
points(c.0.y~c.0.x, col = (gen+1), data = c.mat, pch = 12, cex =2)
legend('topleft', pch = 21, col = 1:(nt+1), legend = c("absent", 1:nt), title = "Time step")
