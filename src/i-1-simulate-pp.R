# Script to simulate data for proposed point-process model

# some useful functions
# Function to compute pairwise distances between two sets of points (in n-dimensions)
pairwise_distances <- function(A, B) {
  # Expand A and B to compute the Euclidean distances
  sqrt(outer(rowSums(A^2), rowSums(B^2), `+`) - 2 * tcrossprod(A, B))
}

# function to tell us if a point is inside or outside a polygon.
# Uses the ray casting algorithm, heading right from the test point.
inside <- function (pt, pg){
  edge.in <- function(pt, pg){ # Is this an edge that we will pass through?
    pg2 <- rbind(pg[-1,], pg[1,])
    b <- (pg2[, "y"] - pg[, "y"])/(pg2[, "x"] - pg[, "x"])  # slope
    a <- pg[, "y"] - b*pg[, "x"] # intercept
    above <- pt[, "y"] > (a+b*pt[, "x"]) # exclude cases on the line
    (above & b > 0) | (!above & b < 0) | (is.infinite(b) & pg[, "x"] > pt[, "x"])
  }
  ei <- edge.in(pt, pg)
  test.y <- pt[,"y"] # y-value of point to be tested
  pg.y <- pg[, "y"] # y-values of polygon vertices
  t1 <- test.y <= pg.y # test is less than a vertex's y value
  t2 <- test.y >= pg.y # test is greater than a vertex's y value
  t2 <- c(t2[-1], t2[1]) # align adjacent vertices
  nc <- sum(((t2+t1)!=1) & ei) # crossings where this is true
  (nc %% 2) == 1
}

# function to extract habitat "raster" to points
raster_to_points <- function(pts, hab.rast, rast.scale){
  # collapse points back to raster indices
  indices <- pts %/% rast.scale + 1
  hab.rast[indices[, 2:1]] # returns a vector of raster values: one for each point
}

# function to calculate probability of at least one being present from distances (from AHM 10.8)
spatial_decay <- function(distance, sigma, lambda.0){
  lambda.ij <- lambda.0*exp(-distance^2/(2*sigma^2)) # poisson encounter rate with distance
  lambda.i <- apply(lambda.ij, 1, sum) # sum over colonies
  prob.pres.i <- 1-exp(-lambda.i) # collapsed to a binary
  prob.pres.i
}

# some parameters
density <- 1e-7 # colonies per m^2
s.rad <- 3 # spot search effective radius
alpha.det <- 0 # intercept for log-odds of detection
beta.det <- 2 # slope of some linear effect on detection
#alpha.extent <- 50 # intercept for change in extent per time period
#beta.extent <- -10 # effect of control on extent growth per time period
#sd.extent <- 30 # standard deviation of random effect on change in extent per time period
sigma <- 500 # parameter affecting spread of the density kernel
lambda.0 <- 300 # expected number of encounters at 0 distance from a colony

# some setup variables
nt <- 1 # number of time steps
x.min <- y.min <- 0 # possible spatial extent of the species across all time, bounding box
x.max <- y.max <- 20000 # in metres

g0.x.min <- g0.y.min <- 7500 # bounding box for possible location of ground zero
g0.x.max <- g0.y.max <- 12500

survey.density <- 30*density # density of survey points

mask.raster <- matrix(c(1,0,1,1), nrow = 2) # make a "raster" denoting available habitat (a triangle in this case)


# place ground zero
g0.x <- 0.5*x.max # location of centre of invasion (putative origin)
g0.y <- 0.5*y.max

# initialisation for dynamic state variables
r0 <- 0.2*max(c(x.max, y.max)) # radius of invasion extent at time = 0, in m
c.n <- density*x.max*y.max # number of "colonies"
c.j0.x <- runif(n = c.n, min = x.min, max = x.max) # place initial colonies at time 0
c.j0.y <- runif(n = c.n, min = y.min, max = y.max)
c.0 <- cbind(c.j0.x, c.j0.y) # matrix of colony locations
z0 <- pairwise_distances(c.0, cbind(g0.x, g0.y)) < r0 # distances from g0 < r0 (i.e. which colonies are real.)
z0 <- z0 * raster_to_points(c.0, mask.raster, rast.scale = x.max/2)
c.0 <- c.0[z0==1,]
s.n <- survey.density*x.max*y.max # number of surveys
s.i0.x <- runif(n = s.n, min = x.min, max = x.max) # place initial surveys at time 0
s.i0.y <- runif(n = s.n, min = y.min, max = y.max)
s.0 <- cbind(s.i0.x, s.i0.y) # matrix of survey locations
s.0 <- s.0[raster_to_points(s.0, mask.raster, rast.scale = x.max/2)==1, ]

# generate data
t.step <- rep(1:nt, each = nrow(s.0)) # vector of time steps
sur.lev.var <- rnorm(length(t.step)) # survey level variable (affects detection)
logit.det <- alpha.det + beta.det*sur.lev.var # linear model on detection
det <- plogis(logit.det) # conditional detection probability
  
d.ij <- pairwise_distances(s.0, c.0) # distances between surveys and colonies

p.pres.i <- spatial_decay(d.ij, sigma, lambda.0)
obs.i <- rbinom(n = length(p.pres.i), size = 1, prob = p.pres.i*det)


sum(obs.i)
plot(s.0[,2]~s.0[,1], col = 1+obs.i)
