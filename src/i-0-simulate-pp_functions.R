# some useful functions for simulating data 

# Function to compute pairwise distances between two sets of points (in n-dimensions)
pairwise_distances <- function(A, B) {
  # Expand A and B to compute the Euclidean distances
  sqrt(outer(rowSums(A^2), rowSums(B^2), `+`) - 2 * tcrossprod(A, B))
}

# function to tell us if a point is inside or outside a polygon.
# Uses the ray casting algorithm, heading right from the test point.
# this function is too much; have switched to raster approach (see raster_to_points)
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
spatial_presence <- function(distance, sigma.u, u.0){
  lambda.ij <- u.0*exp(-distance^2/(2*sigma.u^2)) # poisson encounter rate with distance
  lambda.i <- apply(lambda.ij, 1, sum) # sum over colonies
  prob.pres.i <- 1-exp(-lambda.i) # collapsed to a binary
  prob.pres.i
}

# generates surveys from global variables
generate_surveys <- function(){
  s.i0.x <- runif(n = s.n, min = x.min, max = x.max) # place surveys
  s.i0.y <- runif(n = s.n, min = y.min, max = y.max)
  s.0 <- cbind(s.i0.x, s.i0.y) # matrix of survey locations
  #s.0 <- s.0[raster_to_points(s.0, mask.raster, rast.scale = x.max/2)==1, ] # remove those outside mask
  s.0
}

# generate data
## time step 1
# function that uses global variables to generate data from current s, and c.
generate_data <- function(s.t, c.t, t = 1) {
  sur.lev.var <- rnorm(nrow(s.t)) # survey level variable (affects detection)
  logit.det <- alpha.det + beta.det*sur.lev.var # linear model on detection
  det <- plogis(logit.det) # conditional detection probability
  
  d.ij <- pairwise_distances(s.t, c.t) # distances between surveys and colonies
  
  p.pres.i <- spatial_presence(d.ij, sigma.u, u.0)
  obs.i <- rbinom(n = length(p.pres.i), size = 1, prob = p.pres.i*det) # generate observations
  cbind(time = rep(t, length(p.pres.i)), x = s.t[,1], y = s.t[,2], det.var = sur.lev.var, obs = obs.i) # output
}

# Grows the population according to lambda.t and sigma.d
pp_growth <- function(c.t, lambda.t, sigma.d){
  n.0 <- nrow(c.t)
  n.1 <- rpois(n.0, lambda.t) # realised number of daughter colonies
  inds <- rep(1:n.0, times = n.1)
  c.1 <- c.t[inds,] # new matrix of colonies, pre dispersal
  n.1 <- nrow(c.1)
  c.1[,1] <- rnorm(n = n.1, mean = c.1[,1], sd = sigma.d) # disperse with gaussian kernel
  c.1[,2] <- rnorm(n = n.1, mean = c.1[,2], sd = sigma.d)
  #z0 <- raster_to_points(c.1, mask.raster, rast.scale = x.max/2) #
  #c.1[z0==1,] # remove those outside mask
  c.1
}
