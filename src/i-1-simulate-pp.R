# Script to simulate data for proposed point-process model

# some parameters
density <- 1e-7 # colonies per m^2
#s.rad <- 3 # spot search effective radius
alpha.det <- 0 # intercept for log-odds of detection
beta.det <- 2 # slope of some linear effect on detection
#alpha.extent <- 50 # intercept for change in extent per time period
#beta.extent <- -10 # effect of control on extent growth per time period
#sd.extent <- 30 # standard deviation of random effect on change in extent per time period
k <- 500 # parameter affecting spread of the density kernel

# some setup variables
nt <- 1 # number of time steps
x.min <- y.min <- 0 # possible spatial extent of the species across all time, bounding box
x.max <- y.max <- 20000 # in metres

g0.x.min <- g0.y.min <- 7500 # bounding box for possible location of ground zero
g0.x.max <- g0.y.max <- 12500

survey.density <- 2*density # density of survey points

# place ground zero
g0.x <- 0.5*x.max # location of centre of invasion (putative origin)
g0.y <- 0.5*y.max

# Function to compute pairwise distances between two sets of points (in n-dimensions)
pairwise_distances <- function(A, B) {
  # Expand A and B to compute the Euclidean distances
  sqrt(outer(rowSums(A^2), rowSums(B^2), `+`) - 2 * tcrossprod(A, B))
}

# initialisation for dynamic state variables
r0 <- 0.2*max(c(x.max, y.max)) # radius of invasion extent at time = 0, in m
c.n <- density*x.max*y.max # number of colonies
c.j0.x <- runif(n = c.n, min = x.min, max = x.max) # place initial colonies at time 0
c.j0.y <- runif(n = c.n, min = y.min, max = y.max)
c.0 <- cbind(c.j0.x, c.j0.y) # matrix of colony locations
z0 <- pairwise_distances(c.0, cbind(g0.x, g0.y)) < r0 # distances from g0 < r0 (i.e. which colonies are real.)
s.n <- survey.density*x.max*y.max # number of surveys
s.i0.x <- runif(n = s.n, min = x.min, max = x.max) # place initial surveys at time 0
s.i0.y <- runif(n = s.n, min = y.min, max = y.max)
s.0 <- cbind(s.i0.x, s.i0.y) # matrix of survey locations

# generate data
t.step <- rep(1:nt, each = s.n) # vector of time steps
sur.lev.var <- rnorm(length(t.step)) # survey level variable (affects detection)
logit.det <- alpha.det + beta.det*sur.lev.var # linear model on detection
det <- plogis(logit.det) # conditional detection probability
  
d.ij <- pairwise_distances(s.0, c.0[z0,]) # distances between surveys and colonies

  # function to calculate probability from distances (placeholder for now)
  spatial_decay <- function(distance, k){
    exp(-distance/k)#/(2*pi*distance) # probability of an individual being at the point
  }

p.pres.ij <- spatial_decay(d.ij, k = k)
p.pres.comp <- 1-p.pres.ij # complement of prob presence
p0.i <- apply(p.pres.comp, MARGIN = 1, FUN = prod) # probability of not present from any source
obs.i <- rbinom(n = length(p0.i), size = 1, prob = (1-p0.i)*det)


sum(obs.i)
plot(s.i0.y~s.i0.x, col = 1+obs.i)
