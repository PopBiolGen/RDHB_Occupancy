# Simulate some data to test g-dynamic-spatial-occupancy.R
source("src/g-1-dynamic-spatial-occupancy_functions.R")

# define dimensions for data objects (as per spOccupancy data setup)
TT <- 3 # number of primary time periods
JJ <- 10  # number of sites
KK <- 10 # maximum number of replicates at a site

# some parameters
k <- 0.7
col.int <- -3
col.space <- 1 
ext.int <- -3
ext.b <- 0.2
init.occ <- 1

p.list <- list(k = k, col.b = rbind(col.int, col.space), ext.b = rbind(ext.int, ext.b))

# matrix for occupancy at each primary time step
occ.mat <- matrix(0, ncol = TT+1, nrow = JJ)
occ.mat[init.occ, 1] <- 1 # initialise with a single colonised cell

# distance matrix assuming 1-D habitat in which each sampling point is 1 unit apart
dist.mat <- as.matrix(dist(1:JJ, diag = TRUE, upper = TRUE))

ext.vars <- cbind(rep(1, JJ), rnorm(JJ)) # some site-level covariates of colonisation and extinction
col.vars <- as.matrix(rep(1, JJ), ncol = 1)

dyn.occ(occ.mat, nsteps = TT, dm = dist.mat, ext.vars = ext.vars, col.vars = col.vars, pars = p.list)
