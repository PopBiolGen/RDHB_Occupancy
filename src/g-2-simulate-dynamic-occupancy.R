# Simulate some data to test g-dynamic-spatial-occupancy.R
source("src/g-1-dynamic-spatial-occupancy_functions.R")

# define dimensions for data objects (as per spOccupancy data setup)
TT <- 5 # number of primary time periods
JJ <- 10  # number of sites
KK <- 10 # maximum number of replicates at a site

# some parameters
k <- 0.2
col.int <- -10
col.space <- 10
ext.int <- -7
ext.b <- 0.2
det.int <- 0
det.b1 <- 1
det.b2 <- 0.2
init.occ <- 1
p.obs <- 0.5 # average proportion of maximum observations per primary period

p.list <- list(k = k, 
               col.b = rbind(col.int, col.space), 
               ext.b = rbind(ext.int, ext.b), 
               det.b = rbind(det.int, det.b1, det.b2))

# matrix for occupancy at each primary time step
occ.mat <- matrix(0, ncol = TT+1, nrow = JJ)
occ.mat[init.occ, 1] <- 1 # initialise with a single colonised cell

# distance matrix assuming 1-D habitat in which each sampling point is 1 unit apart
dist.mat <- as.matrix(dist(1:JJ, diag = TRUE, upper = TRUE))

ext.vars <- cbind(rep(1, JJ), rnorm(JJ)) # some site-level covariates of colonisation and extinction
col.vars <- as.matrix(rep(1, JJ), ncol = 1)

o.state <- dyn.occ(occ.mat, nsteps = TT, dm = dist.mat, ext.vars = ext.vars, col.vars = col.vars, pars = p.list)

# make an array to take observations
obs.array <- array(NA, dim = c(JJ, TT, KK))
obs.vars <- list(x1 = array(rnorm(JJ*TT*KK), dim = c(JJ, TT, KK)), x2 = array(rnorm(JJ*TT*KK), dim = c(JJ, TT, KK)))

for (jj in 1:JJ){
  for (tt in 1:TT){
    occ <- rbinom(1, 1, o.state[jj, tt+1]) #occupied, or not?
    k.realised <- rbinom(1, KK, p.obs) # random number of samples in this time.site
    if (k.realised == 0) next
    det.vars <- do.call("cbind", lapply(obs.vars, FUN = function(x){x[jj, tt,]}))
    det.vars <- det.vars[1:k.realised, ]
    det.vars <- cbind(rep(1, nrow(det.vars)), det.vars)
    pred.det <- det.vars %*% p.list$det.b # predictor of log-odds of detection
    prob.det <- plogis(pred.det)
    obs.array[jj, tt, 1:k.realised] <- rbinom(n = k.realised, size = 1, prob = occ*prob.det)
  }
}

