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
            "sigma.u")

# mcmc settings
nb <- 5000
ni <- 2000
nc <- 3

# the model
a <- jags.model(file = "src/model-files/pp-static-JAGS.txt", 
                data = data.list, 
                inits = init.list,
                n.chains = nc)

update(a, n.iter = nb) # burn in
b<-coda.samples(a, 
                variable.names = params, 
                n.iter = ni, 
                thin = 5)

gelman.diag(b)


summary(b)

temp <- MCMCchains(b)
pdf(file = "out/pairs-parameters.pdf")
  pairs(temp)
dev.off()
