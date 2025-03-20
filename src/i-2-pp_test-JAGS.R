# A JAGS model to estimate parameters of the pp model.
# testing the model on simulated data

library(rjags)
library(MCMCvis)

# simulate the data
source("src/i-1-simulate-pp.R")

# Data
max.c <- 50 # maximum colonies (for data augmentation)
n.s <- tapply(sim.dat[,"time"], sim.dat[,"time"], length) # n surveys in each time period

data.list <- list(
  nt = max(sim.dat[, "time"]), # number of time steps, T in model description
  s.0 = sim.dat[, c("x", "y")], # survey site locations (matrix of X and Y coordinates)
  sur.lev.var = sim.dat[, "det.var"], # detection covariate (vector of length II)
  t.id = sim.dat[, "time"], # time step
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
  alpha.det = 0,
  beta.det = 1,
  Z = matrix(1, nrow = max.c, ncol = max(sim.dat[, "time"])),
  v = rep(1, nrow(sim.dat)),
  sigma.u = 500,
  sigma.d = 900)

# parameters to monitor
params <- c("psi",
            "alpha.det",
            "beta.det",
            "sigma.u",
            "sigma.d",
            "lambda")

# mcmc settings
nb <- 5000
ni <- 2000
nc <- 3

# the model
a <- jags.model(file = "src/model-files/pp-JAGS.txt", 
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
