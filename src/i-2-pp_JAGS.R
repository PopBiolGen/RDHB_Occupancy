# A JAGS model to estimate parameters of the pp model.

library(rjags)
library(MCMCvis)

# simulate the data
source("src/i-1-simulate-pp.R")

# Data
max.c <- 25 # maximum colonies (for data augmentation)
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
  g0.x.min = 7500, # bounding box for possible location of ground zero
  g0.y.min = 7500,
  g0.x.max = 12500,
  g0.y.max = 12500,
  II = nrow(sim.dat), # total number of surveys
  M = max.c # maximum number of colonies (data-augmentation approach)
)

# initials
init.list <- list(
  g0 = matrix(c(8000, 9000), nrow = 1),
  alpha.det = 0,
  beta.det = 1,
  r0 = 5000,
  Z = matrix(1, nrow = max.c, ncol = max(sim.dat[, "time"])),
  v = rep(1, nrow(sim.dat)),
  sigma.u = 1000,
  sigma.d = 1000)

# parameters to monitor
params <- c("psi",
            "r0",
            "alpha.det",
            "beta.det",
            "sigma.u",
            "sigma.d",
            "lambda",
            "g0")

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
