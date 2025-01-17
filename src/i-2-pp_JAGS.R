# A JAGS model to estimate parameters of the pp model.

library(rjags)
library(MCMCvis)

# simulate the data
source("src/i-1-simulate-pp.R")

# Data
max.c <- ceiling(2*sum(z0)) # maximum colonies (for data augmentation)

data.list <- list(
  s.0 = s.0, # survey site locations (matrix of X and Y coordinates)
  sur.lev.var = sur.lev.var, # detection covariate (vector of length II)
  obs.i = obs.i,
  #nt = 1, # number of time steps
  lambda.0 = 300, # as data for now
  x.min = 0,
  y.min = 0, # possible spatial extent of the species across all time, bounding box
  x.max = 20000,
  y.max = 20000, # in metres
  g0.x.min = 7500, # bounding box for possible location of ground zero
  g0.y.min = 7500,
  g0.x.max = 12500,
  g0.y.max = 12500,
  II = length(s.i0.x), # number of sites
  JJ = max.c # maximum number of colonies (data-augmentation approach)
)

# initials
init.list <- list(
  g0 = matrix(c(8000, 9000), nrow = 1),
  alpha.det = 0,
  beta.det = 1,
  r0 = 5000,
  da = rep(1,max.c),
  pres = rep(1, length(obs.i)),
  sigma = 1000)

# parameters to monitor
params <- c("psi",
            "r0",
            "alpha.det",
            "beta.det",
            "sigma",
            "g0")

# mcmc settings
nb <- 5000
ni <- 10000
nc = 3

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
pairs(temp)
