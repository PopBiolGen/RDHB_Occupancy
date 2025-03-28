# To fit static point process model to response data
source("src/j-pp-static-get-data.R")

##### Organise into data for pp-static model #####
max.c <- 20
# data list for JAGS
data.list <- list(
  fd = fd$distance, # foraging distances
  n.fd = nrow(fd), # number of foraging distance data points
  s.0 = as.matrix(df.mr[, c("X", "Y")]), # survey site locations (matrix of X and Y coordinates)
  sur.lev.var.1 = df.mr$flowering, # detection covariate 1
  sur.lev.var.2 = df.mr$water, # detection covariate 2
  obs.i = df.mr$presence, # presence/absence observations
  x.min = min(df.mr$X),
  y.min = min(df.mr$Y), # possible spatial extent of the species across all time, bounding box
  x.max = max(df.mr$X),
  y.max = max(df.mr$Y), # in metres
  II = nrow(df.mr), # total number of surveys
  M = max.c # maximum number of colonies (data-augmentation approach)
)

##### Organise other pieces for JAGS #####

# initials
init.list <- list(
  psi = 0.1,
  sigma.det = 100,
  alpha.u = 0,
  beta.1.u = 0,
  beta.2.u = 0,
  c.j0 = matrix(c(seq(min(df.mr$X), max(df.mr$X), length.out = max.c),
           seq(min(df.mr$Y), max(df.mr$Y), length.out = max.c)),
           nrow = max.c),
  Z = rep(1, data.list$M))

# parameters to monitor
params <- c("psi",
            "alpha.u",
            "beta.1.u",
            "beta.2.u",
            "sigma.det",
            "loc")

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

gelman.diag(b, multivariate = FALSE)


summary(b)

save(b, file = "out/temp-coda.RData")