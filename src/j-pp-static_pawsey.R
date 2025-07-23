
# Survey dataset = 773 days
# Analyse 90-day windows
# Windows starting every 30 days
# = 26 30-day periods

rm(list=ls()) # Clear workspace
args <- commandArgs(trailingOnly = TRUE) # Create command line for interacting with job_submission.slurm script
iter <- as.numeric(args[1]) # iter corresponds to array number (iteration) of job

# To fit static point process model to response data
source("src/j-pp-static-get-data_pawsey.R")

##### Organise into data for pp-static model #####
max.c <- 100
# data list for JAGS
data.list <- list(
  fd = fd$distance, # foraging distances
  n.fd = nrow(fd), # number of foraging distance data points
  s.0 = as.matrix(df.mr[, c("X", "Y")]), # survey site locations (matrix of X and Y coordinates)
  sur.lev.var.1 = df.mr$flowering, # detection covariate 1
  sur.lev.var.2 = df.mr$water, # detection covariate 2
  obs.i = df.mr$presence, # presence/absence observations
  x.min = round(min(df.mr$X)),
  y.min = round(min(df.mr$Y)), # possible spatial extent of the species across all time, bounding box
  x.max = round(max(df.mr$X)),
  y.max = round(max(df.mr$Y)), # in metres
  II = nrow(df.mr), # total number of surveys
  M = max.c, # maximum number of colonies (data-augmentation approach)

    m.shore = m.shore, # matrix of shoreline map 
  
  shore.x.min = as.numeric(colnames(m.shore)[1]), # Min X cell coord in m.shore
  shore.y.min = as.numeric(rownames(m.shore)[1]) # Min Y cell coord

#  n.x = c(rep((1/length(round(min(df.mr$X)):round(max(df.mr$X)))), times=length(round(min(df.mr$X)):round(max(df.mr$X))))),
#  n.y = c(rep((1/length(round(min(df.mr$Y)):round(max(df.mr$Y)))), times=length(round(min(df.mr$Y)):round(max(df.mr$Y)))))
#index
#y.index = as.integer((trunc(c.j0[1,2]/100)*100 - shore.y.min) / 100 + 1),
#x.index = as.integer((trunc(c.j0[1,2]/100)*100 - shore.x.min) / 100 + 1)
)

##### Organise other pieces for JAGS #####

# initials
init.list <- list(
  psi = 0.1,
  sigma.det = 100,
  alpha.u = 0,
  beta.1.u = 0,
  beta.2.u = 0,
  c.j0 = round(matrix(c(seq(min(df.mr$X), max(df.mr$X), length.out = max.c),
           seq(min(df.mr$Y), max(df.mr$Y), length.out = max.c)),
           nrow = max.c)),
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
a <- jags.model(file = "src/model-files/pp-static-JAGS_ocean.txt", 
                data = data.list, 
                inits = init.list,
                n.chains = nc)

update(a, n.iter = nb) # burn in
b<-coda.samples(a, 
                variable.names = params, 
                n.iter = ni, 
                thin = 5)

# gelman.diag(b, multivariate = FALSE)


# summary(b)

save(b, 
     file = sprintf("out/temp-coda-start_%s.RData", iter))


source("src/figures/j-pp-static-figures_pawsey.R")
