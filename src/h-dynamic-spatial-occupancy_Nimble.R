# model to run dynamic spatial occupancy model on RDHB data - Using Nimble

library(nimble)

# get the data and organise it
source("src/h-1-dynamic-spatial-occupancy_data-organisation.R")

# arrange the data for Nimble
# given a list, output named objects for JAGS
assign_list <- function(ls, prefix){
  for (ll in 1:length(ls)){
    assign(paste0(prefix, ".x", ll), ls[[ll]], envir = .GlobalEnv)
  }
}

# assign out variables with standard names
assign_list(det.var, "det")
assign_list(ext.var, "ext")

# takes an observation level array (spOccupancy format) and returns it as a vector
obs_level_to_vector <- function(obs.level.array){
  as.vector(obs.level.array)
}

# takes a list of observation level arrays (spOccupancy format) and returns them as a dataframe
obs_level_to_df <- function(obs.level.array.list, na.rm = TRUE){
  site.time.rep <- dim(obs.level.array.list[[1]]) # get dimensions
  var.list <- lapply(obs.level.array.list, FUN = obs_level_to_vector) # throw everything into a dataframe
  indices <- expand.grid(1:site.time.rep[1], 1:site.time.rep[2], 1:site.time.rep[3]) # make indices for site, time, rep
  site.index <- indices[,1] # make a site index
  time.index <- indices[,2] # make a time index
  var.df <- as.data.frame(var.list)
  var.df <- data.frame(site = site.index, time = time.index, var.df)
  if (na.rm) var.df <- na.omit(var.df)
  var.df
}

# take obs level stuff in spOccupancy format and send it back to simple vector format
obs.lev.list <- list(obs = y, det.x1 = det.x1,
                     det.x2 = det.x2,
                     det.x3 = det.x3,
                     det.x4 = det.x4)
obs.df <- obs_level_to_df(obs.lev.list)


# organise JAGS code into Nimble Code
dso.code <- nimbleCode(
  {
    # priors
    rho.int ~ dnorm(0, 1/2) # initial occupancy, slightly regularising
    #rho.b ~ dunif(-4, 0)
    col.int ~ dnorm(-5, 1/4)
    #col.b ~ dnorm(0, 1/9) # colonisation coefficients (intercept and immigration)
    ext.int ~ dnorm(-5, 1/9)
    ext.b ~ dnorm(0, 1/4) # extinction coefficients (intercept and control effort)
    det.int ~ dnorm(0, 1/4)
    det.b1 ~ dnorm(0, 1/4)
    det.b2 ~ dnorm(0, 1/4)
    det.b3 ~ dnorm(0, 1/4)
    det.b4 ~ dnorm(0, 1/4)
    k ~ dunif(0.1, 1) # rate of spatial decay
    
    ### process model
    ## initial occupied cells (immediately prior to observations)
    for (jj in 1:JJ){
      logit(rho[jj]) <- rho.int #+ rho.b*init.dist[jj]
      occ[jj, 1] ~ dbern(rho[jj])
    }
    
    # spatial decay
    for (jj in 1:JJ){
      for (ii in 1:JJ){
        decay.jj.ii[jj, ii] <- exp(-k*dist.mat[jj,ii])
      }
    }
    
    ## state variable over time
    # dynamics prior to and after observation
    for (tt in 2:(TT+1)){
      for (jj in 1:JJ){
        sp.jj.tt[jj, tt] <- sum(decay.jj.ii[jj, 1:JJ]*occ[1:JJ, tt-1]) # spatial colonisation effect
        logit(p.col[jj, tt]) <- col.int + sp.jj.tt[jj, tt]# *col.b
        logit(p.ext[jj, tt]) <- ext.int + ext.b*ext.x1[jj, tt-1]
        ext.col[jj, tt] <- (1-occ[jj,tt-1])*p.col[jj, tt] + 
          occ[jj, tt-1]*(1-p.ext[jj, tt])
        occ[jj,tt] ~ dbern(ext.col[jj, tt])
      }
    }
    
    # observation model
    # only calculate likelihoods where there is data
    logit(p.obs[1:n.obs]) <- det.int + 
        det.b1*det.x1[1:n.obs] + 
        det.b2*det.x2[1:n.obs] +
        det.b3*det.x3[1:n.obs] +
        det.b4*det.x4[1:n.obs]
    for (nn in 1:n.obs){
      obs[nn] ~ dbern(p.obs[nn]*occ[site[nn], time[nn]])
    }
    
    
    # derived variables
    # occupancy over time
    for (tt in 1:TT){
      o.t[tt] <- sum(occ[1:JJ,tt])
    }
  }  
)


# Data list
constant.list <- list(init.dist = init.dist,
                  JJ = JJ,
                  TT = TT,
                  dist.mat = dist.mat,
                  n.obs = nrow(obs.df),
                  site = obs.df$site,
                  time = obs.df$time)
data.list <- list(obs = obs.df$obs,
                  ext.x1 = ext.x1,
                  det.x1 = obs.df$det.x1,
                  det.x2 = obs.df$det.x2,
                  det.x3 = obs.df$det.x3,
                  det.x4 = obs.df$det.x4)

# initials
occ.init <- apply(y, MARGIN = c(1, 2), FUN = sum, na.rm = TRUE) > 0
occ.init <- occ.init+0

init.list <- list(occ = cbind(rep(1, JJ),occ.init))


# parameters to monitor
params <- c("rho.int",
                   #"rho.b",
                   "det.int", 
                   "det.b1", 
                   "det.b2",
                   "det.b3",
                   "det.b4",
                   "col.int",
                   #"col.b",
                   "ext.int",
                   "ext.b",
                   "k",
                   "o.t")

# mcmc settings
nb <- 10000
ni <- 12000
nc = 3

# the model
a.n <- nimbleMCMC(code = dso.code, 
                init = init.list, 
                monitors = params, 
                constants = constant.list,
                data = data.list,
                niter = ni, 
                nburnin = nb, 
                nchains = nc,
                check = FALSE,
                samplesAsCodaMCMC = TRUE)

coda::gelman.diag(a.n)
summary(a.n)

save(a.n, file = "out/dynamic-spatial-occupancy_RDHB_Nimble_Coda.RData")



