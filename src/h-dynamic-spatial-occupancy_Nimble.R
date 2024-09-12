# model to run dynamic spatial occupancy model on RDHB data - Using Nimble

library(nimble)

# get the data and organise it
source("src/h-1-dynamic-spatial-occupancy_data-organisation.R")

# organise JAGS code into Nimble Code
dso.code <- nimbleCode(
  {
    # priors
    rho.int ~ dnorm(0, 1/4) # initial occupancy
    rho.b ~ dnorm(0, 1/4)
    col.int ~ dnorm(-5, 1/4)
    col.b ~ dnorm(0, 1/36) # colonisation coefficients (intercept and immigration)
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
      logit(rho[jj]) <- rho.int + rho.b*init.dist[jj]
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
        logit(p.col[jj, tt]) <- col.int + col.b*sp.jj.tt[jj, tt]
        logit(p.ext[jj, tt]) <- ext.int + ext.b*ext.x1[jj, tt-1]
        ext.col[jj, tt] <- (1-occ[jj,tt-1])*p.col[jj, tt] + 
          occ[jj, tt-1]*(1-p.ext[jj, tt])
        occ[jj,tt] ~ dbern(ext.col[jj, tt])
      }
    }
    
    # observation model
    for (tt in 1:TT){
      for (jj in 1:JJ){
        for (kk in 1:KK){
          logit(p.obs[jj, tt, kk]) <- det.int + 
            det.b1*det.x1[jj, tt, kk] + 
            det.b2*det.x2[jj, tt, kk] +
            det.b3*det.x3[jj, tt, kk] +
            det.b4*det.x4[jj, tt, kk]
          obs[jj, tt, kk] ~ dbern(p.obs[jj, tt, kk]*occ[jj, tt+1]*step(y.real[jj, tt, kk]))
        }
      }
    }
    
    # derived variables
    # occupancy over time
    for (tt in 1:TT){
      o.t[tt] <- sum(occ[1:JJ,tt])
    }
  }  
)

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


# Data list
data.list <- list(init.dist = init.dist,
                  ext.x1 = ext.x1,
                  det.x1 = det.x1,
                  det.x2 = det.x2,
                  det.x3 = det.x3,
                  det.x4 = det.x4,
                  obs = y,
                  y.real = (y.real - 1), # setup for use with JAGS step function
                  JJ = JJ,
                  TT = TT,
                  KK = KK,
                  dist.mat = dist.mat)

# initials
occ.init <- apply(y, MARGIN = c(1, 2), FUN = sum, na.rm = TRUE) > 0
occ.init <- occ.init+0

init.list <- list(occ = cbind(rep(1, JJ),occ.init))


# parameters to monitor
params <- c("rho.int",
                   "rho.b",
                   "det.int", 
                   "det.b1", 
                   "det.b2",
                   "det.b3",
                   "det.b4",
                   "col.int",
                   "col.b",
                   "ext.int",
                   "ext.b",
                   "k",
                   "o.t")

# mcmc settings
nb <- 3000
ni <- 1000
nc = 3

# the model
a.n <- nimbleMCMC(code = dso.code, 
                init = init.list, 
                monitors = params, 
                constants = data.list, 
                niter = ni, 
                nburnin = nb, 
                nchains = nc)

coda::gelman.diag(a.n)
summary(a.n)

save(a.n, file = "out/dynamic-spatial-occupancy_RDHB_Nimble_Coda.RData")



