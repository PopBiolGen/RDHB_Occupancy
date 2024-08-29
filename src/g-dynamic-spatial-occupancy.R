# Here we fit a dynamic occupancy model in which the colonisation probability
# is a function of number of occupancies in surrounding cells (i.e., spatially explicit)

# use data structures of SpOccupancy, which are considerably nicer than those of unmarked.
#library(greta.dynamics)

# simulate some data
source("src/g-2-simulate-dynamic-occupancy.R")

###### Priors ######
col.int <- normal(mean = -5, sd = 2)
col.b <- normal(mean = 0, sd = 6) # colonisation coefficients (intercept and immigration)
ext.int <- normal(mean = -5, sd = 2)
ext.b <- normal(mean = 0, sd = 2) # extinction coefficients (intercept and control effort)
det.int <- normal(mean = 0, sd = 2)
det.b1 <- normal(mean = 0, sd = 2)
det.b2 <- normal(mean = 0, sd = 2)
k <- uniform(0.1, 1) # rate of spatial decay

# parameter list
p.list  <- list(col.b = c(col.int, col.b),
                ext.b = c(ext.int, ext.b),
                det.b = c(det.int, det.b1, det.b2),
                k = k)


##### Data #####
# Set up a 1D distance matrix
dist.mat <- as_data(dist.mat)

# initial cell
init.cell <- 1 #initially colonised cell
init.state <- rep(0, JJ)
init.state[init.cell] <- 1
init.state <- as_data(init.state)

# variables for extinction colonisation, and detection
ext.vars <- as_data(ext.vars)
col.vars <- as_data(col.vars)
det.vars <- as_data(det.vars)

# observed presence/absence at time within primary period
n.obs.jj.tt <- apply(obs.array, MARGIN = c(1, 2), FUN = function(x){sum(!is.na(x))})
n.obs.total <- sum(n.obs.jj.tt)
obs.array[is.na(obs.array)] <- 9 # dummy number to replace NAs
obs.array <- as_data(obs.array)

##### Dynamic model of occupancy over time #####
occ.state <- dyn.occ(init.site = init.cell, 
                     nsteps = TT, 
                     dm = dist.mat, 
                     ext.vars = ext.vars, 
                     col.vars = col.vars, 
                     pars = p.list, 
                     sims = FALSE)

##### calculate likelihoods #####
lhood <- zeros(n.obs.total) # to take likelihoods
ii <- 1 # to index lhood
for (jj in 1:JJ){
  for (tt in 1:TT){
    os <- occ.state[jj, tt] # Probability occupied
    n.obs <- n.obs.jj.tt[jj, tt] # number of actual observations
    obs.jj.tt <- obs.array[jj, tt, 1:n.obs] # vector of observations
    dim(obs.jj.tt) <- n.obs
    if (n.obs == 0) next
    det.vars <- do.call("cbind", lapply(obs.vars, FUN = function(x){x[jj, tt,]}))
    det.vars <- det.vars[1:n.obs, ,drop = FALSE]
    det.vars <- cbind(rep(1, nrow(det.vars)), det.vars)
    pred.det <- det.vars %*% p.list$det.b # predictor of log-odds of detection
    prob.det <- ilogit(pred.det)
    ind.var <- !(1 %in% obs.jj.tt) # indicator variable - no detections in time period
    lhood[ii:(ii+n.obs-1)] <- os*(prob.det*obs.jj.tt + (1-prob.det)*(1-obs.jj.tt)) + (1-os)*ind.var # insert likelihood formula
    ii <- ii+n.obs
  }
}

neg.llhood <- -log(lhood)

zt <- rep(0, length(neg.llhood))

distribution(zt) <- poisson(neg.llhood) 

##### Some plots and diagnostics #####
o.mod <- model(occ.state, col.int, col.b, ext.int, ext.b, det.int, det.b1, det.b2, k)

#plot(o.mod)

draws <- mcmc(model = o.mod, n_cores = 2)

summary(draws)