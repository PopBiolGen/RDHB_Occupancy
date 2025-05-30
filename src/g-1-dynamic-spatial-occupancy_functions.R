##### functions #####
# spatial decay in occupancy with distance, d
decay <- function(d, k){
  exp(-k*d)
}

# mean occupancy arising from sum across all sites with weighting for distance
mean.occ <- function(occ, k, dm){
  if (length(occ) != ncol(dm)) stop("Error in length of occupancy vector passed")
  wght <- decay(dm, k)
  wght <- sweep(wght, MARGIN = 2, STATS = occ, FUN = "*")
  ld <- apply(X = wght, MARGIN = 1, FUN = "sum")
  ld
}

# change occupancy probability according to parameters
update.occ.prob <- function(occ.prob, dm, ext.vars, col.vars, pars, sims){
  col.vars <- cbind(col.vars, mean.occ(occ.prob, pars$k, dm)) # get weighted occupancy
  if (nrow(pars$col.b) != ncol(col.vars)) stop("Dimension mismatch between variables and parameters")
  if (nrow(pars$ext.b) != ncol(ext.vars)) stop("Dimension mismatch between variables and parameters")
  pred.col <- col.vars %*% pars$col.b # predictor of log-odds of colonisation
  pred.ext <- ext.vars %*% pars$ext.b # predictor of log-odds of extinction
  if (!sims) {
    (1-occ.prob)*ilogit(pred.col) + occ.prob*(1-ilogit(pred.ext)) # marginalised occupancy probability
  } else {
    (1-occ.prob)*plogis(pred.col) + occ.prob*(1-plogis(pred.ext))
  }
}

# run over t time primary steps
# occ.prob is JJ * (TT+1) matrix with first column giving initial occupancy at each site.
dyn.occ <- function(init.site, nsteps, dm, ext.vars, col.vars, pars, sims = TRUE){
  if (sims){
    occ.prob <- matrix(0, nrow = nrow(col.vars), ncol = nsteps+1) # make occupancy matrix
  } else {
    occ.prob <- zeros(nrow(col.vars), nsteps+1) # greta version of same  
  }
  occ.prob[init.site, 1] <- 1 # initialise
  for (tt in 1:nsteps){
    occ.prob[,tt+1] <- update.occ.prob(occ.prob = occ.prob[,tt], dm, ext.vars, col.vars, pars, sims)
  }
  occ.prob[, -1] # remove initial column to return JJ * TT matrix
}

# function to return detection probability (P(det | occupied))
det.prob <- function(obs.vars, pars){
  if (nrow(pars$det) != ncol(obs.vars)) stop("Dimension mismatch between variables and parameters")
  pred.det <- obs.vars %*% pars$det.b
  ilogit(pred.det) # detection probability from covariates
}