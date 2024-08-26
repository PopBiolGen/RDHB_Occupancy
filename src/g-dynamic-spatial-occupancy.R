# Here we fit a dynamic occupancy model in which the colonisation probability
# is a function of number of occupancies in surrounding cells (i.e., spatially explicit)

# use data structures of SpOccupancy, which are considerably nicer than those of unmarked.

library(greta)
#library(greta.dynamics)

###### Priors ######
col.b <- normal(mean = 0, sd = 2, dim = c(2,1)) # colonisation coefficients (intercept and immigration)
ext.b <- normal(mean = 0, sd = 2, dim = c(2,1)) # extinction coefficients (intercept and control effort)

k <- uniform(0, 1) # rate of spatial decay

# parameter list
p.list  <- list(col.b = col.b,
                ext.b = ext.b,
                k = k)

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

# change occupancy state according to parameters
update.occ.state <- function(state, dm, ext.vars, col.vars, pars){
  col.vars <- cbind(col.vars, mean.occ(state, pars$k, dm)) # get weighted occupancy
  if (nrow(pars$col.b) != ncol(col.vars)) stop("Dimension mismatch between variables and parameters")
  if (nrow(pars$ext.b) != ncol(ext.vars)) stop("Dimension mismatch between variables and parameters")
  pred.col <- col.vars %*% pars$col.b # predictor of log-odds of colonisation
  pred.ext <- ext.vars %*% pars$ext.b # predictor of log-odds of extinction
  colonised <- bernoulli(prob = plogis(pred.col), dim = length(pred.col)) # colonised? define as random draw from a binomial
  extirpated <- bernoulli(prob = plogis(pred.ext), dim = length(pred.ext)) # extirpated? define as random draw from a binomial
  (1-state)*colonised - state*extirpated # updated state
}

# run over t time steps
dyn.occ <- function(init.state, nsteps, dm, ext.vars, col.vars, pars){
  occ.state <- greta_array(dim = c(nrow(ext.vars), nsteps + 1))
  occ.state[, 1] <- init.state
  for (tt in 1:nsteps){
    occ.state[,tt+1] <- update.occ.state(state = occ.state[,tt], dm, ext.vars, col.vars, pars)
  }
  occ.state
}

##### Data #####
# Set up a 1D distance matrix
nsites <- 3
dm <- as.matrix(dist(cbind(rep(0, nsites), 1:nsites), diag = TRUE, upper = TRUE))
dm <- as_data(dm)

# initial cell
init.cell <- 1 #initially colonised cell
init.state <- rep(0, nsites)
init.state[init.cell] <- 1
init.state <- as_data(init.state)

# variables for extinction and colonisation
ev <- data.frame(intercept = rep(1, nsites),
                 effort = rep(0, nsites)) |>
  as_data()
cv <- data.frame(intercept = rep(1, nsites)) |>
  as_data()

# n time steps
nt <- 3
# some made-up data
occ.dat <- matrix(c(1, 0, 0, 1, 1, 0, 1, 1, 1), ncol = nt) |> as_data()

occ.state <- dyn.occ(init.state = c(1, 0, 0), nsteps = nt, dm = dm, ext.vars = ev, col.vars = cv, pars = p.list)

distribution(occ.dat) <- bernoulli(occ.state[,-1]) 

o.mod <- model(occ.state)


