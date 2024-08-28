# Here we fit a dynamic occupancy model in which the colonisation probability
# is a function of number of occupancies in surrounding cells (i.e., spatially explicit)

# use data structures of SpOccupancy, which are considerably nicer than those of unmarked.
#library(greta.dynamics)

###### Priors ######
col.b <- normal(mean = 0, sd = 2, dim = c(2,1)) # colonisation coefficients (intercept and immigration)
ext.b <- normal(mean = 0, sd = 2, dim = c(2,1)) # extinction coefficients (intercept and control effort)

k <- uniform(0.5, 1) # rate of spatial decay

# parameter list
p.list  <- list(col.b = col.b,
                ext.b = ext.b,
                k = k)


##### Data #####
# Set up a 1D distance matrix
nsites <- 10
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
occ.dat <- matrix(sample(0:1, size = nt*nsites, replace = TRUE), ncol = nt) |> as_data()

##### Dynamic model of occupancy over time #####
occ.state <- dyn.occ(init.prob = c(1, rep(0, nsites-1)), nsteps = nt, dm = dm, ext.vars = ev, col.vars = cv, pars = p.list)

distribution(occ.dat) <- bernoulli(occ.state[,-1]) 

##### Some plots and diagnostics #####
o.mod <- model(occ.state)

plot(o.mod)

sims <- calculate(occ.state, nsim = 5)
sims$occ.state[1,,]

## Currently we have a working model in which occupancy state is observed directly.  Next layer is to add detection within primary session

