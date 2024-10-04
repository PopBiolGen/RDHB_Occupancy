# Here we fit a dynamic occupancy model in which the colonisation probability
# is a function of number of occupancies in surrounding cells (i.e., spatially explicit)

# use data structures of SpOccupancy, which are considerably nicer than those of unmarked.

# here we try out JAGS..
library(rjags)

# simulate some data
source("src/g-2-simulate-dynamic-occupancy.R")

ext.var <- ext.vars[,2]
det.x1 <- obs.vars[[1]]
det.x2 <- obs.vars[[2]]


##### Data #####
# initial cell
init.cell <- 1 #initially colonised cell
init.state <- rep(0, JJ)
init.state[init.cell] <- 1


# observed presence/absence at time within primary period
n.obs.jj.tt <- apply(obs.array, MARGIN = c(1, 2), FUN = function(x){sum(!is.na(x))})
n.obs.total <- sum(n.obs.jj.tt)

obs.array.zero.pad <- obs.array
obs.array.real <- !is.na(obs.array.zero.pad) 
obs.array.zero.pad[!obs.array.real] <- 0 # add fake observations to sidestep the NA issue
obs.array.real <- obs.array.real - 1 # setup for use with JAGS step function

# these zeros are obviated in the JAGS script


# Data list
data.list <- list(ext.var = ext.var,
                 det.x1 = det.x1,
                 det.x2 = det.x2,
                 init.state = init.state,
                 obs = obs.array.zero.pad,
                 obs.real = obs.array.real,
                 JJ = JJ,
                 TT = TT,
                 KK = KK,
                 dist.mat = dist.mat)

occ.init <- apply(obs.array, MARGIN = c(1, 2), FUN = sum, na.rm = TRUE) > 0
occ.init <- cbind(rep(NA, JJ), occ.init+0)

init.list <- list(occ = occ.init,
                  det.int = 0,
                  det.b1 = 0,
                  det.b2 = 0,
                  col.int = -5,
                  col.b = 7,
                  ext.int = -5,
                  ext.b = 0)

# JAGS model
a <- jags.model(file = "src/model-files/invasion-occ_JAGS.txt", 
                data = data.list, 
                inits = init.list,
                n.chains = 3)
update(a, n.iter = 5000) # burn in
b<-coda.samples(a, 
                variable.names = c("det.int", 
                                   "det.b1", 
                                   "det.b2",
                                   "col.int",
                                   "col.b",
                                   "ext.int",
                                   "ext.b",
                                   "k"), 
                n.iter = 10000, 
                thin = 5)
gelman.diag(b)

summary(b)

save(b, "out/simtest.RData")
