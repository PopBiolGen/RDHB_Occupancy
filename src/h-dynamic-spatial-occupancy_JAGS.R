# model to run dynamic spatial occupancy model on RDHB data

library(rjags)

# get the data and organise it
source("src/h-1-dynamic-spatial-occupancy_data-organisation.R")

# arrange the data for JAGS
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

# the model
a <- jags.model(file = "src/model-files/invasion-occ_RDHB_JAGS.txt", 
                data = data.list, 
                inits = init.list,
                n.chains = 3)

update(a, n.iter = 5000) # burn in
b<-coda.samples(a, 
                variable.names = c("rho.int",
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
                                   "o.t"), 
                n.iter = 10000, 
                thin = 5)

gelman.diag(b)


summary(b)

save(a, b, file = "out/dynamic-spatial-occupancy_RDHB_Coda.RData")
