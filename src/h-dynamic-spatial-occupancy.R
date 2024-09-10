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
data.list <- list(ext.x1 = ext.x1,
                  det.x1 = det.x1,
                  det.x2 = det.x2,
                  det.x3 = det.x3,
                  det.x4 = det.x4,
                  n.obs.jj.tt = n.obs.jj.tt,
                  obs = y,
                  JJ = JJ,
                  TT = TT,
                  dist.mat = dist.mat)

# initials
occ.init <- apply(y, MARGIN = c(1, 2), FUN = sum, na.rm = TRUE) > 0
occ.init <- occ.init+0

init.list <- list(y = occ.init)

