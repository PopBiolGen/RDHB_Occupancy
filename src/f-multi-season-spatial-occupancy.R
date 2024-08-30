# script to run multi-season occupancy model on RDHB data with spatial and temporal autocorrelation
# Occupancy is aggregated to grid cells, time is aggregated to (approximately) months

# Get the data
source("src/b-data-organisation.R")

library(spOccupancy)
# https://doserlab.com/files/spoccupancy-web/articles/spacetimemodelshtml#data-structure-and-example-data-set

##### Organise data #####
# data is a list containing data necessary for model fitting. Valid tags are y, occ.covs, det.covs, coords, and grid.index.

# aggregate data (time and space), drop unsampled grid cells
agg_data <- df %>% aggregate_data()
agg_data$df_grid <- filter(agg_data$df_grid, !is.na(mean.prop)) 

# drop geometry
agg_data.ng <- lapply(agg_data, st_drop_geometry)

# select data to use
data_select <- select(agg_data.ng$df, cell.id, time.step, date.time, presence, water, dist.0, hive.removed) %>%
  mutate(time.step2 = time.step^2) %>% # this is a fudge compared to a proper 1st- or 2nd-order Fourier function; exploration
  arrange(cell.id, date.time) %>%
  group_by(cell.id, time.step) %>% # 
  mutate(obs = paste0("obs_", row_number())) %>% # J <- length(unique(data_select$obs))
  ungroup() %>%
  select(time.step, cell.id, date.time, obs, presence, water, dist.0, hive.removed) 

# define dimensions for data objects
TT <- length(unique(data_select$time.step)) # number of primary time periods
JJ <- length(unique(data_select$cell.id)) # number of sites
KK <- length(unique(data_select$obs)) # maximum number of replicates at a site

# y is a three-dimensional array with first dimension equal to the number of sites (J), 
# second dimension equal to the maximum number of primary time periods (i.e., years or seasons), and 
# third dimension equal to the maximum number of replicates at a given site.
y <- array(dim = c(JJ, TT, KK))

pa <- data_select %>%
  select(cell.id, time.step, presence)

cell <- unique(data_select$cell.id)
 for (jj in 1:JJ){ #for each cell
   temp <- filter(pa, cell.id == cell[jj])
   for (tt in 1:TT){ # for each primary time period
     temp.vec <- temp$presence[temp$time.step == tt]
     if (length(temp.vec) == 0) next
     y[jj, tt, 1:length(temp.vec)] <- temp.vec
   }
 }

rm(pa)

# occ.covs is a list of variables included in the occurrence portion of the model. 
# Each list element is a different occurrence covariate, which can be site level or site/primary time period level. 
# Site-level covariates are specified as a vector of length J while 
# site/primary time period level covariates are specified as a matrix 
#with rows corresponding to sites and columns correspond to primary time periods.

oc <- data_select %>%
  select(cell.id, time.step, dist.0) %>%
  group_by(cell.id, time.step) %>%
  summarise(mean.dist = mean(dist.0, na.rm = TRUE)) %>%
  ungroup()


var.vec <- c("mean.dist")
occ.var <- vector(mode = "list", length = length(var.vec))
names(occ.var) <- var.vec

for (vv in var.vec){ # for each variable
  occ.covs.i <- matrix(nrow = JJ, ncol = TT) # empty matrix to take one covariate
  for (jj in 1:JJ){ # for each cell
    temp <- filter(oc, cell.id == cell[jj]) %>%
      select(time.step, {vv})
    if (length(temp) == 0) next
    #occ.covs.i[jj, temp$time.step] <- as.vector(temp[,vv])
  }
  occ.var[[vv]] <- occ.covs.i
}


# Similarly, det.covs is a list of variables included in the detection portion of the model, 
# with each list element corresponding to an individual variable. 
# In addition to site-level and/or site/primary time period-level, 
# detection covariates can also be observational-level. Observation-level covariates are specified 
# as a three-dimensional array with first dimension corresponding to sites, 
# second dimension corresponding to primary time period, and third dimension corresponding to replicate.

dc <- data_select %>%
  select(cell.id, time.step, date.time, water) %>%
  mutate(doy = yday(date.time), # day of year in radians
         doy.rad = doy/365*2*pi,
         sin.doy = sin(doy.rad), # sine and cosine for mean date in a cell/time.step
         cos.doy = cos(doy.rad))


var.vec <- c("sin.doy", "cos.doy", "water")
det.var <- vector(mode = "list", length = length(var.vec))
names(det.var) <- var.vec

for (vv in var.vec){ # for each variable
  da <- array(dim = c(JJ, TT, KK)) # empty array to take one covariate
  for (jj in 1:JJ){ # for each cell
    temp <- filter(dc, cell.id == cell[jj]) %>%
      select(time.step, {vv})
    for (tt in 1:TT){ # each primary time period
      temp.vec <- temp[[vv]][temp$time.step == tt]
      if (length(temp.vec) == 0) next
      da[jj, tt, 1:length(temp.vec)] <- temp.vec
    }
  }
  det.var[[vv]] <- da
}


# coords is a matrix of the observation coordinates used to estimate the spatial random effect for each site. 
# coords has two columns for the easting and northing coordinate, respectively. 
# Typically, each site in the data set will have it's own coordinate, such that coords is a 
# J×2 matrix and grid.index should not be specified. 
# If you desire to estimate spatial random effects at some larger spatial level, 
# e.g., if points fall within grid cells and you want to estimate a spatial random effect for each grid cell instead of each point,
# coords can be specified as the coordinate for each grid cell. 
# In such a case, grid.index is an indexing vector of length J, 
# where each value of grid.index indicates the corresponding row in coords that the given site corresponds to. 
# Note that spOccupancy assumes coordinates are specified in a projected coordinate system.
coords <- select(agg_data$df_grid, geometry, cell.id) %>%
  st_transform(crs = 3577) %>% # transform to albers
  arrange(cell.id) %>%
  mutate(centroid = st_centroid(geometry),
         easting = st_coordinates(centroid)[,1],
         northing = st_coordinates(centroid)[,2]) %>%
  select(easting, northing) %>%
  st_drop_geometry() %>%
  as.matrix()
  

occ.data <- list(y = y, occ.covs = occ.var, det.covs = det.var, coords = coords)


ms.so.fit <- stPGOcc(occ.formula = ~1, 
                     det.formula = ~1 + water + sin.doy + cos.doy, 
                     data = occ.data, 
                     cov.model = "exponential", 
                     NNGP = TRUE, 
                     n.neighbors = 8, 
                     n.batch = 40, 
                     batch.length = 500,
                     ar1 = TRUE,
                     n.chains = 3)

summary(ms.so.fit)
save(ms.so.fit, file = "out/f-ms-so_fit.RData")
