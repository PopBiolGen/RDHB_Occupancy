# Get the data
current.data <- "RDHBSurveillance_2024-09-09.xlsx"
source("src/b-data-organisation.R")

##### Organise data #####
# data is a list containing data necessary for model fitting. Valid tags are y, occ.covs, det.covs, coords, and grid.index.

# aggregate data (time and space), drop unsampled grid cells
agg_data <- df %>% aggregate_data()
agg_data$df <- filter(agg_data$df, !is.na(agg_data$df$cell.id))
agg_data$df_grid <- filter(agg_data$df_grid, agg_data$df_grid$cell.id %in% agg_data$df$cell.id) 

# drop geometry
agg_data.ng <- lapply(agg_data, st_drop_geometry)

# select data to use
data_select <- select(agg_data.ng$df, 
                      cell.id, 
                      time.step, 
                      date.time, 
                      presence, 
                      water,
                      flowering,
                      dist.0, 
                      hive.removed) %>%
  arrange(cell.id, date.time) %>%
  group_by(cell.id, time.step) %>% # 
  mutate(obs = paste0("obs_", row_number())) %>% # J <- length(unique(data_select$obs))
  ungroup() %>%
  select(time.step, cell.id, date.time, obs, presence, water, flowering, dist.0, hive.removed) 

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

# observed presence/absence at time within primary period (for indexing in JAGS)
n.obs.jj.tt <- apply(y, MARGIN = c(1, 2), FUN = function(x){sum(!is.na(x))})

# col.var is a J x nVars matrix of variables included in the colonisation portion of the model. 
# Each list element is a different colonisation covariate, which is recorded at site level. 
# Site-level covariates are specified as a vector of length J

# col.var is currently not used

# ext.var is a list of variables included in the extinction portion of the model. 
# Each list element is a different extinction covariate, which is at site-primary time period level. 
# site/primary time period level covariates are specified as a matrix 
#with rows corresponding to sites and columns correspond to primary time periods.
ev <- data_select %>%
  group_by(cell.id, time.step) %>%
  summarise(hive.removed = sum(hive.removed, na.rm = TRUE)) %>%
  ungroup()


# makes an empty list given a set of variable names
make_var_list <- function(varnames){
  vl <- vector(mode = "list", length = length(varnames))
  names(vl) <- varnames
  vl
}

ext.var <- make_var_list("hive.removed")

# organises variables names in v.list into a list each with JJ x TT x KK matrix of data
extract_vars <- function(df, v.list, obs.level = FALSE){
  for (vv in names(v.list)){
    if (obs.level){
      da <- array(dim = c(JJ, TT, KK)) # empty array to take one covariate
    } else {
      da <- array(dim = c(JJ, TT)) # empty array to take one covariate
    }
    for (jj in 1:JJ){ # for each cell
      temp <- filter(df, cell.id == cell[jj]) %>%
        select(time.step, {vv})
      for (tt in 1:TT){ # each primary time period
        temp.vec <- temp[[vv]][temp$time.step == tt]
        if (length(temp.vec) == 0) next
        if (obs.level){
          da[jj, tt, 1:length(temp.vec)] <- temp.vec
        } else {
          da[jj, tt] <- temp.vec
        }
        
      }
    }
    v.list[[vv]] <- da
  }
  v.list
}

ext.var <- extract_vars(ev, ext.var)
rm(ev)

# Similarly, det.covs is a list of variables included in the detection portion of the model, 
# with each list element corresponding to an individual variable. 
# In addition to site-level and/or site/primary time period-level, 
# detection covariates can also be observational-level. Observation-level covariates are specified 
# as a three-dimensional array with first dimension corresponding to sites, 
# second dimension corresponding to primary time period, and third dimension corresponding to replicate.

dc <- data_select %>%
  select(cell.id, time.step, date.time, water, flowering) %>%
  mutate(doy = yday(date.time), # day of year in radians
         doy.rad = doy/365*2*pi,
         sin.doy = sin(doy.rad), # sine and cosine for mean date in a cell/time.step
         cos.doy = cos(doy.rad))

det.var <- make_var_list(c("sin.doy", "cos.doy", "water", "flowering")) 

det.var <- extract_vars(dc, det.var, obs.level = TRUE)


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

