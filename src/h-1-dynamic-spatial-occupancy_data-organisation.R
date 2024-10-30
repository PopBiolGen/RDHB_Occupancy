# Get the data
current.data <- "RDHBSurveillance_2024-10-09.xlsx"
# these two variables to control data amounts during model development.  To be removed for full fit.
max.obs.jj.tt <- 30 # maximum number of (non hive removal) observations per site/time to allow in the data
max.tt <- 20
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
  mutate(obs.num = row_number(), obs = paste0("obs_", obs.num)) %>% # J <- length(unique(data_select$obs))
  ungroup() %>%
  filter(time.step < max.tt) %>%
  filter(obs.num < max.obs.jj.tt) %>% # set a maximum number of observations per site.time
  select(time.step, cell.id, date.time, obs, presence, water, flowering, dist.0, hive.removed) 

# define dimensions for data objects
TT <- length(unique(data_select$time.step)) # number of primary time periods
JJ <- length(unique(data_select$cell.id)) # number of sites
KK <- length(unique(data_select$obs)) # maximum number of replicates at a site

init.dist <- tapply(data_select$dist.0, INDEX = data_select$cell.id, FUN = mean)

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

# number of observations for each tt, jj (used to trim loops in Nimble)
n.obs.jj.tt <- apply(y, MARGIN = c(1, 2), FUN = function(x){sum(!is.na(x))})

# for handling NAs in the data
y.zero.pad <- y
y.real <- !is.na(y.zero.pad) 
y.zero.pad[!y.real] <- 0 # add fake observations to sidestep the NA issue (these are switched off in JAGS)

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

# organises variables names in v.list into a list each with JJ x TT (optional x KK) array of data
extract_vars <- function(df, v.list, obs.level = FALSE, fill.zeros = FALSE){
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
        if (length(temp.vec) == 0) {
          if (fill.zeros) {
            temp.vec <- 0
          } else next
        }
        if (obs.level){
          if (fill.zeros){
            da[jj, tt, ] <- c(temp.vec, rep(0, KK-length(temp.vec))) # pad with zeros
          } else {
            da[jj, tt, 1:length(temp.vec)] <- temp.vec  # else leave NAs
          }
          
        } else {
          da[jj, tt] <- temp.vec
        }
        
      }
    }
    v.list[[vv]] <- da
  }
  v.list
}

ext.var <- extract_vars(ev, ext.var, fill.zeros = TRUE) # zeros here are actually true zeros
rm(ev)

# Similarly, det.vars is a list of variables included in the detection portion of the model, 
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

det.var <- extract_vars(dc, det.var, obs.level = TRUE, fill.zeros = TRUE) # zeros here are bypassed with step in JAGS


# coords to get dist.mat
coords <- select(agg_data$df_grid, geometry, cell.id) %>%
  st_transform(crs = 3577) %>% # transform to albers
  arrange(cell.id) 

dist.mat <- as.matrix(dist(coords, diag = TRUE, upper = TRUE))
