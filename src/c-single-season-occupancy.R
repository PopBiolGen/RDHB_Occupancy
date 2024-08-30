# script to run single-season occupancy model on RDHB data

# Get the data
source("src/b-data-organisation.R")

# additional libraries required
library(unmarked)

# choose only the latest time interval to assume a single-season occupancy, drop unsampled grid cells
agg_data <- df %>% temporal_aggregation() %>%
          filter(time.step == max(time.step, na.rm = TRUE)) %>%
          aggregate_data()
agg_data$df_grid <- filter(agg_data$df_grid, !is.na(mean.prop)) 

# make a map
z <- map_point_grid(agg_data$df, agg_data$df_grid, summ.col = mean.prop)
z

# drop geometry
agg_data <- lapply(agg_data, st_drop_geometry)


##### make occu inputs ####
# select data to use
data_select <- select(agg_data$df, cell.id, date.time, pres, hour, hour2, water, dist_0, hive.removed) %>%
              arrange(cell.id, date.time) %>%
              group_by(cell.id) %>%
              mutate(obs = paste0("obs_", row_number())) %>%
              ungroup() %>%
              select(-date.time) 

# make a site by time observation matrix
obs_matrix <- data_select %>%
              select(cell.id, obs, pres) %>%
              tidyr::pivot_wider(names_from = obs, values_from = pres) %>%
              select(-cell.id) %>%
              st_drop_geometry() %>%
              as.matrix()

# make a site by n_covariates dataframe
site_covs <- data_select %>%
              select(cell.id, dist_0, pres, hive.removed) %>% # site covariates
              group_by(cell.id) %>%
              summarise(mean.dist = mean(dist_0, na.rm = TRUE),
                        mean.prop = mean(pres, na.rm = TRUE),
                        n.hive.removed = sum(hive.removed),
                        mean.dist2 = mean.dist^2) 
# append a distance from grid cell with highest prevalence
#max_prevalence <- filter(site_covs, mean.prop == max(mean.prop))
#site_covs$dist_prev <- st_distance(site_covs, max_prevalence) 

#site_covs <- st_drop_geometry(site_covs) %>%
#              as.data.frame()

# make observation-level covariate list
# to do this
make_obs_covs_list <- function(df, cov.names = c("hour", "hour2", "water")){
  # to make site by time dataframe for a given covariate 
  ind_obs_cov <- function(vec){
    tibble(cell.id = df$cell.id, obs = df$obs, vec = vec) %>%
    tidyr::pivot_wider(names_from = obs, values_from = vec) %>%
    select(-cell.id)
  }
  # applied to all covariates
  out <- lapply(df[, cov.names], ind_obs_cov)
  names(out) <- cov.names
  out
}

obs_covs <- data_select %>%
            st_drop_geometry() %>%
            select(cell.id, obs, hour, hour2, water) %>%
            make_obs_covs_list()

# create an unmarked frame
umf <- unmarkedFrameOccu(y = obs_matrix, siteCovs = site_covs, obsCovs = obs_covs)

fit <- occu(~ 1 + water + hour + hour2
            ~ 1 + mean.dist, 
            data = umf)
