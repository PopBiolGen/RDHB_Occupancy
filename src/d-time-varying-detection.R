# script to run occupancy model on RDHB data allowing time-varying detection
# Assume occupancy varies by grid cell but is is fixed for entire period of data collection
# Model allows detection to vary across "seasons" as well as within day
# In this case, "seasons" equate to around 2 months of data.

# Get the data
source("src/b-data-organisation.R")

# additional libraries required
library(unmarked)

# aggregate data (time and space), drop unsampled grid cells
agg_data <- df %>% aggregate_data()
agg_data$df_grid <- filter(agg_data$df_grid, !is.na(mean.prop)) 

# drop geometry
agg_data.ng <- lapply(agg_data, st_drop_geometry)


##### make occu inputs ####
# select data to use
data_select <- select(agg_data.ng$df, cell.id, time.step, date.time, presence, hour, hour2, water, dist.0, hive.removed) %>%
              mutate(time.step2 = time.step^2) %>% # this is a fudge compared to a proper 1st- or 2nd-order Fourier function; exploration
              arrange(cell.id, date.time) %>%
              group_by(cell.id) %>%
              mutate(obs = paste0("obs_", row_number())) %>%
              ungroup() %>%
              select(-date.time) 

# make a site by time observation matrix
obs_matrix <- data_select %>%
              select(cell.id, presence, obs) %>% 
              group_by(cell.id) %>%
              tidyr::pivot_wider(names_from = obs, values_from = presence) %>%
              ungroup() %>%
              select(-cell.id) %>%
              as.matrix()

# make a site by n_covariates dataframe
site_covs <- data_select %>%
              select(cell.id, dist.0, presence, hive.removed) %>% # site covariates
              group_by(cell.id) %>%
              summarise(mean.dist = mean(dist.0, na.rm = TRUE),
                        mean.prop = mean(presence, na.rm = TRUE),
                        n.hive.removed = sum(hive.removed),
                        mean.dist2 = mean.dist^2) %>%
              ungroup() %>%
              mutate(cell.id = factor(cell.id)) # make cell.id a factor
            

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
            select(cell.id, obs, time.step, time.step2, hour, hour2, water) %>%
            make_obs_covs_list(cov.names = c("time.step", "time.step2", "hour", "hour2", "water"))

# create an unmarked frame
umf.tvd <- unmarkedFrameOccu(y = obs_matrix, siteCovs = site_covs, obsCovs = obs_covs)

# fit a model
fit.tvd <- occu(~ 1 + time.step + time.step2 + water + hour + hour2
            ~ 1 + mean.dist, 
            data = umf.tvd)
summary(fit.tvd)