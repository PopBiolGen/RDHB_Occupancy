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
data_select <- select(agg_data$df, cell.id, date.time, presence, water, dist.0, hive.removed) %>%
              arrange(cell.id, date.time) %>%
              group_by(cell.id) %>%
              mutate(obs = paste0("obs_", row_number())) %>%
              ungroup() %>%
              select(-date.time) 

# make a site by time observation matrix
obs_matrix <- data_select %>%
              select(cell.id, obs, presence) %>%
              tidyr::pivot_wider(names_from = obs, values_from = presence) %>%
              select(-cell.id) %>%
              st_drop_geometry() %>%
              as.matrix()

# make a site by n_covariates dataframe
site_covs <- data_select %>%
              select(cell.id, dist.0, presence, hive.removed) %>% # site covariates
              group_by(cell.id) %>%
              summarise(mean.dist = mean(dist.0, na.rm = TRUE),
                        mean.prop = mean(presence, na.rm = TRUE),
                        n.hive.removed = sum(hive.removed),
                        mean.dist2 = mean.dist^2) 

# make observation-level covariate list
# to do this
make_obs_covs_list <- function(df, cov.names = c("water")){
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
            select(cell.id, obs, water) %>%
            make_obs_covs_list()

# create an unmarked frame
umf <- unmarkedFrameOccu(y = obs_matrix, siteCovs = site_covs, obsCovs = obs_covs)

fit <- occu(~ 1 + water
            ~ 1 + mean.dist, 
            data = umf)
summary(fit)
