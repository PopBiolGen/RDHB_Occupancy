# script to run static occupancy model on RDHB data

# Get the data
source("src/b-data-aggregation.R")

# additional libraries required
library(unmarked)

# choose only the latest time interval to assume a static occupancy, drop unsampled grid cells, drop geometry
df_t_x <- df_grid %>%
          filter(!st_is_empty(geometry)) %>%
          st_drop_geometry() %>%
          filter(time.step == max(time.step, na.rm = TRUE))
          

##### make occu inputs ####
# select data to use
data_select <- select(df_t_x, cell.id, date.time, pres, hour, hour2, water, dist_0) %>%
              arrange(cell.id, date.time) %>%
              mutate(date = paste0("obs_", row_number())) %>%
              select(-date.time) 

# make a site by time observation matrix
obs_matrix <- data_select %>%
              select(cell.id, date, pres) %>%
              tidyr::pivot_wider(names_from = date, values_from = pres) %>%
              select(-cell.id) %>%
              as.matrix()

# make a site by n_covriates dataframe
site_covs <- data_select %>%
              select(cell.id, dist_0) %>% # site covariates
              group_by(cell.id) %>%
              summarise(mean.dist = mean(dist_0, na.rm = TRUE)) %>%
              as.data.frame()

# make observation-level covariate list
# to do this
make_obs_covs_list <- function(df, cov.names = c("hour", "hour2", "water")){
  # to make site by time dataframe for a given covariate 
  ind_obs_cov <- function(vec){
    tibble(cell.id = df$cell.id, date = df$date, vec = vec) %>%
    tidyr::pivot_wider(names_from = date, values_from = vec) %>%
    select(-cell.id)
  }
  # applied to all covariates
  out <- lapply(df[, cov.names], ind_obs_cov)
  names(out) <- cov.names
  out
}

obs_covs <- data_select %>%
            select(cell.id, date, hour, hour2, water) %>%
            make_obs_covs_list()

# create an unmarked frame
umf <- unmarkedFrameOccu(y = obs_matrix, siteCovs = site_covs, obsCovs = obs_covs)

fit <- occu(~ 1 + water + hour + hour2
            ~ 1 + mean.dist, 
            data = umf)
