# script to run dynamic occupancy model on RDHB data allowing time-varying occupancy and detection
# Occupancy is aggregated to grid cells, time is aggregated to (approximately) months

# Get the data
source("src/b-data-organisation.R")

# additional libraries required
library(unmarked)

n.times <- 10 # number of time slices to break the data into

# aggregate data (time and space), drop unsampled grid cells
agg_data <- df %>% aggregate_data(n.periods = n.times)
agg_data$df_grid <- filter(agg_data$df_grid, !is.na(mean.prop)) 

# drop geometry
agg_data.ng <- lapply(agg_data, st_drop_geometry)


##### make occu inputs ####
# select data to use and toss into multiframe
data_select <- select(agg_data.ng$df, cell.id, time.step, date.time, pres, water, dist_0, hive.removed) %>%
  mutate(time.step2 = time.step^2) %>% # this is a fudge compared to a proper 1st- or 2nd-order Fourier function; exploration
  arrange(cell.id, date.time) %>%
  group_by(cell.id, time.step) %>% # T <- max(data_select$time.step)
  mutate(obs = paste0("obs_", row_number())) %>% # J <- length(unique(data_select$obs))
  ungroup() %>%
  select(time.step, cell.id, date.time, obs, pres, water, dist_0, hive.removed) 

# make a site by n_covariates dataframe (M rows) M <- nrow(site_covs)
site_covs <- data_select %>%
  select(cell.id, dist_0, hive.removed) %>% # site covariates
  group_by(cell.id) %>%
  summarise(mean.dist = mean(dist_0, na.rm = TRUE), # in kms
            n.hive.removed = sum(hive.removed),
            n.records = n()) %>%
  ungroup() %>%
  mutate(mean.dist2 = mean.dist^2,
         cell.id = factor(cell.id)) # make cell.id a factor

# make a site by time observation matrix (M x TJ matrix, columns in time.step-major observation-minor order)
obs_matrix <- data_select %>%
  select(cell.id, time.step, pres, obs) %>% 
  mutate(time.step = paste0("time_", time.step)) %>%
  tidyr::pivot_wider(names_from = time.step, values_from = pres) %>%
  tidyr::pivot_wider(names_from = obs, values_from = contains("time_")) %>%
  select(-cell.id) %>%
  as.matrix()

# re-order columns... could almost certainly be done more cleanly than this!
order.var <- colnames(obs_matrix)
order.var <- gsub("time_", "", x = order.var)
order.var <- gsub("obs_", "", x = order.var)
order.var <- stringr::str_split(order.var, pattern = "_") %>% data.frame() %>% t()
order.var <- data.frame(time.step = as.numeric(order.var[,1]), obs = as.numeric(order.var[,2]))
ord.vec <- order(order.var$time.step, order.var$obs)
obs_matrix <- obs_matrix[,ord.vec] #now in time-step, observation order
rm(ord.vec, order.var)


# make observation-level covariate list (MTJ rows - site-time.step-observation order)
obs_covs <- data_select %>%
  select(cell.id, time.step, obs, water) %>%
  mutate(obs = factor(obs)) %>% # what follows is an awful hack to get all factor levels into the frame...
  tidyr::pivot_wider(names_from = obs, values_from = water) %>%
  tidyr::pivot_longer(cols = contains("obs_"), cols_vary = "fastest", names_to = "obs", values_to = "water") %>%
  mutate(cell.id = factor(paste0("cell_", cell.id))) %>%
  tidyr::pivot_wider(names_from = cell.id, values_from = water) %>%
  tidyr::pivot_longer(cols = contains("cell_"), cols_vary = "fastest", names_to = "cell.id", values_to = "water") %>%
  mutate(cell.id = as.numeric(gsub("\\D", "", cell.id)), obs = as.numeric(gsub("\\D", "", obs))) %>%
  arrange(cell.id, time.step, obs) %>%
  select(cell.id, time.step, obs, water)

# create an unmarked frame
umf.do <- unmarkedMultFrame(y = obs_matrix, siteCovs = site_covs, obsCovs = obs_covs, numPrimary = n.times)

# fit a dynamic occupancy model
start.values <- c(0, 0, 0, 0, -2, 0)

fit.do <- colext(psiformula = ~1 + mean.dist, 
                 gammaformula = ~1, 
                 epsilonformula = ~1, 
                 pformula = ~1 + water, 
                 data = umf.do, 
                 starts = start.values,
                 method = "BFGS")
summary(fit.do)