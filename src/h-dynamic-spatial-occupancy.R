# model to run dynamic spatial occupancy model on RDHB data

library(rjags)

# get the data and organise it
source("src/h-1-dynamic-spatial-occupancy_data-organisation.R")

# arrange the data for JAGS
