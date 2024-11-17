# script to run multi-season occupancy model on RDHB data with spatial and temporal autocorrelation
# Occupancy is aggregated to grid cells, time is aggregated to (approximately) months


library(spOccupancy)
# https://doserlab.com/files/spoccupancy-web/articles/spacetimemodelshtml#data-structure-and-example-data-set

# get the data and organise it
source("src/f-1-multi-season-spatial-occupancy_data-organisation.R")

# fit the model
ms.so.fit <- stPGOcc(occ.formula = ~1, 
                     det.formula = ~1 + water + flowering + sin.doy + cos.doy, 
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