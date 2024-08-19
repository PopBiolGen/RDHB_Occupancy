# to visualise results from the multi-season-spatial-occupancy model
# If using an already fitted model..
source("src/b-data-organisation.R")
# aggregate data (time and space), drop unsampled grid cells
n.times <- 10
agg_data <- df %>% aggregate_data(n.periods = n.times)
agg_data$df_grid <- filter(agg_data$df_grid, !is.na(mean.prop)) 
load(file = "out/f-ms-so_fit.RData")

# else
# source("src/f-multi-season-spatial-occupancy.R")

library(gganimate)
library(spOccupancy)


##### Make an animated map of occupancy over time #####
# estimated occupancy at each site and time.step
est.psi <- apply(ms.so.fit$psi.samples, c(2, 3), mean) 
# join occupancy estimates onto grid
grid.dat <- cbind(agg_data$df_grid, est.psi)

gd.long <- grid.dat %>%
  tidyr::pivot_longer(cols = starts_with("X"), names_to = "time.step", values_to = "Occupancy") %>%
  mutate(time.step = as.numeric(gsub("X", "", time.step)))

# grab a shoreline
shoreline <- st_read(file.path(Sys.getenv("DATA_PATH"), "Spatial/Igis_aus_outline/Australia_boundary.shp"))

# set bounding box to extent of data
bbox <- st_bbox(gd.long)

p <- ggplot() +
  geom_sf(data = shoreline, fill = "lightblue", color = "blue") +
  geom_sf(data = gd.long, aes(fill = Occupancy)) +
  coord_sf(xlim = c(bbox["xmin"], bbox["xmax"]), ylim = c(bbox["ymin"], bbox["ymax"])) +  # Apply bounding box
  scale_fill_viridis_c(option = "plasma") +
  theme_minimal() +
  labs(title = "Time step: {closest_state}")

anim <- p +
  transition_states(time.step, transition_length = 2, state_length = 1) +
  ease_aes("cubic-in-out")

animate(anim, nframes = max(gd.long$time.step), fps = 1)


##### Make map of occupancy in last time step #####
p <- ggplot() +
  geom_sf(data = shoreline, fill = "lightblue", color = "blue") +
  geom_sf(data = filter(gd.long, time.step == max(time.step)), aes(fill = Occupancy)) +
  coord_sf(xlim = c(bbox["xmin"], bbox["xmax"]), ylim = c(bbox["ymin"], bbox["ymax"])) +  # Apply bounding box
  scale_fill_viridis_c(option = "plasma") +
  theme_minimal() +
  labs(title = paste0("Time step: ", max(gd.long$time.step)))

p
ggsave(filename = "out/multi-season-spatial-occupancy.png")


##### Examine how occupancy changes over time #####
est.psi.time <- apply(est.psi, MARGIN = 2, mean)
est.psi.time.se <- apply(est.psi, MARGIN = 2, sd)/sqrt(nrow(est.psi)) # loses parameter uncertainty

p <- ggplot() +
  geom_point(aes(x = 1:n.times, y = est.psi.time)) +
  geom_errorbar(aes(x = 1:n.times, 
                    ymin = est.psi.time-2*est.psi.time.se,
                    ymax = est.psi.time+2*est.psi.time.se),
                width = 0) +
  labs(x = "Time step", y = "Mean occupancy") +
  ylim(0, 0.3) +
  theme_minimal()
p
ggsave(filename = "out/multi-season-spatial-occupancy-mean-occ-over-time.png")


##### Look at how detection changes throughout the year #####
doy <- 1:365
doy.rad <- doy/365*2*pi
 
# function to return predictions for each day of year given posterior sample of parameters
det.pred.fun <- function(samp.vec, water = 0){
  pred <- samp.vec["(Intercept)"] + water*samp.vec["water"] + sin(doy.rad)*samp.vec["sin.doy"] + cos(doy.rad)*samp.vec["cos.doy"]
  plogis(pred) # back onto the probability scale
}

# without water
det.fun.samps <- matrix(ncol = length(doy), nrow = nrow(ms.so.fit$alpha.samples))
for (ss in 1:nrow(ms.so.fit$alpha.samples)){
  det.fun.samps[ss,] <- det.pred.fun(ms.so.fit$alpha.samples[ss,])
}
det.fun.mean <- apply(det.fun.samps, 2, mean)
det.fun.sd <- apply(det.fun.samps, 2, sd)

# with water
det.fun.samps <- matrix(ncol = length(doy), nrow = nrow(ms.so.fit$alpha.samples))
for (ss in 1:nrow(ms.so.fit$alpha.samples)){
  det.fun.samps[ss,] <- det.pred.fun(ms.so.fit$alpha.samples[ss,], water = 1)
}
det.fun.mean.water <- apply(det.fun.samps, 2, mean)
det.fun.water.sd <- apply(det.fun.samps, 2, sd)

p <- ggplot() +
  geom_errorbar(aes(x = doy, 
                    ymin = det.fun.mean-2*det.fun.sd,
                    ymax = det.fun.mean+2*det.fun.sd),
                width = 0,
                col = "lightgrey") +
  geom_line(aes(x = doy, y = det.fun.mean, col = "No water")) +
  geom_errorbar(aes(x = doy, 
                    ymin = det.fun.mean.water-2*det.fun.water.sd,
                    ymax = det.fun.mean.water+2*det.fun.water.sd),
                width = 0,
                col = "lightgrey") +
  geom_line(aes(x = doy, y = det.fun.mean.water, col = "Water")) +
  labs(x = "Day of year", y = "Detection probability") +
  scale_fill_manual('Legend Title', values=c('Water', 'No water')) +
  ylim(0, 0.4) +
  theme_minimal()
p
ggsave(filename = "out/multi-season-spatial-occupancy-detection-over-time.png")

rm(det.fun.samps, det.fun.mean, det.fun.mean.water, det.fun.sd, det.fun.water.sd)

##### Look at how correlation in occupancy changes with distance #####
x <- seq(0, 10000, 10) # vector of distances, in km
cor.fun.samps <- matrix(ncol = length(x), nrow = nrow(ms.so.fit$theta.samples))
for (ss in 1:nrow(ms.so.fit$theta.samples)){
  cor.fun.samps[ss,] <- exp(-ms.so.fit$theta.samples[ss, "phi"]*x)
}
cor.fun.mean <- apply(cor.fun.samps, 2, mean)
cor.fun.sd <- apply(cor.fun.samps, 2, sd)

p <- ggplot() +
  geom_errorbar(aes(x = x, 
                    ymin = cor.fun.mean-2*cor.fun.sd,
                    ymax = cor.fun.mean+2*cor.fun.sd),
                width = 0,
                col = "lightgrey") +
  geom_line(aes(x = x, y = cor.fun.mean)) +
  labs(x = "Distance (m)", y = "Correlation between sites") +
  theme_minimal()
p
ggsave(filename = "out/multi-season-spatial-occupancy-correlation-over-distance.png")
