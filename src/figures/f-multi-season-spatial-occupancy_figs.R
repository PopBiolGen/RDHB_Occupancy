# to visualise results from the multi-season-spatial-occupancy model
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
  theme_bw() +
  labs(title = "Time step: {closest_state}")

anim <- p +
  transition_states(time.step, transition_length = 2, state_length = 1) +
  ease_aes("cubic-in-out")

#animate(anim, nframes = max(gd.long$time.step), fps = 1)


##### Make map of occupancy in last time step #####
p <- ggplot() +
  geom_sf(data = shoreline, fill = "lightblue", color = "blue") +
  geom_sf(data = filter(gd.long, time.step == max(time.step)), aes(fill = Occupancy)) +
  coord_sf(xlim = c(bbox["xmin"], bbox["xmax"]), ylim = c(bbox["ymin"], bbox["ymax"])) +  # Apply bounding box
  scale_fill_viridis_c(option = "plasma") +
  theme_bw() +
  labs(title = paste0("Month of response: ", max(gd.long$time.step)))



ggsave(filename = "out/multi-season-spatial-occupancy.png")


##### Examine how occupancy changes over time #####
est.psi.time <- apply(est.psi, MARGIN = 2, mean)
est.psi.time.se <- apply(est.psi, MARGIN = 2, sd)/sqrt(nrow(est.psi)) # loses parameter uncertainty

p <- ggplot() +
  geom_point(aes(x = 1:TT, y = est.psi.time)) +
  geom_errorbar(aes(x = 1:TT, 
                    ymin = est.psi.time-2*est.psi.time.se,
                    ymax = est.psi.time+2*est.psi.time.se),
                width = 0) +
  labs(x = "Months of response", y = "Mean occupancy") +
  ylim(0, 0.3) +
  theme_bw()

ggsave(filename = "out/multi-season-spatial-occupancy-mean-occ-over-time.png")


##### Look at how detection changes throughout the year #####
doy <- 1:365
doy.rad <- doy/365*2*pi
 
# function to return predictions for each day of year given posterior sample of parameters
det.pred.fun <- function(samp.vec, water = 0, flower = 0){
  pred <- samp.vec["(Intercept)"] + 
    water*samp.vec["water"] + 
    flower*samp.vec["flowering"] +
    sin(doy.rad)*samp.vec["sin.doy"] + 
    cos(doy.rad)*samp.vec["cos.doy"]
  plogis(pred) # back onto the probability scale
}

# function to generate mean and se of detection for each day of the year
det.generator <- function(doy = 1:365, fit = ms.so.fit, water = 0, flower = 0) {
  det.fun.samps <- matrix(ncol = length(doy), nrow = nrow(fit$alpha.samples))
  for (ss in 1:nrow(fit$alpha.samples)){
    det.fun.samps[ss,] <- det.pred.fun(fit$alpha.samples[ss,], water = water, flower = flower)
  }
  out.mean <- apply(det.fun.samps, 2, mean)
  out.sd <- apply(det.fun.samps, 2, sd)
  water <- rep(water, length(out.mean))
  flower <- rep(flower, length(out.mean))
  data.frame(water = water, flower = flower, doy = doy, mean = out.mean, se = out.sd)
}


# No flowering
# without water
nw.nf <- det.generator(water = 0, flower = 0)
# with water
w.nf <- det.generator(water = 1, flower = 0)

# With flowering
# without water
nw.f <- det.generator(water = 0, flower = 1)
# with water
w.f <- det.generator(water = 1, flower = 1)

plot.det <- rbind(nw.nf, w.nf, nw.f, w.f)
rm(nw.nf, w.nf, nw.f, w.f)

p <- ggplot(data = plot.det, aes(x = doy)) +
  geom_errorbar(aes(ymin = mean-2*se,
                    ymax = mean+2*se),
                width = 0,
                col = "lightgrey") +
  geom_line(aes(y = mean)) +
  labs(x = "Day of year", y = "Detection probability") +
  scale_fill_manual('Legend Title', values=c('Water', 'No water')) +
  facet_grid(rows = vars(water), 
             cols = vars(flower),
             labeller = labeller(flower = c(`0` = "No flowers", `1` = "Flowers"),
                                 water = c(`0` = "No water", `1` = "Water"))) +
  #ylim(0, 0.4) +
  theme_bw()
p
ggsave(filename = "out/multi-season-spatial-occupancy-detection-over-time.png")

rm(plot.det)

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
  theme_bw()
p
ggsave(filename = "out/multi-season-spatial-occupancy-correlation-over-distance.png")
