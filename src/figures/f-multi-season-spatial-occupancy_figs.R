# to visualise results from the multi-season-spatial-occupancy model
# source("src/f-multi-season-spatial-occupancy.R")

library(gganimate)
library(spOccupancy)

load(file = "out/f-ms-so_fit.RData")

est.psi <- apply(ms.so.fit$psi.samples, c(2, 3), mean) # estimated occupancy at each time.step

# join occupancy estimates onto grid
grid.dat <- cbind(agg_data$df_grid, est.psi)


##### Make an animated map of occupancy over time #####
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

#animate(anim, nframes = max(gd.long$time.step), fps = 1)


##### Make map of occupancy in last time step #####
p <- ggplot() +
  geom_sf(data = shoreline, fill = "lightblue", color = "blue") +
  geom_sf(data = filter(gd.long, time.step == max(time.step)), aes(fill = Occupancy)) +
  coord_sf(xlim = c(bbox["xmin"], bbox["xmax"]), ylim = c(bbox["ymin"], bbox["ymax"])) +  # Apply bounding box
  scale_fill_viridis_c(option = "plasma") +
  theme_minimal() +
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
  theme_minimal()

ggsave(filename = "out/multi-season-spatial-occupancy-mean-occ-over-time.png")
