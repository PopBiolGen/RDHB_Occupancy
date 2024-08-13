# to visualise results from the dynamic occupancy model
source("src/e-dynamic-occupancy.R")
library(gganimate)

# look at predicted occupancy per time step 
occ.dyn <- fit.do@smoothed[2, ,] #returns a time x site smoothed estimate of occupancy probability
occ.dyn.mean <- fit.do@smoothed.mean[2,] # returns mean probability of occupancy across all sites over time

# join occupancy estimates onto grid
grid.dat <- cbind(agg_data$df_grid, t(occ.dyn))

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