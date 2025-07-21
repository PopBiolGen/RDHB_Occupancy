
# Understanding plots

# Possible colony locations??
ggplot(point.data, 
       aes(x=x, y=y))+
  geom_point()+
  coord_cartesian(xlim = c(474000, 476000),
                  ylim = c(7720000, 7722000))
# Smoothing over points to give density estimates?
ggplot(density_df, 
       aes(x=x, y=y, fill=density))+
  geom_tile()+
  coord_cartesian(xlim = c(474000, 476000),
                  ylim = c(7720000, 7722000))
# Turn into raster, and add contours
ggplot(density_df, 
       aes(x=x, y=y)) +
  geom_raster(aes(fill = density), interpolate = TRUE) +
  geom_contour(aes(z = density))

# Overlay density with map raster, survey point data, and estimated colony locations
ggplot(density_df, aes(x=x, y=y)) +
  geom_sf(data = shoreline, fill = "lightblue", color = "blue", inherit.aes = FALSE) + # Draw coastline
  coord_sf(xlim = range(density_df$x), ylim = range(density_df$y)) +  # Apply bounding box
  geom_raster(aes(fill = density), # use geom_raster when every cell has data (otherwise use geom_tile)
              interpolate = TRUE) + # interpolate smooths between cells
  geom_contour(aes(z = density), color = "black", alpha = 0.5) +
  scale_fill_viridis_c(alpha = 0.4) +
  geom_point(data = df.mr[df.mr$presence==0,], 
             aes(x = X, y = Y), 
             colour = "black",
             inherit.aes = FALSE,
             alpha = 0.1) +
  geom_point(data = df.mr[df.mr$presence==1,], 
             aes(x = X, y = Y), 
             colour = "red",
             inherit.aes = FALSE) +
  theme_minimal() +
  labs(title = "2D Kernel Density Estimation",
       x = "X Coordinate",
       y = "Y Coordinate",
       fill = "Density")+
  geom_point(data=point.data, # Check if colonies being dropped in ocean
             aes(x=x, y=y), alpha=0.1, col="navy")


# Save density as raster
plot(r)

### Ben's script: ####
# Generate contour lines
clines <- contourLines(x = unique(density_df$x), 
                       y = unique(density_df$y), 
                       z = matrix(density_df$density, 
                                  nrow = length(unique(density_df$y)), 
                                  byrow = TRUE))
# Convert the contourLines list to an sf object
# Convert to sf with labels
contour_sf <- do.call(rbind, lapply(seq_along(clines), function(i) {
  cl <- clines[[i]]
  coords <- cbind(cl$x, cl$y)
  sf::st_sf(
    level = cl$level,                         # Contour value for labeling
    label = paste0("Level: ", cl$level),      # Optional: formatted label string
    geometry = sf::st_sfc(sf::st_linestring(coords)),
    crs = 32750  # or your projected CRS
  )
}))

#####
ggplot(contour_sf) +
  geom_sf()

# How to extract value at ANY point? 

row.n <- 125
m.point <- matrix(c(df.mr$X[row.n], df.mr$Y[row.n]), ncol=2)
# Extracting values from raster?
extract(x = r, 
        y = m.point)

# Seems to work now, even if colony coords are NOT the exact coords used in df.mr
