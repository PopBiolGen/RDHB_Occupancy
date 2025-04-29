##### make some plots #####
library(ggplot2)
source("src/a-setup.R")
source("src/j-pp-static-get-data.R")
load(file = "out/temp-coda.RData")

# grab a shoreline
shoreline <- st_read(file.path(Sys.getenv("DATA_PATH"), 
                               "Spatial/Igis_aus_outline/Australia_boundary.shp")) |>
             st_transform(crs = 32750) # UTM 50S)

### Non-loc parameters
temp <- MCMCchains(b, params = c("psi",
                                 "alpha.u",
                                 "beta.1.u",
                                 "beta.2.u",
                                 "sigma.det"))
pdf(file = "out/static-pairs-parameters.pdf")
pairs(temp)
dev.off()

### density map
temp <- MCMCchains(b, params = c("loc"))
x <- as.vector(temp[,1:(ncol(temp)/2)])
y <- as.vector(temp[,(ncol(temp)/2+1):ncol(temp)])
point.data <- data.frame(x = x, y = y) |> subset(x!=0)
# Estimate 2D kernel density
density_est <- MASS::kde2d(point.data$x, point.data$y, n = 100, lims = c(min(df.mr$X), max(df.mr$X), min(df.mr$Y), max(df.mr$Y)))  # 100x100 grid

# Convert to a data frame for ggplot (note the transpose for graphing)
density_df <- data.frame(
  x = rep(density_est$x, each = length(density_est$x)),
  y = rep(density_est$y, times = length(density_est$y)),
  density = as.vector(t(density_est$z))
)

# Plot the density with ggplot2

ggplot(density_df, aes(x=x, y=y)) +
  geom_sf(data = shoreline, fill = "lightblue", color = "blue", inherit.aes = FALSE) +
  coord_sf(xlim = range(density_df$x), ylim = range(density_df$y)) +  # Apply bounding box
  geom_raster(aes(fill = density), interpolate = TRUE) +
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
       fill = "Density")
ggsave("out/static-density-map.pdf")

# export density raster and contours into a geopackage
made_date <- today()
library(terra)
# Density raster
# Create an empty raster with the correct dimensions
r <- rast(
  nrows = length(density_est$x), 
  ncols = length(density_est$y),
  xmin = min(density_est$x), xmax = max(density_est$x),
  ymin = min(density_est$y), ymax = max(density_est$y)
)
crs(r) <- "EPSG:32750" # set CRS
# Assign density values to the raster (note transpose and flip)
values(r) <- t(density_est$z[,nrow(density_est$z):1]) 
# write it out
rast.fname <- paste0("out/density_raster", made_date, ".tif")
writeRaster(r, 
            filename = rast.fname, 
            overwrite = TRUE)

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

# Some metadata
contour_sf$created_by <- "Ben Phillips"
contour_sf$created_on <- Sys.Date()
contour_sf$description <- "Contour lines generated from 2D density surface predicting the locality of hives."

# Path to GDB (will be created if needed)
gpkg_path <- "out/output_contours.gpkg"

# Save as a feature class inside the GDB
st_write(contour_sf, dsn = gpkg_path, layer = "density_contours", driver = "GPKG", delete_layer = TRUE)

# List of files to include in the zip
files_to_zip <- c("out/output_contours.gpkg", rast.fname, "out/static-density-map.pdf")

# Create the zip archive
zip(zipfile = paste0("out/density_outputs", made_date, ".zip"), files = files_to_zip)