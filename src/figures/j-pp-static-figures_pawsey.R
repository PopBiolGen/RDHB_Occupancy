
# Remove before using Pawsey:
#iter <- 3

#install.packages('MCMCvis', 
#                 repos = repo,
#                 lib = mylib,
#                 dependencies = TRUE)
#library(MCMCvis)

library(rjags)
library(dplyr)


##### make some plots #####
library(ggplot2)
source("src/a-setup.R")

# Use pawsey version! (subsets surveillance data by date)
source("src/j-pp-static-get-data_pawsey.R")

load(file = paste("out/temp-coda-start_", 
                  iter, 
                  ".RData", sep="")) # Load data

# grab a shoreline
# Slightly different directory -- make consistent!
# Already crs transformed
shoreline <- st_read(file.path(Sys.getenv("DATA_PATH"), 
                               "RDHB/Spatial/Australia_boundary.shp")) |>
  st_transform(crs = 32750)


### Non-loc parameters
#temp <- MCMCchains(b, params = c("psi",
#                                 "alpha.u",
#                                 "beta.1.u",
#                                 "beta.2.u",
#                                 "sigma.det"))
#pdf(file = "out/static-pairs-parameters.pdf")
#pairs(temp)
#dev.off()

### density map
#temp <- MCMCchains(b, params = c("loc")) # Produces matrix of MCMC output, where each row is an interation...
  # So should be the coords of each predicted colony per iter, that best fits with data...

# *MCMCvis not installing on Pawsey, but this gets the same thing:
temp <- as.data.frame(as.matrix(b)) # as.matrix produces the same as above, just prints other parameters as well (alpha, beta.1, beta.2, psi, sigma.det)
temp <- temp %>% # Clunky, but go to df and back to use dplyr to select only loc columns
  select(contains("loc", ignore.case = F))
temp <- as.matrix(temp)

x <- as.vector(temp[,1:(ncol(temp)/2)]) # First half of columns are long values
y <- as.vector(temp[,(ncol(temp)/2+1):ncol(temp)]) # Second half are lat values 
point.data <- data.frame(x = x, y = y) |> subset(x!=0) # Grid of all locations... 

#... of inferred location od colonies across all iterations of model
# I think it looks weird because so many are dropped in the ocean (where there are no surveys to show that of course they aren't there)
# Need to ignore all oceans points when running model...

# Estimate 2D kernel density
density_est <- MASS::kde2d(point.data$x, point.data$y, n = 100, lims = c(min(df.mr$X), max(df.mr$X), min(df.mr$Y), max(df.mr$Y)))  # 100x100 grid

# ^ What is this doing?? Smoothing a kernel over all those points to give a general density of points? 


# Convert to a data frame for ggplot (note the transpose for graphing)
density_df <- data.frame(
  x = rep(density_est$x, each = length(density_est$x)),
  y = rep(density_est$y, times = length(density_est$y)),
  density = as.vector(t(density_est$z))
)

# Plot the density with ggplot2

dens_plot <- ggplot(density_df, aes(x=x, y=y)) +
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
       fill = "Density")


ggsave(dens_plot, 
     file = sprintf("out/figs/static-density-map-start_%s.pdf", iter))

# export density raster and contours into a geopackage
#made_date <- today()
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
#rast.fname <- paste0("out/density_raster", made_date, ".tif")
# Save density raster
writeRaster(r, 
            filename = paste0("out/figs/density_raster_start_", iter, ".tif"))


##### Save as Geopackage ####
## See Ben's figure script