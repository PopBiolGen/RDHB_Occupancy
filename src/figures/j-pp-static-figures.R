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
                                 "alpha.sig",
                                 "beta.1.sig",
                                 "beta.2.sig",
                                 "u.0"))
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
# Convert to a data frame for ggplot
density_df <- data.frame(
  x = rep(density_est$x, each = length(density_est$y)),
  y = rep(density_est$y, times = length(density_est$x)),
  density = as.vector(t(density_est$z))
)

# Plot the density with ggplot2
# set bounding box to extent of data
bbox <- st_bbox(density_df)

ggplot(density_df, aes(x=x, y=y)) +
  geom_sf(data = shoreline, fill = "lightblue", color = "blue", inherit.aes = FALSE) +
  coord_sf(xlim = range(density_df$x), ylim = range(density_df$y)) +  # Apply bounding box
  geom_raster(aes(fill = density), interpolate = TRUE) +
  geom_contour(aes(z = density), color = "black", alpha = 0.5) +
  scale_fill_viridis_c(alpha = 0.4) +
  geom_point(data = df.mr[df.mr$presence==0,], 
             aes(x = X, y = Y), 
             colour = "black",
             inherit.aes = FALSE) +
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
