# read in setup stuff
source("src/a-setup.R")

# get survey data (contains points, polylines and polygons)
df <- st_read(file.path(data_dir, "RDHBSurveillance_05062024.csv"), 
              options = c("GEOM_POSSIBLE_NAMES=Spatial"),
              crs = 4326)

# remove empty geometries
df <- df %>%
  filter(!st_is_empty(Spatial))
  
# identify first point in the remaining geometries
first_points <- st_sfc(lapply(df$Spatial, extract_first_point))

# replace the complex multipoint geometry collections with a single point
df$Spatial <- first_points

# trim down to only data where we have a point record
df <- subset(df,!st_is_empty(df$Spatial))

# make a grid and spatial join to point data 
df_grid <- spatial_aggregation(df)