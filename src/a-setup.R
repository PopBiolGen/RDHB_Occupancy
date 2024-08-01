# A file containing useful bits and pieces that are used by all downstream scripts
library(dplyr)
library(lubridate)
library(sf)

# Define local directory containing data files
data_dir <- file.path(Sys.getenv("DATA_PATH"), "RDHB")
# Note: usethis::edit_r_environ()
# and set DATA_PATH="your/local/path/to/data/directory"
# and then restart R

##### Custom functions #####

# To extract the first point from a geometry
# Used here to simplify the mess of geometries in the dataset to just the first POINT in each geometry
extract_first_point <- function(geometry) {
  if (inherits(geometry, "sfc_GEOMETRYCOLLECTION") || inherits(geometry, "GEOMETRYCOLLECTION")) {
    # Extract the first POINT from the geometry collection
    for (geom in geometry) {
      if (inherits(geom, "sfc_POINT") || inherits(geom, "POINT")) {
        return(geom)
      }
    }
  } else if (inherits(geometry, "sfc_POINT") || inherits(geometry, "POINT")) {
    # If it's already a POINT, return it
    return(geometry)
  }
  return(st_point())  # Return an empty point if no point found
}

# To find neighbours for each cell
find_neighbours <- function(cell, grid) {
  neighbours <- st_touches(cell, grid, sparse = FALSE)
  neighbour_ids <- which(neighbours)
  return(neighbour_ids)
}

# To undertake spatial aggregation
# takes sf dataframe of point data
# makes a grid and spatial joins
# identifies neighbours for each grid cell
spatial_aggregation <- function(sf.df, cell.size = 0.1){
  # make a grid polygon using bbox of sf.df
  grid <- st_make_grid(sf.df, cellsize = cell.size, square = TRUE) 
  # Convert grid to sf dataframe object and give id numbers
  sf_grid <- st_sf(geometry = st_sfc(grid), cell.id = 1:length(grid)) 
  # find neighbours for each cell
  sf_grid$neighbours <- lapply(1:nrow(sf_grid), function(i) {
    find_neighbours(sf_grid[i, ], grid)
  })
  # spatial join
  grid_d <- st_join(sf_grid, sf.df, join = st_intersects)
}

# To undertake temporal aggregation
# take sf points dataframe
temporal_aggregation <- function(sf.df, n.chunks = 6){
  sf.df <- sf.df %>%
    mutate(time.step = as.numeric(cut(date_time, breaks = n.chunks)))
  sf.df
}
