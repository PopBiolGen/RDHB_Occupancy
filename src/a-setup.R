# A file containing useful bits and pieces that are used by all downstream scripts
library(dplyr)
library(lubridate)
library(sf)
library(leaflet)
library(ggplot2)
library(readxl)


# Define local directory containing data files
data_dir <- file.path(Sys.getenv("DATA_PATH"), "RDHB")
# Note: usethis::edit_r_environ()
# and set DATA_PATH="your/local/path/to/data/directory"
# and then restart R

##### Custom functions #####

# To aggregate data in space and time
aggregate_data <- function(df, cell.size = 0.005){
  # make a grid and spatial join to point data 
  df_grid <- spatial_aggregation(df, cell.size)
  # place cell.id onto points data
  cells <- select(df_grid, ID, cell.id) %>% st_drop_geometry()
  df <- left_join(df, cells)
  # make grid summaries
  df_grid <- make_grid_summary(df_grid)
  # make time aggregations
  if (!("time.step" %in% names(df))){
    df <- temporal_aggregation(df)  
  }
  return(list(df = df, df_grid = df_grid))
}

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

# Makes grid-level summary after spatial join of grid onto point-level data
# returns grid dataframe with one row per grid cell
make_grid_summary <- function(df){
  get_first <- function(x){x[1]}
  grid_summ <- df %>%
    group_by(cell.id) %>%
    summarise(mean.dist = mean(dist.0),
              mean.prop = mean(presence),
              n.hive.removed = sum(hive.removed))
  # find neighbours for each cell
  grid_summ$neighbours <- lapply(1:nrow(grid_summ), function(i) {
    find_neighbours(grid_summ[i, ], grid_summ)
  })
  grid_summ
}

# To map points given points and grid_dataframe
map_point_grid <- function(df, df_grid, summ.col = mean.prop){
  # extract lat and long
  df <- cbind(df, st_coordinates(df)) %>% rename(long = X, lat = Y)
  # make plots of the data and grid
  plot_grid <- df_grid %>%
    filter(!st_is_empty(geometry)) %>%
    mutate(plot.col = ifelse(is.na({{summ.col}}), 0, {{summ.col}})) %>%
    mutate(plot.col = plot.col/max(plot.col))
  
  map <- leaflet() %>%
    addTiles() %>%
    
    ## Add detection points
    addCircleMarkers(data = df[df$presence == 0, ], ~long, ~lat, radius = 3, color = "blue", fillColor = "blue", fillOpacity = 1, group = "Absent") %>%
    addCircleMarkers(data = df[df$presence == 1, ], ~long, ~lat, radius = 3, color = "red", fillColor = "red", fillOpacity = 1, group = "Present") %>%
    
    
    addPolygons(data = plot_grid,
                color = "blue",          # Color of the polygon borders
                weight = 2,              # Weight of the polygon borders
                fillColor = "blue",      # Fill color of the polygons
                fillOpacity = plot_grid$plot.col,       # Opacity of the fill color
                # fillOpacity = results_pol$prob,       # Opacity of the fill color - note too dark to be meaningful
                label = ~cell.id,             # Labels for the polygons
                group = "Summary")   %>%          
    
    addScaleBar(position = "bottomleft")%>%
    fitBounds(116.72, -20.82, 116.76, -20.42)%>%  # Set the bounding box
    
    addLayersControl(
      overlayGroups = c("Present", "Absent", 
                        "Summary"),
      options = layersControlOptions(collapsed = FALSE)) 
  
  map
  
}


# To undertake spatial aggregation
# takes sf dataframe of point data
# makes a grid and spatial joins
# identifies neighbours for each grid cell
spatial_aggregation <- function(sf.df, cell.size){
  # make a grid polygon using bbox of sf.df
  grid <- st_make_grid(sf.df, cellsize = cell.size, square = TRUE) 
  # Convert grid to sf dataframe object and give id numbers
  sf_grid <- st_sf(geometry = st_sfc(grid), cell.id = 1:length(grid)) 
  # spatial join
  grid_d <- st_join(sf_grid, sf.df, join = st_intersects)
  grid_d
}

# To undertake temporal aggregation
# take sf points dataframe
# aggregates to year-month
temporal_aggregation <- function(sf.df){
  sf.df <- sf.df %>%
    mutate(ym = format(date.time, "%Y-%m"),
      time.step = as.numeric(as.factor(ym)))
  sf.df
}
