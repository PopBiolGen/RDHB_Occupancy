# read in setup stuff
source("src/a-setup.R")

# get survey data (contains points, polylines and polygons)
df <- st_read(file.path(data_dir, "RDHBSurveillance_05062024.csv"), 
              options = c("GEOM_POSSIBLE_NAMES=Spatial"),
              crs = 4326)

# remove empty geometries and define date/time, presence/absence, select only the useful columns
df <- df %>%
  rename(geometry = Spatial) %>%
  filter(!st_is_empty(geometry)) %>%
  mutate(date_time = dmy_hm(dateOfActivity, tz = "Australia/Perth"), # date/time/hour
         hour = hour(date_time)) %>%
  mutate(pres = ifelse(grepl("Red", SpeciesObservedInFieldTXT), 1, 0)) %>% # present/absent data
  select(ID, Title, date_time, hour, pres, Notes, HostLookup, SpeciesCount, # only take columns we need
           HostOther, HostFlowering, CaseLink, SmallGridID, SurveillanceActivity, 
           ColonyNumber, CaseTXT,
           SpeciesActivity, 
           SmallGridID.ID, 
           DistanceFromRoad_meters, SlopeOrientation, 
           DistanceFromKnownForage_meters, AltitudeInMeters, 
           geometry)
  
# identify first point in the remaining geometries
first_points <- st_sfc(lapply(df$geometry, extract_first_point))

# replace the complex multipoint geometry collections with a single point
df$geometry <- first_points
rm(first_points)

# trim down to only data where we have a point record
df <- subset(df,!st_is_empty(df$geometry))

# filter out spatial / temporal outliers (data entry or other issues?)
df <- df %>%
  filter(year(date_time) != 2005) %>% #remove surprising year
  mutate(lat = unlist(lapply(geometry,function(x){x[2]})), # extract lat and long for ease of access
         long = unlist(lapply(geometry,function(x){x[1]}))) %>%
  filter(lat > -20.8) # remove surprising points a long way south

# make a grid and spatial join to point data 
df_grid <- spatial_aggregation(df)

# make time aggregations
df <- temporal_aggregation(df)
