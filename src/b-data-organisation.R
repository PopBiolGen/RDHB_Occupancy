# read in setup stuff
source("src/a-setup.R")

# get survey data (contains points, polylines and polygons)
coltypes <- c("numeric", "numeric", "numeric", "text", 
              "numeric", "numeric", "numeric", "numeric", "text", "numeric", "numeric", 
              "numeric", "numeric", "numeric", "numeric", "date", "date", "text", 
              "text", "numeric", "numeric", "numeric", "numeric", "numeric", 
              "text", "text", "text", "text", "text", "numeric", "text", 
              "numeric", "text", "text", "numeric", "numeric", "numeric", "text", 
              "text", "numeric", "numeric", "numeric", "numeric", "text", "numeric", "text", "numeric", 
              "text", "text", "text", "text", "text", "numeric", 
              "numeric", "numeric", "text", "text", "text", "text", "text", "text", "text", 
              "text", "text")

df <- read_xlsx(path = file.path(data_dir, "RDHBSurveillance_2024-08-12.xlsx"), 
                sheet = "Export",
                col_types = coltypes)

### to do, cast to sf ###

# remove empty geometries and define date/time, presence/absence, select only the useful columns
df <- df %>%
  rename(geometry = Spatial) %>%
  filter(!st_is_empty(geometry)) %>%
  mutate(date.time = dmy_hm(dateOfActivity, tz = "Australia/Perth"), # date/time/hour
         hour = hour(date.time)) %>%
  mutate(pres = ifelse(grepl("Red", SpeciesObservedInFieldTXT), 1, 0)) %>% # present/absent data
  select(ID, Title, date.time, hour, pres, Notes, HostLookup, SpeciesCount, # only take columns we need
           HostOther, HostFlowering, CaseLink, SmallGridID, SurveillanceActivity, 
           ActivityTXT, ColonyNumber, CaseTXT,
           SpeciesActivity, 
           SmallGridID.ID, 
           DistanceFromRoad_meters, SlopeOrientation, 
           DistanceFromKnownForage_meters, AltitudeInMeters, 
           geometry)
  
# identify first point in the remaining geometries
first_points <- st_sfc(lapply(df$geometry, extract_first_point))

# replace the complex multipoint geometry collections with a single point
df$geometry <- first_points
st_crs(df) <- 4326
rm(first_points)

# trim down to only data where we have a point record
df <- subset(df,!st_is_empty(df$geometry))

# filter out spatial / temporal outliers (data entry or other issues?)
df <- df %>%
  filter(year(date.time) != 2005) %>% #remove surprising year
  mutate(lat = unlist(lapply(geometry,function(x){x[2]})), # extract lat and long for ease of access
         long = unlist(lapply(geometry,function(x){x[1]}))) %>%
  filter(lat > -20.7) # remove surprising points a long way south

# calculate location of earliest record and distance from there to all other records
# we want distance in metres, so first cast to Australian Albers (CRS = 3577)
df_albers <- select(df, geometry, pres, date.time) %>%
              st_transform(crs = 3577)
earliest_record <- filter(df_albers, pres==1) %>%
                    filter(date.time == min(date.time))
df$dist_0 <- as.numeric(st_distance(earliest_record, df_albers))
rm(df_albers)

# make other useful covariates
df <- mutate(df, time_0 = (date.time-earliest_record$date.time)/(60*60*24), # time since incursion detected
             water = ifelse(grepl("water", Notes, ignore.case = TRUE) | # water around?
                              grepl("water", HostOther, ignore.case = TRUE), 1, 0),
             flowering = ifelse(HostFlowering %in% c("2;#Partially Flowering", "3;#Fully Flowering"), 1, 0), # flowering host?
             food.water = ifelse(grepl(1, flowering, ignore.case = TRUE) | # food or water
                                   grepl(1, water, ignore.case = TRUE), 1, 0),
             hive.removed = ifelse(grepl("Colony found", SurveillanceActivity, ignore.case = TRUE) | # hive removed?
                                     grepl("Colony found", ActivityTXT, ignore.case = TRUE), 1, 0),
             hour2 = hour^2) %>%
  select(ID,
         date.time, # grab useful stuff, ditch the rest
         hour, 
         hour2,
         pres, 
         lat, 
         long, 
         dist_0, 
         time_0, 
         water, 
         flowering, 
         food.water,
         hive.removed,
         geometry)
