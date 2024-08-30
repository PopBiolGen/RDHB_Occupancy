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

# extract the colony data
cny.df <- df %>%
  filter(!is.na(ColonyNumber) & ColonyNumber != 0 ) %>%
  select(date = dateOfActivity,
         lat = Lat,
         long = Long,
         contains("Comb"), 
         contains("Cell"), 
         contains("Colony"), 
         DiagnosticNotes) %>%
  st_as_sf(coords = c("long", "lat"), crs = 4326)


# tidy up, select only the useful columns, cast to sf
df <- df %>%
  filter(is.na(ColonyNumber) | ColonyNumber == 0 ) %>% # remove colony data
  mutate(presence = ifelse(grepl("Red", SpeciesObservedValue), 1, 0),
         hive.removed = ifelse(grepl("Colony found", SurveillanceActivityValue, ignore.case = TRUE), 1, 0)) %>% # present/absent data
  select(date.time = dateOfActivity, 
         lat = Lat, 
         long = Long, 
         dist.forage.m = DistanceFromKnownForage_meters,
         dist.water.m = DistanceFromKnownWaterSource_met,
         dist.road.m = DistanceFromRoad_meters,
         flight.line.deg = FlightLineDegrees,
         host.flower = HostFloweringValue,
         slope.orientation = SlopeOrientation,
         temperatue = Temperature,
         hive.removed,
         presence, 
         abundance = AbundanceValue) %>%
  filter(!is.na(lat)) %>%
  filter(lat > -20.7) %>% # remove surprising points a long way south
  filter(year(date.time) != 2005) %>% #remove surprising year
  st_as_sf(coords = c("long", "lat"), crs = 4326) # cast to a spatial object

# calculate location of earliest record and distance from there to all other records
earliest_record <- filter(df, presence == 1) %>%
                    filter(date.time == min(date.time))
df$dist_0 <- as.numeric(st_distance(earliest_record, df))/1000 # in kms

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
#Tamrika Testing if this gets through to Ben