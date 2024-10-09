# read in setup stuff
source("src/a-setup.R")

# get survey data (contains points only)
coltypes <- c("numeric", "numeric", "numeric", "text", 
              "numeric", "numeric", "numeric", "numeric", "text", "numeric", "numeric", 
              "numeric", "numeric", "numeric", "numeric", "date", "date", "text", 
              "text", "numeric", "numeric", "numeric", "numeric", "numeric", 
              "text", "text", "text", "text", "text", "numeric", "text", 
              "numeric", "text", "text", "numeric", "numeric", "numeric", "text", 
              "text", "numeric", "numeric", "numeric", "numeric", "text", "numeric", "text", "numeric", 
              "text", "text", "text", "text", "text", "numeric", 
              "numeric", "numeric", "text", "text", "text", "text", "text", "text", "text", 
              "text", "text", "text", "numeric", "text")

df <- read_xlsx(path = file.path(data_dir, current.data), 
                sheet = "Export",
                col_types = coltypes)

df <- mutate(df, ID = row_number()) # add an ID column

# extract the colony data
cny.df <- df %>%
  filter(!is.na(ColonyNumber) & ColonyNumber != 0)

# get first record for each colony for binding back to df
cny.detected <- cny.df %>% 
  group_by(ColonyNumber) %>%
  arrange(ColonyNumber, dateOfActivityTreatment) %>%
  filter(row_number()==1) %>%
  mutate(SurveillanceActivityValue = rep("Colony found", n()))

# organise colony data as its own entity
cny.df <- cny.df %>%
  select(ID,
         date = dateOfActivity,
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
  rbind(cny.detected) %>% # put single detection of each colony back in
  mutate(presence = ifelse(grepl("Red", SpeciesObservedValue), 1, 0), # present/absent data
         hive.removed = ifelse(grepl("Colony found", SurveillanceActivityValue, ignore.case = TRUE), 1, 0),
         hour = hour(dateOfActivity),
         hour2 = hour^2) %>% 
  select(ID,
         date.time = dateOfActivity, 
         lat = Lat, 
         long = Long, 
         hour, 
         hour2,
         dist.forage.m = DistanceFromKnownForage_meters,
         dist.water.m = DistanceFromKnownWaterSource_met,
         dist.road.m = DistanceFromRoad_meters,
         flight.line.deg = FlightLineDegrees,
         host.flower = HostFloweringValue,
         slope.orientation = SlopeOrientation,
         temperatue = Temperature,
         Notes,
         WaterForagingValue,
         HazardTypeValue,
         HostOther,
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
df$dist.0 <- as.numeric(st_distance(earliest_record, df))/1000 # in kms

# make other useful covariates
df <- mutate(df, 
             time.0 = (date.time-earliest_record$date.time)/(60*60*24), # time since incursion detected
             water = ifelse(grepl("water", Notes, ignore.case = TRUE) | 
                              grepl("water", WaterForagingValue, ignore.case = TRUE) |
                              grepl("water", HazardTypeValue, ignore.case = TRUE) |
                              grepl("water", HostOther, ignore.case = TRUE), 1, 0),
             flowering = ifelse(grepl("Fully Flowering", host.flower) | 
                                  grepl("Partially Flowering", host.flower), 1, 0)) %>% #water around?
  select(-Notes, -WaterForagingValue,
         -HazardTypeValue,
         -HostOther, -host.flower)

rm(cny.detected, earliest_record, coltypes)
