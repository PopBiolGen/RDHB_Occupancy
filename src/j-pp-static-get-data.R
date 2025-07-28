library(rjags)
library(MCMCvis)

##### Get the data #####
# Surveillance data
current.data <- "RDHBSurveillance_2025-07-21.xlsx"
source("src/b-data-organisation.R")
# foraging data
fd <- read_xlsx(path = file.path(data_dir, "Abrol-foraging-data.xlsx")) |> 
  mutate(n.round = round(n.bees)) |>
  tidyr::uncount(n.round)


##### Filter and get coordinates #####
up.to.date <- as.Date("2025-07-21") # analyse data in the 90 days up to this date
# most recent (last 3 months) data: ".mr"
df.mr <- df |> 
  filter(date.time <= up.to.date & date.time > (up.to.date-90)) 
# get coordinates in UTM 50S
df.utm <- df.mr |>
  st_transform(crs = 32750) |> # UTM 50S
  st_coordinates()
# cbind to filtered dataframe
df.mr <- cbind(df.mr, df.utm) |> 
  st_drop_geometry()
rm(df.utm)