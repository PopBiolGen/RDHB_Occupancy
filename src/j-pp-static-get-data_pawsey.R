
repo <- "https://cran.ms.unimelb.edu.au/" # Set mirror to download packages.
mylib <- "/software/projects/pawsey1163/acoates/setonix/2024.05/r/4.3"

#install.packages(c('MCMCvis', ## Difficulty installing MCMCvis, but might not need it now..?
#                   'rstan',
#                   'StanHeaders',
#                   'loo',
#                   'posterior'), 
#                 repos = repo,
#                 lib = mylib)

library(rjags)
library(readr)
#library(MCMCvis)

##### Get the data #####
# Surveillance data
current.data <- "RDHBSurveillance_2025-05-02.xlsx"
source("src/b-data-organisation.R")
# foraging data
fd <- read_xlsx(path = file.path(data_dir, "Abrol-foraging-data.xlsx")) |> 
  mutate(n.round = round(n.bees)) |>
  tidyr::uncount(n.round)

## Shoreline map matrix - coords are 1 for land, 0 for ocean
m.shore <- as.matrix(read.csv('src/figures/matrix.shoreline.csv',
                              header = T,
                              row.names = 1,
                              check.names = F)) # Make sure disable otherwise puts X in front of colnames

##### Filter and get coordinates #####

first.date <- as.Date(min(df$date.time)) # Earliest survey date

start.date <- first.date + ((iter - 1)*30) # From this date (start every 30 days)
up.to.date <- start.date + 90 # To this date (run for 90 days)
#up.to.date <- as.Date("2025-04-08") # analyse data in the 90 days up to this date

# most recent (last 3 months) data: ".mr"
df.mr <- df |> 
  filter(date.time <= up.to.date & date.time > start.date) 
# get coordinates in UTM 50S
df.utm <- df.mr |>
  st_transform(crs = 32750) |> # UTM 50S
  st_coordinates()
# cbind to filtered dataframe
df.mr <- cbind(df.mr, df.utm) |> 
  st_drop_geometry()
rm(df.utm)