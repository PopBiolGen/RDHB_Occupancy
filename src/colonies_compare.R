#iter<-3
#start.date <- first.date + ((iter - 1)*30) # From this date
#up.to.date <- start.date + 90 # To this date
#start.date <- first.date + ((iter - 1)*90) # From this date
#up.to.date <- start.date + 90 # To this date

current.data <- "RDHBSurveillance_2025-05-02.xlsx"
source("src/b-data-organisation.R")

# 'b-data-org..' creates a separate df of detected colonies 
View(cny.df)

library(lubridate)




cny.df <-  cny.df %>% 
#  mutate(date = dmy_hms(date))
mutate(date = as.Date(date, format = "%d.%m.%Y"))

df <- df %>% 
  #  mutate(date = dmy_hms(date))
  mutate(date = as.Date(date.time, format = "%d.%m.%Y"))

ggplot()+
  geom_point(data=df, 
             aes(x=date, y=presence))+
  geom_point(data=cny.df, 
             aes(x=date, y=2),
             col="red", size=2)+
  scale_x_date(date_breaks = "3 month",
               date_minor_breaks = "1 month")+
  theme(panel.grid.major = element_line(size = 2))

#

### Look at colonies detected in the 1 month AFTER the 3-month survey period...

#e.g.
iter <- 3

first.date <- as.Date(min(df$date.time)) # Earliest survey date
start.date <- first.date + ((iter - 1)*30) # From this date (start every 30 days)
up.to.date <- start.date + 90 # To this date (run for 90 days)
#up.to.date <- as.Date("2025-04-08") # analyse data in the 90 days up to this date
# most recent (last 3 months) data: ".mr"

#### Subset survey data: ###
df.mr <- df |> 
  filter(date.time <= up.to.date & date.time > start.date) 
# get coordinates in UTM 50S
df.utm <- df.mr |>
  st_transform(crs = 32750) |> # UTM 50S
  st_coordinates()
# cbind to filtered dataframe
df.mr <- cbind(df.mr, df.utm) |> 
  st_drop_geometry()

c(min(df.mr$date), max(df.mr$date))

### Subset colony data ###
cny.mr <- cny.df |> 
  filter(date >= # Can play around with time-frame here...
        # up.to.date & # Only look at colonies found AFTER survey window
           up.to.date - 30 &  # Include last month of survey window?
         #  up.to.date - 60 &  # Include last 2 months of survey window?
           date < up.to.date + 30) 
cny.utm <- cny.mr |>
  st_transform(crs = 32750) |> # UTM 50S
  st_coordinates()
# cbind to filtered dataframe
cny.mr <- cbind(cny.mr, cny.utm) |> 
  st_drop_geometry()
# remove duplicated colonies
cny.mr <- cny.mr[!duplicated(cny.mr[, c("X","Y")]),]
  


### Run figures script to get dens_plot & r for iter

ggplot(density_df, aes(x=x, y=y)) +
  geom_sf(data = shoreline, fill = "lightblue", color = "blue", inherit.aes = FALSE) + # Draw coastline
  coord_sf(xlim = range(density_df$x), ylim = range(density_df$y)) +  # Apply bounding box
  geom_raster(aes(fill = density), # use geom_raster when every cell has data (otherwise use geom_tile)
              interpolate = TRUE) + # interpolate smooths between cells
  geom_contour(aes(z = density), color = "black", alpha = 0.5) +
  scale_fill_viridis_c(alpha = 0.4) +
#  geom_point(data = df.mr[df.mr$presence==0,], 
#             aes(x = X, y = Y), 
#             colour = "black",
#             inherit.aes = FALSE,
#             alpha = 0.1) +
  geom_point(data = cny.mr, 
             aes(x = X, y = Y), 
             colour = "red",
             inherit.aes = FALSE) +
  theme_minimal() +
  labs(title = "2D Kernel Density Estimation",
       x = "X Coordinate",
       y = "Y Coordinate",
       fill = "Density")


# Extract density values from raster for colony locations


m.point <- matrix(c(cny.mr$X, cny.mr$Y), ncol=2)
# Extracting values from raster?
cny.dens <- as.matrix(extract(x = r, 
                      y = m.point),
                      ncol=1)
cny.mr <- cbind(cny.mr, cny.dens)
colnames(cny.mr)[ncol(cny.mr)] <- 'dens'

summary(cny.mr$dens)
sum(cny.mr$dens)

# Compare against randomly dropped locations?

n.cny <- nrow(cny.mr)

r.cny <- matrix(nrow = n.cny, ncol=3)

x.min <- round(min(df.mr$X))
y.min <- round(min(df.mr$Y)) #
x.max <- round(max(df.mr$X))
y.max <- round(max(df.mr$Y)) #

shore.x.min <- as.numeric(colnames(m.shore)[1]) # Min X cell coord in m.shore
shore.y.min <- as.numeric(rownames(m.shore)[1]) # Min Y cell coord

r.cny[,1] <- runif(n.cny, x.min, x.max)
r.cny[,2] <- runif(n.cny, y.min, y.max)


#shore.r <- rast(shoreline)
#crs(shore.r) <- "EPSG:32750"
#extract(x = shore.r, 
#        y = r.cny[2,])


for(j in 1:n.cny){

 r.cny[j, 3] <- m.shore[trunc((((r.cny[j, 2]/100*100) - shore.y.min) / 100) + 1), # As in JAGS script,
                        trunc((((r.cny[j, 1]/100*100) - shore.x.min) / 100) + 1)] # Match coord against matrix of landscape (1 = land, 0 = ocean)
}

while(sum(r.cny[, 3]) != n.cny){ # As long as there are 0s (colsum != n.cny)
  
  for(j in 1:n.cny){
    
    if(r.cny[j, 3] == 0) { # Redraw coords for those rows
    
      r.cny[j, 1] <- runif(1, x.min, x.max)
      r.cny[j, 2] <- runif(1, y.min, y.max)
      
      r.cny[j, 3] <- m.shore[trunc((((r.cny[j, 2]/100*100) - shore.y.min) / 100) + 1), # As in JAGS script,
                           trunc((((r.cny[j, 1]/100*100) - shore.x.min) / 100) + 1)] # Match coord against matrix of landscape (1 = land, 0 = ocean)
    }
  }
}
  

r.dens <- as.matrix(extract(x = r, # raster of densities 
                              y = r.cny[,c(1,2)]),
                      ncol=1)
r.cny <- cbind(r.cny, r.dens)
colnames(r.cny)[ncol(r.cny)] <- 'dens'
r.cny <- as.data.frame(r.cny)

# Random
summary(r.cny$dens)
sum(r.cny$dens) 

# vs. Actual colony locs
summary(cny.mr$dens)
sum(cny.mr$dens)  





