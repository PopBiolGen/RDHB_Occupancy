
# From: https://koordinates.com/layer/113994-australian-bureau-of-statistics-2021-australia-boundary/
file.name <- "RDHB/Spatial/aus_abs_3857.shp" #

shoreline <- st_read(file.path(Sys.getenv("DATA_PATH"), 
                               file.name),
                     crs = 'EPSG:3857')

ggplot() +
  geom_sf(data = shoreline, 
          fill = "lightblue", color = "blue", inherit.aes = FALSE)

shoreline <- st_transform(shoreline, crs = 32750) # Change CRS to match df

ggplot() +
  geom_sf(data = shoreline, 
          fill = "lightblue", color = "blue", inherit.aes = FALSE)


ggplot(density_df, aes(x=x, y=y)) +
  geom_sf(data = shoreline, fill = "lightblue", color = "blue", inherit.aes = FALSE) +
  coord_sf(xlim = range(density_df$x), ylim = range(density_df$y)) +  # Apply bounding box
  geom_raster(aes(fill = density), interpolate = TRUE) +
  geom_contour(aes(z = density), color = "black", alpha = 0.5) +
  scale_fill_viridis_c(alpha = 0.4) +
  geom_point(data = df.mr[df.mr$presence==0,], 
             aes(x = X, y = Y), 
             colour = "black",
             inherit.aes = FALSE,
             alpha = 0.1) +
  geom_point(data = df.mr[df.mr$presence==1,], 
             aes(x = X, y = Y), 
             colour = "red",
             inherit.aes = FALSE) +
  theme_minimal() +
  labs(title = "2D Kernel Density Estimation",
       x = "X Coordinate",
       y = "Y Coordinate",
       fill = "Density")

# Save
st_write(shoreline, "Australia_boundary.shp")

#################

file.name <- "RDHB/Spatial/aus_abs_3857.shp" #
shoreline <- st_read(file.path(Sys.getenv("DATA_PATH"), 
                               file.name),
                     crs = 'EPSG:3857')
# Australia
library(ozmaps)
sf_oz <- subset(ozmap("country")) # All Australia

# All Aus
ggplot() +
  geom_sf(data = shoreline, fill = "lightblue", color = "blue", inherit.aes = FALSE) +
  geom_sf(data = sf_oz, fill = "pink")


df$geometry <- st_transform(df$geometry, crs = crs(shoreline))

coords <- st_coordinates(df$geometry)
summary(coords[,'X'])
summary(coords[,'Y'])

# study area
ggplot() +
  geom_sf(data = sf_oz, fill = "pink")+
  geom_sf(data = shoreline, fill = NA, color = "blue", inherit.aes = FALSE) +
  coord_sf(ylim = c(min(coords[,'Y']), min(coords[,'Y'])),
                  xlim = c(min(coords[,'X']), min(coords[,'X'])))

#pnt <- df[1,"geometry"]
#pnt <- st_transform(pnt, crs = crs(shoreline))
#

ggplot()+
  geom_sf(data = shoreline, fill = "lightblue", color = "blue", 
          inherit.aes = FALSE)+
  
  geom_point(aes(x=coords[1,'X'], y=coords[1,'Y']), col="red", size=4) +
#  geom_point(data=pnt, aes(x=st_coordinates(pnt)[1],
 #                          y=st_coordinates(pnt)[2]))+
coord_sf(xlim = c(min(coords[,'X']),
                  max(coords[,'X'])),
         ylim = c(min(coords[,'Y']),
                 max(coords[,'Y'])))
 # coord_sf(xlim = c(min(st_coordinates(df$geometry)[1]),
  #                  max(st_coordinates(df$geometry)[1])),
   #        ylim = c(min(st_coordinates(df$geometry)[1]),
    #                max(st_coordinates(df$geometry)[1])))
                    
                      
#Land
(as.numeric(st_intersects(df$geometry[1], shoreline))) # Can also use this against raster??

#Ocean
ocean <- st_as_sf(expand.grid(116.74, -20.56), coords=1:2, # Convert coords to sf object
                        crs=st_crs(sf_oz)) # Coordinate reference system
ocean <- st_transform(ocean, crs = crs(shoreline))

ggplot()+
  geom_sf(data = shoreline, fill = "lightblue", color = "blue", 
          inherit.aes = FALSE)+
  
  geom_point(aes(x=coords[1,'X'], y=coords[1,'Y']), col="red", size=4) +
  
  geom_point(aes(x=st_coordinates(ocean)[1], y=st_coordinates(ocean)[2]), col="navy", size=4) +
  
  coord_sf(xlim = c(min(coords[,'X']),
                    max(coords[,'X'])),
           ylim = c(min(coords[,'Y']),
                    max(coords[,'Y'])))

(as.numeric(st_intersects(ocean, shoreline)))

### THIS WORKS! Just need to compare dropped colony against shoreline (which should already be in the same CRS)

# JAGS

df.mr$X[1]
df.mr$Y[1]

pnt <- expand.grid(df.mr$X[1], df.mr$Y[1] + 2000)

pnt <- st_as_sf(pnt, coords=1:2, # Convert coords to sf object
                crs=st_crs(shoreline)) # Coordinate reference system
pnt



coords <- st_transform(df$geometry, crs = crs(shoreline))
coords <- st_coordinates(coords)

ggplot()+
  geom_sf(data = shoreline, fill = "lightblue", color = "blue", 
          inherit.aes = FALSE) +
  
  geom_point(aes(x=st_coordinates(pnt)[1], 
                 y=st_coordinates(pnt)[2]), 
             col="navy", size=4) +
  
  coord_sf(xlim = c(min(coords[,'X']),
                    max(coords[,'X'])),
           ylim = c(min(coords[,'Y']),
                    max(coords[,'Y'])))

(as.numeric(st_intersects(pnt, shoreline)))


### JAGS ####

# Need to create matrix of 0s and 1s (ocean/land)
# Feed into JAGS
# - determine which cell random colony belongs to,
# - multiply z[i] by that value
# - *Make this* z2[i] which is the variable that turns colonies on/off



shoreline <- st_read(file.path(Sys.getenv("DATA_PATH"), 
                               "RDHB/Spatial/Australia_boundary.shp"),
                     crs = 32750)

# CHECK - what is resolution of study area??

df.sa <- df |> # Use complete dataset (ALL survey points)
  st_transform(crs = 32750) |> # UTM 50S
  st_coordinates()
df.sa <- as.data.frame(df.sa)

# Check resolution
library(geosphere)
coords <- as.matrix(st_coordinates(df$geometry))
distm(coords[1,], coords[2,])
# 63.5m distance between points, using coords

c1 <- df[1,] |>
  st_transform(crs = 32750) |> # UTM 50S
  st_coordinates()
c2 <- df[2,] |>
  st_transform(crs = 32750) |> # UTM 50S
  st_coordinates()
sqrt((c1[,'X'] - c2[,'X'])^2 + (c1[,'Y'] - c2[,'Y'])^2)
# 63.5 distance beteen points using CRS = 32750
# Confirm points = metres

max(df.sa$X) - min(df.sa$X)
max(df.sa$Y) - min(df.sa$Y)
# ~15km x 15km study area

# Use 100m resolution? Better than nothing, and saves memory!

ceiling(max(df.sa$X)) - floor(min(df.sa$X))

library(plyr)
library(dplyr)
df.xmin <- round_any(min(df.sa$X), 100, f=floor)
df.xmax <- round_any(max(df.sa$X), 100, f=ceiling)
df.ymin <- round_any(min(df.sa$Y), 100, f=floor)
df.ymax <- round_any(max(df.sa$Y), 100, f=ceiling)
x.dim <- length(seq(df.xmin, df.xmax, by = 100))
y.dim <- length(seq(df.ymin, df.ymax, by = 100))

m.shore <- matrix(0, ncol = x.dim, nrow=y.dim)
colnames(m.shore) <- c(seq(df.xmin, df.xmax, by = 100))
rownames(m.shore) <- c(seq(df.ymin, df.ymax, by = 100))

##### TEST ####
for (yi in 1:nrow(m.shore)){ 
  for (xi in 1:2){ # TEST with 2 col
    
y.pnt <- as.numeric(rownames(m.shore)[yi])
x.pnt <- as.numeric(colnames(m.shore)[xi])

pnt <- st_as_sf(expand.grid(x.pnt, y.pnt), # Make sure order is right! (x, y)
                coords=1:2, 
                crs=st_crs(shoreline)) 

m.shore[yi, xi] <- (as.numeric(st_intersects(pnt, shoreline)))

  }
}

test.map <- data.frame(cell = c(m.shore[,2]),
                       x = c(rep(as.numeric(colnames(m.shore)[2]), times=nrow(m.shore))),
                       y = as.numeric(rownames(m.shore)))
# Check matches
ggplot()+
  geom_sf(data = shoreline, fill = "lightblue", color = "blue", 
          inherit.aes = FALSE) +
  geom_point(data=test.map,
             aes(x=x, 
                 y=y, 
             col=cell), size=2) +
  coord_sf(xlim = c(df.xmin, df.xmax),
           ylim = c(df.ymin, df.ymax))
# Works
#####

for (yi in 1:nrow(m.shore)){ 
  for (xi in 1:ncol(m.shore)){ # TEST with 2 col
    
    y.pnt <- as.numeric(rownames(m.shore)[yi])
    x.pnt <- as.numeric(colnames(m.shore)[xi])
    
    pnt <- st_as_sf(expand.grid(x.pnt, y.pnt), # Make sure order is right! (x, y)
                    coords=1:2, 
                    crs=st_crs(shoreline)) 
    
    m.shore[yi, xi] <- (as.numeric(st_intersects(pnt, shoreline)))
    
  }
}

m.shore[is.na(m.shore)] <- 0

# The matrix, when viewed as a map, is actually upside down...
# But this doesn't seem to affect the coords or results...
write.csv(m.shore, 'src/figures/matrix.shoreline.csv')


View(m.shore)

### Next input nearest coord into Jags cj

m.shore <- as.matrix(read.csv('src/figures/matrix.shoreline.csv',
                              header = T,
                              row.names = 1,
                              check.names = F)) # Make sure disable otherwise puts X in front of colnames

# test
cj <- c(477000, 7724600)
library(plyr)
library(dplyr)
m.shore[as.character(round_any(cj[2], 100)), # y value as row
        as.character(round_any(cj[1], 100))] # x value as col

# divide by cell size
# trunc()
round(477699/100)*100
#trunc(477699/100)*100
# 

cj <- c(467525, 7713002)


# Need to create matrix of 0s and 1s (ocean/land)
# Feed into JAGS
# - determine which cell random colony belongs to,
# - multiply z[i] by that value
# - *Make this* z2[i] which is the variable that turns colonies on/off


m.shore[as.character(round(cj[2]/100)*100), # y value as row
        as.character(round(cj[1]/100)*100)] # x value as col




m.shore[((trunc(7713025.79/100)*100) - shore.y.min) / 100 + 1,
        ((trunc(467205.725/100)*100) - shore.x.min) / 100 + 1]


cols <- c(as.numeric(colnames(m.shore)))
(cols[146] - cols[1]) / 100 + 1


rows <- c(as.numeric(rownames(m.shore)))
(rows[150] - rows[1]) / 100 + 1



length(cols)


cols <- matrix(0, ncol = ncol(m.shore), nrow=2)
cols[1,] <- as.numeric(colnames(m.shore))
cols[2,] <- c(1:ncol(cols))

library(reshape2)
melt(m.shore)

m.shore <- melt(m.shore)

c.j0[jj, 1] ~ dunif(x.min, x.max)



ggplot()+
  geom_sf(data = shoreline, fill = "lightblue", color = "blue", 
          inherit.aes = FALSE) +
  geom_point(aes(x=475400, y=7718500), col="yellow", size=3)+
  geom_point(aes(x=477000, y=7724600), col="red", size=3)+
  
  coord_sf(xlim = c(df.xmin, df.xmax),
           ylim = c(df.ymin, df.ymax))

# Old idea
for (jj in 1:M){ 
  sf_c <- st_as_sf(c.j0, coords=1:2, crs=st_crs(shoreline))
  land[jj] <- ifelse((as.numeric(st_intersects(sf_c[jj,], shoreline))) == 1,
                     1,
                     0)
}


