
#iter <- 


library(lubridate)

# Get datasets
source("src/j-pp-static-get-data_pawsey.R")

# df.mr = survey data (for 3 month time window)
# cny.df = separate df of detected colonies 

cny.df <-  cny.df %>% # Reformat dates
mutate(date = as.Date(date, format = "%d.%m.%Y"))

# Surveys vs colony detections through time
#ggplot()+
#  geom_point(data=df, 
#             aes(x=date, y=presence))+
#  geom_point(data=cny.df, 
#             aes(x=date, y=2),
#             col="red", size=2)+
#  scale_x_date(date_breaks = "3 month",
#               date_minor_breaks = "1 month")+
#  theme(panel.grid.major = element_line(size = 2))

### Subset colony data by time window ###
# Which time frame to use???
cny.mr <- cny.df |> 
  filter(date >= # Can play around with time-frame here...
         up.to.date & # Only look at colonies found 1 month AFTER survey window
        # up.to.date - 30 &  # Include last month of survey window?
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
  
### Repeat this from figure script ###
# Load predictions for matching iter
load(file = paste("out/temp-coda-start_", 
                  iter, 
                  ".RData", sep="")) 
# load shoreline
shoreline <- st_read(file.path(Sys.getenv("DATA_PATH"), 
                               "RDHB/Spatial/Australia_boundary.shp")) |>
  st_transform(crs = 32750)

temp <- as.data.frame(as.matrix(b)) # as.matrix produces the same as above, just prints other parameters as well (alpha, beta.1, beta.2, psi, sigma.det)
temp <- temp %>% # Clunky, but go to df and back to use dplyr to select only loc columns
  select(contains("loc", ignore.case = F))
temp <- as.matrix(temp)

x <- as.vector(temp[,1:(ncol(temp)/2)]) # First half of columns are long values
y <- as.vector(temp[,(ncol(temp)/2+1):ncol(temp)]) # Second half are lat values 
point.data <- data.frame(x = x, y = y) |> subset(x!=0) # Grid of all locations... 

density_est <- MASS::kde2d(point.data$x, point.data$y, n = 100, lims = c(min(df.mr$X), max(df.mr$X), min(df.mr$Y), max(df.mr$Y)))  # 100x100 grid

density_df <- data.frame(
  x = rep(density_est$x, each = length(density_est$x)),
  y = rep(density_est$y, times = length(density_est$y)),
  density = as.vector(t(density_est$z))
)


# Load raster

r <- rast(paste("out/figs/rasters/density_raster_start_", 
                iter, 
                ".tif", sep=""))

#######

# Plot colonies against prediction map
colony_plot <- ggplot(density_df, aes(x=x, y=y)) +
  geom_sf(data = shoreline, # Add coastline
          fill = "lightblue", color = "blue", inherit.aes = FALSE) +
  coord_sf(xlim = range(density_df$x), 
           ylim = range(density_df$y)) +  # Apply bounding box
  geom_raster(aes(fill = density), # Plot density
              interpolate = TRUE) + # interpolate smooths between cells
  geom_contour(aes(z = density), # Add contours
               color = "black", alpha = 0.5) +
  scale_fill_viridis_c(alpha = 0.4) +

  geom_point(data = cny.mr, # Plot colony locations
             aes(x = X, y = Y), 
             colour = "red",
             inherit.aes = FALSE) +
  theme_minimal() +
  labs(title = "2D Kernel Density Estimation",
       x = "X Coordinate",
       y = "Y Coordinate",
       fill = "Density")

ggsave(colony_plot, 
       file = sprintf("out/figs/colonies/colonies-vs-density-map-iter_%s.pdf", iter))

# Extract density values from raster for colony locations

m.point <- matrix(c(cny.mr$X, cny.mr$Y), ncol=2)
# Extracting values from raster?
cny.dens <- as.matrix(extract(x = r, 
                      y = m.point),
                      ncol=1)
cny.mr <- cbind(cny.mr, cny.dens)
colnames(cny.mr)[ncol(cny.mr)] <- 'dens'
cny.mr <- cny.mr[,c('X','Y','dens')]

# Save density per colony
#write.csv(cny.mr,
#          file = sprintf("out/colonies/colonies-dens-iter_%s.csv", iter))


### SUMMARISE DENSITIES, COMPARE WITH SURVEYS, RANDOM ###

# ALSO calculate density for each survey (in df.mr)
# * plot hist of all survey loc densities, and compare against density for each colony
# GET QUANTILE SCORE (q_c) OF EACH COLONY AGAINST ALL SURVEY DENSITIES
# EXAMINE DISTRIBUTION OF q_c

m.point.df <- matrix(c(df.mr$X, df.mr$Y), ncol=2)
# Extracting values from raster?
df.dens <- as.matrix(extract(x = r, 
                              y = m.point.df),
                      ncol=1)
df.mr <- cbind(df.mr, df.dens)
colnames(df.mr)[ncol(df.mr)] <- 'dens'

# Maybe this?
ecdf(df.mr$dens)(cny.mr$dens[1])
# ecdf = Empirical Cumulative Distribution Function
# 1st term is vector of values (densities across surveys)
#2nd term is a value you want to compare against distribution (density for 1 colony loc)

cny.mr$qc <- NA

for(n in 1:nrow(cny.mr)){
  
  cny.mr$qc[n] <- ecdf(df.mr$dens)(cny.mr$dens[n])
#  cny.mr$qc[n] <- ecdf(df.pos$dens)(cny.mr$dens[n])
}

# USE ALL SURVEYS? OR ONLY THOSE WITH POSITIVE DETECTION?
df.pos <- subset(df.mr, presence==1)
hist(df.pos$dens)

ggplot()+
  geom_boxplot(data=df.mr, 
          aes(y=dens))+
  geom_boxplot(data=df.pos, 
               aes(y=dens),
               col="blue", width=0.5)+
  geom_boxplot(data=cny.mr, 
               aes(y=dens),
               col="red", width=0.25)


# What to save ??

cny.mr$data <- 'colonies'
df.mr$data <- 'surveys'

df.all <- rbind(cny.mr[,c('X','Y','dens','data')], # Include qc??
                df.mr[,c('X','Y','dens','data')])

write.csv(df.all,
          file = sprintf("out/dens-colonies-vs-surveys-iter_%s.csv", iter))


# SUM DENSITIES ACROSS COLONIES
# Create summary matrix to fill
cny_summary <- matrix(NA, 
                      ncol=3, nrow=2)
cny_summary[1,1] <- sum(cny.mr$dens) # Sum of densities of colony locs
colnames(cny_summary) <- c('mean_sum', # Sum of densities (average sum across iterations for random points)
                           'sd_mean', # SD of average mean sums
                           'qc') # quantile score for summed colonies vs. distribution of 100 random sums
rownames(cny_summary) <- c(paste('iter_',iter, sep=""), 
                           paste('random_',iter, sep=""))

# Compare against same number of RANDOMLY dropped locations

n.cny <- nrow(cny.mr) # number actual colonies
x.min <- round(min(df.mr$X)) # min & max coords
y.min <- round(min(df.mr$Y)) 
x.max <- round(max(df.mr$X))
y.max <- round(max(df.mr$Y)) 
shore.x.min <- as.numeric(colnames(m.shore)[1]) # Min X cell coord in m.shore
shore.y.min <- as.numeric(rownames(m.shore)[1]) # Min Y cell coord

r.cny_sums <- c(rep(0, times=100)) # empty matrix to fill with summed densities

for(i in 1:100) { # Repeat the following random process 100x

r.cny <- matrix(nrow = n.cny, ncol=3)
r.cny[,1] <- runif(n.cny, x.min, x.max) # Randomly drop colonies in x and y coordinates (same number of colonies as cny.mr)
r.cny[,2] <- runif(n.cny, y.min, y.max)

for(j in 1:n.cny){ # Determine if random points fall on land or in sea

 r.cny[j, 3] <- m.shore[trunc((((r.cny[j, 2]/100*100) - shore.y.min) / 100) + 1), # As in JAGS script,
                        trunc((((r.cny[j, 1]/100*100) - shore.x.min) / 100) + 1)] # Match coord against matrix of landscape (1 = land, 0 = ocean)
}

while(sum(r.cny[, 3]) != n.cny){ # As long as there are 0s (colsum != n.cny) (ie some points in ocean)...
  
  for(j in 1:n.cny){
    
    if(r.cny[j, 3] == 0) { # Redraw coords for those rows
    
      r.cny[j, 1] <- runif(1, x.min, x.max)
      r.cny[j, 2] <- runif(1, y.min, y.max)
      
      r.cny[j, 3] <- m.shore[trunc((((r.cny[j, 2]/100*100) - shore.y.min) / 100) + 1), # As in JAGS script,
                           trunc((((r.cny[j, 1]/100*100) - shore.x.min) / 100) + 1)] # Match coord against matrix of landscape (1 = land, 0 = ocean)
    }
  }
}
  
r.dens <- as.matrix(extract(x = r, # Then extract prob densities from raster for those random points 
                              y = r.cny[,c(1,2)]),
                      ncol=1)
r.cny <- cbind(r.cny, r.dens)
colnames(r.cny)[ncol(r.cny)] <- 'dens'

r.cny_sums[i] <- sum(r.cny[,'dens']) # Sum densities at put in ith row

}

# Look at distribution of summed densities
hist(r.cny_sums)

# Calculate mean and SD of these summed densities over 100 iterations
cny_summary[2,] <- c(mean(r.cny_sums), # Average summed density of random points (over all iterations)
                     sd(r.cny_sums), NA) # SD of this distribution

# What is quantile for summed densities of colony locs (compared to random loc?)
cny_summary[1,3] <- ecdf(r.cny_sums)(cny_summary[1,1])


write.csv(cny_summary,
          file = sprintf("out/colonies/colonies-vs-random-sums-iter_%s.csv", 
                         iter))


