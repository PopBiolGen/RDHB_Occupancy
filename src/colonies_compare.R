
#### Comparing model predictions (estimated colony density) with READ DATA (actual colony locations)

#iter <- 

#### Get data ####

library(terra)
#library(lubridate)

# Get datasets
source("src/j-pp-static-get-data_pawsey.R")

# df.mr = survey data (for 3 month time window)
# cny.df = separate df of detected colonies 

cny.df <-  cny.df %>% # Reformat dates
mutate(date = as.Date(date, format = "%d.%m.%Y"))

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

# Using full df bounding box
density_est <- MASS::kde2d(point.data$x, point.data$y, n = 100, 
                           lims = c(min(df$X), max(df$X), min(df$Y), max(df$Y)))  # 100x100 grid

density_df <- data.frame(
  x = rep(density_est$x, each = length(density_est$x)),
  y = rep(density_est$y, times = length(density_est$y)),
  density = as.vector(t(density_est$z))
)


# Load raster

r <- rast(paste("out/figs/rasters/density_raster_start_", 
                iter, 
                ".tif", sep=""))

#### Plot colonies against predictions ####

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


#### SUMMARISE DENSITIES, COMPARE WITH SURVEYS, RANDOM ####

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
df.mr$qc <- NA

for(n in 1:nrow(cny.mr)){
  
  cny.mr$qc[n] <- ecdf(df.mr$dens)(cny.mr$dens[n])
#  cny.mr$qc[n] <- ecdf(df.pos$dens)(cny.mr$dens[n])
}

#ggplot(data=df.mr)+
#  geom_histogram(aes(x=dens))+
#  geom_vline(xintercept = c(cny.mr$dens), col="red")


# What to save ??

cny.mr$data <- 'colonies'
df.mr$data <- 'surveys'

df.all <- rbind(cny.mr[,c('X','Y','dens','qc','data')], # Include qc??
                df.mr[,c('X','Y','dens','qc','data')])

write.csv(df.all,
          file = sprintf("out/dens-colonies-vs-surveys-iter_%s.csv", iter))

#### Compare against random locations ####

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
#hist(r.cny_sums)

# Calculate mean and SD of these summed densities over 100 iterations
cny_summary[2,] <- c(mean(r.cny_sums), # Average summed density of random points (over all iterations)
                     sd(r.cny_sums), NA) # SD of this distribution

# What is quantile for summed densities of colony locs (compared to random loc?)
cny_summary[1,3] <- ecdf(r.cny_sums)(cny_summary[1,1])


write.csv(cny_summary,
          file = sprintf("out/colonies/colonies-vs-random-sums-iter_%s.csv", 
                         iter))


### ALT - plot area under cumulative distribution'


#### POSTERIOR PREDICTIVE CHECKS ####

post.draws <- data.frame(x = x, y = y) # Locations of 100 colonies, for EACH of 1200 draws (not sure why 1200??)
post.draws$nd <- rep(c(1:1200), each=100)
post.draws <- post.draws |> subset(x!=0)


# Option 1: keep using density.

# Estimate density for each draw separately
# E.g to TEST
point.data <- subset(post.draws, nd==5)
density_est_v2 <- MASS::kde2d(point.data$x, point.data$y, n = 100, lims = c(min(df.mr$X), max(df.mr$X), min(df.mr$Y), max(df.mr$Y)))  # 100x100 grid
density_df_v2 <- data.frame(
  x = rep(density_est_v2$x, each = length(density_est_v2$x)),
  y = rep(density_est_v2$y, times = length(density_est_v2$y)),
  density = as.vector(t(density_est_v2$z))
)
# Where did model drop colonies IN ONE ITERATION, compared to actual locations colonies found 
ggplot(density_df_v2, aes(x=x, y=y)) +
  geom_sf(data = shoreline, # Add coastline
          fill = "lightblue", color = "blue", inherit.aes = FALSE) +
  coord_sf(xlim = range(density_df_v2$x), 
           ylim = range(density_df_v2$y)) +  # Apply bounding box
  geom_raster(aes(fill = density), # Plot density
              interpolate = TRUE) + # interpolate smooths between cells
  geom_contour(aes(z = density), # Add contours
               color = "black", alpha = 0.5) +
  scale_fill_viridis_c(alpha = 0.4) +
  
  geom_point(data = point.data, # Plot colony locations
             aes(x = x, y = y), 
             colour = "black",
             inherit.aes = FALSE) +
  
  geom_point(data = cny.mr, # Plot colony locations
             aes(x = X, y = Y), 
             colour = "red",
             inherit.aes = FALSE) +
  theme_minimal()

# An option would be to then simulate NEW, RANDOMLY DETERMINED colony locations, 
# * where the probability of being pulled is dependent on the DENSITY from one draw *
# Do this a bunch of times, for a bunch of draws
# Match these locations to new/previous/total (?) density values, and compare against actual locations
# BUT ALL THIS SEEMS TO BE GOING BACK AND FORTH BETWEEN STEPS CALCULATING DENSITY
# Better to just cut out density entirely, and instead directly use data given!

# Option 2: 
# Work directly with hive location draws
# Treat each posterior draw’s hive locations as a simulated replicate dataset.
# Compare summary statistics (e.g., number of hives near observed hives, clustering indices, distance to nearest neighbour, etc.) between the simulated datasets and the observed dataset.
# This skips the KDE step entirely and keeps everything in the point-process framework


hives_obs <- cny.mr[,c(1:2)]
colnames(hives_obs) <- c('x','y')
# Split post.draws into a list of data frames (each element of list corresponding to a draw)
hives_sims <- split(post.draws, f = post.draws$nd)

# From ChatGPT
# What this does
# Defines a window (study area).
# Converts observed and simulated hive sets into ppp objects (spatstat point patterns).
# Computes simple summary statistics (examples: count, mean nearest-neighbour distance, Clark-Evans index).
# Builds distributions of these stats across posterior draws.
# Compares the observed statistic against the simulated distribution.
## -------------------------
## Step 1: Define the observation window
## -------------------------
# We'll take the convex hull or bounding box of all observed + simulated points
all_x <- c(hives_obs$x, unlist(lapply(hives_sims, function(df) df$x)))
all_y <- c(hives_obs$y, unlist(lapply(hives_sims, function(df) df$y)))
obs_window <- spatstat.geom::owin(xrange = range(all_x), yrange = range(all_y))

## -------------------------
## Step 2: Convert to point pattern objects
## -------------------------
pp_obs <- spatstat.geom::ppp(hives_obs$x, hives_obs$y, window = obs_window)
pp_sims <- lapply(hives_sims, function(df) {
  spatstat.geom::ppp(df$x, df$y, window = obs_window)
})

## -------------------------
## Step 3: Choose test statistics
## -------------------------
# Example 1: Number of points
stat_count <- function(pp) { pp$n }

# Example 2: Mean nearest-neighbour distance
stat_nnd   <- function(pp) { mean(nndist(pp)) }

# Example 3: Clark-Evans aggregation index (ratio of observed NN distance to CSR expectation)
stat_clark_evans <- function(pp) {
  mean(nndist(pp)) / (0.5 / sqrt(pp$n / area.owin(pp$window)))
}

## -------------------------
## Step 4: Compute stats for observed + simulated
## -------------------------
obs_stats <- c(
  count = stat_count(pp_obs),
  nnd   = stat_nnd(pp_obs),
  clark = stat_clark_evans(pp_obs)
)

sim_stats <- data.frame(
  count = sapply(pp_sims, stat_count),
  nnd   = sapply(pp_sims, stat_nnd),
  clark = sapply(pp_sims, stat_clark_evans)
)

## -------------------------
## Step 5: Compare distributions
## -------------------------
par(mfrow = c(1,3))
for (stat in names(obs_stats)) {
  hist(sim_stats[[stat]], main = paste("PPC for", stat),
       xlab = stat, breaks = 20, col = "lightgray", border = "white")
  abline(v = obs_stats[[stat]], col = "red", lwd = 2)
}

## Bayesian p-values
bayes_p <- colMeans(sweep(sim_stats, 2, obs_stats, FUN = ">="))
print(bayes_p)

## -------------------------
## Optional: ggplot version
## -------------------------
sim_stats_long <- tidyr::pivot_longer(sim_stats, everything(),
                                      names_to = "stat", values_to = "value")

ggplot(sim_stats_long, aes(x = value)) +
  geom_histogram(bins = 20, fill = "grey70") +
  facet_wrap(~stat, scales = "free") +
  geom_vline(data = data.frame(stat = names(obs_stats), value = obs_stats),
             aes(xintercept = value), colour = "red") +
  theme_minimal() +
  labs(title = "Posterior Predictive Checks on Hive Point Patterns")


#The Bayesian p-value (pB) is just a tail probability:
#𝑝𝐵= Pr(𝑇sim ≥𝑇obs  ∣ posterior)
# Near 0.5: The observed statistic is right in the middle of the simulated distribution → model is consistent with the data (at least for that statistic).
# Near 0 or 1: The observed statistic is extreme compared to simulations → the model is not capturing that aspect of the data.

#Important nuance:
#Bayesian p-values are not like frequentist p-values. They don’t measure “probability the null is true.”
# Instead, they measure how surprising the observed data is, if the model were true. So they’re mainly a diagnostic for model fit, not a formal hypothesis test.


### The problem with this approach is there isn't an inherent summary statistic that I wish to measure!
# Number of colonies or nearest distance to neighbour isn't relevant for my draws!
# I just want to know if the predicted locations (over all draws) are close (by what metric?) to the actual locations?
# Measuring something like density (calculated by averaging over all probably 

#########

### ANOTHER OPTION ####
# Randomly drop survey points in study area -> probability of drop SCALED BY predicted density
# i.e. simulate surveys that are informed by model output
# Repeat many iterations
# * How often do these intersect with ACTUAL colony locations?? *
  



##############################################
# colonies over FULL 4 month window
cny.full <- cny.df |> 
  filter(date >= 
           start.date & # Start of model window
           date < up.to.date + 30) # Up to 1 month AFTER model window
cny.utm <- cny.full |>
  st_transform(crs = 32750) |> # UTM 50S
  st_coordinates()
# cbind to filtered dataframe
cny.full <- cbind(cny.full, cny.utm) |> 
  st_drop_geometry()
# remove duplicated colonies
cny.full <- cny.full[!duplicated(cny.full[, c("X","Y")]),]

colony_plot + geom_point(data = cny.full, # Plot colony locations
                         aes(x = X, y = Y), 
                         colour = "red",
                         inherit.aes = FALSE)

#####
xy <- point.data
nbins <- 40
x.bin <- seq(floor(min(xy[,1])), ceiling(max(xy[,1])), length=nbins)
y.bin <- seq(floor(min(xy[,2])), ceiling(max(xy[,2])), length=nbins)

freq <-  as.data.frame(table(findInterval(xy[,1], x.bin),findInterval(xy[,2], y.bin)))
freq[,1] <- as.numeric(freq[,1])
freq[,2] <- as.numeric(freq[,2])

freq2D <- diag(nbins)*0
freq2D[cbind(freq[,1], freq[,2])] <- freq[,3]

par(mfrow=c(1,2))
image(x.bin, y.bin, freq2D, col=topo.colors(max(freq2D)))
contour(x.bin, y.bin, freq2D, add=TRUE, col=rgb(1,1,1,.7))

palette(rainbow(max(freq2D)))
cols <- (freq2D[-1,-1] + freq2D[-1,-(nbins-1)] + freq2D[-(nbins-1),-(nbins-1)] + freq2D[-(nbins-1),-1])/4
persp(freq2D, col=cols)
