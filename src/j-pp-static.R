# To fit static point process model to response data
library(rjags)
library(MCMCvis)

##### Get the data #####
# Surveillance data
current.data <- "RDHBSurveillance_2024-10-09.xlsx"
source("src/b-data-organisation.R")
# foraging data
fd <- read_xlsx(path = file.path(data_dir, "Abrol-foraging-data.xlsx")) |> 
  mutate(n.round = round(n.bees)) |>
  tidyr::uncount(n.round)


##### Filter and get coordinates #####

# most recent (last 6 months) data: ".mr"
df.mr <- df |> 
  filter(time.0 >= (max(time.0)-60)) 
# get coordinates in UTM 50S
df.utm <- df.mr |>
  st_transform(crs = 32750) |> # UTM 50S
  st_coordinates()
# cbind to filtered dataframe
df.mr <- cbind(df.mr, df.utm) |> 
  st_drop_geometry()
rm(df.utm)

##### Organise into data for pp-static model #####
max.c <- 20
# data list for JAGS
data.list <- list(
  fd = fd$distance, # foraging distances
  n.fd = nrow(fd), # number of foraging distance data points
  s.0 = as.matrix(df.mr[, c("X", "Y")]), # survey site locations (matrix of X and Y coordinates)
  sur.lev.var.1 = df.mr$flowering, # detection covariate 1
  sur.lev.var.2 = df.mr$water, # detection covariate 2
  obs.i = df.mr$presence, # presence/absence observations
  x.min = min(df.mr$X),
  y.min = min(df.mr$Y), # possible spatial extent of the species across all time, bounding box
  x.max = max(df.mr$X),
  y.max = max(df.mr$Y), # in metres
  II = nrow(df.mr), # total number of surveys
  M = max.c # maximum number of colonies (data-augmentation approach)
)

##### Organise other pieces for JAGS #####

# initials
init.list <- list(
  psi = 0.1,
  alpha.sig = log(100),
  beta.1.sig = log(1),
  beta.2.sig = log(1),
  l.u.0 = log(50),
  c.j0 = matrix(c(seq(min(df.mr$X), max(df.mr$X), length.out = max.c),
           seq(min(df.mr$Y), max(df.mr$Y), length.out = max.c)),
           nrow = max.c),
  Z = rep(1, data.list$M))

# parameters to monitor
params <- c("psi",
            "alpha.sig",
            "beta.1.sig",
            "beta.2.sig",
            "u.0",
            "loc")

# mcmc settings
nb <- 5000
ni <- 2000
nc <- 3

# the model
a <- jags.model(file = "src/model-files/pp-static-JAGS.txt", 
                data = data.list, 
                inits = init.list,
                n.chains = nc)

update(a, n.iter = nb) # burn in
b<-coda.samples(a, 
                variable.names = params, 
                n.iter = ni, 
                thin = 5)

gelman.diag(b, multivariate = FALSE)


summary(b)

save(b, file = "out/temp-coda.RData")
##### make some plots #####


### Non-loc parameters
temp <- MCMCchains(b, params = c("psi",
                                 "alpha.sig",
                                 "beta.1.sig",
                                 "beta.2.sig",
                                 "u.0"))
pdf(file = "out/static-pairs-parameters.pdf")
pairs(temp)
dev.off()

### density map
temp <- MCMCchains(b, params = c("loc"))
x <- as.vector(temp[,1:(ncol(temp)/2)])
y <- as.vector(temp[,(ncol(temp)/2+1):ncol(temp)])
point.data <- data.frame(x = x, y = y) |> subset(x!=0)
# Estimate 2D kernel density
density_est <- MASS::kde2d(point.data$x, point.data$y, n = 100, lims = c(min(df.mr$X), max(df.mr$X), min(df.mr$Y), max(df.mr$Y)))  # 100x100 grid
# Convert to a data frame for ggplot
density_df <- data.frame(
  x = rep(density_est$x, each = length(density_est$y)),
  y = rep(density_est$y, times = length(density_est$x)),
  density = as.vector(t(density_est$z))
)

# Plot the density with ggplot2
library(ggplot2)
ggplot(density_df, aes(x=x, y=y)) +
  geom_raster(aes(fill = density), interpolate = TRUE) +
  geom_contour(aes(z = density), color = "black", alpha = 0.5) +
  scale_fill_viridis_c(alpha = 0.4) +
  geom_point(data = df.mr, 
             aes(x = X, y = Y), 
             colour = (df.mr$presence+1),
             inherit.aes = FALSE) +
  theme_minimal() +
  labs(title = "2D Kernel Density Estimation",
       x = "X Coordinate",
       y = "Y Coordinate",
       fill = "Density")
ggsave("out/static-density.map.pdf")
