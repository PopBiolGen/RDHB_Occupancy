# script for making figures from coda output of dynamic spatial occupancy model

library(rjags)
library(ggplot2)

load("out/dynamic-spatial-occupancy_RDHB_Nimble_Coda.RData")
gelman.diag(a.n)
core.pars <- as.matrix(a.n)
core.pars <- core.pars[, !grepl(pattern = "o.t", colnames(core.pars))]
pairs(core.pars)
densplot(a.n)

##### Dispersal plot #####

# function to calculate the dispersal kernel at each distance, with uncertainty from posterior samples
calculate_kernel <- function(coda.obj, x = seq(0, 9, 0.1)){
  kernel_func <- function(k, x){exp(-k*x)} # kernel function
  samples <- as.matrix(coda.obj)[, "k"] # samples of the parameter
  kern.samps.x <- sapply(X = samples, FUN = kernel_func, x = x) # estimates at each x for each sample
  summ.samps <- apply(kern.samps.x, 1, quantile, probs = c(0.025, 0.25, 0.5, 0.75, 0.975))
  data.frame(x = x, t(summ.samps))
}

# plot the kernel
c.k <- calculate_kernel(a.n)
p <- ggplot(data = c.k, aes(x = x)) +
  geom_errorbar(aes(ymin = X2.5.,
                    ymax = X97.5.),
                width = 0,
                col = "lightgrey") +
  geom_errorbar(aes(ymin = X25.,
                    ymax = X75.),
                width = 0,
                col = "darkgrey",
                linewidth = 3) +
  geom_line(aes(x = x, y = X50.)) +
  labs(x = "Distance (km)", y = "Dispersal density") +
  theme_bw()
p
ggsave(filename = "out/dso-dispersal.png")

##### Detection plots #####

# function to return predictions for each day of year given a sample from posterior of parameters
det.pred.fun <- function(samp.vec, water, flower){
  doy <- 1:365
  doy.rad <- doy/365*2*pi
  pred <- samp.vec["det.int"] + 
    water*samp.vec["det.b3"] +
    flower*samp.vec["det.b4"] +
    sin(doy.rad)*samp.vec["det.b1"] + 
    cos(doy.rad)*samp.vec["det.b2"]
  plogis(pred) # back onto the probability scale
}

# return summaries of predictions for detection 
calculate_detection <- function(coda.obj, water = 0, flower = 0){
  samples <- as.matrix(coda.obj)[, grepl(pattern = "det", colnames(coda.obj[[1]]))] # samples of the parameters
  kern.samps.det <- apply(samples, 1, det.pred.fun, water = water, flower = flower)
  summ.samps.det <- apply(kern.samps.det, 1, quantile, probs = c(0.025, 0.25, 0.5, 0.75, 0.975))
  data.frame(doy = 1:365, water = rep(water, 365), flower = rep(flower, 365), t(summ.samps.det))
}

wat.flow.comb <- expand.grid(water = 0:1, flower = 0:1)
plot.det <- c()
for (cc in 1:nrow(wat.flow.comb)){
  plot.det <- rbind(plot.det, calculate_detection(a.n, water = wat.flow.comb[cc, "water"], flower = wat.flow.comb[cc, "flower"]))
}

p <- ggplot(data = plot.det, aes(x = doy)) +
  geom_errorbar(aes(ymin = X2.5.,
                    ymax = X97.5.),
                width = 0,
                col = "lightgrey") +
  geom_errorbar(aes(ymin = X25.,
                    ymax = X75.),
                width = 0,
                col = "darkgrey",
                linewidth = 3) +
  geom_line(aes(y = X50.)) +
  labs(x = "Day of year", y = "Detection probability") +
  scale_fill_manual('Legend Title', values=c('Water', 'No water')) +
  facet_grid(rows = vars(water), 
             cols = vars(flower),
             labeller = labeller(flower = c(`0` = "No flowers", `1` = "Flowers"),
                                 water = c(`0` = "No water", `1` = "Water"))) +
  ylim(0, 1) +
  theme_bw()
p
ggsave(filename = "out/dso-detection-over-time.png")



##### Occupancy over time #####

occ.post <- as.matrix(a.n)[, grepl("o.t", colnames(a.n[[1]]))]
occ.post <- t(apply(occ.post, 2, quantile, probs = c(0.025, 0.25, 0.5, 0.75, 0.975)))
occ.post <- data.frame(time = 1:nrow(occ.post), occ.post)

p <- ggplot(data = occ.post, aes(x = time)) +
  geom_errorbar(aes(ymin = X2.5.,
                    ymax = X97.5.),
                width = 0,
                col = "lightgrey") +
  geom_errorbar(aes(ymin = X25.,
                    ymax = X75.),
                width = 0,
                col = "darkgrey",
                linewidth = 3) +
  geom_point(aes(y = X50.)) +
  labs(x = "Month of response", y = "Number of occupied grid cells") +
  theme_bw()
p

ggsave(filename = "out/dso-mean-occ-over-time.png")


