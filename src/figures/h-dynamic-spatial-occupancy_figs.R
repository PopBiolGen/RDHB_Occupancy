# script for making figures from coda output of dynamic spatial occupancy model

library(rjags)
library(ggplot2)

load("out/dynamic-spatial-occupancy_RDHB_Nimble_Coda.RData")

densplot(a.n)

# function to calculate the dispersal kernel at each distance, with uncertainty from posterior samples
calculate_kernel <- function(coda.obj, x = seq(0, 5, 0.1)){
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
