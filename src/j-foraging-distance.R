## ---------------------------
##
## Script name: j-foraging-distance.R
##
## Purpose of script: Fit a couple of potential kernels to the foraging distance data extracted from
  # graph in Abrol 1988, "Foraging range of subtropical bees..."
##
## Author: Ben Phillips
##
## Date Created: 2025-03-20
##
##
## ---------------------------
##
## Notes:
##
## --------------------------
## load up the packages we will need 
library(readxl)
library(MASS)
library(tidyr)
## ---------------------------

## load up our functions into memory
source("src/a-setup.R")
## ---------------------------


# read data digitised from Abrol 1988 figure: "foraging-data.jpg
fd <- read_xlsx(path = file.path(data_dir, "Abrol-foraging-data.xlsx"))

fd.long <- fd |> mutate(n.round = round(n.bees)) |>
            uncount(n.round) # expand to one distance row per observation

k1 <- function(x, sd){dnorm(x, mean = 0, sd = sd)} # normal with mean of 0
fit.k1 <- fitdistr(fd.long$distance, densfun = k1, start = list(sd = 35))
k1.line <- data.frame(x = 1:750, 
                 y = k1(1:750, 
                           sd = fit.k1$estimate["sd"])*50*2) # correct for dx and reflection

k2 <- function(x, s, df){dt(x/s, df)/s} # non-standard t distribution 
fit.k2 <- fitdistr(fd.long$distance, densfun = k2, start = list(s = 35, df = 40))
k2.line <- data.frame(x = 1:750, 
                      y = k2(1:750, 
                             s = fit.k2$estimate["s"],
                             df = fit.k2$estimate["df"])*50*2) # correct for dx and reflection

ggplot() +
  geom_col(data = fd, aes(x = distance, y = n.bees/sum(n.bees)), alpha = 0.6) +
  geom_line(data = k1.line, aes(x = x, y = y)) +
  geom_line(data = k2.line, aes(x = x, y = y), linetype = "dashed") +
  xlab("Distance") +
  ylab("Density of observations") +
  theme_minimal()

ggsave(filename = "out/kernel-fit.png")
# so t-distribution doesn't give us much change.  Stick with normal.
