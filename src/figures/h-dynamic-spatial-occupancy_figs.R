# script for making figures from coda output of dynamic spatial occupancy model

library(rjags)

load("out/dynamic-spatial-occupancy_RDHB_Nimble_Coda.RData")

densplot(a.n)
