# to visualise results from the time-varying-detection model
source("src/d-time-varying-detection.R")

###### make a map of the entire dataset ######
z <- map_point_grid(agg_data$df, agg_data$df_grid, summ.col = mean.prop)
z

###### How does detection vary across the time steps? ######
# get summary variable to plot
ts.plot.df <- select(agg_data.ng$df, date.time, time.step) %>%
              group_by(time.step) %>%
              summarise(mid.date = median(date.time)) %>%
              ungroup() %>%
              mutate(time.step2 = time.step^2)
# get model coefficients for detection
tvd.ests <- coef(fit.tvd)
det.ests <- tvd.ests[grepl("p\\(", names(tvd.ests))]
# data to predict onto
n <- nrow(ts.plot.df)
pred.mat <- cbind(int = rep(1, n), ts = ts.plot.df$time.step, ts2 = ts.plot.df$time.step2, h20 = rep(0, n), hr = rep(12, n), hr2 = rep(12^2, n))
preds <- pred.mat %*% det.ests # linear predictor
p <- exp(preds)/(1+exp(preds)) # on the probability scale

pdf(file = "out/detection-date-w-static-occupancy.pdf")
plot(p~ts.plot.df$mid.date, 
     ylab = "Probability of detection (assuming static occupancy)", 
     xlab = "Date", 
     bty = "l")
dev.off()

###### How does detection vary across the day? ######
n <- length(6:18)
pred.mat <- cbind(int = rep(1, n), ts = rep(3, n), ts2 = rep(3^2, n), h20 = rep(0, n), hr = 6:18, hr2 = (6:18)^2)
preds <- pred.mat %*% det.ests # linear predictor
p <- exp(preds)/(1+exp(preds)) # on the probability scale

pdf(file = "out/detection-hour-w-static-occupancy.pdf")
plot(p~pred.mat[,"hr"], 
     xlab = "Hour of the day", 
     ylab = "Detection probability",
     bty = "l")
dev.off()

###### How does occupancy vary through space? ######
psi.ests <- tvd.ests[grepl("psi\\(", names(tvd.ests))]

x.plot.df <- select(agg_data$df_grid, cell.id, mean.dist) %>%
            mutate(pocc.logit = psi.ests[1]+psi.ests[2]*mean.dist,
                   pocc = exp(pocc.logit)/(1+exp(pocc.logit)))
c.pal <- colorNumeric(palette = colorRamp(c("yellow", "red")), domain = x.plot.df$pocc)

map <- leaflet() %>%
  addTiles() %>%
  addPolygons(data = x.plot.df,
              color = "blue",          # Color of the polygon borders
              weight = 2,              # Weight of the polygon borders
              fillColor = c.pal(x.plot.df$pocc),      # Fill color of the polygons
              fillOpacity = 0.5,       # Opacity of the fill color
              # fillOpacity = results_pol$prob,       # Opacity of the fill color - note too dark to be meaningful
              label = ~cell.id,             # Labels for the polygons
              group = "Probability of occupancy")   %>%          
  
  addScaleBar(position = "bottomleft")%>%
  fitBounds(116.72, -20.82, 116.76, -20.42)%>%  # Set the bounding box
  
  addLegend(
    pal = c.pal,
    values = x.plot.df$pocc,
    title = "Values",
    position = "bottomright") %>%

  addLayersControl(
    overlayGroups = c("Probability of occupancy"),
    options = layersControlOptions(collapsed = FALSE))

map
