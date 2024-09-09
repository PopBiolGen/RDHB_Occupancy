# One ring to rule them all...

# filename of current data set
current.data <- "RDHBSurveillance_2024-09-09.xlsx"

# setup and load these data
source("src/b-data-organisation.R")

# update CPUE figures
source("src/figures/CpUE-figure.R")

# run multi-season spatial occupancy
source("src/f-multi-season-spatial-occupancy.R")

# update occupancy figures
source("src/figures/f-multi-season-spatial-occupancy_figs.R")

# update the report
rmarkdown::render("ms/RDHB_report.Rmd", 
                  output_file = file.path(paste0("RDHB_report_", Sys.Date(), ".html")))
