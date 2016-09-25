### Deploys Shiny App to Cloud ###
### ROD Real-Time Distribution ###

# set working directory to app location
appDir<-"Y:/PICCC_analysis/ROD_distribution/rod_app/"
setwd(appDir)

# load necessary packages
library(rsconnect)
library(shiny)

# deploy shiny app
deployApp()
