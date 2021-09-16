# How to launch the ROD Portal app via shinyapp.io

###########################
##### GETTING STARTED #####
###########################

# shinyapps.io account information
# login: lkaiser15@gmail.com
# password: usgsrodapp
# list shinyapps.io accounts: accounts(server = NULL)

# install R packages for shiny apps (if not already done)
# install.packages('rsconnect')
library('rsconnect')

# authorize shinyapps.io account 
rsconnect::setAccountInfo(name='rapid-ohia-death',
                          token='F071E993F006968F36A036F3B11E0E9F',
                          secret='bEYvdfSw68acrukdasZZwUFFhJhdhmFqJD+7PX6V')

# # set drive
# drive<-"E:/"
# # set root directory 
# rootDir<-paste0(drive, "Dropbox/PIERC/ROD_updates/")  
# setwd(rootDir)

############################
##### UPDATED APP DATA #####
############################

# RUN EACH MANUALLY FOR DATA QUALITY CHECK #

# # 1. update rod data format
# source("1_rod_data-UPDATE.R")
# 
# # 2. extract variables at rod sites
# source("2_data_extract-UPDATE.R")
# source("2a_all_data_extract-UPDATE.R")
# 
# # 3.thin rod data
# source("3_thin_rod-UPDATE.R")
# 
# # 4. run biomod2 sdm rod model
# source("4_rod_model-UPDATE.R")

# # list files to copy over to app portal
# app_files<-list.files("data/data_to_use/")[which(list.files("data/data_to_use/") != 
#                                                    "previous_versions" & 
#                                                  list.files("data/data_to_use/") !=
#                                                    "ROD_README_data_extract.txt")]
# # copy updated data to app portal folder
# file.copy(paste0("data/data_to_use/", app_files), "rod_portal/",
#           overwrite = TRUE, recursive = FALSE, copy.mode = TRUE)

######################
##### UPDATE APP #####
######################

# load necessary packages
library(shiny)
library(rsconnect)

# # for version control of packages with new errors
# sessionInfo()
# remove.packages("ggplot2")
# require(devtools)
# install_version("ggplot2", version = "2.2.1", repos = "http://cran.us.r-project.org")

# deploy shiny app
# rsconnect::deployApp('C:/Users/lkaiser/Dropbox/PIERC/ROD_updates/', 
#                      account = 'rapid-ohia-death')

# rsconnect::deployApp('G:/Dropbox/PIERC/ROD_distribution/rod_portal', 
#                      account = 'rapid-ohia-death')

# --- OR --- #

# set drive
drive<-"E:/"
# set working directory to app location
appDir<-paste0(drive, "Dropbox/PIERC/ROD_updates/rod_portal/")
setwd(appDir)

# deploy shiny app
deployApp(account = 'rapid-ohia-death')

##### END ROD PORTAL UPDATE #####
# app last successfully deployed 09/15/2021
# https://rapid-ohia-death.shinyapps.io/rod_portal/
# password for data sensitivity: usgsrodapp 
# new url at: hawaiirodresearch.org