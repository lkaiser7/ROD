### UI for ROD Shiny App ###

##################
##### SET UP #####
##################

# load necessary packages
library(shiny)

# set working directory to app location
appDir<-"Y:/PICCC_analysis/ROD_distribution/rod_app/"
setwd(appDir)

# load master data set
rod<-read.csv("rod_dist_vars.csv", header = TRUE)

##############
##### UI #####
##############

# run shiny UI with tabs and sidebar layout pages
shinyUI<-navbarPage(inverse = TRUE, 
                    tags$div(tags$b("Rapid Ohia Death (ROD) on Hawai'i Island")),
                    tabsetPanel(position = "above", type = "tabs", 
                                
                                # map of rod distribution data
                                tabPanel("ROD Distribution", 
                                         sidebarLayout(
                                           sidebarPanel(
                                             checkboxGroupInput("type1", 
                                                                "Select Ceratocystis fimbriata strain:",
                                                                choices = c("A", "B", "BOTH" = "Yes"), 
                                                                selected = c("A", "B", "Yes"))
                                           ),
                                           mainPanel(plotOutput("RODmap")))
                                ),
                                # plots of rod bioclimatic variables
                                tabPanel("Climatic Associations",
                                         sidebarLayout(
                                           sidebarPanel(
                                             checkboxGroupInput("type2", 
                                                                "Select Ceratocystis fimbriata strain:",
                                                                choices = c("A", "B", "BOTH" = "Yes"), 
                                                                selected = c("A", "B", "Yes")),
                                             radioButtons("bioclim", "Select Variable:",
                                                          choices = c("Annual Mean Temperature" = "bio1",
                                                                      "Annual Precipitation" = "bio12",
                                                                      "Temperature Seasonality" = "bio4",
                                                                      "Precipitation Seasonality" = "bio15",
                                                                      "Temperature Annual Range" = "bio7",
                                                                      "Mean Diurnal Range" = "bio2"))                      
                                           ),
                                           mainPanel(plotOutput("bc_plots")))
                                ),
                                # plots of other climate factors
                                tabPanel("Vegetation Characteristics",
                                         sidebarLayout(
                                           sidebarPanel(
                                             checkboxGroupInput("type3", 
                                                                "Select Ceratocystis fimbriata strain:",
                                                                choices = c("A", "B", "BOTH" = "Yes"), 
                                                                selected = c("A", "B", "Yes")),
                                             radioButtons("ecovar", "Select Data to Display:",
                                                          choices = c("Bioregions" = "bioreg",
                                                                      # "Vegetation Classes" = "COMM",
                                                                      # "Biome Type" = "BIOME",
                                                                      "Land Cover Type" = "higap"))
                                           ),
                                           mainPanel(plotOutput("bv_plots")))
                                ),
                                tabPanel("Geologic Classifications",
                                         sidebarLayout(
                                           sidebarPanel(
                                             checkboxGroupInput("type4", 
                                                                "Select Ceratocystis fimbriata strain:",
                                                                choices = c("A", "B", "BOTH" = "Yes"), 
                                                                selected = c("A", "B", "Yes")),
                                             radioButtons("geol", "Select Data to Display:",
                                                          choices = c("Age Range" = "geol_age",
                                                                      "Lithology" = "geol_type"))
                                           ),
                                           mainPanel(plotOutput("geol_plots")))
                                )
                    ) # end tabsetPanel
                    #), tabPanel("Road Maps")
                    # end shiny UI with navbarPage
)
### END UI ###