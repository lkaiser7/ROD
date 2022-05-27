### UI for ROD Shiny App ###

##################
##### SET UP #####
##################

# load necessary packages
library(shiny)
library(shinyjs)

# set working directory to app location
appDir<-"C:/Users/Lauren/Dropbox/NPS_drive_backup/PIERC/ROD_distribution/rod_portal/"  # laptop
appDir<-"E:/Dropbox/PIERC/ROD_updates/rod_portal/"  # hard drive
# DO NOT RUN setwd() ON SHINY SERVER (DOES NOT EXIST): setwd(appDir)

#----- LOADING -----#

##############
##### UI #####
##############

# run shiny UI with tabs and sidebar layout pages
shinyUI<-navbarPage(inverse = TRUE, 
                    "Rapid Ohia Death (ROD) by Ceratocystis on Hawai'i Island",
                    tabsetPanel(type = "tabs", #position = "above", 
                                
                                # Model Projections
                                tabPanel("Threat",
                                         sidebarLayout(
                                           sidebarPanel(em(code("Please allow up to 1 minute for app to load")), 
                                                        br(), br(),
                                                        radioButtons("sp", "Select Ceratocystis Species",
                                                                     choices = c("C. lukuohia", "C. huliohia")),
                                                        br(),
                                                        em("April 2022 Update: The data and analysis shown here have inherent caveats and limitations 
                                                        that are well described in the associated manuscript:"), 
                                                        uiOutput("fem_doi")
                                                        ),
                                           mainPanel(tags$h4(uiOutput("sdm_info")), 
                                                     plotOutput("sre_image")))
                                           ),
                                
                                # Location Data
                                tabPanel("Distribution", 
                                         sidebarLayout(
                                           sidebarPanel(em(code("Please allow up to 1 minute for app to load")), 
                                                        br(), br(),
                                            passwordInput("password", 
                                                          "For data sensitivity purposes, please enter the password to 
                                                          view the interactive map and precise location data used in this analysis:"),
                                            checkboxGroupInput("type1", 
                                                              "Select Ceratocystis Species:",
                                                              choices = c("C. lukuohia" = "C. lukuohia", 
                                                                          "C. huliohia" = "C. huliohia"), 
                                                              selected = c("C. lukuohia")), #, "C. huliohia")),
                                             radioButtons("neg", "Add Locations of Negative Samples To Map?",
                                                          choices = c("NO", "YES")),
                                            br(),
                                            em("NOTE: The data and analysis shown here have inherent caveats and limitations 
                                                that are well described in the associated manuscript:"),
                                            uiOutput("fem_doi0")
                                           ),
                                           mainPanel(
                                             uiOutput("rod_url"),
                                             #uiOutput("map_url"),
                                             imageOutput("ctahr_img"),
                                             br(), br(), br(), br(), br(), br(), br(),
                                             br(), br(), br(), br(), br(), br(), br(),
                                             htmlOutput("ctahr_map"),
                                             plotOutput("RODmap")))
                                ),
                                                                
                                # plot spread of ROD over time
                                tabPanel("Spread",
                                         sidebarLayout(
                                           sidebarPanel(em(code("Please allow up to 1 minute for app to load")), 
                                                        br(), br(),
                                                        passwordInput("password1", 
                                                                      "For data sensitivity purposes, please enter the password to 
                                                                       view the interactive map and precise location data used in this analysis:"),
                                                        checkboxGroupInput("type",
                                                                           "Select Year of Sample Site:",
                                                                           choices = c("2014", "2015", "2016", "2017", "2018", "2019", "2020"),
                                                                           selected = c("2014", "2015", "2016", "2017", "2018", "2019", "2020")),
                                                        br(),
                                                        em("NOTE: The data and analysis shown here have inherent caveats and limitations 
                                                                that are well described in the associated manuscript:"),
                                                        uiOutput("fem_doi1")
                                           ),
                                           mainPanel(plotOutput("TimeMap")))
                                ),
                                
                                # plots of rod bioclimatic variables
                                tabPanel("Climatic Associations",
                                         sidebarLayout(
                                           sidebarPanel(em(code("Please allow up to 1 minute for app to load")), 
                                                        br(), br(),
                                             radioButtons("type2", "Select Ceratocystis Species:",
                                                           choices = c("C. lukuohia" = "C. lukuohia", "C. huliohia" = "C. huliohia")),
                                             radioButtons("bioclim", "Select Variable:",
                                                          choices = c("Annual Mean Temperature (C)" = "bio1",
                                                                      "Temperature Seasonality (C)" = "bio4",
                                                                      "Min. Temperature of Coldest Month (C)" = "bio6",
                                                                      "Temperature Annual Range (C)" = "bio7",
                                                                      "Mean Temperature of Driest Quarter (C)" = "bio9",
                                                                      "Mean Temperature of Warmest Quarter (C)" = "bio10",
                                                                      "Mean Temperature of Coldest Quarter (C)"  = "bio11",
                                                                      "Annual Precipitation (mm)" = "bio12",
                                                                      "Precipitation of Driest Month (mm)" = "bio14",
                                                                      "Precipitation Seasonality (mm)" = "bio15",
                                                                      "Precipitation of Wettest Quarter (mm)" = "bio16")),
                                             br(),
                                             passwordInput("password2", 
                                                           "For data sensitivity purposes, please enter the password to 
                                                            view the interactive map and precise location data used in this analysis:"),
                                             em("NOTE: The data and analysis shown here have inherent caveats and limitations 
                                                     that are well described in the associated manuscript:"),
                                             uiOutput("fem_doi2")
                                           ),
                                           mainPanel(tags$h4(uiOutput("bc_info")),
                                          fluidRow(column(12, tags$div(imageOutput("bc_plots"),
                                                                 tags$br(), tags$hr()))),
                                           fluidRow(column(12, tags$div(imageOutput("bio_img"),
                                                                 tags$br(), tags$hr()))),
                                           fluidRow(column(12, align = "left", tags$div(plotOutput("bc_map"))))
                                           ))
                                ),
                                
                                # plots of other climate factors
                                tabPanel("Vegetation Characteristics",
                                         sidebarLayout(
                                           sidebarPanel(em(code("Please allow up to 1 minute for app to load")), 
                                                        br(), br(),
                                             radioButtons("type3", "Select Ceratocystis Species:",
                                                          choices = c("C. lukuohia" = "C. lukuohia", "C. huliohia" = "C. huliohia")),
                                             radioButtons("ecovar", "Select Data to Display:",
                                                          choices = c("Bioregions" = "BIOREG",
                                                                      "Land Cover Type" = "HIGAP")),
                                             br(),
                                             p(),
                                             passwordInput("password3", 
                                                           "For data sensitivity purposes, please enter the password to 
                                                            view the interactive map and precise location data used in this analysis:"),
                                             em("NOTE: The data and analysis shown here have inherent caveats and limitations 
                                                     that are well described in the associated manuscript:"),
                                             uiOutput("fem_doi3")
                                           ),
                                           mainPanel(tags$h4(uiOutput("vb_info")),
                                                     fluidRow(column(12, tags$div(imageOutput("vb_img"),
                                                                                  tags$br(), tags$hr()))),
                                                     fluidRow(column(12, tags$div(plotOutput("bv_plot"))))
                                                     ))
                                ),
                                
                                tabPanel("Geologic Classifications",
                                         sidebarLayout(
                                           sidebarPanel(em(code("Please allow up to 1 minute for app to load")), 
                                                        br(), br(),
                                             radioButtons("type4", "Select Ceratocystis Species:",
                                                          choices = c("C. lukuohia" = "C. lukuohia", "C. huliohia" = "C. huliohia")),
                                             radioButtons("geol", "Select Data to Display:",
                                                          choices = c("Age Range" = "AGE_GROUP",
                                                                      "Lithology" = "NAME")),
                                             br(),
                                             passwordInput("password4", 
                                                           "For data sensitivity purposes, please enter the password to 
                                                            view the interactive map and precise location data used in this analysis:"),
                                             em("NOTE: The data and analysis shown here have inherent caveats and limitations 
                                                     that are well described in the associated manuscript:"),
                                             uiOutput("fem_doi4")
                                           ),
                                           mainPanel(tags$h4(uiOutput("geol_info")),
                                                     fluidRow(column(12, tags$div(imageOutput("geol_img"),
                                                                                  tags$br(), tags$hr()))),
                                                     fluidRow(column(12, tags$div(plotOutput("geol_plots"))))
                                                     ))
                                ),
                                
                                # plots accessibility to roads
                                tabPanel("Road Accessibility",
                                         sidebarLayout(
                                           sidebarPanel(em(code("Please allow up to 1 minute for app to load")), 
                                                        br(), br(),
                                             checkboxGroupInput("type5", "Select Ceratocystis Species:",
                                                                choices = c("C. lukuohia" = "C. lukuohia", "C. huliohia" = "C. huliohia"),   
                                                                selected = c("C. lukuohia")),
                                             br(),
                                             passwordInput("password5", 
                                                           "For data sensitivity purposes, please enter the password to 
                                                            view the interactive map and precise location data used in this analysis:"),
                                             em("NOTE: The data and analysis shown here have inherent caveats and limitations 
                                                     that are well described in the associated manuscript:"),
                                             uiOutput("fem_doi5")
                                           ),
                                           mainPanel(plotOutput("road_plots")))
                                )
                    ) # end tabsetPanel
)
### END UI ###
