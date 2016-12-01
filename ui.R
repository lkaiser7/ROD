### UI for ROD Shiny App ###

##################
##### SET UP #####
##################

# load necessary packages
library(shiny)
library(shinyjs)

# set working directory to app location
appDir<-"Y:/PICCC_analysis/ROD_distribution/rod_app/"
# setwd(appDir)

# load master data set
rod<-read.csv("rod_dist_vars.csv", header = TRUE)

#----- LOADING -----#

appCSS<-"
#loading-content {
  position: absolute;
  background: #000000;
  opacity: 0.9;
  z-index: 100;
  left: 0;
  right: 0;
  height: 100%;
  text-align: center;
  color: #FFFFFF;
}
"

##############
##### UI #####
##############

# run shiny UI with tabs and sidebar layout pages
shinyUI<-navbarPage(inverse = TRUE, 
                    "Rapid Ohia Death (ROD) on Hawai'i Island",
                    tabsetPanel(position = "above", type = "tabs", 
                                
                                # useShinyjs(),
                                # inlineCSS(appCSS),
                                
                                # loading messge
                                # div(id = "loading-content", h2("Loading...")),
                                # hidden content while loading
                                # hidden(div(id = "app-content", 
                                
                                # home page tab
                                # tabPanel("Home",
                                #          sidebarLayout(
                                #            sidebarPanel(
                                #              helpText("Rapid ??Ohi??a Death (ROD) (Ceratocystis fimbriata) is a fungal disease",
                                #                       "that is currently attaching and killing ??ohi??a (Metrosideros polymorpha)",
                                #                       "on the Big Island of Hawai??i. First identified in 2014, ROD has killed",
                                #                       "hundreds of thousands of trees across 34,000 acres. This app is designed",
                                #                       "to help analyze the most up-to-date ROD data to better address this issue.",
                                #                       "Please allow time for the tabs to load for the first time once selected.",
                                #                       "Photos show impacts of the disease and ??ohi??a remains of infected trees",
                                #                       "Photos JB Friday")
                                #              ),
                                #          mainPanel(plotOutput("ROD")))
                                # ),
                                
                                # map of rod distribution data
                                tabPanel("ROD Distribution", 
                                         sidebarLayout(
                                           sidebarPanel(
                                             checkboxGroupInput("type1", 
                                                                "Select Ceratocystis fimbriata strain:",
                                                                choices = c("A", "B", "Undefined" = "Undefined"), 
                                                                selected = c("A", "B", "Undefined"))
                                           ),
                                           mainPanel(plotOutput("RODmap")))
                                ),
                                
                                # plots of rod climate envelopes
                                tabPanel("Threat of ROD",
                                         sidebarLayout(
                                           sidebarPanel(
                                             radioButtons("sres", "Select Ceratocystis fimbriata strain:",
                                                          choices = c("All Strains" = "all_strain",
                                                                      "A" = "a_strain",
                                                                      "B" = "b_strain",
                                                                      "Undefined" = "u_strain")),
                                             radioButtons("dofaw", "Add DOFAW Aerial Surveys To Map?",
                                                          choices = c("NO", "YES")) 
                                           ),
                                           mainPanel(plotOutput("sre_plots")))
                                ),
                                                                
                                # plots of rod bioclimatic variables
                                tabPanel("Climatic Associations",
                                         sidebarLayout(
                                           sidebarPanel(
                                             checkboxGroupInput("type2", 
                                                                "Select Ceratocystis fimbriata strain:",
                                                                choices = c("A", "B", "Undefined" = "Undefined"), 
                                                                selected = c("A", "B", "Undefined")),
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
                                                                choices = c("A", "B", "Undefined" = "Undefined"), 
                                                                selected = c("A", "B", "Undefined")),
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
                                                                choices = c("A", "B", "Undefined" = "Undefined"), 
                                                                selected = c("A", "B", "Undefined")),
                                             radioButtons("geol", "Select Data to Display:",
                                                          choices = c("Age Range" = "geol_age",
                                                                      "Lithology" = "geol_type"))
                                           ),
                                           mainPanel(plotOutput("geol_plots")))
                                ),
                                
                                # plots accessibility to roads
                                tabPanel("Road Accessibility",
                                         sidebarLayout(
                                           sidebarPanel(
                                             checkboxGroupInput("type5", 
                                                                "Select Ceratocystis fimbriata strain:",
                                                                choices = c("A", "B", "Undefined" = "Undefined"), 
                                                                selected = c("A", "B", "Undefined"))
                                           ),
                                           mainPanel(plotOutput("road_plots")))
                                )
                                
                                # )) # end hidden content
                    ) # end tabsetPanel
                    #), tabPanel("Road Maps")
                    # end shiny UI with navbarPage
)
### END UI ###