
# ROD
Rapid Ohia Death R Shiny Portal

http://hawaiirodresearch.org/ 

UPDATED data extraction and analysis for shiny app

### DATA PROCESSING ###
"1_rod_data-UPDATE.R"

formats ROD data collected from field surveys uploaded to ArcGIS Online

monthly data packages received from Brian Tucker (<bjtucker@hawaii.edu>)

"2_HI_data_extract-UPDATE.R"
collects variables at ROD locations for Hawaii Island

(ie. bioclimatic variables, vegetation, geology, etc.)

"2_all_data_extract-UPDATE.R"

collects variables at ROD locations for all islands

"3_this_rod-UPDATE.R"

thins ROD collection points for prevalence points

selects unique samples to avoid oversampling error

"4_rod_model-UPDATE.R"

ensemble modeling approach for threat levels of ROD

uses BIOMOD2 package to project ROD distribution 


### SHINY APP ###
"ui.R" - user inputs selects that set up server outputs

"server.R" - server function to produce app responses to user inputs
