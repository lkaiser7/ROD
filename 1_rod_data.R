### ROD Data Formatting and Collection ###
### raw rod in situ data to be updated ###

##################
##### SET UP #####
##################

# set root directory 
rootDir<-"Y:/PICCC_analysis/ROD_distribution/"
setwd(rootDir)

# load necessary packages
library(rgdal)

# path to location data
locDir<-paste0(rootDir, "data/ROD_locations/")
# path to environmental data
envDir<-paste0(rootDir, "data/environmental/")
# path to processed and saved output files
outDir<-paste0("data/data_to_use/")

# set coordinate system to be used with data as needed
LatLon<-'+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0'
utmSys<-'+proj=utm +zone=4 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0'

################################
##### RAW IN SITU ROD FILE #####
################################

#----- NOTE (as of 09/30/2016) -----#
# change the 'rod_file' object to match name of new rod data file
# make sure rod dates in data are manually converted from excel date format
# this data and app builds on an original file with the following format:
# Latitude,N,16,6 | Longitude,N,16,6 | Positive | Date_Sampl,N,16,6

# store file name of in situ rod data to be update
rod_file<-"original_raw_rod"
# load raw rod file to be processed
raw_rod<-read.csv(paste0(locDir, rod_file, ".csv"), header = T, stringsAsFactors = F)
# rename column headers
names(raw_rod)<-c("LAT", "LON", "POS_NEG", "DATE")
# check the final data set for formatting - 262 observations as of 09/20/2016
head(raw_rod)

####################################
##### FORMAT ORIGINAL ROD DATA #####
####################################

# remove points where no rod was detected (POS_NEG = 0) - 88 observations as of 09/20/2016
only_rod<-raw_rod[which(raw_rod$POS_NEG != "No"), ]
# remove points with missing Lat/Lon data - 86 observations as of 09/20/2016
only_rod<-only_rod[which(only_rod$LAT != 0 & only_rod$LON != 0), ]

# rename AC observations as Type A
only_rod$POS_NEG[only_rod$POS_NEG == "AC"]<-"A"
# separate observations classified as AC & B
ab<-only_rod[which(only_rod$POS_NEG == "AC & B"), ]  # 3 observations
# rename separate AC & B types as B
ab$POS_NEG[ab$POS_NEG == "AC & B"]<-"B"
# rename remaining AC & B types as A
only_rod$POS_NEG[only_rod$POS_NEG == "AC & B"]<-"A"
# add separated Type B values
only_rod<-rbind(only_rod, ab)
# change Yes values to Undefined
only_rod$POS_NEG[only_rod$POS_NEG == "Yes"]<-"Undefined"
# 89 observations total as of 09/20/2016 - 24 A, 15 B, and 50 Undefined

# convert Lat/Lon to utm for data extraction
utm_coords<-project(as.matrix(only_rod[ , c("LON", "LAT")]), utmSys)

# create final data frame for rod occurrences 
rod_occ<-data.frame(only_rod$DATE, only_rod$POS_NEG, only_rod$LON, only_rod$LAT, utm_coords)
# rename column headers
names(rod_occ)<-c("DATE", "POS_NEG", "LON", "LAT", "X_UTM", "Y_UTM")

# add jitter to rod points
rod_jitter<-data.frame(jitter(rod_occ$LON, amount = 0), jitter(rod_occ$LAT, amount = 0),
                       jitter(rod_occ$X_UTM, amount = 0), jitter(rod_occ$Y_UTM, amount = 0))
# rename column headers
names(rod_jitter)<-c("j_LON", "j_LAT", "j_X_UTM", "j_Y_UTM")
# add jitter to the rest of the rod data
rod_occ<-data.frame(rod_occ, rod_jitter)

# save final output of formatted data in output directory
write.csv(rod_occ, paste0(outDir, "ROD_OCC.csv"))

########################################
### END ROD LOCATION DATA FORMATTING ###
########################################