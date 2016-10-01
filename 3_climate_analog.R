### ROD Climate Envelope Model and Roads ###
### Mapping climate distribution of ROD  ###

##################
##### SET UP #####
##################

# load necessary packages
library(raster)
library(rgdal)
library(biomod2)

# set root directory 
rootDir<-"Y:/PICCC_analysis/ROD_distribution/"
setwd(rootDir)

# path to location data
locDir<-paste0(rootDir, "data/ROD_locations/")
# path to environmental data
envDir<-paste0(rootDir, "data/environmental/")
# path to processed and saved output files
outDir<-paste0("data/data_to_use/")

# set coordinate system to be used with data as needed
LatLon<-'+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0'
utmSys<-'+proj=utm +zone=4 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0'

# load Hawaii Island coasts shapefile
hi_coast<-readOGR(paste0(envDir, "hawaii_coast"), "bi_coast")
# load archipelago coast shapefile
all_coast<-readOGR(paste0(envDir, "hawaii_coast"), "Main_Hawaiian_Islands_simple3")

#################################
##### LOAD AND EXTRACT DATA #####
#################################

#----- FORMAT ROD SPECIES DATA -----#

# load rod data formatted from script '1_rod_data.R'
rod<-read.csv(paste0(outDir, "ROD_OCC.csv"), header = T)
# check the loaded data set - 89 observations as of 09/20/2016
head(rod)

# separate species occurrences by ROD type
a_rod<-rod[which(rod$POS_NEG == "A"), ]    # 10 observations 
b_rod<-rod[which(rod$POS_NEG == "B"), ]    # 15 observations
y_rod<-rod[which(rod$POS_NEG == "Yes"), ]  # 50 observations 

# keep only lat/lon coordinates and add 1 for occurrences in PA column
all_rod<-data.frame(rod$LON, rod$LAT)
all_rod$pa<-rep(1, length(all_rod[,1]))
names(all_rod)<-c("x", "y", "PA")

a_rod<-data.frame(a_rod$LON, a_rod$LAT)
a_rod$pa<-rep(1, length(a_rod[,1]))
names(a_rod)<-c("x", "y", "PA")

b_rod<-data.frame(b_rod$LON, b_rod$LAT)
b_rod$pa<-rep(1, length(b_rod[,1]))
names(b_rod)<-c("x", "y", "PA")

y_rod<-data.frame(y_rod$LON, y_rod$LAT)
y_rod$pa<-rep(1, length(y_rod[,1]))
names(y_rod)<-c("x", "y", "PA")

#----- BIOCLIMATIC VARIABLES -----#

# list bioclimatic variables
bioclims<-list.files(paste0(envDir, "hrcm_bios/current_500m/"))

# load bio1 and bio12 for modeling
bio1<-raster(paste0(envDir, "hrcm_bios/current_500m/", bioclims[1]))
bio12<-raster(paste0(envDir, "hrcm_bios/current_500m/", bioclims[4]))
# stack predictors together
bio_stack<-stack(bio1, bio12)
# check loaded predictors
bio_stack

#----- PSEUDO-ABSENCES FOR PRESENCE ONLY DATA -----#

# reclassify raster values of predictors 
mySREresp<-reclassify(subset(bio_stack, 1, drop = TRUE), c(-Inf, Inf, 0))
# save for each rod type
a_SREresp<-mySREresp
b_SREresp<-mySREresp
y_SREresp<-mySREresp

# assign shared points between bioclims and presence points to 1
mySREresp[cellFromXY(mySREresp, all_rod[,1:2])]<-1
a_SREresp[cellFromXY(a_SREresp, a_rod[,1:2])]<-1
b_SREresp[cellFromXY(b_SREresp, b_rod[,1:2])]<-1
y_SREresp[cellFromXY(y_SREresp, y_rod[,1:2])]<-1

# create negative rasters for absences (opposite of presence points)
neg<-mySREresp == 0
a_neg<-a_SREresp == 0
b_neg<-b_SREresp == 0
y_neg<-y_SREresp == 0

# assign number of pseudo-absence points to be selected
n_PA_pts = 1000

# create matrix of potential pseudo-absence candidate points
potPAs<-rasterToPoints(neg, fun = function(x){x==1})
a_potPAs<-rasterToPoints(a_neg, fun = function(x){x==1})
b_potPAs<-rasterToPoints(b_neg, fun = function(x){x==1})
y_potPAs<-rasterToPoints(y_neg, fun = function(x){x==1})

# extract lat/lon data for potential pseud0-absences
potPAs<-as.data.frame(potPAs[,1:2])
a_potPAs<-as.data.frame(a_potPAs[,1:2])
b_potPAs<-as.data.frame(b_potPAs[,1:2])
y_potPAs<-as.data.frame(y_potPAs[,1:2])

# add NA for pseudo-absences in PA column
potPAs<-cbind(potPAs, PA = rep('NA', dim(potPAs)[1], 1))
a_potPAs<-cbind(a_potPAs, PA = rep('NA', dim(a_potPAs)[1], 1))
b_potPAs<-cbind(b_potPAs, PA = rep('NA', dim(b_potPAs)[1], 1))
y_potPAs<-cbind(y_potPAs, PA = rep('NA', dim(y_potPAs)[1], 1))

# merge pseudo-absences with occurrence points for rod data
all_data<-data.frame(rbind(all_rod, potPAs))
a_data<-data.frame(rbind(a_rod, a_potPAs))
b_data<-data.frame(rbind(b_rod, b_potPAs))
y_data<-data.frame(rbind(y_rod, y_potPAs))

#----- DATA FOR MODEL FITTING -----#

# extract bioclimatic variables for species data
all_bc<-extract(bio_stack, all_data[,1:2], cellnumbers = TRUE)
a_bc<-extract(bio_stack, a_data[,1:2], cellnumbers = TRUE)
b_bc<-extract(bio_stack, b_data[,1:2], cellnumbers = TRUE)
y_bc<-extract(bio_stack, y_data[,1:2], cellnumbers = TRUE)

# create data frame with species points and bioclimatic variables
all_bc_data<-data.frame(cbind(all_data, all_bc))
a_bc_data<-data.frame(cbind(a_data, a_bc))
b_bc_data<-data.frame(cbind(b_data, b_bc))
y_bc_data<-data.frame(cbind(y_data, y_bc))

# find duplicate cells
all_dups<-duplicated(all_bc_data$cells)
a_dups<-duplicated(a_bc_data$cells)
b_dups<-duplicated(b_bc_data$cells)
y_dups<-duplicated(y_bc_data$cells)

# remove data in duplicated cells - PA should be sorted in order of 1, 0, NA
all_bc_data<-all_bc_data[!all_dups, -4]
a_bc_data<-a_bc_data[!a_dups, -4]
b_bc_data<-b_bc_data[!b_dups, -4]
y_bc_data<-y_bc_data[!y_dups, -4]

#----- SURFACE RANGE ENVELOPE -----#

# create function to run model for each strain of rod
rod_sre<-function(rod_type, strain_data){
  
  # create vector of PA points
  myResp = data.frame(strain_data[,3])
  # reclassify 'NA' values as NA
  myResp[myResp == 'NA']<-NA
  # coerce presence-absence points as numeric vector
  myResp<-as.numeric(myResp[,1])
  
  # select xy data coordinates for presences only
  myRespXY = strain_data[which(myResp == 1), c("x", "y")]
  
  # select bioclimatic values only
  myExpl<-bio_stack
  
  # build raster layer based on bioclimatic variables
  myResp<-reclassify(subset(myExpl, 1, drop = TRUE), c(-Inf, Inf, 0))
  # add presence values to reclassified raster
  myResp[cellFromXY(myResp, myRespXY)]<-1
  
  # run surface range envelope
  sre_map<-sre(Response = myResp, Explanatory = myExpl, NewData = myExpl, Quant = 0.0001)
  # save sre raster output
  writeRaster(sre_map, paste0(outDir, rod_type, "_current_CE.tif"), 
              format = "GTiff", overwrite = TRUE)
  
  # save mapped output as jpeg for viewing 
  jpeg(paste0(outDir, rod_type, "_current_CE.jpg"), 
       width = 10, height = 10, units = "in",
       pointsize = 12, quality = 90, bg = "white", res = 300)
  # plot rod analog climates
  plot(sre_map, legend = FALSE, main = paste0(rod_type, " surface range envelope"),
       xlab = "Longitude", ylab = "Latitude")
  legend("topright", legend = "ROD Climatic Range",
         pch = 15, col = "darkgreen")    
  # save image file
  dev.off()
  
  # final function output
  return(sre_map)
}

# run for function testing and debugging
# rod_type = "a"
# strain_data = a_bc_data
# rm(rod_type, strain_data)

# run sre function for each strain of rod
all_sre<-rod_sre("all_strains", all_bc_data)
a_sre<-rod_sre("a_strain", a_bc_data)
b_sre<-rod_sre("b_strain", b_bc_data)
y_sre<-rod_sre("y_strain", y_bc_data)

# create mask of sre points
sre_mask<-all_sre
sre_mask[sre_mask == 0]<-1
# create archipelago-wide mask for background plotting
island_mask<-as.data.frame(as(sre_mask, "SpatialPixelsDataFrame"))
# save archpelago-wide extracted ohia distribution
write.csv(island_mask, paste0(outDir, "all_islands_mask.csv"))

# stack all sre rasters
sre_stack<-stack(all_sre, a_sre, b_sre, y_sre)
names(sre_stack)<-c("all_strain", "a_strain", "b_strain", "y_strain")

#----- MERGE WITH OHIA DISTRIBUTION -----#

# load raster of ohia distribution
ohia_map<-raster(paste0(envDir, "HIGAP_500m_ohia.tif"))
# extract all values for the entire archipelago
all_ohia_pts<-as.data.frame(as(ohia_map, "SpatialPixelsDataFrame"))
# save archpelago-wide extracted ohia distribution
write.csv(all_ohia_pts, paste0(outDir, "all_islands_ohia.csv"))

# mosaic sre raster stack with ohia distribution
sre_ohia<-mosaic(sre_stack, ohia_map, fun = sum)
# 2 = rod & ohia, 1 = rod | ohia, 0 = NA
sre_ohia

# rename stack layers
names(sre_ohia)<-c("all_strain", "a_strain", "b_strain", "y_strain")

# extract raster values from sres
sre_ohia_pts<-as.data.frame(as(sre_ohia, "SpatialPixelsDataFrame"))
# save extracted sre points for ggplots in shiny app
write.csv(sre_ohia_pts, paste0(outDir, "all_islands_rod_current_CE.csv"))

#----- DOFAW POLYGONS -----#

# load DOFAW survey shapefile layer of potential rod areas
dofaw_survey<-readOGR(paste0(envDir, "DOFAW"), "potential_rod_areas")
head(dofaw_survey)

# convert to lat/lon projection
dofaw_shp<-spTransform(dofaw_survey, CRS(LatLon))

# create an identifiying reference column for spatial data
dofaw_shp@data$id<-rownames(dofaw_shp@data)
# create data frame from spatial object
dofaw_pts<-fortify(dofaw_shp, region = "id")
# merge fortified object with spatial data for geologic variables
all_dofaw<-merge(dofaw_pts, dofaw_shp@data, by = "id")
head(all_dofaw)

# plot polygons
# ggplot(all_dofaw, aes(x = long, y = lat, group = group)) + 
#   geom_polygon(fill = NA) + geom_path(color = "white") + coord_fixed()

# create table of rank for potential of rod
dofaw_table<-data.frame(table(all_dofaw$NAME, useNA = "ifany"))
# add numbers to match rank order from lightest to very severe 
dofaw_table$Count<-c(5:7,4, 8, 11, 8, 12, 10, 13:15, 1:3, 16:17)
# match ranked numbers to table
all_dofaw$name<-dofaw_table$Count[match(all_dofaw$NAME, dofaw_table$Var1)]
# replace 17 category with NAs
all_dofaw$name[all_dofaw$name == 17]<-NA

# save final DOFAW output file of survey data
write.csv(all_dofaw, paste0(outDir, "all_islands_dofaw_survey.csv"))

####################################
##### END OF CLIMATE ENVELOPES #####
####################################

# plot polygons with colorw ranking severity of rod potential
# ggplot(all_dofaw, aes(x = long, y = lat, group = group, color = name)) + 
  # geom_polygon(fill = NA) + geom_path() +
  # scale_color_gradient(low = "yellow", high = "red", breaks = c(1, 6, 11, 16),
  #                     labels = c("Light < 10%", "Moderate 10-30%", 
  #                                "Severe 30-50%", "Very Severe > 50%"), 
  #                     na.value = "black") +
  # coord_fixed() + theme_gray()

#----- TEST PLOTS -----#

# ROD POTENTIAL MAP - ADD TO SHINY APP AS NEW TAB
# ggplot(sre_ohia_pts) + 
#   geom_tile(data = island_mask, aes(x = x, y = y, fill = layer), colour = "gray") +
#   geom_tile(data = all_ohia_pts, aes(x = x, y = y, fill = HIGAP_500m_ohia), colour = "darkgreen") +
#   geom_raster(data = subset(sre_ohia_pts, sre_ohia_pts$a_strain == 2), 
#               aes(x = x, y = y, fill = a_strain)) + 
#   scale_fill_continuous(breaks = 2, 
#                         guide = guide_legend(title = "Suitable Climatic\nRange for ROD", 
#                                              title.position = "right", label = F)) +
#   geom_path(data = all_dofaw, aes(x = long, y = lat, group = group, color = name)) +
#   scale_color_gradient(low = "yellow", high = "red", breaks = c(1, 6, 11, 16),
#                        labels = c("Light < 10%", "Moderate 10-30%", 
#                                   "Severe 30-50%", "Very Severe > 50%"), 
#                        na.value = "black") +
#   labs(x = "Longitude", y = "Latitude", colour = "ROD Potential") +
#   ggtitle("Potential ROD Sites") + coord_fixed() + theme_gray()
