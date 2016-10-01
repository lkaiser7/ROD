### ROD Data Formatting and Collection ###
### Formatting data variables for ROD ###

##################
##### SET UP #####
##################

# load necessary packages
library(rgdal)
library(raster)
library(ggplot2)
library(maptools)

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
# convert Hawaii Island coastlines to utm coordinates
hi_utm<-spTransform(hi_coast, CRS(utmSys))

#################################
##### LOAD AND EXTRACT DATA #####
#################################

#----- FORMATTED ROD DATA -----#

# load rod data formatted from script '1_rod_data.R'
rod<-read.csv(paste0(outDir, "ROD_OCC.csv"), header = T)
# check the loaded data set - 89 observations as of 09/20/2016
head(rod)

# save rod object for final master rod data set
rod_master<-rod

#----- OHIA DISTRIBUTION -----#

# load ohia distribution raster
ohia_rast<-raster(paste0(envDir, "HIGAP_500m_ohia.tif"))
# check loaded ohia map
ohia_rast
# ohia archipelago-wide by Lat/Lon at 500m resolution

# crop to Hawaii Island only
hi_ohia<-crop(ohia_rast, hi_coast)
# extract all values for Hawaii Island
all_ohia_spdf<-as.data.frame(as(hi_ohia, "SpatialPixelsDataFrame"))
# save island-wide extracted ohia distribution
write.csv(all_ohia_spdf, paste0(outDir, "all_hi_ohia.csv"))

# convert Lat/Lon to utm for data extraction
utm_ohia<-project(as.matrix(all_ohia_spdf[ , c("x", "y")]), utmSys)
# save ohia object for final master ohia distribution data set
ohia_vars<-data.frame(all_ohia_spdf, utm_ohia)
# rename column headers
names(ohia_vars)[4:5]<-c("x_utm", "y_utm")

# extract values from Hawaii Island at rod points
rod_ohia<-extract(hi_ohia, cbind(rod$LON, rod$LAT))
# add extracted rod ohia values to master rod data set
rod_master$OHIA<-rod_ohia

#----- BIOCLIMATIC VARIABLES -----#

# list bioclimatic variables
bioclims<-list.files(paste0(envDir, "hrcm_bios/current_500m/"))

# load bioclimatic variable
predictors<-raster(paste0(envDir, "hrcm_bios/current_500m/", bioclims[1]))
# create raster stack of bioclimatic variables 
for(bc in 2:length(bioclims)){
  # load single variables as temporar rasters
  temp_bc<-raster(paste0(envDir, "hrcm_bios/current_500m/", bioclims[bc]))
  # add raster layers to stack
  predictors<-addLayer(predictors, temp_bc)
}
# remove temporary objects
rm(bc, temp_bc)

# check loaded predictors raster stack
predictors
# 19 variables archipelago-wide by Lat/Lon at 500m resolution 

# crop to Hawaii Island only
hi_preds<-crop(predictors, hi_coast)
# extract all values for Hawaii Island
all_bc_spdf<-as.data.frame(as(hi_preds, "SpatialPixelsDataFrame"))
# save island-wide extracted bioclimatic variables
write.csv(all_bc_spdf, paste0(outDir, "all_hi_2000_bioclims.csv"))

# extract bioclims from Hawaii across ohia distribution
ohia_bc<-extract(hi_preds, cbind(ohia_vars$x, ohia_vars$y))
# add extracted bioclims across ohia distribution to master ohia data set 
ohia_vars<-data.frame(ohia_vars, ohia_bc)

# extract bioclims from Hawaii Island at rod points
rod_bc<-extract(hi_preds, cbind(rod$LON, rod$LAT))
# add extracted rod bioclim values to master rod data set
rod_master<-data.frame(rod_master, rod_bc)

#----- VEGETATION -----#

# load vegetation HIGAP raster
veg_rast<-raster(paste0(envDir, "HIGAP.tif"))
# check loaded landcover map
veg_rast
# HIGAP archipelago-wide by UTM at 30x30 resolution

# crop to Hawaii Island only
hi_veg<-crop(veg_rast, hi_utm)

# reduce resolution size of raster by aggregating pixels to 120x120 resolution
hi_veg<-aggregate(hi_veg, fact = 4, fun = modal)
# reproject vegetation raster to lat/long
hi_veg<-projectRaster(hi_veg, crs = LatLon, method = "ngb")
# extract all values for Hawaii Island
all_veg_spdf<-as.data.frame(as(hi_veg, "SpatialPixelsDataFrame"))
# save island-wide extracted vegetation classifications
write.csv(all_veg_spdf, paste0(outDir, "all_hi_veg.csv"))

# load classifcation table to match number values to categories
veg_class<-read.csv(paste0(envDir,"HIGAP Revised - Land Cover Classes FINAL for TABLES.csv"))

# extract vegetation classes from Hawaii across ohia distribution
ohia_veg<-extract(hi_veg, cbind(ohia_vars$x, ohia_vars$y))
# add extracted vegetation classes across ohia distribution to master ohia data set 
ohia_vars$higap<-ohia_veg
# match classification to HIGAP value
ohia_vars$HIGAP<-veg_class$Detailed.Map.Unit..Revised.HI.GAP.[match(ohia_veg, veg_class$Value)]

# extract vegetation classes from Hawaii Island at rod points
rod_veg<-extract(hi_veg, cbind(rod$LON, rod$LAT))
# add extracted rod vegetation values to master rod data set
rod_master$higap<-rod_veg
# match classification to HIGAP value
rod_master$HIGAP<-veg_class$Detailed.Map.Unit..Revised.HI.GAP.[match(rod_veg, veg_class$Value)]

#----- BIOREGIONS -----#

# load bioregions raster
br_rast<-raster(paste0(envDir, "bioregions.tif"))
# check loaded bioregions map
br_rast
# bioregions archipelago-wide by UTM at 100x100 resolution

# crop to Hawaii Island only
hi_br<-crop(br_rast, hi_utm)

# reduce resolution size of raster by aggregating pixels to 120x120 resolution
hi_br<-aggregate(hi_br, fact = 2, fun = modal)
# reproject bioregions raster to lat/long
hi_bioreg<-projectRaster(hi_br, crs = LatLon, method = "ngb")
# extract all values for Hawaii Island
all_br_spdf<-as.data.frame(as(hi_bioreg, "SpatialPixelsDataFrame"))
# save island-wide extracted bioregions
write.csv(all_br_spdf, paste0(outDir, "all_hi_bioreg.csv"))

# load classifcation table to match number values to categories
bioregs<-read.csv(paste0(envDir,"bioregion_classes.csv"))

# extract bioregions from Hawaii across ohia distribution
ohia_br<-extract(hi_bioreg, cbind(ohia_vars$x, ohia_vars$y))
# add extracted bioregions across ohia distribution to master ohia data set 
ohia_vars$bioreg<-ohia_br
# match classification to bioregions FID value
ohia_vars$BIOREG<-bioregs$VOLCANO[match(ohia_br, bioregs$FID)]

# extract bioregions from Hawaii Island at rod points
rod_br<-extract(hi_bioreg, cbind(rod$LON, rod$LAT))
# add extracted rod bioregions to master rod data set
rod_master$bioreg<-rod_br
# match classification to bioregions FID value
rod_master$BIOREG<-bioregs$VOLCANO[match(rod_br, bioregs$FID)]

#----- GEOLOGY -----# 

# load geology shapefile data 
geol_shp<-readOGR(paste0(envDir, "hawaii_geology"), "Haw_St_geo_20070426_region")
# check loaded geology shapefile
geol_shp
# 15 geology layers archipelago-wide by UTM at no resolution

# crop to Hawaii Island only
hi_geol<-geol_shp[which(geol_shp$ISLAND == "Hawaii"), ]
# convert to lat/lon projection
hi_geol<-spTransform(hi_geol, CRS(LatLon))

# create an identifiying reference column for spatial data
# hi_geol@data$id<-rownames(hi_geol@data)
# create data frame from spatial object
# geol_pts<-fortify(hi_geol, region = "id")
# merge fortified object with spatial data for geologic variables
# all_geol_df<-merge(geol_pts, hi_geol@data, by = "id")

# convert the shapefile to a raster based on a standardised background raster
for(var in 1:length(names(hi_geol))){  # set var = 1 for debugging
  # get variable name
  var_nm<-names(hi_geol)[var]
  # rasterize variable
  geo_rast<-rasterize(hi_geol, hi_ohia, var_nm)
  
  # stack rasters
  if(var == 1){
    # create raster stack
    geo_stack<-geo_rast
  }else{
    geo_stack<-stack(geo_stack, geo_rast)    
  }  
}
# rename raster layers
names(geo_stack)<-names(hi_geol)

# extract all values for Hawaii Island
all_geol_spdf<-as.data.frame(as(geo_stack, "SpatialPixelsDataFrame"))

# create classifcation table to match number values to categories from Table 6
age_table<-data.frame(0:14)
age_table$AGE_GROUP<-c("multiple", "0-200 yr", "200-750 yr", "750-1500 yr", "1500-3000 yr",
                   "3000-5000 yr", "5000-10000 yr", "10-30 ka", "30-50 ka", "50-140 ka",
                   "140-780 ka", "780-1000 ka", "1-2 ma", "2-4 ma", "4-6 ma")
names(age_table)<-c("age_group", "AGE_GROUP")
# add number for classification number based on name type
name_table<-data.frame(sort(unique(all_geol_spdf$NAME)), sort(unique(hi_geol$NAME)))
names(name_table)<-c("name", "NAME")

# add classification column based on tables
all_geol_spdf$age_group<-age_table$AGE_GROUP[match(all_geol_spdf$AGE_GROUP, age_table$age_group)]
all_geol_spdf$name<-name_table$NAME[match(all_geol_spdf$NAME, name_table$name)]

# save island-wide extracted geological variables
write.csv(all_geol_spdf, paste0(outDir, "all_hi_geol.csv"))

# extract geological variables from Hawaii across ohia distribution
ohia_geol<-extract(hi_geol, cbind(ohia_vars$x, ohia_vars$y))
# add desired geology columns to master ohia data set
ohia_vars$age_group<-ohia_geol$AGE_GROUP
ohia_vars$AGE_GROUP<-age_table$AGE_GROUP[match(ohia_geol$AGE_GROUP, age_table$age_group)]
ohia_vars$name<-name_table$name[match(ohia_geol$NAME, name_table$NAME)]
ohia_vars$NAME<-ohia_geol$NAME

# extract geological variables from Hawaii Island at rod points
rod_geol<-extract(hi_geol, cbind(rod$LON, rod$LAT))
# add desired geology columns to master rod data set
rod_master$age_group<-rod_geol$AGE_GROUP
rod_master$AGE_GROUP<-age_table$AGE_GROUP[match(rod_geol$AGE_GROUP, age_table$age_group)]
rod_master$name<-name_table$name[match(rod_geol$NAME, name_table$NAME)]
rod_master$NAME<-rod_geol$NAME

#----- ROAD PROXIMITY -----#

# load road data 
road_map<-raster(paste0(envDir, "all_HI_roads_distance_wgs84_resampled_aligned.tif"))
# distances from roads archipelago-wide by Lat/Lon

# extract distance values from Hawaii Island at rod points
rod_roads<-extract(road_map, cbind(rod$LON, rod$LAT))
# add extracted road distances to master rod data set
rod_master$roads<-rod_roads*103.597  # convert degree to distance (old value = 110.70428 km)

#----- SAVE MASTER FILES -----#

# save master files with extracted data for shiny app
write.csv(ohia_vars, paste0(outDir, "ohia_dist_vars.csv"))  # ohia master file 
write.csv(rod_master, paste0(outDir, "rod_dist_vars.csv"))  # rod master file

######################################
##### END OF VARIABLE EXTRACTION #####
######################################