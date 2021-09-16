### ROD Data Analysis for Publication ###
### Analyzing select variables of ROD ###
### Use presence/absence point ratios ###

##################
##### SET UP #####
##################

# load necessary packages
library(rgdal)
library(raster)
library(ggplot2)
library(rgeos)

# set drive
drive<-"E:/"
# set root directory 
rootDir<-paste0(drive, "Dropbox/PIERC/ROD_updates/")  
setwd(rootDir)

# path to location data
locDir<-paste0(rootDir, "data/ROD_locations/")
# path to environmental data
envDir<-paste0(rootDir, "data/environmental/")
# path to processed and saved output files
dataDir<-paste0("data/data_to_use/")

# set coordinate system to be used with data as needed
LatLon<-'+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0'
utmSys<-'+proj=utm +zone=4 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0'

# load Hawaii Island coasts shapefile
hi_coast<-readOGR(paste0(envDir, "hawaii_coast"), "bi_coast")
# convert Hawaii Island coastlines to utm coordinates
hi_utm<-spTransform(hi_coast, CRS(utmSys))

# load archipelago coasts shapefile
all_coast<-readOGR(paste0(envDir, "hawaii_coast"), "Main_Hawaiian_Islands_simple3")
all_utm<-spTransform(all_coast, CRS(utmSys))

#################################
##### LOAD AND EXTRACT DATA #####
#################################

#----- FORMATTED ROD DATA -----#

# load rod data formatted from script '1_rod_data.R'
rod<-read.csv(paste0(dataDir, "ALL_ROD.csv"), header = T)
head(rod)

# save rod object for final master rod data set
rod_master<-rod

#----- OHIA DISTRIBUTION -----#

# load ohia distribution raster archipelago-wide by Lat/Lon at 500m resolution
ohia_rast<-raster(paste0(envDir, "HIGAP_500m_ohia.tif"))
# crop to Hawaii Island only
# hi_ohia<-crop(ohia_rast, hi_coast)
hi_ohia<-ohia_rast

#----- BIOCLIMATIC VARIABLES -----#

# list bioclimatic variables
bioclims<-list.files(paste0(envDir, "hrcm_bios/current_500m/"))

# load 19 bioclimatic variables archipelago-wide by Lat/Lon at 500m resolution
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
# crop to Hawaii Island only
# hi_preds<-crop(predictors, hi_coast)
hi_preds<-predictors

# extract bioclims from Hawaii Island at rod points
rod_bc<-extract(hi_preds, cbind(rod$X_LON, rod$Y_LAT))
# add extracted rod bioclim values to master rod data set
rod_master<-data.frame(rod_master, rod_bc)

#----- VEGETATION -----#

# load vegetation HIGAP raster archipelago-wide by UTM at 30x30 resolution
veg_rast<-raster(paste0(envDir, "HIGAP.tif"))
# crop to Hawaii Island only
# hi_veg<-crop(veg_rast, hi_utm)
hi_veg<-veg_rast

# reduce resolution size of raster by aggregating pixels to 120x120 resolution
hi_veg<-aggregate(hi_veg, fact = 4, fun = modal)
# reproject vegetation raster to lat/long
hi_veg<-projectRaster(hi_veg, crs = LatLon, method = "ngb")

# load classifcation table to match number values to categories
veg_class<-read.csv(paste0(envDir,"HIGAP Revised - Land Cover Classes FINAL for TABLES.csv"))
# add matching color for each categorical value
veg_class$COLOR<-rev(topo.colors(length(unique(veg_class$Value))))

# extract vegetation classes from Hawaii Island at rod points
rod_veg<-extract(hi_veg, cbind(rod$X_LON, rod$Y_LAT)) # if not converted
rod_veg<-extract(hi_veg, cbind(rod$X_UTM, rod$Y_UTM))
# add extracted rod vegetation values to master rod data set
rod_master$higap<-rod_veg
# match classification to HIGAP value
rod_master$HIGAP<-veg_class$Detailed.Map.Unit..Revised.HI.GAP.[match(rod_veg, veg_class$Value)]
# match color to classification value
rod_master$HIGAP_COLOR<-veg_class$COLOR[match(rod_veg, veg_class$Value)]

#----- BIOREGIONS -----#

# load bioregions raster archipelago-wide by UTM at 100x100 resolution
br_rast<-raster(paste0(envDir, "bioregions.tif"))
# crop to Hawaii Island only
# hi_br<-crop(br_rast, hi_utm)
hi_br<-br_rast

# reduce resolution size of raster by aggregating pixels to 120x120 resolution
hi_br<-aggregate(hi_br, fact = 2, fun = modal)
# reproject bioregions raster to lat/long
hi_bioreg<-projectRaster(hi_br, crs = LatLon, method = "ngb")

# load classifcation table to match number values to categories
bioregs<-read.csv(paste0(envDir,"bioregion_classes.csv"))
# add matching color for each categorical value
bioregs$COLOR<-rev(topo.colors(length(unique(bioregs$FID))))

# extract bioregions from Hawaii Island at rod points
rod_br<-extract(hi_bioreg, cbind(rod$X_LON, rod$Y_LAT))
# add extracted rod bioregions to master rod data set
rod_master$bioreg<-rod_br
# match classification to bioregions FID value
rod_master$BIOREG<-bioregs$VOLCANO[match(rod_br, bioregs$FID)]
# match color to classification value
rod_master$BIOREG_COLOR<-bioregs$COLOR[match(rod_br, bioregs$FID)]

#----- GEOLOGY -----# 

# load geology shapefile data 15 geology layers archipelago-wide by UTM at no resolution
geol_shp<-readOGR(paste0(envDir, "hawaii_geology"), "Haw_St_geo_20070426_region")
# crop to Hawaii Island only
# hi_geol<-geol_shp[which(geol_shp$ISLAND == "Hawaii"), ]
hi_geol<-geol_shp
# convert to lat/lon projection
hi_geol<-spTransform(hi_geol, CRS(LatLon))

# convert the shapefile to a raster based on a standardised background raster
for(var in 1:length(names(hi_geol))){  # set var = 1 for debugging
  # get variable name
  var_nm<-names(hi_geol)[var]
  hi_geol[[var]]<-as.factor(hi_geol[[var]])
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
age_table$COLOR<-rev(cm.colors(length(unique(age_table$AGE_GROUP))))

# add number for classification number based on name type
name_table<-data.frame(sort(unique(all_geol_spdf$NAME)), sort(unique(hi_geol$NAME))[c(-19,-24,-34,-38,-44,-49)])
names(name_table)<-c("name", "NAME")
name_table$COLOR<-rev(cm.colors(length(unique(name_table$NAME))))

# extract geological variables from Hawaii Island at rod points
rod_geol<-extract(hi_geol, cbind(rod$X_LON, rod$Y_LAT))
# add desired geology columns to master rod data set
rod_master$age_group<-rod_geol$AGE_GROUP
rod_master$AGE_GROUP<-age_table$AGE_GROUP[match(rod_geol$AGE_GROUP, age_table$age_group)]
# match color to classification value
rod_master$AGE_COLOR<-age_table$COLOR[match(rod_geol$AGE_GROUP, age_table$age_group)]
# add desired geology columns to master rod data set
rod_master$name<-name_table$name[match(rod_geol$NAME, name_table$NAME)]
rod_master$NAME<-rod_geol$NAME
# match color to classification value
rod_master$NAME_COLOR<-name_table$COLOR[match(rod_geol$NAME, name_table$NAME)]

#----- ROAD PROXIMITY -----#

# load distances from roads archipelago-wide by Lat/Lon
road_map<-raster(paste0(envDir, "all_HI_roads_distance_wgs84_resampled_aligned.tif"))
# crop to Hawaii Island only
# hi_roads<-crop(road_map, hi_coast)
hi_roads<-road_map

# remove background values
# coast_rast<-rasterize(hi_coast, hi_roads)
coast_rast<-rasterize(all_coast, hi_roads)
hi_rds<-hi_roads*coast_rast
# reduce resolution size of raster by aggregating pixels
hi_roads<-aggregate(hi_rds, fact = 2, fun = modal)

# extract distance values from Hawaii Island at rod points
rod_roads<-extract(road_map, cbind(rod$X_LON, rod$Y_LAT))
# add extracted road distances to master rod data set
rod_master$ROADS<-rod_roads*103.597  # convert degree to distance (old value = 110.70428 km)

#----- FENCED AREAS -----#

# load island-wide fences shapefile
all_fence<-readOGR(paste0(envDir, "fenced_units"), "AllUngulateUnit_2018")
# convert coordinate system from utm to lat/lon
all_fence<-spTransform(all_fence, CRS(LatLon))
# crop to Hawaii Island only
# hi_fence<-all_fence[which(all_fence$Island == "Hawaii"), ]
hi_fence<-all_fence
# head(hi_fence@data)
# plot(hi_fence)
# sum(as.numeric(data.frame(table(hi_fence@data$Shape__Are))$Var1))

# check fence status
table(hi_fence@data$FenceStat)
# keep fenced areas and remove unfenced areas
hi_fence<-hi_fence[which(hi_fence@data$FenceStat == "Compliant" | hi_fence@data$FenceStat == "Fenced"), ]

# convert the shapefile to a raster based on a standardised background ohia raster
for(var in 1:length(names(hi_fence))){  # set var = 1 for debugging
  # get variable name
  var_nm<-names(hi_fence)[var]
  hi_fence[[var]]<-as.factor(hi_fence[[var]])
  # rasterize variable
  unit_rast<-rasterize(hi_fence, hi_ohia, var_nm)
  
  # stack rasters
  if(var == 1){
    # create raster stack
    fence_stack<-unit_rast
  }else{
    fence_stack<-stack(fence_stack, unit_rast)    
  }  
}
# rename raster layers
names(fence_stack)<-names(hi_fence)
# save raster tiff data file
# SEE: rod_fenced_units.R SCRIPT
# writeRaster(fence_stack$GlobalID, format="GTiff", overwrite = T, paste0(envDir, "fenced_units.tif"))

# extract all values for fenced areas
hi_fence_spdf<-as.data.frame(as(fence_stack, "SpatialPixelsDataFrame"))
head(hi_fence_spdf)

# extract fence values from Hawaii Island at rod points
rod_fences<-extract(fence_stack$GlobalID, cbind(rod$X_LON, rod$Y_LAT))
# add extracted fence column to master rod data set
rod_master$fenced<-rod_fences
# add label column for fenced and unfenced units
rod_master$FENCED<-rod_master$fenced
rod_master$FENCED[which(rod_master$FENCED > 0)]<-"Fenced"
rod_master$FENCED[which(is.na(rod_master$FENCED))]<-"Not Fenced"

#----- SAVE MASTER FILES -----#

# save master files with extracted data for all islands
write.csv(rod_master, paste0(dataDir, "all_islands_dist_vars.csv"))  # rod master file

# remove points beyond Hawaii Island
rod_master<-rod_master[complete.cases(rod_master$bio1),]
# plot(rod_master$LON, rod_master$LAT)

# save master files with extracted data
write.csv(rod_master, paste0(dataDir, "all_dist_vars.csv"))  # rod master file

######################################
##### END OF VARIABLE EXTRACTION #####
######################################

# load rod master file
# rod_master<-read.csv(paste0(dataDir, "all_dist_vars.csv"), header = T)
# rod_master<-rod_master[,-1]