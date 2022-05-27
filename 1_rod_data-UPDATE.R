### ROD Data Formatting and Collection ###
### raw rod in situ data to be updated ###

##################
##### SET UP #####
##################

# load necessary packages
library(rgdal)
library(raster)
library(dplyr)

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
outDir<-paste0("data/data_to_use/")

# set coordinate system to be used with data as needed
LatLon<-'+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0'
utmSys4<-'+proj=utm +zone=4 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0'
# utmSys5<-'+proj=utm +zone=5 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0'
ArcSys<-'+init=epsg:3857'  # web standard for ArcGIS Online data

# load statewide coast shapefile
arch_coast<-readOGR(paste0(envDir, "hawaii_coast"), "Main_Hawaiian_Islands_simple3")
# convert coastlines to utm coordinates
arch_utm<-spTransform(arch_coast, CRS(utmSys4))
# convert coastlines to arcgis online standard coordinates
arch_arc<-spTransform(arch_coast, CRS(ArcSys))

# load Hawaii Island coast shapefile
hi_coast<-readOGR(paste0(envDir, "hawaii_coast"), "bi_coast")
# convert Hawaii Island coast to utm coordinates
hi_utm<-spTransform(hi_coast, CRS(utmSys4))

# get extent for Hawaii Island only
hi_ext<-extent(hi_coast)

################################
##### RAW IN SITU ROD FILE #####
################################

#----- UPDATE (07/25/2019) -----#
# updated manually with arcgis online data

#----- UPDATE (03/02/2020) ------#
# updated with arcgis ROD geodatabase online
# Brian Tucker <bjtucker@hawaii.edu>
# rod_gdb<-"ROD_Database_20200302.gdb"

#----- UPDATE (04/01/2020) -----#
# updated with latest arcgis ROD geodatabase
# Brian Tucker <bjtucker@hawaii.edu>
# rod_gdb<-"ROD_Database_20200331.gdb"

#----- UPDATE (06/01/2020) -----#
# updated with latest arcgis ROD geodatabase
# Brian Tucker <bjtucker@hawaii.edu>
# rod_gdb<-"ROD_Database_20200604.gdb"
# still does not include Kauai private landowner data

# #----- UPDATE (10/02/2020) -----#
# # updated with latest arcgis ROD geodatabase
# # Brian Tucker <bjtucker@hawaii.edu>
# rod_gdb<-"ROD_Database_20201001.gdb"

# #----- UPDATE (01/07/2021) -----#
# # updated with latest arcgis ROD geodatabase
# # Brian Tucker <bjtucker@hawaii.edu>
# rod_gdb<-"ROD_Database_20210106.gdb"

# #----- UPDATE (05/05/2021) -----#
# # updated with latest arcgis ROD geodatabase
# # Brian Tucker <bjtucker@hawaii.edu>
# rod_gdb<-"ROD_Database_20210503.gdb"

# #----- UPDATE (09/01/2021) -----#
# # updated with latest arcgis ROD geodatabase
# # Brian Tucker <bjtucker@hawaii.edu>
# rod_gdb<-"ROD_Database_20210901.gdb"

# #----- UPDATE (12/01/2021) -----#
# # updated with latest arcgis ROD geodatabase
# # Brian Tucker <bjtucker@hawaii.edu>
# rod_gdb<-"ROD_Database_20211201.gdb"

#----- UPDATE (04/11/2022) -----#
# updated with latest arcgis ROD geodatabase
# Brian Tucker <bjtucker@hawaii.edu>
rod_gdb<-"ROD_Database_20220411.gdb"
# NOTE: returns Field_Form_v7 version of database

# locate arcgis geodatabase files
gdbDir<-paste0(locDir, rod_gdb)
# load data
gdb_data<-readOGR(dsn = gdbDir)
# summary(gdb_data); plot(gdb_data)

# get arcgis dataframe
gdb_df<-gdb_data@data
# head(gdb_df); names(gdb_df)

# get gdb projection
gdb_crs<-gdb_data@proj4string
# get arcgis coordinates
gdb_xy<-gdb_data@coords

# create spatial points object
gdb_sp<-SpatialPointsDataFrame(gdb_xy, data = gdb_df, proj4string = CRS(as.character(gdb_crs)))
# convert arcgis system to lat/lon
gdb_ll<-spTransform(gdb_sp, CRS(LatLon))
# convert arcgis system to utm
gdb_utm<-spTransform(gdb_sp, CRS(utmSys4))

# plot(gdb_sp@coords, pch = 20, main = "ArcGIS"); plot(arch_arc, add = T)
# plot(gdb_ll@coords, pch = 20, main = "Lat/Lon", col = "blue"); plot(arch_coast, add = T)
# plot(gdb_utm@coords, pch = 20, main = "UTM", col = "red"); plot(arch_utm, add = T)

###########################
##### FORMAT ROD DATA #####
###########################

# create final dataframe for rod data
rod<-data.frame(X = gdb_sp@coords[,1], Y = gdb_sp@coords[,2],
                X_LON = gdb_ll@coords[,1], Y_LAT = gdb_ll@coords[,2],
                X_UTM = gdb_utm@coords[,1], Y_UTM = gdb_utm@coords[,2])

# add date column
rod$DATE<-as.Date(gdb_df$sample_date, format = "%Y%m%d")
# add day
rod$DAY<-format(rod$DATE, "%d")
# add month
rod$MONTH<-format(rod$DATE, "%m")
# add year
rod$YEAR<-format(rod$DATE, "%Y")

# check column names
if(names(gdb_df)[1] == "globalid"){
  # select data columns to keep 
  rod_update<-subset(gdb_df, select = c(globalid, sample_id, location_area, elevation_m,
                                        tree_diameter, tree_height, tree_symptom, disease_sign,
                                        staining_yes_no, ceratocystis_species))
}else{ # for Field_Form_v7 version
  # select data columns to keep 
  rod_update<-subset(gdb_df, select = c(sample_id_short, sample_id, location_area, elevation_m,
                                        tree_diameter, tree_height, tree_symptom, disease_sign,
                                        staining_yes_no, ceratocystis_species))
}
# head(rod_update)

# add selected columns to final dataset
rod<-cbind(rod, rod_update)

# add duplicate ceratocystis strains
ab_rod<-rod[which(rod$ceratocystis_species == "A,B"),]
# create separate records for A & B strains
a_rod<-ab_rod; b_rod<-ab_rod
# rename strains
a_rod$ceratocystis_species<-"A"; b_rod$ceratocystis_species<-"B"
# add strains to final dataset
rod<-rbind(rod, a_rod, b_rod)
# remove duplicate A,B strains
rod<-rod[which(rod$ceratocystis_species != "A,B"),]

# format positive (1) or negative (0) detections
rod$PN_NUM<-as.character(rod$ceratocystis_species)
rod$PN_NUM[which(rod$PN_NUM == "A" | rod$PN_NUM == "B")]<-1
rod$PN_NUM[which(rod$PN_NUM == "ND")]<-0

# format strain names as scientific names
rod$SP_NM<-as.character(rod$ceratocystis_species)
rod$SP_NM[which(rod$SP_NM == "A")]<-"C. lukuohia"
rod$SP_NM[which(rod$SP_NM == "B")]<-"C. huliohia"
rod$SP_NM[which(rod$SP_NM == "ND")]<-"Not Detected"
# set factor levels
rod$SP_NM<-factor(rod$SP_NM, levels = c("C. lukuohia", "C. huliohia", "Not Detected"))

# table(rod$ceratocystis_species)
# table(rod$PN_NUM)
# table(rod$SP_NM)

# remove points with bad coordinates
rod<-rod[which(rod$Y_LAT > 15),]

# check final dataset
head(rod); tail(rod)
plot(rod$X_LON, rod$Y_LAT, pch = 20)

#####################
##### SAVE DATA #####
#####################

# order points by date
rod<-arrange(rod, rod$DATE)

# identify duplicate records for at the same location on a single day per strain
rod_dups<-duplicated(cbind(rod$DATE, rod$X_LON, rod$Y_LAT, rod$SP_NM))
table(rod_dups)
# remove duplicate entries
rod<-rod[which(duplicated(cbind(rod$DATE, rod$X_LON, rod$Y_LAT, rod$SP_NM)) == F),]
# rod<-all_rod[which(duplicated(all_rod$Sample_ID) == F),]

# # plot points
# plot(arch_coast)
# points(rod$X_LON, rod$Y_LAT)
# points(rod$X_LON[which(rod$PN_NUM == 1)], rod$Y_LAT[which(rod$PN_NUM == 1)], col = 'yellow')
# points(rod$X_LON[which(rod$SP_NM == "C. huliohia")], rod$Y_LAT[which(rod$SP_NM == "C. huliohia")], col = 'blue')
# points(rod$X_LON[which(rod$SP_NM == "C. lukuohia")], rod$Y_LAT[which(rod$SP_NM == "C. lukuohia")], col = 'red')

# save final output of all formatted rod data in output directory
write.csv(rod, paste0(outDir, "ALL_ROD.csv"))
# write.csv(rod, paste0(outDir, "NEW_ROD.csv"))

# #----- POSITIVE POINTS ONLY -----#

# # select points where rod was detected (PN_NUM = 1)
# pos_rod<-rod[which(rod$PN_NUM == 1), ]
# table(pos_rod$SP_NM, useNA = "ifany")

# ### data duplicate check ###
# # identify duplicate records for at the same location on a single day per strain
# pos_rod$dups<-duplicated(cbind(pos_rod$DATE, pos_rod$X_LON, pos_rod$Y_LAT, pos_rod$SP_NM))
# table(pos_rod$dups)  # 0 duplicates found (previously removed)

# # plot points
# plot(arch_coast, main = "ROD Lab & Field Data")
# points(rod$X_LON, rod$Y_LAT)
# points(pos_rod$X_LON, pos_rod$Y_LAT, col = 'red')

# # save final output of only positive rod data in output directory (Hawaii Island only)
# write.csv(pos_rod, paste0(outDir, "ROD_OCC.csv"))

########################################
### END ROD LOCATION DATA FORMATTING ###
########################################