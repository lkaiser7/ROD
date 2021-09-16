### ROD Distribution Models for Publication ###
### Ensemble modeling using BIOMOD2 package ###

##################
##### SET UP #####
##################

# load necessary packages
library(rgdal)
library(biomod2)
library(raster)
library(ggplot2)
# library(dismo)

# set raster options
rasterOptions(todisk = TRUE)

# set drive
drive<-"E:/"
# set root directory 
rootDir<-paste0(drive, "Dropbox/PIERC/ROD_updates/")  
setwd(rootDir)

# path to formatted data
dataDir<-paste0(rootDir, "data/data_to_use/")
# path to HRCM bioclimatic variables
climDir<-paste0(rootDir, "data/environmental/hrcm_bios/")
# path to processed and saved output files
outDir<-paste0(rootDir, "model_output/")

#--- RASTER FUNCTIONS ---#

# function to save binary rasters
binary_raster = function(bin_rast, rast_nm){
  
  # tiff image file
  tiff(paste0("rod_portal/www/", rast_nm, "_image.tif"), 
       #tiff(paste0(rast_nm, "_image.tif"), 
       res = 300, units = "in", pointsize = 12, width = 10, height = 8, compression = "lzw")
  plot(bin_rast, legend = F)
  legend("bottomleft",'ROD Binary Occurrences', pch = 15, col = 'forestgreen', bty = 'n')
  dev.off()
  
  # tiff data file
  writeRaster(bin_rast, paste0("rod_portal/www/", rast_nm, ".tif"), format="GTiff", overwrite = TRUE)
  #writeRaster(bin_rast, paste0(rast_nm, ".tif"), format="GTiff", overwrite = TRUE)
}

# function to save suitability rasters
suit_raster = function(suit_rast, rast_nm){
  # tiff image file
  tiff(paste0("rod_portal/www/", rast_nm, "_image.tif"), 
       #tiff(paste0(rast_nm, "_image.tif"), 
       res = 300, units = "in", pointsize = 12, width = 10, height = 8, compression = "lzw")
  col_ramp<-colorRampPalette(c('gray96', 'darkgreen'))
  plot(suit_rast, col = col_ramp(n = 99), legend = TRUE, legend.width = 1, legend.shrink = 0.75,
       breaks = seq(minValue(suit_rast), maxValue(suit_rast), length.out = 100),
       legend.args = list(text = "", side = 4, font = 2, line = 2.5, cex = 0.8),
       axis.args = list(at = seq(minValue(suit_rast), maxValue(suit_rast), (maxValue(suit_rast)-minValue(suit_rast))/10),
                        labels = seq(minValue(suit_rast), maxValue(suit_rast), (maxValue(suit_rast)-minValue(suit_rast))/10)))
  plot(hi_state, add = T)
  dev.off()
  
  # tiff data file
  writeRaster(suit_rast, paste0("rod_portal/www/", rast_nm, ".tif"), format="GTiff", overwrite = TRUE)
  #writeRaster(suit_rast, paste0(rast_nm, ".tif"), format="GTiff", overwrite = TRUE)
}

# set palette for color ramp
col_ramp<-colorRampPalette(c('gray96', 'darkolivegreen', 'yellow', 'orange', 'red'))

#####################
##### LOAD DATA #####
#####################

# load Hawaiian Archipelago shapefile
hi_state<-readOGR(paste0(rootDir, "data/environmental/hawaii_coast"), "Main_Hawaiian_Islands_simple3")
hi_spdf<-fortify(hi_state)

# load Hawaii Island coasts shapefile
hi_coast<-readOGR(paste0(rootDir, "data/environmental/hawaii_coast"), "bi_coast")
hi_df<-fortify(hi_coast)
hi_extent<-extent(hi_coast)

# load ohia distribution raster
ohia_rast<-raster(paste0(rootDir, "data/environmental/HIGAP_500m_ohia.tif"))
# extract all ohia values
ohia_spdf<-as.data.frame(as(ohia_rast, "SpatialPixelsDataFrame"))

# create ohia mask
ohia_mask<-ohia_rast
ohia_mask[ohia_mask >= 1]<-1
ohia_mask[ohia_mask < 1]<-NA

# list selected bioclimatic predictor files
bc_vars<-c("bio6.tif", "bio12.tif")

# load predictors
current_bc<-stack(paste0(climDir, "current_500m/", bc_vars))
future_bc<-stack(paste0(climDir, "future_500m/", bc_vars))

# create blank (0) raster layer across the Hawaiian Islands
hi_mask<-reclassify(subset(current_bc, 1, drp = T), c(-Inf, Inf, 0))
# crop to Hawaii Island only
# hi_mask<-crop(hi_mask, extent(hi_coast))
# crop to Hawaiian Islands
hi_mask<-crop(hi_mask, extent(hi_state))
# reduce resolution size of raster by aggregating pixels
hi_mask<-aggregate(hi_mask, fact = 2, fun = mean)
# plot(hi_mask, legend = F)

# list species choices
sp_choices = c("CL", "CH")
# list species names
sp_names = c("C. lukuohia", "C. huliohia")
# list species files names
sp_files = c("C_lukuohia", "C_huliohia")

for(s in 1:length(sp_choices)){  # s = 1
  # select species 
  sp = sp_choices[s]
  sp_nm = sp_names[s]
  sp_f = sp_files[s]
  # table(rod_data$SP_NM, useNA = 'ifany')
  
  # load all rod data
  rod_data<-read.csv(paste0(dataDir, "ALL_ROD.csv"), header = T)
  # format date column (sorts data by date)
  rod_data$DATE<-as.Date(rod_data$DATE)
  
  # subset negative rod samples
  rod_neg<-rod_data[which(rod_data$SP_NM != "C. lukuohia" & rod_data$SP_NM != "C. huliohia"),]
  
  # subset data to focus on Strain A of ROD only
  rod_data<-rod_data[which(rod_data$SP_NM == sp_nm),]
  # dim(rod_data)
  
  # extract cell values for rod xy coordinates
  rod_data$CELLS<-cellFromXY(hi_mask, cbind(rod_data$X_LON, rod_data$Y_LAT))
  rod_neg$CELLS<-cellFromXY(hi_mask, cbind(rod_neg$X_LON, rod_neg$Y_LAT))
  # head(rod_data)
 
  # identify duplicate cells
  rod_data$DUPS<-duplicated(rod_data$CELLS)
  rod_neg$DUPS<-duplicated(rod_neg$CELLS)
  # table(rod_data$DUPS); table(rod_neg$DUPS)
  
  # remove duplicated rod samples per cell
  rod_data<-rod_data[which(rod_data$DUPS == FALSE),]  
  rod_neg<-rod_neg[which(rod_neg$DUPS == FALSE),]  
  # dim(rod_data); dim(rod_neg)
  
  # create data frame for modeling purposes
  rod_pts<-data.frame(x = rod_data$X_LON, y = rod_data$Y_LAT)
  neg_pts<-data.frame(x = rod_neg$X_LON, y = rod_neg$Y_LAT)
  
  #######################
  ##### FORMAT DATA #####
  #######################
  
  # add presence/absence column to rod occurrences
  rod_pts$pa<-rep(1, dim(rod_pts)[1])
  neg_pts$pa<-rep(0, dim(neg_pts)[1])
  
  # store number of presence points
  occ_count<-dim(rod_pts)[1]
  abs_count<-dim(neg_pts)[1]
  
  # create blank (0) raster layer across the Hawaiian Islands
  pred_occ<-reclassify(subset(current_bc, 1, drp = T), c(-Inf, Inf, 0))
  # crop to Hawaii Island only
  # pred_occ<-crop(pred_occ, extent(hi_coast))
  # crop to Hawaiian Islands
  pred_occ<-crop(pred_occ, extent(hi_state))
  # assign presence of rod occurrences to blank raster
  pred_occ[cellFromXY(pred_occ, rod_pts[,1:2])]<-1
  
  # create raster layer of absence points (where rod does not occur)
  pred_na<-pred_occ == 0
  # convert raster to points to get potential absence locations
  pa_pts<-as.data.frame(rasterToPoints(pred_na, fun = function(x){x == 1}))[-3]
  # add column of NA values for absences
  pa_pts$pa<-NA  #creates incomplete cases if looking for duplicates
  # set select number of absence points to 1000
  pa_count<-1000
  
  # merge real rod occurrence data with created candidate pseudo-absence points
  rod_model_data<-rbind(rod_pts, neg_pts, pa_pts)
  # head(rod_model_data); table(rod_model_data$pa, useNA = "ifany"); dim(rod_model_data)[1]
  
  # create matrix with cell numbers of real data from selected bioclim variables
  bc_data<-extract(current_bc, rod_model_data[, 1:2], cellnumbers = T) 
  # combine rod data with extracted bioclim values
  rod_bc_data<-cbind(rod_model_data, bc_data)
  # remove incomplete cases, if any
  rod_bc_data<-rod_bc_data[complete.cases(cbind(rod_bc_data$bio1, rod_bc_data$bio12)),]
  
  # identify duplicates within the cell column
  dup_data<-duplicated(rod_bc_data$cells) 
  # remove duplicate cells and drop cell column  
  rod_bc_data<-rod_bc_data[!dup_data, -4] 
  # head(rod_bc_data); tail(rod_bc_data); table(rod_bc_data$pa, useNA = "ifany")
  # remove points from other islands
  # rod_bc_data<-rod_bc_data[which(rod_bc_data$x >= hi_extent@xmin),]
  
  # change directory path to save model output files and access .jar file
  setwd(outDir)
  # save final data set used for modeling
  write.csv(rod_bc_data, paste0(sp, "_ModelPoints.csv"))
  
  # store modeling start time
  start_time<-Sys.time()
  
  ### BIOMOD Fitting ###
  
  # load data for proper BIOMOD2 formatting
  myBiomodData<-BIOMOD_FormatingData(
    resp.name = sp_f,  #species name
    resp.var = data.frame(rod_bc_data$pa),  #pres/abs/pa points
    expl.var = data.frame(rod_bc_data[,4:dim(rod_bc_data)[2]]),  #bioclim values
    resp.xy = data.frame(rod_bc_data$x, rod_bc_data$y),  #xy coordinates
    PA.nb.rep = 10,  #number of PAs selections
    PA.nb.absences = 1000,  #number of PAs to select
    PA.strategy = "random",  #how to select PAs
    PA.dist.min = 0.1886026)  #minimum distance to presences (calculated 20km)
  
  # set selected modeling techniques 
  # myBiomodOption<-BIOMOD_ModelingOptions(
  #   # MAXENT: Maximum Entropy
  #   MAXENT.Phillips = list(maximumiterations = 100, visible = FALSE, linear = TRUE, 
  #                          quadratic = TRUE, product = TRUE, threshold = TRUE, hinge = TRUE, 
  #                          lq2lqptthreshold = 80, l2lqthreshold = 10, hingethreshold = 15, 
  #                          beta_threshold = -1, beta_categorical = -1, beta_lqp = -1, 
  #                          beta_hinge = -1, defaultprevalence = 0.5))
  myBiomodOption<-BIOMOD_ModelingOptions(
    # GAM: 
    GAM = list(maximumiterations = 100, visible = FALSE, linear = TRUE, 
               quadratic = TRUE, product = TRUE, threshold = TRUE, hinge = TRUE, 
               lq2lqptthreshold = 80, l2lqthreshold = 10, hingethreshold = 15, 
               beta_threshold = -1, beta_categorical = -1, beta_lqp = -1, 
               beta_hinge = -1, defaultprevalence = 0.5))
  
  # create ensemble model from formatting data and modeling options
  myBiomodModelOut<-BIOMOD_Modeling(
    data = myBiomodData,  #formatted biomod data
    #models = 'MAXENT.Phillips',  #select model types to run
    models = 'GAM',  #select model types to run
    models.options = myBiomodOption,  #biomod options object
    NbRunEval = 10,  #number of evaluation runs
    DataSplit = 80,  #amount of data to use for training
    Yweights = NULL,  #response points weights
    Prevalence = NULL, #used to build 'weighted response weights' 
    VarImport = 4,  #permuations to estimate variable importance*** 10
    do.full.models = TRUE,  #calibrate and evaluate to all data
    # models.eval.meth = c("ROC", "TSS"),  #evaluation metrics
    models.eval.meth = c("ROC"),  #evaluation metrics
    SaveObj = TRUE)  #save output object
  
  # MODELING SUMMARY: return summary of BIOMOD2 model fitting
  myBiomodModelOut
  # save BIOMOD model fitting output from workspace
  save("myBiomodModelOut", file = paste0(sp, "_ModelFit.RData"))
  
  # return output model evaluation metrics results
  write.csv(t(data.frame(get_evaluations(myBiomodModelOut))), paste0(sp, "_model_evaluations.csv"))
  # get variable importance of selected bioclim variables
  write.csv(t(data.frame(get_variables_importance(myBiomodModelOut))), paste0(sp, "_variable_importance.csv"))
  
  ### BIOMOD Ensembles ###
  
  # list models computed from model fitting to use to build ensemble 
  remaining_models<-myBiomodModelOut@models.computed
  
  # combine models and make ensemble predictions from model fitting 
  myBiomodEM<-BIOMOD_EnsembleModeling( 
    modeling.output = myBiomodModelOut,  #BIOMOD.models.out from model fitting
    chosen.models = remaining_models,  #vector of model runs to use
    em.by = 'all',  #how models will be combined
    # eval.metric = c("ROC", "TSS"), #evaluation metrics to build ensemble
    eval.metric = c("ROC"), #evaluation metrics to build ensemble
    # ERROR: match number of metrics eval.metric.quality.threshold = rep(0.5, 2),  #threshold to exclude models per eval stats
    eval.metric.quality.threshold = rep(0.5),  #threshold to exclude models per eval stats
    prob.mean = TRUE,  #estimate mean probabilities 
    prob.cv = TRUE,  #estimate coefficient of variation
    prob.ci = TRUE,  #estimate confidence interval of prob.mean
    prob.ci.alpha = 0.05,  #signficance level for estimating confidence interval
    prob.median = TRUE,  #estimate median
    committee.averaging = TRUE,  #estimate committee averaging
    prob.mean.weight = TRUE,  #estimate weighted sums
    prob.mean.weight.decay = 'proportional') #define relative importance of weights
  
  # MODELING SUMMARY: return summary of BIOMOD2 ensemble models
  myBiomodEM
  # save BIOMOD ensemble model output from workspace
  save("myBiomodEM", "myBiomodModelOut", "remaining_models", file = paste0(sp, "_EnsembleFit.RData"))
  
  # get evaluation scores and statistics
  write.csv(t(data.frame(get_evaluations(myBiomodEM))), paste0(sp, "_ensemble_evaluations.csv"))
  
  ### Model Projections ###
  
  # run model projections archipelago-wide
  myBiomodProj<-BIOMOD_Projection(
    modeling.output = myBiomodModelOut,  #BIOMOD.models.out from model fitting
    new.env = current_bc,  #new explanatory variables to project model
    proj.name = "current",  #projection directory
    selected.models = remaining_models,  #which models to use for projections
    # binary.meth = c("ROC", "TSS"),  #evaluation method statistics 
    binary.meth = c("ROC"),  #evaluation method statistics 
    compress = 'xz',  #compression format of files on hard drive
    build.clamping.mask = FALSE,  #if clamping mask should be saved or not
    keep.in.memory = TRUE) #if clamping mask should be saved to hard disk or not
  
  # save projection R workspace environement 
  save("myBiomodProj", file = paste0(sp, "_ModelProjection.RData"))
  
  # run ensemble projections for species
  myBiomodEF <- BIOMOD_EnsembleForecasting(
    projection.output = myBiomodProj,  #BIOMOD.projection.out from projections
    total.consensus = TRUE,  #mean of all combined model projections
    EM.output = myBiomodEM,  #BIOMOD.EnsembleModeling.out from ensemble modeling
    # ERROR: binary.meth = c("ROC", "TSS"),  #evaluation method statistics 
    binary.meth = c("ROC"),  #evaluation method statistics 
    keep.in.memory = TRUE)  #if output should be saved to hard disk or not
  
  # save ensemble projections R workspace environement 
  save("myBiomodProj", "myBiomodEF", file = paste0(sp, "_EmsembleForecast.RData"))
  
  # store modeling end time
  end_time<-Sys.time()
  # calculate processing time
  end_time - start_time
  
  ################################
  ##### STORE RASTER OUTPUTS #####
  ################################
  
  # set folder directory
  if(s == 1){
    spDir<-"C.lukuohia"
  }else{
    spDir<-"C.huliohia"  
  }
  
  # load rod ensemble grid raster as stack for different bands
  cf_ensemble<-stack(paste0(spDir, "/proj_current/proj_current_", spDir, "_ensemble.grd"))
  # names(cf_ensemble)
  # plot(cf_ensemble)
  # cf_ensemble@layers
  
  # raster layer band with weighted ROC and TSS means (wmean) from ensemble models
  roc_band<-which(names(cf_ensemble) == paste0(spDir, "_EMwmeanByROC_mergedAlgo_mergedRun_mergedData"))
  # store roc suitability raster layer
  roc_suit<-raster(cf_ensemble, layer = roc_band)/1000
  # load evaluation ensemble rasters for roc and tss (28 layers each)
  roc_rast<-stack(paste0(spDir, "/proj_current/proj_current_", spDir, "_ensemble_ROCbin.grd"))
  
  # store roc binary raster layer
  roc_bin<-raster(roc_rast, layer = roc_band)
  
  ###################################
  ##### RASTER OUTPUT FUNCTIONS #####
  ###################################
  
  # change directory path to save raster output files
  setwd(rootDir)
  
  ### Save Rasters ###
  
  # save binary rasters
  binary_raster(roc_bin, paste0(sp, "_Binary_current_ROCwmean"))
  # save suitability rasters
  suit_raster(roc_suit, paste0(sp, "_Suitability_current_ROCwmean"))
  # save clipped suitability rasters
  suit_raster(roc_bin*roc_suit, paste0(sp, "_ClippedSuitability_current_ROCwmean"))
  
  ### MODEL FIGURES ###

  # plot rod suitability (Figure 4)
  # tiff(paste0("rod_portal/www/Fig4_", sp, "_suitability_rod_gam_current500m_image.tif"), 
  #      res = 300, units = "in", pointsize = 12, width = 10, height = 8, compression = "lzw")
  png(paste0("rod_portal/www/Fig4_", sp, "_suitability_rod_gam_current500m_image.png"), 
       res = 300, units = "in", pointsize = 12, width = 10, height = 8)
  plot(roc_suit, col = col_ramp(n = 99)) 
  plot(hi_state, add = T)
  mtext("High", side = 4, at = 21.6, las = 1, adj = -0.65)
  mtext("Low", side = 4, at = 19.55, las = 1, adj = -0.75)
  dev.off()

  # crop rod model to ohia distribution only
  crop_rod<-roc_suit + ohia_mask
  crop_rod[crop_rod < 1]<-NA
  
  # plot rod model within ohia (Figure S2)
  # tiff(paste0("rod_portal/www/FigS2_", sp, "_ohia_rod_gam_current500m_image.tif"), 
  #      res = 300, units = "in", pointsize = 12, width = 10, height = 8, compression = "lzw")
  png(paste0("rod_portal/www/FigS2_", sp, "_ohia_rod_gam_current500m_image.png"), 
       res = 300, units = "in", pointsize = 12, width = 10, height = 8)
  plot(crop_rod-1, col = col_ramp(n = 99), main = "Suitability within Ohia Distribution")
  #legend("bottomleft", legend = "Ohia Distribution", pch = 15, col = 'gray86', bty = 'n', pt.cex = 2)
  #legend("bottomleft", legend = "", pch = 0, bty = 'n', pt.cex = 2)
  plot(hi_state, add = T)
  mtext("High", side = 4, at = 21.6, las = 1, adj = -0.65)
  mtext("Low", side = 4, at = 19.55, las = 1, adj = -0.75)
  dev.off()
  
}

# RUN IF NEEDED ONLY #
# not yet edited to new data or individual cerartocystis species

#-------------------#
### PLANT OVERLAP ###
#-------------------#

# # load ohia distribution raster
# ohia_rast<-raster(paste0(rootDir, "data/environmental/HIGAP_500m_ohia.tif"))
# # crop to Hawaii Island only
# hi_ohia<-crop(ohia_rast, hi_coast)
# 
# # group tiff files
# plant_tiff<-list.files("data/environmental/plant_present_dist/", pattern = ".tif$")
# # plants 186:223 do not have data and result in error as cannot create a RasterLayer object
# plant_tiff<-plant_tiff[-186:-223]
# length(plant_tiff)
# 
# # load previously assembled formatting summary file (SEE: rod_plant_overlap.R)
# plant_df<-read.csv(paste0(rootDir, "data/environmental/HI_Plant_Present_Dist_List.csv"), header = T)
# head(plant_df)
# 
# # list Hawaii Island plant files
# bi_plants<-list.files(paste0(rootDir, "data/environmental/hawaii_plants/"))
# 
# # add columns for plant cell sums
# plant_df$HI_ARCH_CELLS<-NA
# plant_df$BI_ONLY_CELLS<-NA
# 
# # load Strain A ROD BIOMOD2 suitability distribution results
# rod_model<-raster(paste0("rod_portal/www/ClippedSuitability_current_ROCwmean.tif"))
# # crop to Hawaii Island
# bi_rod_model<-crop(rod_model, hi_extent)
# # plot(bi_rod_model, col = rev(heat.colors(30))); plot(hi_coast, add = T)
# 
# # load Strain A ROD BIOMOD2 binary distribution results
# rod_bin<-raster(paste0("rod_portal/www/Binary_current_ROCwmean.tif"))
# # crop to Hawaii Island
# bi_rod_bin<-crop(rod_bin, hi_extent)
# 
# # create data frame for cell sums
# cell_sum<-data.frame(species = "ohia/rod", 
#                      ohia_dist = cellStats(ohia_rast, sum), rod_model = cellStats(rod_model, sum),
#                      rod_bin = cellStats(rod_bin, sum), bi_ohia = cellStats(hi_ohia, sum), 
#                      bi_rod = cellStats(bi_rod_model, sum), bi_bin = cellStats(bi_rod_bin, sum))
# # create data frame for percent overlap
# pct_overlap<-data.frame(species = "ohia/rod", cell_sum[1,-1]/cell_sum[1,-1]*100)
# 
# # calculate suitability scores for ROD
# summary(as.data.frame(as(rod_model, "SpatialPixelsDataFrame"))[1])
# rod_suit<-summary(bi_rod_model@data@values) # use to create data frame only!
# bi_rod_suit<-summary(bi_rod_model@data@values)
# 
# # set path to save plant overlap outputs
# plantPath<-paste0(rootDir, "data/environmental/")
# 
# # start
# start_time<-Sys.time()
# 
# for(m in 1:length(plant_tiff)){  # set m = 3 for debugging
#   # store species name
#   plant<-unlist(strsplit(plant_tiff[m], split = " "))[-1]
#   plant_file<-paste0(plant[1], "_", plant[2])
#   plant_nm<-paste0(plant[1], "_", file_path_sans_ext(plant[2]))
#   
#   # load archipelago wide plant raster
#   arch_plant<-raster(paste0(rootDir, "data/environmental/hi_arch_plants/", plant_file))
#   
#   # resample plant extent
#   arch_plant<-resample(arch_plant, rod_model, "bilinear")
#   # calculate sum of all plant cells
#   plant_sum<-cellStats(arch_plant, stat = sum)
#   # add to plant data frame
#   plant_df$HI_ARCH_CELLS[m]<-plant_sum
#   
#   # mask plant across the Hawaiian Archipelago
#   ohia_mask<-mask(ohia_rast, arch_plant)
#   rod_mask<-mask(rod_model, arch_plant)
#   bin_mask<-mask(rod_bin, arch_plant)
#   
#   # save raster masks
#   writeRaster(ohia_mask, "GTiff", overwrite = T,
#               file = paste0(plantPath, "hi_arch_plants/", plant_nm, "_ohia_mask.tif"))
#   writeRaster(rod_mask, "GTiff", overwrite = T,
#               file = paste0(plantPath, "hi_arch_plants/", plant_nm, "_rod_mask.tif"))
#   writeRaster(bin_mask, "GTiff", overwrite = T,
#               file = paste0(plantPath, "hi_arch_plants/", plant_nm, "_bin_mask.tif"))
# 
#   # save tiff image file
#   tiff(paste0(plantPath, "hi_arch_plants/", plant_nm, "_masks_images.tif"),
#        res = 300, units = "in", pointsize = 12, width = 10, height = 8, compression = "lzw")
#   # plot image for visual output
#   par(mfrow = c(2,2))
#   plot(arch_plant, main = paste0(plant_nm, " Distribution"), legend = F)
#   plot(hi_state, add = T)
#   plot(ohia_mask, main = "Within Ohia Distribution", legend = F)
#   plot(hi_state, add = T)
#   plot(rod_mask, main = "Projected ROD Suitability")
#   plot(hi_state, add = T)
#   plot(bin_mask, main = "Within ROD Binary Distribution", legend = F)
#   plot(hi_state, add = T)
#   dev.off()
#   
#   # calculate the overlap sum of cells 
#   ohia_sum<-cellStats(ohia_mask, stat = sum)
#   rod_sum<-cellStats(rod_mask, stat = sum)
#   bin_sum<-cellStats(bin_mask, stat = sum)
#   # create data frame of sums
#   arch_sum<-data.frame(species = plant_nm, ohia_dist = ohia_sum, rod_model = rod_sum, rod_bin = bin_sum)
#   
#   # calculate suitability score summary values
#   plant_suit<-summary(rod_mask@data@values)
#   # add to data frame for final output
#   rod_suit<-rbind(rod_suit, plant_suit)
#   
#   # check if plant is present on Hawaii Island
#   if(plant_df$MIN[m] == "1"){
#     # load Hawaii Island plant raster
#     hi_plant<-raster(paste0(rootDir, "data/environmental/hawaii_plants/", plant_file))
#     # resample plant extents
#     bi_plant<-resample(hi_plant, bi_rod_model, "bilinear")
#     
#     # calculate sum of plant cells
#     bi_plant_sum<-cellStats(bi_plant, stat = sum)
#     # add to plant data frame
#     plant_df$BI_ONLY_CELLS[m]<-bi_plant_sum
#     
#     # mask plant across the Hawaii Island
#     bi_ohia_mask<-mask(hi_ohia, bi_plant)
#     bi_rod_mask<-mask(bi_rod_model, bi_plant)
#     bi_bin_mask<-mask(bi_rod_bin, bi_plant)
#     
#     # save raster masks
#     writeRaster(bi_ohia_mask, "GTiff", overwrite = T,
#                 file = paste0(plantPath, "hawaii_plants/", plant_nm, "_ohia_mask.tif"))
#     writeRaster(bi_rod_mask, "GTiff", overwrite = T,
#                 file = paste0(plantPath, "hawaii_plants/", plant_nm, "_rod_mask.tif"))
#     writeRaster(bi_bin_mask, "GTiff", overwrite = T,
#                 file = paste0(plantPath, "hawaii_plants/", plant_nm, "_bin_mask.tif"))
# 
#     # save tiff image file
#     tiff(paste0(plantPath, "hawaii_plants/", plant_nm, "_masks_images.tif"),
#          res = 300, units = "in", pointsize = 12, width = 10, height = 8, compression = "lzw")
#     # plot image for visual output
#     par(mfrow = c(2,2))
#     plot(bi_plant, main = paste0(plant_nm, " Distribution"), legend = F)
#     plot(hi_coast, add = T)
#     plot(bi_ohia_mask, main = "Within Ohia Distribution", legend = F)
#     plot(hi_coast, add = T)
#     plot(bi_rod_mask, main = "Projected ROD Suitability")
#     plot(hi_coast, add = T)
#     plot(bi_bin_mask, main = "Within ROD Binary Distribution", legend = F)
#     plot(hi_coast, add = T)
#     dev.off()
#     
#     # calculate the overlap sum of cells 
#     bi_sum<-data.frame(bi_ohia = cellStats(bi_ohia_mask, stat = sum),
#                        bi_rod = cellStats(bi_rod_mask, stat = sum), 
#                        bi_bin = cellStats(bi_bin_mask, stat = sum))
#     
#     # calculate suitability score summary values
#     bi_suit<-summary(bi_rod_mask@data@values)
#     # add to data frame for final output
#     bi_rod_suit<-rbind(bi_rod_suit, bi_suit)
#     
#   }else{
#     bi_sum<-data.frame(bi_ohia = NA, bi_rod = NA, bi_bin = NA)
#   }
#   
#   # combine sum data frames
#   sum_df<-cbind(arch_sum, bi_sum)
#   # add data frame to final output
#   cell_sum<-rbind(cell_sum, sum_df)
#   
#   # calculate overlap of cells
#   plant_overlap<-data.frame(species = plant_nm, sum_df[,-1]/cell_sum[1,-1]*100)
#   # add data frame to final output
#   pct_overlap<-rbind(pct_overlap, plant_overlap)
# }
# 
# # end
# end_time<-Sys.time()
# end_time - start_time
# 
# # save final outputs 
# # write.csv(data.frame(Species = cell_sum$species, rod_suit),
# #           #write.csv(data.frame(Species = c("bi_rod", plant_df$SPECIES), rod_suit),
# #           file = paste0(plantPath, "HI_Arch_ROD_Suitability_Scores.csv"))
# # write.csv(data.frame(Species = cell_sum$species[which(cell_sum$bi_rod != "NA")], bi_rod_suit),
# #           #write.csv(data.frame(Species = c("bi_rod", file_path_sans_ext(bi_plants)), bi_rod_suit),
# #           file = paste0(plantPath, "HawaiiIsland_ROD_Suitability_Scores.csv"))
# write.csv(cell_sum, file = paste0(plantPath, "Plant_Overlap_Cell_Values.csv"))
# write.csv(pct_overlap, file = paste0(plantPath, "Plant_Overlap_Percent.csv"))
# write.csv(plant_df, file = paste0(plantPath, "Plant_Present_Dist_Cell_Values.csv"))
# 
# #----------------------#
# ### OUTPUT SUMMARIES ###
# #----------------------#
# 
# # load final datasets
# # rod_suit<-read.csv("THIN_summary_analysis/plant_overlap/HI_Arch_ROD_Suitability_Scores.csv", header = T)
# # bi_rod_suit<-read.csv("THIN_summary_analysis/plant_overlap/HawaiiIsland_ROD_Suitability_Scores.csv", header = T)
# cell_sum<-read.csv("data/environmental/Plant_Overlap_Cell_Values.csv", header = T)
# pct_overlap<-read.csv("data/environmental/Plant_Overlap_Percent.csv", header = T)
# plant_df<-read.csv("data/environmental/Plant_Present_Dist_Cell_Values.csv", header = T)
# 
# # recalculate overlap of plant distribution per rod binary model over plant area
# arch_overlap<-cell_sum$rod_bin[-1]/plant_df$HI_ARCH_CELLS*100
# bi_overlap<-cell_sum$bi_bin[-1]/plant_df$BI_ONLY_CELLS*100
# plant_overlap<-data.frame(Scientific.Name = str_replace(cell_sum$species[-1], "_", " "), 
#                           df_species = plant_df$SPECIES,
#                           arch_overlap, bi_overlap)
# write.csv(plant_overlap, file = "data/environmental/Binary_Plant_Overlap_Percent.csv")
# head(plant_overlap)
# 
# # combine and merge final data set
# names(cell_sum)[2]<-"SPECIES"
# final_overlap<-merge(plant_df[,-1:-2], cell_sum[-1,-1], by = "SPECIES")
# names(plant_overlap)[2]<-"SPECIES"
# final_overlap<-merge(final_overlap, plant_overlap[,-1])
# head(final_overlap)
# # save final output file
# write.csv(final_overlap, file = "data/environmental/FINAL_Plant_Overlap.csv")
# 
# # show archipelago-wide overlap with species (@ 100% overlap)
# dim(plant_overlap[which(plant_overlap$arch_overlap == 100),])[1]
# 
# # select species on Hawaii Island only from percentage of overlap
# bi_pct_overlap<-plant_overlap[which(plant_overlap$bi_overlap != "NA"),]
# head(bi_pct_overlap)
# # Big Island overlap with species (@ 90% overlap)
# dim(bi_pct_overlap[which(bi_pct_overlap$bi_overlap > 90),])[1]
# 
# #-----------------#
# ### T&E SPECIES ###
# #-----------------#
# # final_overlap<-read.csv("THIN_summary_analysis/plant_overlap/FINAL_Plant_Overlap.csv", header = T)
# 
# # load threatened and endangered FWS plant species list
# # ECOS/Listed Plants https://ecos.fws.gov/ecp/ 
# fws_te<-read.csv("data/environmental/FWS_TE_ListedPlants.csv", header = T)
# head(fws_te); dim(fws_te)
# 
# # merge with plant overlap values to see which are most T&E
# fws_overlap<-merge(plant_overlap, fws_te, "Scientific.Name")
# head(fws_overlap); dim(fws_overlap)
# table(fws_overlap$Federal.Listing.Status)
# # save table output of TE overlap species
# write.csv(fws_overlap, "data/environmental/FWS_TE_OverlapPlants.csv")
# # fws_overlap<-read.csv("data/environmental/FWS_TE_OverlapPlants.csv", header = T)
# 
# # show 100% overlap with T&E species archipelago-wide
# dim(fws_overlap[which(fws_overlap$arch_overlap == 100),])[1]
# # show 100% overlap with T&E species on Hawaii Island only
# fws_overlap[which(fws_overlap$bi_overlap > 98),]
# 
# ### FIGURES ###
# 
# # show overlap with T&E species on Hawaii Island only 
# dim(fws_overlap[which(fws_overlap$bi_overlap > 0),])
# ggplot(fws_overlap[which(fws_overlap$bi_overlap > 0),]) +
#   geom_bar(aes(x = reorder(Scientific.Name, -bi_overlap), y = bi_overlap), stat = "identity") + 
#   labs(x = "Species", y = "Percent Range Overlap (%)") +
#   geom_text(aes(label = paste0(round(bi_overlap, 1), "%"), 
#                 x = Scientific.Name, y = bi_overlap + 3), size = 3, color = "black") +
#   ggtitle(expression(paste("FWS T&E Species Overlap with ", italic("C. lukuohia"), " on Hawai'i Island"))) +
#   coord_flip() + theme_classic() + theme(plot.title = element_text(hjust = 0.5))
# # save ggplot output
# ggsave(paste0("rod_portal/www/Fig5_HI_TE_Overlap_ROD.png"), width = 9, height = 8)

