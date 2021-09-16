### Thinned Prevalence Graphs ###

##################
##### SET UP #####
##################

# load necessary packages
library(raster)
library(plyr)
library(zoo)
library(ggplot2)
# library(rgdal)
# library(maptools)
# library(tools)
# library(stringr)
# library(data.table)
# library(geosphere)

# set drive
drive<-"E:/"
# set root directory 
rootDir<-paste0(drive, "Dropbox/PIERC/ROD_updates/")  
setwd(rootDir)

######################
##### COPY FILES #####
######################

# path to extracted data
extractDir<-paste0(rootDir, "data/data_to_use/")
# path to formatted data
dataDir<-paste0(rootDir, "rod_portal/")

# list of files to copy
data_list<-list.files(extractDir, pattern = ".csv")
# tif_list<-list.files(extractDir, pattern = ".tif")

# copy new extracted data files to portal for summary analysis
file.copy(paste0(extractDir, data_list), dataDir, overwrite = T)
# file.copy(paste0(extractDir, tif_list), dataDir, overwrite = T)

# path to processed and saved output files
# outDir<-paste0(rootDir, "THIN_summary_analysis/")

#####################
##### LOAD DATA #####
#####################

# load rod data 
ROD_data<-read.csv("rod_portal/rod_dist_vars.csv")
head(ROD_data)
# format date column
ROD_data$DATE<-as.Date(ROD_data$DATE, format = "%Y-%m-%d")

# remove incomplete data to crop to hawaii island only
ROD_data<-ROD_data[which(is.na(ROD_data$bio1) == F), ]
# any(sapply(ROD_data, is.infinite))
ROD_data<-ROD_data[which(is.na(ROD_data$SP_NM) == F), ]

# set parameters for graph outputs
min_N_for_comparison = 10  # drop classes with less than 10 samples
n_bins = 8  # bins for bioclim prevalence gradient graphs 
n_bins_combo_all = 10  # bins for double bioclim gradient vars

# thinning data options
thinned = T  # use thinned data set of all unique samples
# symptomatic type options
only_symptomatic = F  # only look at positive samples that have ROD-like symptoms
# create directory for thinned analysis outputs
dir.create("rod_portal/www/", showWarnings = F)

# set string to use in output legends if thinned
if(thinned){
  thin_legend_str = "unique "
}else{
  thin_legend_str = ""
}

# load ohia data
ohia<-read.csv(paste0(dataDir, "ohia_dist_vars.csv"), header = T)

# select variables and do statistical tests
vars_to_plot<-c("bio1", "bio12", "bio4", "bio7", "bio15", "bio6", 
                "bio9", "bio10", "bio11", "bio14", "bio16")

# list bioclimatic variable names
bc_names<-c("Mean annual temperature", "Mean annual precipitation", "Temperature seasonality", 
            "Temperature annual range", "Precipitation seasonality", "Minimum T at coldest month", 
            "Mean T of driest quarter", "Mean T of warmest quarter", "Mean T of coldest quarter",
            "Precipitation of driest month", "Precipitation of wettest quarter")

# set names for bioclimatic columns 
variable_cols = c("bio1", "bio12", "bio4", "bio7", "bio15", "bio6", "bio9", "bio10", "bio11",  
                  "bio14", "bio16", "ROADS", "HIGAP",  "BIOREG", "AGE_GROUP", "NAME")
# set names for bioclimatic variables
variable_cols_names = c("Mean annual temperature", "Mean annual precipitation", "Temperature seasonality", 
                        "Temperature annual range", "Precipitation seasonality", "Minimum T at coldest month", 
                        "Mean T of driest quarter", "Mean T of warmest quarter", "Mean T of coldest quarter",
                        "Precipitation of driest month", "Precipitation of wettest quarter", "Distance to roads", 
                        "Vegetation class", "Biogeographical region", "Substrate age", "Substrate type")

########################
##### RUN ANALYSIS #####
########################

#--- RUNNING MEANS FOR BIOCLIMS ---#

# format data per species
ROD_data_CL<-ROD_data[which(ROD_data$SP_NM != "C. huliohia"),]
ROD_data_CH<-ROD_data[which(ROD_data$SP_NM != "C. lukuohia"),]
# arrange data in chronological order
ROD_data_CL<-arrange(ROD_data_CL, ROD_data_CL$DATE)
ROD_data_CH<-arrange(ROD_data_CH, ROD_data_CH$DATE)

#--- THINNING ANALYSIS ---#

# function to summarize data by categories and tabulate across land cover classes
# SEE: http://www.cookbook-r.com/Manipulating_data/Summarizing_data/
summarySE<-function(data = NULL, measurevar, groupvars = NULL, na.rm = FALSE,
                    conf.interval = 0.95, .drop = TRUE){
  # length which can handle NA's: if na.rm == T, don't count them
  length2<-function (x, na.rm = FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  # Summary for each group's data frame, return a vector with N, mean, and sd
  datac<-ddply(data, groupvars, .drop = .drop,
               .fun = function(xx, col){
                 c(N = length2(xx[[col]], na.rm = na.rm),
                   mean = mean   (xx[[col]], na.rm = na.rm),
                   sd = sd     (xx[[col]], na.rm = na.rm),
                   min = min   (xx[[col]], na.rm = na.rm),
                   max = max   (xx[[col]], na.rm = na.rm),
                   sum = sum   (xx[[col]], na.rm = na.rm)
                 )
               },
               measurevar
  )
  # Rename the "mean" column    
  datac<-rename(datac, c("mean" = measurevar))
  # calculate standard error of the mean
  datac$se<-datac$sd / sqrt(datac$N)
  # calculate data range
  datac$range <- datac$max-datac$min
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is 0.95, use 0.975 (above/below), and use df=N-1
  ciMult<-qt(conf.interval/2 + 0.5, datac$N - 1)
  datac$ci<-datac$se * ciMult
  
  # function output
  return(datac)
}

# set strains of ceratocystis levels
ROD_data$PosVsNeg<-NA

# list species choices
sp_choices = c("CL", "CH")
# select species 
sp = "CL"  # select "CL" or "CH" 

for(s in 1:length(sp_choices)){  # s = 1
  # select species 
  sp = sp_choices[s]  # select "CL" or "CH" 
  
  # name species of strain to be used
  proj_str = paste0(sp, "_all")
  
  # for species
  if (sp == "CL"){
    ROD_data[ROD_data$SP_NM == "C. lukuohia","PosVsNeg"]<-"Ceratocystis positive"
    ROD_data[ROD_data$SP_NM %in% c("C. huliohia", "Not Detected"), "PosVsNeg"]<-"Negative"
    # check PosVsNeg column
    # table(ROD_data$PosVsNeg, useNA = "ifany")
    ROD_data_filtered<-ROD_data_CL
    # table(ROD_data_filtered$SP_NM, useNA = 'ifany')
    
    # subset A strain of CF only
    a_rod<-ROD_data_filtered[which(ROD_data_filtered$SP_NM == "C. lukuohia"), ]
  }else{
    ROD_data[ROD_data$SP_NM == "C. huliohia", "PosVsNeg"]<-"Ceratocystis positive"
    ROD_data[ROD_data$SP_NM %in% c("C. lukuohia", "Not Detected"), "PosVsNeg"]<-"Negative"
    # check PosVsNeg column
    # table(ROD_data$PosVsNeg, useNA = "ifany")
    ROD_data_filtered<-ROD_data_CH
    # table(ROD_data_filtered$SP_NM, useNA = 'ifany')
    
    # subset A strain of CF only
    a_rod<-ROD_data_filtered[which(ROD_data_filtered$SP_NM == "C. huliohia"), ]
  }

  # subset samples not detected
  nd_rod<-ROD_data_filtered[which(ROD_data_filtered$SP_NM == "Not Detected"), ]
  
  # loop to calculate running summary of statistics for bioclimatic variables
  for(var_to_plot in vars_to_plot){  # set var_to_plot=vars_to_plot[1] for debugging
    # set bio name
    bio_var = var_to_plot
    bio_nm = data.frame(vars_to_plot, bc_names)$bc_names[which(vars_to_plot == var_to_plot)]
    
    # select bio column
    rod_bc = a_rod[,var_to_plot]
    # select ohia distibution
    ohia_bc<-ohia[,var_to_plot]
    
    # calculate ohia min, mean, and max
    ohia_min<-min(ohia_bc, na.rm = T)
    ohia_mean<-mean(ohia_bc, na.rm = T)
    ohia_max<-max(ohia_bc, na.rm = T)
    
    # plot running means
    mean_plot<-ggplot() + 
      # OHIA
      geom_hline(aes(yintercept = ohia_min, linetype = 'minimum'), color = 'darkblue', size = 1) +
      geom_hline(aes(yintercept = ohia_mean, linetype = 'mean'), color = 'ivory4', size = 1) +
      geom_hline(aes(yintercept = ohia_max, linetype = 'maximum'), color = 'darkred', size = 1) +
      scale_linetype_manual(name = 'Ohia Range', values = rep(1, 3), 
                            guide = guide_legend(override.aes = list(color = c('darkred', 'ivory4', 'darkblue')), order = 1)) + 
      
      # NEGATIVES
      geom_point(data = nd_rod, aes(DATE, nd_rod[,var_to_plot], shape = SP_NM), color = 'gray') +
      scale_shape_manual(name = 'Ceratocystis', values = c(1), labels = c('Not Detected')) +
      # POSITIVES
      geom_point(data = a_rod, aes(DATE,  rod_bc, color = 'black')) +
      geom_line(data = a_rod, aes(DATE,  cumsum(rod_bc)/seq(along = as.matrix(rod_bc)), color = 'blue'), size = 1, linetype = 2) +
      geom_line(data = a_rod, aes(DATE,  cummax(rod_bc), color = 'orange'), size = 1) +
      geom_line(data = a_rod, aes(DATE,  cummin(rod_bc), color = 'red'), size = 1) + 
      scale_color_manual(name = 'C. lukuohia', values = c('black', 'orange', 'red', 'blue'), 
                         labels = c('variable value', 'running mean', 'cumulative max', 'cumulative min'),
                         guide = guide_legend(override.aes = list(shape = c(19, NA, NA, NA), linetype = c('blank', 'dotted', 'solid', 'solid')))) +
      # plot settings
      labs(title = bio_nm, x = NULL, y = "Bioclimatic Variable") + 
      theme_bw() + theme(plot.title = element_text(hjust = 0.5)) +
      guides(color = guide_legend(order = 2), shape = guide_legend(order = 3))
    
    # save plot
    # ggsave(filename = paste0(dataDir, "www/", proj_str, "_", var_to_plot, "_means.tiff"), mean_plot, 
             # compress = "lzw", height = 6, width = 8, dpi = 300)
    ggsave(filename = paste0(dataDir, "www/", proj_str, "_", var_to_plot, "_means.png"), mean_plot, 
           height = 6, width = 8, dpi = 300)
  }
  
  # create strings for graph outputs
  graph_string = paste0("Negative proportion of ", thin_legend_str, "samples")
  graph_string2 = paste0("Prevalence in ", thin_legend_str, "samples")
  # remove remaining NA values
  ROD_data<-ROD_data[!is.na(ROD_data$PosVsNeg),]
  
  # thin samples
  if(thinned){
    # set thinning string
    thin_str = "rod_portal/www/"
    # load raster mask
    raster_500m<-raster("data/environmental/hrcm_bios/current_500m/bio1.tif")
    # extract cell values for rod xy coordinates
    ROD_data$CELLS<-cellFromXY(raster_500m, cbind(ROD_data$X_LON, ROD_data$Y_LAT))
    # select unique independent variables
    ind_vars<-c("age_group", "higap", "CELLS", "NAME", "bioreg")
    # apply unique variable to each rod sample
    ROD_data$env_ID<-apply(ROD_data[,ind_vars], 1, paste, collapse = ".")
    # select unique variables only
    unique_envs<-unique(ROD_data$env_ID)
    length(unique_envs)
    # create thinned data column to keep
    ROD_data$thin<- T
    # create temporary ID column
    ROD_data$tmpRowId<-c(1:dim(ROD_data)[1])
    # select unique environmental variables to keep
    unique_env = unique_envs[4]
    # set counter
    i = 1
    for(unique_env in unique_envs){
      tmp_df<-ROD_data[ROD_data$env_ID %in% unique_env, ]
      # select positive ceratocystis samples first
      if(any(tmp_df$PosVsNeg %in% "Ceratocystis positive")){
        row_to_pick_tmpDF<-which(tmp_df$PosVsNeg == "Ceratocystis positive")[1]
      }else{
        # select next sample available
        row_to_pick_tmpDF<-1
      }
      # select row to remove
      row_to_pick<-tmp_df[row_to_pick_tmpDF, "tmpRowId"]
      ROD_data[ROD_data$tmpRowId == row_to_pick, "thin"]<-F
      cat("doing ", i, "found row ", row_to_pick, " to keep \n")
      # set counter for next sample
      i=i+1
    }
    # remove duplicate rows to thin samples
    ROD_data<-ROD_data[ROD_data$thin == F,]
  }else{
    # set string for all data
    thin_str = ""
  }
  # remove remaining NA values
  ROD_data<-ROD_data[!is.na(ROD_data$PosVsNeg),]
  # count samples per category
  total_samples<-dim(ROD_data)[1]
  positive_samples<-length(ROD_data$PosVsNeg[ROD_data$PosVsNeg == "Ceratocystis positive"])
  negative_samples<-total_samples - positive_samples
  # print sample counts
  cat("##########################", "\n", total_samples, " total samples ", " with ", 
      positive_samples, " positive samples (",  positive_samples/total_samples, ") \n", 
      "##########################", "\n")
  
  #--- PREVALENCE OVER TIME ---#
  
  # create temporal rod dataset
  ROD_data_time<-ROD_data
  # create year-month column
  ROD_data_time$Yr_Mo<-format(ROD_data_time$DATE, "%Y-%m")
  
  # create positive numeric column
  ROD_data_time$pos<-0
  # set positive samples to 1, negatives to 0
  ROD_data_time[ROD_data_time$PosVsNeg %in% "Ceratocystis positive", "pos"]<-1
  ROD_data_time$neg<-0
  ROD_data_time[ROD_data_time$PN_NUM %in% "0", "neg"]<-1
  
  # run summary function
  pos_df<-summarySE(ROD_data_time, measurevar = "pos", groupvars = c("Yr_Mo"))
  neg_df<-summarySE(ROD_data_time, measurevar="neg", groupvars=c("Yr_Mo"))
  # subset data
  pos_df<-pos_df[,c("Yr_Mo", "sum")]
  neg_df<-neg_df[,c("Yr_Mo", "sum")]
  
  # merged data to calculate prevalence
  df_prevalence_by_time<-merge(pos_df, neg_df, by = "Yr_Mo", all = T)
  names(df_prevalence_by_time) = c("Yr_Mo", "Positive", "Not_detected")
  # calculate prevalence
  df_prevalence_by_time$prevalence<-df_prevalence_by_time$Positive/(df_prevalence_by_time$Positive + df_prevalence_by_time$Not_detected)
  df_prevalence_by_time$N<-(df_prevalence_by_time$Positive + df_prevalence_by_time$Not_detected)
  df_prevalence_by_time<-na.omit(df_prevalence_by_time)
  df_prevalence_by_time<-df_prevalence_by_time[,c("Yr_Mo", "N", "prevalence")]
  
  # # create zoo series data frame
  # z<-read.zoo(df_prevalence_by_time, FUN = as.yearmon)
  # # plot number of data samples collected over time
  # filename<-paste0(thin_str, proj_str, "_prevalence_over_time.tiff")
  # tiff(filename, width = 1500, height = 1200, units = "px", pointsize = 12, compression ="lzw")
  # plot(z, xlab = "Date", main = "")
  # dev.off()
  
  # check global prevalence (should match the numbers above)
  ROD_data_filtered<-ROD_data[!is.na(ROD_data$PosVsNeg),]
  # total_samples<-dim(ROD_data_filtered)[1]
  # jnk<-dim(ROD_data[ROD_data$PosVsNeg == "Ceratocystis positive",])[1]
  # jnk/total_samples
  
  #--- BIOCLIMS ---#
  
  # plot continuous variables
  # var_to_plot<-vars_to_plot[4]
  for(var_to_plot in vars_to_plot){
    # run summary function
    summary_DF<-summarySE(ROD_data_filtered, measurevar = var_to_plot, groupvars = c("PosVsNeg"))
    # calculate confidence intervals
    summary_DF$low_ci<-summary_DF[,var_to_plot] - summary_DF$ci
    summary_DF$high_ci<-summary_DF[,var_to_plot] + summary_DF$ci
    # calculate standard deviations
    summary_DF$low_sd<-summary_DF[,var_to_plot] - summary_DF$sd
    summary_DF$high_sd<-summary_DF[,var_to_plot] + summary_DF$sd
    
    # t-test calculations
    if(only_symptomatic){
      chars<-capture.output(print(t.test(ROD_data_filtered[ROD_data_filtered$PosVsNeg == "Ceratocystis positive", var_to_plot], 
                                         ROD_data_filtered[ROD_data_filtered$PosVsNeg == "Symptomatic only",var_to_plot])))
    }else{
      chars<-capture.output(print(t.test(ROD_data_filtered[ROD_data_filtered$PosVsNeg == "Ceratocystis positive", var_to_plot], 
                                         ROD_data_filtered[ROD_data_filtered$PosVsNeg == "Negative",var_to_plot])))
    }
    # save t-test output
    # writeLines(chars, con = file(paste0(rootDir, thin_str, proj_str, "_", var_to_plot, "_ttest.txt")))
    
    # plots with confidence intervals
    p_CI<-ggplot(data = summary_DF, aes_string(x = "PosVsNeg", y = var_to_plot)) +
      geom_bar(stat="identity") +
      geom_errorbar(aes(ymin = low_ci, ymax = high_ci), width = 0.2, position = position_dodge(0.9))
    p_CI<-p_CI + geom_text(aes(label = N), vjust = 1.6, color="white", size = 3.5) + theme_minimal()
    p_CI<-p_CI + theme(axis.title = element_blank())
    
    # # plot standard deviations
    # p_SD<-ggplot(data = summary_DF, aes_string(x = "PosVsNeg", y = var_to_plot)) +
    #   geom_bar(stat = "identity") + ylab(variable_cols_names[which(variable_cols == var_to_plot)]) +
    #   geom_errorbar(aes(ymin = low_sd, ymax = high_sd), width = 0.2, position=position_dodge(0.9))
    # save plot
    # ggsave(filename = paste0(thin_str, proj_str, "_", var_to_plot, "_mean_and_sd.tiff"), p_CI, compress="lzw")
    
    # percent of positive samples by bins
    bins<-quantile(x = ROD_data_filtered[,var_to_plot], probs = seq(0, 1, length.out = n_bins))
    bins_midpoints<-(bins[c(2:n_bins)]+bins[c(1:(n_bins-1))])/2
    bin_class<-signif(bins, 3)
    bin_class<-paste0(bin_class[1:(length(bin_class)-1)], "-", bin_class[2:(length(bin_class))])
    ROD_data_filtered_contEnv<-ROD_data_filtered
    ROD_data_filtered_contEnv$contClass<-.bincode(x=ROD_data_filtered_contEnv[,var_to_plot], breaks=bins)
    ROD_data_filtered_contEnv$contClass<-bin_class[ROD_data_filtered_contEnv$contClass]
    ROD_data_filtered_contEnv$count<-1
    
    # run summary function
    summary_DF<-summarySE(ROD_data_filtered_contEnv, measurevar = "count", groupvars = c("PosVsNeg", "contClass"))
    summary_DF<-summary_DF[,c("PosVsNeg", "contClass", "N")]
    summary_DF_pos<-summary_DF[summary_DF$PosVsNeg == "Ceratocystis positive", c("contClass", "N")]
    summary_DF_sym<-summary_DF[!summary_DF$PosVsNeg == "Ceratocystis positive", c("contClass", "N")]
    summary_DF_ratio<-merge(x = summary_DF_pos, y = summary_DF_sym, by = "contClass", all = T)
    summary_DF_ratio[is.na(summary_DF_ratio$N.x),"N.x"]<-0
    summary_DF_ratio[is.na(summary_DF_ratio$N.y),"N.y"]<-0
    summary_DF_ratio$Ratio_symp_only<-summary_DF_ratio$N.x/(summary_DF_ratio$N.x + summary_DF_ratio$N.y)
    summary_DF_ratio$N<-summary_DF_ratio$N.x + summary_DF_ratio$N.y
    summary_DF_ratio$dev_mean<-NA
    summary_DF_ratio$dev_sd<-NA
    summary_DF_ratio$dev_dist<-NA
    summary_DF_ratio$dev_p<-NA
    
    summary_DF_ratio<-summary_DF_ratio[!is.na(summary_DF_ratio$contClass),]
    i = 1
    for(i in c(1:dim(summary_DF_ratio)[1])){
      test_val<-summary_DF_ratio[i, "Ratio_symp_only"]
      summary_DF_ratio_all_else<-summary_DF_ratio[-i,]
      summary_DF_ratio_all_else<-summary_DF_ratio_all_else[summary_DF_ratio_all_else$N>min_N_for_comparison,]
      temp_dev_mean<-mean(summary_DF_ratio_all_else$Ratio_symp_only, na.rm = T)
      temp_dev_sd<-sd(summary_DF_ratio_all_else$Ratio_symp_only, na.rm = T)
      temp_dev_dist<-(test_val-temp_dev_mean)/temp_dev_sd
      temp_dev_p<-pnorm(test_val, mean = temp_dev_mean, sd = temp_dev_sd, lower.tail = F, log.p = FALSE)
      
      summary_DF_ratio$dev_mean[i]=temp_dev_mean
      summary_DF_ratio$dev_sd[i]=temp_dev_sd
      summary_DF_ratio$dev_dist[i]=temp_dev_dist
      summary_DF_ratio$dev_p[i]=1-temp_dev_p
    }
    # format output
    reorder<-order(match(summary_DF_ratio[,"contClass"], bin_class))
    summary_DF_ratio<-summary_DF_ratio[reorder, ]
    summary_DF_ratio$bin_midpoint<-bins_midpoints
    # save continuous variable calculations
    # write.csv(summary_DF_ratio, file = paste0(thin_str, proj_str, "_", var_to_plot, "_continuos_variable_by_bin.csv"), row.names = F)
    
    # save calculated correlations
    chars<-capture.output(print(cor.test(summary_DF_ratio$bin_midpoint, summary_DF_ratio$Ratio_symp_only)))
    # writeLines(chars, con = file(paste0(rootDir, thin_str, proj_str, "_", var_to_plot, "_prevalence_correlation.txt")))
    
    summary_DF_ratio_screened<-summary_DF_ratio[summary_DF_ratio$N > min_N_for_comparison,]
    summary_DF_ratio_screened<-summary_DF_ratio
    summary_DF_ratio_screened[,"contClass"]<-factor(summary_DF_ratio_screened[,"contClass"], levels = bin_class)
    p<-ggplot(data = summary_DF_ratio_screened, aes_string(x = "contClass", y = "Ratio_symp_only")) +
      geom_bar(stat = "identity") + xlab("") + ylab(graph_string2)
    p<-p + geom_text(aes(label = N), vjust = 1.6, color = "white", size = 3.5) + theme_minimal()
    p<-p + theme(axis.text.x = element_text(angle = 90, hjust = 1)) #+ ylim(ylim_vars)
    p<-p + xlab(variable_cols_names[which(variable_cols == var_to_plot)])
    # bio1: p<-p + xlab(expression(paste("Mean annual temperature ", "(", degree, "C)")))
    # bio12: p<-p + xlab(paste(variable_cols_names[which(variable_cols == var_to_plot)], "(mm)"))
    # bio7: p<-p + xlab(expression(paste("Temperature annual range ", "(", degree, "C)")))
    # save continuous variable plots by bins
    # ggsave(filename = paste0(thin_str, proj_str, "_", var_to_plot, "_cont_var_by_bins.tiff"), p, compress = "lzw")
    ggsave(filename = paste0(thin_str, proj_str, "_", var_to_plot, "_cont_var_by_bins.png"), p)
  }
  
  #--- LANDSCAPE CHARACTERISTICS ---#
  
  # for categorical values calculate ratio of symptomatic only to positive
  vars_to_plot2<-c("HIGAP", "BIOREG", "AGE_GROUP", "NAME")
  
  # var_to_plot = vars_to_plot2[1]
  ROD_data_filtered$count = 1
  for (var_to_plot in vars_to_plot2){
    summary_DF<-summarySE(ROD_data_filtered, measurevar = "count", groupvars = c("PosVsNeg", var_to_plot))
    summary_DF<-summary_DF[,c("PosVsNeg", var_to_plot, "N")]
    summary_DF_pos<-summary_DF[summary_DF$PosVsNeg == "Ceratocystis positive", c(var_to_plot, "N")]
    summary_DF_sym<-summary_DF[!summary_DF$PosVsNeg == "Ceratocystis positive", c(var_to_plot, "N")]
    summary_DF_ratio<-merge(x=summary_DF_pos, y=summary_DF_sym, by = var_to_plot, all = T)
    summary_DF_ratio[is.na(summary_DF_ratio$N.x),"N.x"] = 0
    summary_DF_ratio[is.na(summary_DF_ratio$N.y),"N.y"] = 0
    summary_DF_ratio$Ratio_symp_only<-summary_DF_ratio$N.y/(summary_DF_ratio$N.x+summary_DF_ratio$N.y)
    summary_DF_ratio$N<-summary_DF_ratio$N.x + summary_DF_ratio$N.y
    summary_DF_ratio$dev_mean = NA
    summary_DF_ratio$dev_sd = NA
    summary_DF_ratio$dev_dist = NA
    summary_DF_ratio$dev_p = NA
    
    i = 1
    for(i in c(1:dim(summary_DF_ratio)[1])){
      test_val<-summary_DF_ratio[i, "Ratio_symp_only"]
      summary_DF_ratio_all_else<-summary_DF_ratio[-i,]
      summary_DF_ratio_all_else<-summary_DF_ratio_all_else[summary_DF_ratio_all_else$N > min_N_for_comparison,]
      temp_dev_mean<-mean(summary_DF_ratio_all_else$Ratio_symp_only, na.rm = T)
      temp_dev_sd<-sd(summary_DF_ratio_all_else$Ratio_symp_only, na.rm = T)
      temp_dev_dist<-(test_val-temp_dev_mean)/temp_dev_sd
      temp_dev_p<-pnorm(test_val, mean = temp_dev_mean, sd = temp_dev_sd, lower.tail = TRUE, log.p = FALSE)
      
      summary_DF_ratio$dev_mean[i]<-temp_dev_mean
      summary_DF_ratio$dev_sd[i]<-temp_dev_sd
      summary_DF_ratio$dev_dist[i]<-temp_dev_dist
      summary_DF_ratio$dev_p[i]<-1-temp_dev_p
    }
    
    #dput(unique(summary_DF[,var_to_plot]))
    if(var_to_plot=="BIOREG"){factor_order=c("Kilauea", "NE Mauna Loa", "Maunakea", "Kohala", "NW Mauna Loa", "Hualalai", "Kona", "Kau")}
    if(var_to_plot=="AGE_GROUP"){factor_order=c("0-200 yr", "200-750 yr", "750-1500 yr", "1500-3000 yr", "3000-5000 yr", "5000-10000 yr", "10-30 ka", "30-50 ka", "50-140 ka", "140-780 ka", "multiple")}
    if(var_to_plot=="HIGAP"){factor_order=c("Closed ohia wet forest", "Low-stature ohia wet forest", "Closed koa-ohia wet forest", "Open koa-ohia wet forest", "Alien wet forest", "Native wet shrubland", 
                                            "Native wet cliff community", "Mixed native-alien wet shrubs and grass/sedges", "Alien wet shrubland", "Alien wet grassland", "Closed ohia mesic forest", 
                                            "Closed koa-ohia mesic forest", "Open ohia mesic forest", "Open koa-ohia mesic forest", "Alien mesic forest", "Native mesic shrubland", "Alien mesic shrubland",
                                            "Alien mesic grassland", "Open koa-mamane dry forest", "Mixed mamane-naio-native trees dry woodland", "Native dry shrubland","Alien dry forest", "Alien dry grassland", 
                                            "Uluhe ferns and native shrubs", "Very sparse vegetation to unvegetated", "Cultivated agriculture", "Developed open space", "High intensity developed",
                                            "Low intensity developed", "Medium intensity developed", "Plantation forest")}
    if(var_to_plot=="NAME"){factor_order=c("Puna Basalt", "Kau Basalt", "Laupahoehoe Volcanics", "Hamakua Volcanics", "Hawi Volcanics","Hualalai Volcanics", "Ninole Basalt", "Pololu Volcanics", "Caldera wall rocks", "Tephra deposits", "Alluvium")}  
    
    reorder<-order(match(summary_DF_ratio[,var_to_plot],factor_order))
    summary_DF_ratio<-summary_DF_ratio[reorder, ]
    
    # write.csv(summary_DF_ratio, file = paste0(thin_str, proj_str, "_", var_to_plot, ".csv"), row.names = F)
    summary_DF_ratio_screened<-summary_DF_ratio[summary_DF_ratio$N > min_N_for_comparison,]
    summary_DF_ratio_screened[,var_to_plot]<-factor(summary_DF_ratio_screened[,var_to_plot], levels = factor_order)
    p<-ggplot(data = summary_DF_ratio_screened, aes_string(x = var_to_plot, y = "Ratio_symp_only")) +
      geom_bar(stat = "identity") + xlab("") + ylab(graph_string)
    p = p + geom_text(aes(label = N), vjust = 1.6, color = "white", size = 3.5)+ theme_minimal()
    p = p + theme(axis.text.x = element_text(angle = 90, hjust = 1)) #+ ylim(0,1)
    p = p + xlab(variable_cols_names[which(variable_cols == var_to_plot)])
    if(var_to_plot == "HIGAP" | var_to_plot == "NAME"){
      # ggsave(filename = paste0(thin_str, proj_str, "_", var_to_plot, ".tiff"), p, compress = "lzw",
      ggsave(filename = paste0(thin_str, proj_str, "_", var_to_plot, ".png"), p,
             height = 5, width = 6, dpi = 300)
    }else{
      # ggsave(filename = paste0(thin_str, proj_str, "_", var_to_plot, ".tiff"), p, compress = "lzw",
      ggsave(filename = paste0(thin_str, proj_str, "_", var_to_plot, ".png"), p,
             height = 4, width = 6, dpi = 300)
    }
    
  }
  
  #--- POSITIVE RATIOS ---#
  
  # plot the inverse (ratio of positives)
  ROD_data_filtered$count = 1
  for(var_to_plot in vars_to_plot2){
    summary_DF<-summarySE(ROD_data_filtered, measurevar = "count", groupvars = c("PosVsNeg", var_to_plot))
    summary_DF<-summary_DF[,c("PosVsNeg", var_to_plot, "N")]
    summary_DF_pos<-summary_DF[summary_DF$PosVsNeg == "Ceratocystis positive", c(var_to_plot, "N")]
    summary_DF_sym<-summary_DF[!summary_DF$PosVsNeg == "Ceratocystis positive", c(var_to_plot, "N")]
    summary_DF_ratio<-merge(x = summary_DF_pos, y = summary_DF_sym, by = var_to_plot, all = T)
    summary_DF_ratio[is.na(summary_DF_ratio$N.x),"N.x"] = 0
    summary_DF_ratio[is.na(summary_DF_ratio$N.y),"N.y"] = 0
    summary_DF_ratio$Ratio_positives = summary_DF_ratio$N.x/(summary_DF_ratio$N.x + summary_DF_ratio$N.y)
    summary_DF_ratio$N<-summary_DF_ratio$N.x + summary_DF_ratio$N.y
    summary_DF_ratio$dev_mean = NA
    summary_DF_ratio$dev_sd = NA
    summary_DF_ratio$dev_dist = NA
    summary_DF_ratio$dev_p = NA
    
    i = 1
    for(i in c(1:dim(summary_DF_ratio)[1])){
      test_val<-summary_DF_ratio[i, "Ratio_positives"]
      summary_DF_ratio_all_else<-summary_DF_ratio[-i,]
      summary_DF_ratio_all_else<-summary_DF_ratio_all_else[summary_DF_ratio_all_else$N > min_N_for_comparison,]
      temp_dev_mean<-mean(summary_DF_ratio_all_else$Ratio_positives, na.rm=T)
      temp_dev_sd<-sd(summary_DF_ratio_all_else$Ratio_positives, na.rm=T)
      temp_dev_dist<-(test_val-temp_dev_mean)/temp_dev_sd
      temp_dev_p<-pnorm(test_val, mean = temp_dev_mean, sd = temp_dev_sd, lower.tail = F, log.p = FALSE)
      
      summary_DF_ratio$dev_mean[i]<-temp_dev_mean
      summary_DF_ratio$dev_sd[i]<-temp_dev_sd
      summary_DF_ratio$dev_dist[i]<-temp_dev_dist
      summary_DF_ratio$dev_p[i]<-1-temp_dev_p
    }
    
    if(var_to_plot=="BIOREG"){factor_order=c("Kilauea", "NE Mauna Loa", "Maunakea", "Kohala", "NW Mauna Loa", "Hualalai", "Kona", "Kau")}
    if(var_to_plot=="AGE_GROUP"){factor_order=c("0-200 yr", "200-750 yr", "750-1500 yr", "1500-3000 yr", "3000-5000 yr", "5000-10000 yr", "10-30 ka", "30-50 ka", "50-140 ka", "140-780 ka", "multiple")}
    if(var_to_plot=="HIGAP"){factor_order=c("Closed ohia wet forest", "Low-stature ohia wet forest", "Closed koa-ohia wet forest", "Open koa-ohia wet forest", "Alien wet forest", "Native wet shrubland", 
                                            "Native wet cliff community", "Mixed native-alien wet shrubs and grass/sedges", "Alien wet shrubland", "Alien wet grassland", "Closed ohia mesic forest", 
                                            "Closed koa-ohia mesic forest", "Open ohia mesic forest", "Open koa-ohia mesic forest", "Alien mesic forest", "Native mesic shrubland", "Alien mesic shrubland",
                                            "Alien mesic grassland", "Open koa-mamane dry forest", "Mixed mamane-naio-native trees dry woodland", "Native dry shrubland","Alien dry forest", "Alien dry grassland",
                                            "Uluhe ferns and native shrubs", "Very sparse vegetation to unvegetated", "Cultivated agriculture", "Developed open space", "High intensity developed",
                                            "Low intensity developed", "Medium intensity developed", "Plantation forest")}
    if(var_to_plot=="NAME"){factor_order=c("Puna Basalt", "Kau Basalt", "Laupahoehoe Volcanics", "Hamakua Volcanics", "Hawi Volcanics","Hualalai Volcanics", "Ninole Basalt", "Pololu Volcanics", "Caldera wall rocks", "Tephra deposits", "Alluvium")}
    
    reorder<-order(match(summary_DF_ratio[,var_to_plot], factor_order))
    summary_DF_ratio<-summary_DF_ratio[reorder, ]
    
    # write.csv(summary_DF_ratio, file=paste0(thin_str, proj_str, "_", var_to_plot, "_positive_ratio.csv"), row.names = F)
    
    summary_DF_ratio_screened<-summary_DF_ratio[summary_DF_ratio$N > min_N_for_comparison,]
    summary_DF_ratio_screened[,var_to_plot]<-factor(summary_DF_ratio_screened[,var_to_plot], levels = factor_order)
    p<-ggplot(data = summary_DF_ratio_screened, aes_string(x = var_to_plot, y = "Ratio_positives")) +
      geom_bar(stat = "identity") + xlab("") + ylab(graph_string2)
    p = p +geom_text(aes(label = N), vjust = 1.6, color = "white", size = 3.5) + theme_minimal()
    p = p + theme(axis.text.x = element_text(angle = 90, hjust = 1))#+ylim(0,1)
    p = p + xlab(variable_cols_names[which(variable_cols == var_to_plot)])
    if(var_to_plot == "HIGAP" | var_to_plot == "NAME"){
      # ggsave(filename = paste0(thin_str, proj_str, "_", var_to_plot, "_positive_ratio.tiff"), p, compress = "lzw")
      ggsave(filename = paste0(thin_str, proj_str, "_", var_to_plot, "_positive_ratio.png"), p,
             height = 5, width = 6, dpi = 300)
    }else{
      # ggsave(filename = paste0(thin_str, proj_str, "_", var_to_plot, "_positive_ratio.tiff"), p, compress = "lzw")
      ggsave(filename = paste0(thin_str, proj_str, "_", var_to_plot, "_positive_ratio.png"), p,
             height = 4, width = 6, dpi = 300)
    }

  }
  
}

### END THINNED SUMMARY ANALYSIS ###