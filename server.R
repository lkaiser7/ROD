### UI for ROD Shiny App ###

##################
##### SET UP #####
##################

# load necessary packages
library(shiny)
library(rgdal)
library(raster)
library(ggplot2)
library(gridExtra)

# set working directory to app location
appDir<-"E:/Dropbox/PIERC/ROD_updates/rod_portal/"  # hard drive
# DO NOT RUN setwd() ON SHINY SERVER (DOES NOT EXIST): setwd(appDir)

#####################
##### LOAD DATA #####
#####################

# load coastlines
hi_coast<-readOGR("hawaii_coast", "Main_Hawaiian_Islands_simple3")
bi_coast<-readOGR("hawaii_coast", "bi_coast")

# load suitability projection
# suit_raster<-raster("www/CL_Suitability_current_ROCwmean.tif")

# load master rod data set
rod<-read.csv("rod_dist_vars.csv", header = TRUE)
# load all rod sampled site data
all_rod<-read.csv("all_dist_vars.csv", header = TRUE)

# # check data
# table(rod$SP_NM, useNA = 'ifany')
# table(all_rod$SP_NM, useNA = 'ifany')

# remove incomplete cases for hawaii island
rod<-rod[which(is.na(rod$bio1) == F),]
rod<-rod[which(is.na(rod$SP_NM) == F), ]

# load elevation data for background maps
hi_elev<-read.csv("hi_elev_spdf.csv", header = TRUE)

# load ohia data set with extracted variables
hi_ohia<-read.csv("ohia_dist_vars.csv", header = TRUE)

# load data for all hawaii island for ggplot base maps
hi_bc<-read.csv("all_hi_2000_bioclims.csv")  # bioclims
hi_veg<-read.csv("all_hi_veg.csv")           # higap classes/biomes
hi_bioreg<-read.csv("all_hi_bioreg.csv")     # bioregions
hi_geol<-read.csv("all_hi_geol.csv")         # geology 
hi_roads<-read.csv("all_hi_roads.csv")       # roads

##################
##### SERVER #####
##################

# run server to display UI selections 
shinyServer<-function(input, output){ 
  
  #----- THREAT DISTRIBUTION -----#
  
  doi_url<-a("Fostering real-time climate adaptation: Analyzing past, current, and forecast temperature 
             to understand the dynamic risk to Hawaiian honeycreepers from avian malaria",
             href = "https://authors.elsevier.com/a/1ZHQS1L~GwKQnm")
  output$fem_doi<-renderUI({
    tagList(doi_url)
  })
  
  # describe plot image
  output$sdm_info<-renderText({
    "To determine the potential threat of these Ceratocystis species for forests on all of the
     main Hawaiian Islands, we created a species distribution model to identify the areas that 
     each of these species could potentially occupy on other neighboring islands based on current 
     confirmed Ceratocystis locations."
  })
  
  output$sre_image<-renderImage({
    
    if(input$sp == "C. lukuohia"){
      sp = "CL"
    }else{
      sp = "CH"
    }
    
    # store image file name
    # sdm_plot<-normalizePath(file.path('./www', paste0('Fig4_', sp, '_suitability_rod_gam_current500m_image', '.tif')))
    sdm_plot<-normalizePath(file.path('./www', paste0('Fig4_', sp, '_suitability_rod_gam_current500m_image', '.png')))
    list(src = sdm_plot, height = "150%")
  }, deleteFile = FALSE)
  
  # output$sre_plots<-renderPlot({
    # plot model results
    # tiff("www/CL_sdm_suitability_image.tif",
    #      res = 300, units = "in", pointsize = 12, width = 10, height = 8, compression = "lzw")
    # plot(suit_raster, col = rev(heat.colors(50)), 
    #      main = 'Projected Suitability of C. lukuohia')
    # plot(hi_coast, add = T)
    # dev.off()
    
  #----- MAP -----# 
  
  output$fem_doi0<-renderUI({
    tagList(doi_url)
  })
  
  url_org<-a("RapidOhiaDeath.org", href = "https://cms.ctahr.hawaii.edu/rod/")
  output$rod_url<-renderUI({
    tags$h4(tagList("For more information, please visit: ", url_org), style = "color:#dd4814")
  })
  
  output$ctahr_img<-renderImage({
    # store file name as reactive
    ctahr_img<-normalizePath(file.path('./www', paste0('CTAHR_ROD', '.png')))
    list(src = ctahr_img, height = "175%")
  }, deleteFile = FALSE)

  # # display CTAHR map
  # output$map_url<-renderText({
  #   if(input$password == "usgsrodapp"){
  #     return(NULL)
  #   }else{
  #     c('<img src="',
  #       # OLD: "https://gms.ctahr.hawaii.edu/gs/handler/getmedia.ashx?moid=31395&dt=3&g=5", '">')
  #       "https://gms.ctahr.hawaii.edu/gs/handler/getmedia.ashx?moid=67071&dt=3&g=5", '">')
  #   }
  # })
  
  # link to CTAHR map
  ctahr_org<-a("UH Manoa CTAHR", href = "https://cms.ctahr.hawaii.edu/rod/THE-DISEASE/DISTRIBUTION")
  output$ctahr_map<-renderUI({
    if(input$password == "usgsrodapp"){
      return(NULL)
    }else{
      tags$h4(tagList("Map courtesy of ", ctahr_org, "and their partners.
                    While this map shows general locations for confirmed cases of Rapid Ohia Death,
                    it does not necessarily indicate the pattern of disease spread over time."))
    }
  })
  
  # create reactive dataset based on rod strain
  map_data<-reactive({
    rod[rod$SP_NM %in% input$type1, ]
  })
  
  output$RODmap<-renderPlot({
    if(input$password == "usgsrodapp"){
      # create base map from elevation raster file
      base_map<-ggplot(data = hi_elev, aes(x = x, y = y)) +
        geom_raster(aes(fill = hawaii)) +
        scale_fill_gradientn(colours = terrain.colors(8924), guide = F) +
        labs(x = "Longitude", y = "Latitude") + coord_fixed() +
        theme_bw() + theme(plot.title = element_text(hjust = 0.5))
      
      # add confirmed locations of rod detections
      if(input$neg == "NO"){
        # map rod occurrence points by type
        base_map + geom_point(data = map_data(), aes(x = X_LON, y = Y_LAT, shape = factor(SP_NM)),
                              size = 4, alpha = 0.8) + scale_shape_discrete(solid = F, name = "Species") +
          ggtitle("Confirmed Ceratocystis Locations")
      }else{
        # add negative detection sites to rod map
        base_map + geom_point(data = map_data(), aes(x = X_LON, y = Y_LAT, shape = factor(SP_NM)),
                              size = 3, alpha = 0.8) + scale_shape_discrete(solid = F, name = "Species") +
          geom_point(data = rod[which(rod$PN_NUM == 0),], aes(x = X_LON, y = Y_LAT),
                     shape = 4, size = 1, colour = "red") + ggtitle("Sample Sites")
      }
    }
  }, height = 500) #, width = 1000)
  
  #----- TEMPORAL SPREAD -----#

  output$fem_doi1<-renderUI({
    tagList(doi_url)
  })
  
  time_data<-reactive({
    rod[rod$YEAR %in% input$type, ]
  })
  
  output$TimeMap<-renderPlot({
    # spread of rod over time
    time_map<-ggplot(data = hi_elev, aes(x = x, y = y)) +
      geom_raster(aes(fill = hawaii), alpha = 0.8) +
      scale_fill_gradientn(colours = terrain.colors(8924), guide = F) +
      scale_color_manual(values = rev(rainbow(12)[c(-3:-7)]), name = "Year") + 
      coord_fixed() + theme_bw() + theme(plot.title = element_text(hjust = 0.5))
    
    time_plot<-ggplot(data = rod[which(rod$SP_NM != "Not Detected"),]) + 
      geom_bar(aes(YEAR, fill = SP_NM), position = "dodge") +
      labs(x = NULL, y = "Number of Samples") + ggtitle("Samples Collected Per Year") + 
      scale_fill_manual(name = "Species", labels = c("C. huliohia", "C. lukuohia"), 
                         values = c("#CC3300", "#FF9933")) + 
      theme_bw() + theme(plot.title = element_text(hjust = 0.5))
    
    # arrange plots
    if(input$password1 == "usgsrodapp"){
      time_map + 
        geom_point(data = time_data()[which(time_data()$SP_NM != "Not Detected"),], 
        # DATA CHECK: geom_point(data = rod[which(rod$SP_NM != "Not Detected"),], 
                   aes(x = X_LON, y = Y_LAT, colour = factor(YEAR)),
                   shape = 21, fill = NA, size = 4, alpha = 0.5) +
        labs(x = "Longitude", y = "Latitude") + ggtitle("Confirmation of Ceratocystis Over Time")
    }else{
      time_plot
    }
  }, height = 500)

  #----- BIOCLIMS -----# 

  output$fem_doi2<-renderUI({
    tagList(doi_url)
  })
  
  # describe plot image
  output$bc_info<-renderText({
    "For these continuous climatic variables, we used sample dates to characterize expansion of the 
     environmental range of C. lukuohia over time within the context of 'ohi'a forest environmental range.
     Additionally, for each of these variables, we examined pathogen prevalence rates
     (i.e., the number of positive detections versus total number of samples) across environmental gradients."
  })
  
  # bioclim ranges
  output$bc_plots<-renderImage({
    if(input$type2 == "C. lukuohia"){
      select_sp = "CL"
    }else{
      select_sp = "CH"
    }
    
    # store file name as reactive 
    # clim_plot<-normalizePath(file.path('./www', paste0(select_sp, '_all_', input$bioclim, '_means', '.tiff')))
    clim_plot<-normalizePath(file.path('./www', paste0(select_sp, '_all_', input$bioclim, '_means', '.png')))
    list(src = clim_plot, height = "100%")
  }, deleteFile = FALSE)
  
  # bioclim prevalence
  output$bio_img<-renderImage({
    if(input$type2 == "C. lukuohia"){
      select_sp = "CL"
    }else{
      select_sp = "CH"
    }

    # store file name as reactive
    # clim_img<-normalizePath(file.path('./www', paste0(select_sp, '_all_', input$bioclim, '_cont_var_by_bins', '.tiff')))
    clim_img<-normalizePath(file.path('./www', paste0(select_sp, '_all_', input$bioclim, '_cont_var_by_bins', '.png')))
    list(src = clim_img, height = "100%")
  }, deleteFile = FALSE)
  
  # create reactive dataset based on rod strain
  bc_data<-reactive({
    rod[rod$SP_NM %in% input$type2, ]
  })
  
  # find bins for histogram based on ohia distribution
  ohia_bins<-reactive({
    ohia_hist<-hist(hi_ohia[input$bioclim][ , 1], plot = F)
    length(ohia_hist$breaks)
  })

  output$bc_map<-renderPlot({
    # plot bioclim raster
    # PREV: rast_plot<-ggplot(data = hi_bc, aes(x = x, y = y)) +
    #         geom_raster(aes(fill = hi_bc[input$bioclim][ , 1])) +
    rast_plot<-ggplot(data = hi_bc, aes(x = round(x, 4), y = round(y, 4))) +
      geom_tile(aes(fill = hi_bc[input$bioclim][ , 1])) +
      scale_fill_gradientn(name = "Variable", colours = rev(terrain.colors(ohia_bins()))) +
      geom_point(data = bc_data(),
                 aes(x = X_LON, y = Y_LAT, shape = factor(SP_NM)), size = 4, alpha = 0.8) +
      scale_shape_discrete(solid = F, name = "Species") +
      labs(x = "Longitude", y = "Latitude") + coord_fixed() +
      theme_bw() + theme(plot.title = element_text(hjust = 0.5))
    
    if(input$password2 == "usgsrodapp"){
      rast_plot
    }
  })
  
  #----- VEG & BIOREGS -----#  

  output$fem_doi3<-renderUI({
    tagList(doi_url)
  })
  
  # describe plot image
  output$vb_info<-renderText({
    "For each of these variables, we examined pathogen prevalence rates 
    (i.e., the number of positive detections versus total number of samples) 
    across different landscape classes."
  })
  
  # vegetation and bioregion prevalence
  output$vb_img<-renderImage({
    if(input$type3 == "C. lukuohia"){
      select_sp = "CL"
    }else{
      select_sp = "CH"
    }
    
    # store file name as reactive
    # vegbio_img<-normalizePath(file.path('./www', paste0(select_sp, '_all_', input$ecovar, '_positive_ratio.tiff')))
    vegbio_img<-normalizePath(file.path('./www', paste0(select_sp, '_all_', input$ecovar, '_positive_ratio.png')))
    list(src = vegbio_img, height = "100%")
  }, deleteFile = FALSE)
  
  # create reactive dataset based on rod strain
  vc_data<-reactive({
    # all_rod[all_rod$STRAIN %in% input$type3, ]
    rod[rod$SP_NM %in% input$type3, ]
  })
  
  output$bv_plot<-renderPlot({
    # plot vegetation classes raster    
    veg_plot<-ggplot(data = hi_veg, aes(x = x, y = y)) +
      geom_raster(aes(fill = HIGAP)) + 
      scale_fill_gradientn(colours = unique(hi_veg$COLOR), guide = F) +
      # ALT: scale_fill_gradientn(colours = as.character(levels(hi_veg$COLOR)), guide = F) +
      geom_point(data = vc_data(), aes(x = X_LON, y = Y_LAT, shape = factor(SP_NM)),
                 size = 4, alpha = 0.8, fill = "gray") +
      scale_shape_discrete(solid = F, name = "Species") +
      labs(x = "Longitude", y = "Latitude") + coord_fixed() + 
      theme_bw() + theme(plot.title = element_text(hjust = 0.5))
    
    # plot bioregions raster 
    bioreg_plot<-ggplot(data = hi_bioreg, aes(x = x, y = y)) +
      geom_raster(aes(fill = bioregions)) + 
      scale_fill_gradientn(colours = unique(hi_bioreg$COLOR), guide = F) +
      # ALT: scale_fill_gradientn(colours = as.character(levels(hi_bioreg$COLOR)), guide = F) +
      geom_point(data = vc_data(),
                 aes(x = X_LON, y = Y_LAT, shape = factor(SP_NM)),
                 size = 4, alpha = 0.8, fill = "gray") +
      scale_shape_discrete(solid = F, name = "Species") +
      labs(x = "Longitude", y = "Latitude") + coord_fixed() + 
      theme_bw() + theme(plot.title = element_text(hjust = 0.5))    
    
    # arrange plots based on user ui
    if(input$ecovar == "HIGAP"){
      if(input$password3 == "usgsrodapp"){
        veg_plot
      }
    } else {
      if(input$password3 == "usgsrodapp"){
        bioreg_plot 
      }
    }
  })
  
  #----- GEOLOGY -----#  

  output$fem_doi4<-renderUI({
    tagList(doi_url)
  })
  
  # describe plot image
  output$geol_info<-renderText({
    "For each of these variables, we examined pathogen prevalence rates 
    (i.e., the number of positive detections versus total number of samples) 
    across different landscape classes."
  })
  
  # geologic prevalence
  output$geol_img<-renderImage({
    if(input$type4 == "C. lukuohia"){
      select_sp = "CL"
    }else{
      select_sp = "CH"
    }
    
    # store file name as reactive
    # geo_img<-normalizePath(file.path('./www', paste0(select_sp, '_all_', input$geol, '_positive_ratio.tiff')))
    geo_img<-normalizePath(file.path('./www', paste0(select_sp, '_all_', input$geol, '_positive_ratio.png')))
    list(src = geo_img, height = "100%")
  }, deleteFile = FALSE)
  
  # create reactive dataset based on rod strain
  geol_data<-reactive({
    # all_rod[all_rod$STRAIN %in% input$type4, ]
    rod[rod$SP_NM %in% input$type4, ]
  })
  
  output$geol_plots<-renderPlot({
    # plot geologic ages raster    
    age_plot<-ggplot(data = hi_geol, aes(x = x, y = y)) +
      geom_raster(aes(fill = AGE_GROUP)) + 
      #scale_fill_gradientn(colours = cbind(levels(hi_ohia$AGE_GROUP), cm.colors(length(levels(hi_ohia$AGE_GROUP))))[,2], guide = F) +
      # ALT: scale_fill_gradientn(colours = as.character(levels(hi_geol$AGE_COLOR)), guide = F) +
      scale_fill_gradientn(colours = unique(hi_geol$AGE_COLOR), guide = F) +
      geom_point(data = geol_data(),
                 aes(x = X_LON, y = Y_LAT, shape = factor(SP_NM)),
                 size = 4, alpha = 0.8, fill = "gray") +
      scale_shape_discrete(solid = F, name = "Species") +
      labs(x = "Longitude", y = "Latitude") + coord_fixed() + 
      theme_bw() + theme(plot.title = element_text(hjust = 0.5))
    
    # plot geologic names raster 
    type_plot<-ggplot(data = hi_geol, aes(x = x, y = y)) +
      geom_raster(aes(fill = NAME)) + 
      #scale_fill_gradientn(colours = cbind(levels(hi_ohia$NAME), cm.colors(length(levels(hi_ohia$NAME))))[,2], guide = F) +
      # ALT: scale_fill_gradientn(colours = as.character(levels(hi_geol$NAME_COLOR)), guide = F) +
      scale_fill_gradientn(colours = unique(hi_geol$NAME_COLOR), guide = F) +
      geom_point(data = geol_data(),
                 aes(x = X_LON, y = Y_LAT, shape = factor(SP_NM)),
                 size = 4, alpha = 0.8, fill = "gray") +
      scale_shape_discrete(solid = F, name = "Species") +
      labs(x = "Longitude", y = "Latitude") + coord_fixed() + 
      theme_bw() + theme(plot.title = element_text(hjust = 0.5))    
    
    # arrange plots based on user ui
    if(input$geol == "AGE_GROUP"){
      if(input$password4 == "usgsrodapp"){
        age_plot
      }
    } else {
      if(input$password4 == "usgsrodapp"){
        type_plot 
      }
    }
  })
  
  #----- ROAD ACCESIBILITY -----#
  
  output$fem_doi5<-renderUI({
    tagList(doi_url)
  })
  
  # create reactive dataset based on rod strain
  road_data<-reactive({
    rod[rod$SP_NM %in% input$type5, ]
  })
  
  output$road_plots<-renderPlot({
    # plot roads raster    
    road_plot<-ggplot(data = hi_roads, aes(x = x, y = y), na.rm = T) +
      geom_raster(na.rm = T, aes(fill = km_road_dist)) + 
      scale_fill_gradientn(name = "Distance (km)", colours = rev(heat.colors(10))) +
      geom_point(data = road_data(), 
                 aes(x = X_LON, y = Y_LAT, shape = factor(SP_NM)), size = 4, alpha = 0.8) + 
      scale_shape_discrete(solid = F, name = "Species") + ggtitle("Distance From Road") +
      labs(x = "Longitude", y = "Latitude") + coord_fixed() + 
      theme_bw() + theme(plot.title = element_text(hjust = 0.5))
    
    # histograms of road distance data
    ohia_plot<-ggplot(data = hi_ohia) + 
      geom_histogram(na.rm = T, aes(x = ROADS), bins = 18,
                     fill = rev(heat.colors(18)), colour = "grey") +  
      ggtitle("Entire Ohia Distribution") + xlim(0, 8) +
      labs(x = NULL, y = "Frequency") + theme_bw() + theme(plot.title = element_text(hjust = 0.5))   
    
    rod_plot<-ggplot(data = road_data()) + 
      #geom_histogram(na.rm = T, aes(x = road_data()["ROADS"]), bins = 13, # OLD
      geom_histogram(na.rm = T, aes(x = road_data()$ROADS), bins = 13,
                     fill = rev(heat.colors(13)), colour = "grey") +  xlim(0, 8) +
      ggtitle("Ceratocystis Distribution") + labs(x = NULL, y = "Frequency") + 
      theme_bw() + theme(plot.title = element_text(hjust = 0.5))

    # display blank screen if no variable is selected
    if(is.null(input$type5))
      return(road_plot)
    
    # arrange plots
    if(input$password5 == "usgsrodapp"){
      grid.arrange(arrangeGrob(ohia_plot, rod_plot), road_plot, nrow = 2, ncol = 1) 
    }else{
      grid.arrange(arrangeGrob(ohia_plot, rod_plot), nrow = 2, ncol = 1)
    }
  }, height = 850)
  
  # end shiny SERVER
}

# run shiny web app locally: shinyApp(ui = shinyUI, server = shinyServer)