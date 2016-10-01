### UI for ROD Shiny App ###

##################
##### SET UP #####
##################

# load necessary packages
library(shiny)
library(ggmap)
library(rgdal)
library(gridExtra)
#library(cowplot)

# set working directory to app location
appDir<-"Y:/PICCC_analysis/ROD_distribution/rod_app/"
setwd(appDir)

#####################
##### LOAD DATA #####
#####################

# load master rod data set
rod<-read.csv("rod_dist_vars.csv", header = TRUE)

# load ohia data set with extracted variables
hi_ohia<-read.csv("ohia_dist_vars.csv", header = TRUE)

# load data for all hawaii island for ggplot base maps
hi_bc<-read.csv("all_HI_2000_bioclims.csv")  # bioclims
hi_veg<-read.csv("all_HI_veg.csv")           # higap classes/biomes
hi_bioreg<-read.csv("all_HI_bioreg.csv")     # bioregions
hi_geol<-read.csv("all_HI_geol.csv")         # geology 

# load data for archipelago-wide maps
all_islands<-read.csv("all_islands_mask.csv")
all_sre<-read.csv("all_islands_rod_current_CE.csv")
all_dofaw<-read.csv("all_islands_dofaw_survey.csv")

# get Google map of Big Island
BI<-get_googlemap(center = c(lon = -155.5, lat = 19.5), 
                  zoom = 9, maptype = "hybrid") 

##################
##### SERVER #####
##################

# run server to display UI selections 
shinyServer<-function(input, output){ 
  
  #----- MAP -----# 
  
  # create reactive dataset based on rod strain
  map_data<-reactive({
    rod[rod$POS_NEG %in% input$type1, ]
  })
  
  output$RODmap<-renderPlot({
    
    # map rod occurrence points by type
    ggmap(BI) + 
      geom_point(data = map_data(), aes(x = j_LON, y = j_LAT, fill = factor(POS_NEG)),
                 shape = 22, size = 5, colour = "grey") +
      labs(x = "Longitude", y = "Latitude") + ggtitle("Confirmed ROD Locations") +
      scale_fill_discrete(name = "Strain") + theme_gray() 
  }, width = 800, height = 800)
  
  #----- CLIMATE ENVELOPE -----#
  
  # create reactive dataset based on climate envelope of rod strain
  ce_fill<-reactive({
    ce_data<-as.numeric(all_sre[input$sres][ , 1])
    ce_data[which(ce_data != 0)]
  })
  
  output$sre_plots<-renderPlot({
    
    # map threat of rod    
    ce_plot<-ggplot(all_sre) + 
      geom_tile(data = all_islands, aes(x = x, y = y, fill = layer), colour = "gray") + 
      geom_raster(data = all_sre[which(all_sre[input$sres] != 0), ], 
                  aes(x = x, y = y, fill = ce_fill())) + 
      scale_fill_continuous(breaks = 1:3,
                            guide = guide_legend(title = "Ohia Distribution\nSuitable Climatic\nRange for ROD",
                                                 title.position = "right", label = F)) + 
      labs(x = "Longitude", y = "Latitude") +
      ggtitle("Potential ROD Sites within Ohia Distribution") + coord_fixed() + theme_gray()
    
    # map dofaw surveys
    dofaw_plot<-ce_plot + geom_path(data = all_dofaw, aes(x = long, y = lat, group = group, color = name)) +
      scale_color_gradient(low = "yellow", high = "red", breaks = c(1, 6, 11, 16),
                           labels = c("Light < 10%", "Moderate 10-30%", 
                                      "Severe 30-50%", "Very Severe > 50%"), 
                           na.value = "black") +
      labs(x = "Longitude", y = "Latitude", colour = "ROD Potential") +
      ggtitle("Potential ROD Sites") + coord_fixed() + theme_gray()
    
    # rod proximity to roads histogram    
    road_hist<-ggplot(data = rod) + 
      geom_histogram(aes(x = roads, fill = POS_NEG), # position = "dodge",
                     bins = 8, colour = "gray") + 
      ggtitle("ROD Proximity to Roads") + labs(x = "Distance to Road (km)", y = "Count") +
      scale_fill_discrete(name = "Strain") + theme_gray()
    
    # add dofaw survey polygons if selected
    if(input$dofaw == "YES"){
      grid.arrange(dofaw_plot, road_hist, nrow = 2, ncol = 1)
    } else {
      grid.arrange(ce_plot, road_hist, nrow = 2, ncol = 1)
    } 
    
  }, width = 800, height = 800)
  
  
  #----- BIOCLIMS -----# 
  
  # create reactive dataset based on rod strain
  bc_data<-reactive({
    rod[rod$POS_NEG %in% input$type2, ]
  })
  
  # find bins for histogram based on ohia distribution
  ohia_bins<-reactive({
    ohia_hist<-hist(hi_ohia[input$bioclim][ , 1])
    length(ohia_hist$breaks)
  })
  
  # set x limits to make plots comparative
  x_range<-reactive({
    range(hi_bc[input$bioclim])
  })
  
  output$bc_plots<-renderPlot({
    
    # plot bioclim raster    
    rast_plot<-ggplot(data = hi_bc, aes(x = x, y = y)) +
      geom_raster(aes(fill = hi_bc[input$bioclim])) + 
      scale_fill_gradientn(name = "Variable", colours = rev(terrain.colors(ohia_bins()))) +
      geom_point(data = bc_data(), 
                 aes(x = j_LON, y = j_LAT, shape = factor(POS_NEG)), size = 3) + 
      scale_shape_discrete(solid = F, name = "Strain") +
      labs(x = "Longitude", y = "Latitude") + coord_fixed() + theme_gray()
    
    # histograms of bioclimatic variable data
    ohia_plot<-ggplot(data = hi_ohia) + 
      geom_histogram(na.rm = T, aes(x = hi_ohia[input$bioclim]), bins = ohia_bins(),
                     fill = rev(terrain.colors(ohia_bins())), colour = "grey") +  
      ggtitle("Entire Ohia Distribution") + xlim(x_range()) +
      labs(x = NULL, y = "Frequency") + theme_gray()   
    
    rod_plot<-ggplot(data = bc_data()) + 
      geom_histogram(na.rm = T, aes(x = bc_data()[input$bioclim]), bins = ohia_bins(),
                     fill = rev(terrain.colors(ohia_bins())), colour = "grey") +  
      ggtitle("ROD Only Distribution") + xlim(x_range()) +
      labs(x = NULL, y = "Frequency") + theme_gray()  
    
    # arrange plots
    grid.arrange(rast_plot, arrangeGrob(ohia_plot, rod_plot), nrow = 2, ncol = 1)

  }, width = 800, height = 800)
  
  #----- VEG & BIOREGS -----#  
  
  # create reactive dataset based on rod strain
  vc_data<-reactive({
    rod[rod$POS_NEG %in% input$type3, ]
  })
  
  # create reactive data table of vegetation classes
  vc_table<-reactive({
    ohia_tab<-data.frame(table(hi_ohia$HIGAP))
    ohia_tab<-merge(ohia_tab, data.frame(table(vc_data()["HIGAP"])), by = "Var1", all = T)
    names(ohia_tab)<-c("landcover", "all_ohia", "strain")
    ohia_tab[is.na(ohia_tab)]<-0
    ohia_tab
  })
  
  # create reactive data table of bioregions
  br_table<-reactive({
    br_tab<-data.frame(table(hi_ohia$BIOREG))
    br_tab<-merge(br_tab, data.frame(table(vc_data()["BIOREG"])), by = "Var1", all = T)
    names(br_tab)<-c("bioregion", "all_ohia", "strain")
    br_tab[is.na(br_tab)]<-0
    br_tab
  })
  
  output$bv_plots<-renderPlot({
    
    # plot vegetation classes raster    
    veg_plot<-ggplot(data = hi_veg, aes(x = x, y = y)) +
      geom_raster(aes(fill = HIGAP)) + 
      scale_fill_gradientn(colours = rev(topo.colors(32)), guide = F) +
      geom_point(data = vc_data(), aes(x = j_LON, y = j_LAT, shape = factor(POS_NEG)), 
                 size = 3, fill = "gray") + 
      scale_shape_discrete(solid = F, name = "Strain") +
      labs(x = "Longitude", y = "Latitude") + coord_fixed() + theme_gray() 
    
    # plot bioregions raster 
    bioreg_plot<-ggplot(data = hi_bioreg, aes(x = x, y = y)) +
      geom_raster(aes(fill = bioregions)) + 
      scale_fill_gradientn(colours = rev(topo.colors(8)), guide = F) +
      geom_point(data = vc_data(), 
                 aes(x = j_LON, y = j_LAT, shape = factor(POS_NEG)), 
                 size = 3, fill = "gray") + 
      scale_shape_discrete(solid = F, name = "Strain") +
      labs(x = "Longitude", y = "Latitude") + coord_fixed() + theme_gray()    
    
    # barplot of vegetation class data
    ohia_v_plot<-ggplot(data = vc_table(), aes(x = landcover, y = all_ohia)) + 
      geom_bar(stat = "identity", fill = rev(topo.colors(32)), color = "grey") +
      labs(x = "", y = "count") + ggtitle("Entire Ohia Distribution") +
      coord_flip() + theme_gray()
    
    rod_v_plot<-ggplot(data = vc_table(), aes(x = landcover, y = strain)) + 
      geom_bar(stat = "identity", fill = rev(topo.colors(32)), color = "grey") +
      labs(x = "", y = "count") + ggtitle("ROD Only Distribution") +
      coord_flip() + theme_gray()
    
    # barplot of bioregion data
    ohia_b_plot<-ggplot(data = br_table(), aes(x = bioregion, y = all_ohia)) + 
      geom_bar(stat = "identity", fill = rev(topo.colors(8)), color = "grey") +
      labs(x = "", y = "count") + ggtitle("Entire Ohia Distribution") +
      coord_flip() + theme_gray()
    
    rod_b_plot<-ggplot(data = br_table(), aes(x = bioregion, y = strain)) + 
      geom_bar(stat = "identity", fill = rev(topo.colors(8)), color = "grey") +
      labs(x = "", y = "count") + ggtitle("ROD Only Distribution") +
      ylim(0, 45) + coord_flip() + theme_gray()
    
    # arrange plots based on user ui
    if(input$ecovar == "higap"){
      grid.arrange(veg_plot, arrangeGrob(ohia_v_plot, rod_v_plot, nrow = 1, ncol = 2), 
                   nrow = 2, ncol = 1)
    } else {
      grid.arrange(bioreg_plot, arrangeGrob(ohia_b_plot, rod_b_plot), nrow = 2, ncol = 1)
    }
    
  }, width = 800, height = 800)
  
  #----- GEOLOGY -----#  
  
  # create reactive dataset based on rod strain
  geol_data<-reactive({
    rod[rod$POS_NEG %in% input$type4, ]
  })
  
  # create reactive data table of age classes
  age_table<-reactive({
    age_tab<-data.frame(table(hi_ohia$AGE_GROUP))
    age_tab<-merge(age_tab, data.frame(table(geol_data()["AGE_GROUP"])), by = "Var1", all = T)
    names(age_tab)<-c("age", "all_ohia", "strain")
    age_tab[is.na(age_tab)]<-0
    age_tab
  })
  
  # create reactive data table of names
  type_table<-reactive({
    type_tab<-data.frame(table(hi_ohia$NAME))
    type_tab<-merge(type_tab, data.frame(table(geol_data()["NAME"])), by="Var1", all=T)
    names(type_tab)<-c("type", "all_ohia", "strain")
    type_tab[is.na(type_tab)]<-0
    type_tab
  })
  
  output$geol_plots<-renderPlot({
    
    # plot geologic ages raster    
    age_plot<-ggplot(data = hi_geol, aes(x = x, y = y)) +
      geom_raster(aes(fill = AGE_GROUP)) + 
      scale_fill_gradientn(name = input$geol, colours = rev(cm.colors(11)), guide = F) +
      geom_point(data = geol_data(), 
                 aes(x = j_LON, y = j_LAT, shape = factor(POS_NEG)), 
                 size = 3, fill = "gray") + 
      scale_shape_discrete(solid = F, name = "Strain") +
      labs(x = "Longitude", y = "Latitude") + coord_fixed() + theme_gray() 
    
    # plot geologic names raster 
    type_plot<-ggplot(data = hi_geol, aes(x = x, y = y)) +
      geom_raster(aes(fill = NAME)) + 
      scale_fill_gradientn(name = input$geol, colours = rev(cm.colors(18)), guide = F) +
      geom_point(data = geol_data(), 
                 aes(x = j_LON, y = j_LAT, shape = factor(POS_NEG)), 
                 size = 3, fill = "gray") + 
      scale_shape_discrete(solid = F, name = "Strain") +
      labs(x = "Longitude", y = "Latitude") + coord_fixed() + theme_gray()    
    
    # barplot of age data
    ohia_a_plot<-ggplot(data = age_table(), aes(x = age, y = all_ohia)) + 
      geom_bar(stat = "identity", fill = rev(cm.colors(11)), color = "grey") +
      labs(x = "", y = "count") + ggtitle("Entire Ohia Distribution") +
      coord_flip() + theme_gray()
    
    rod_a_plot<-ggplot(data = age_table(), aes(x = age, y = strain)) + 
      geom_bar(stat = "identity", fill = rev(cm.colors(11)), color = "grey") +
      labs(x = "", y = "count") + ggtitle("ROD Only Distribution") +
      ylim(0, 35) + coord_flip() + theme_gray()
    
    # barplot of lithology data
    ohia_t_plot<-ggplot(data = type_table(), aes(x = type, y = all_ohia)) + 
      geom_bar(stat = "identity", fill = rev(cm.colors(11)), color = "grey") +
      labs(x = "", y = "count") + ggtitle("Entire Ohia Distribution") +
      coord_flip() + theme_gray()
    
    rod_t_plot<-ggplot(data = type_table(), aes(x = type, y = strain)) + 
      geom_bar(stat = "identity", fill = rev(cm.colors(11)), color = "grey") +
      labs(x = "", y = "count") + ggtitle("ROD Only Distribution") +
      ylim(0, 45) + coord_flip() + theme_gray()
    
    # arrange plots based on user ui
    if(input$geol == "geol_age"){
      grid.arrange(age_plot, arrangeGrob(ohia_a_plot, rod_a_plot, nrow = 1, ncol = 2), 
                   nrow = 2, ncol = 1)
    } else {
      grid.arrange(type_plot, arrangeGrob(ohia_t_plot, rod_t_plot), nrow = 2, ncol = 1)
    }
    
  }, width = 800, height = 800)
  
  
  # end shiny SERVER
}

# run shiny web application
shinyApp(ui = shinyUI, server = shinyServer)
# deployApp()