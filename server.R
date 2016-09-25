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
  
  # calculate density of selected bioclim - DO NOT NEED ANYMORE 

  
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
      ggtitle("Entire Ohia Distribution") +
      labs(x = NULL, y = "Frequency") + theme_gray()   
    
    rod_plot<-ggplot(data = bc_data()) + 
      geom_histogram(na.rm = T, aes(x = bc_data()[input$bioclim]), bins = ohia_bins(),
                     fill = rev(terrain.colors(ohia_bins())), colour = "grey") +  
      ggtitle("ROD Only Distribution") + ylim(0, 30) +
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
      scale_fill_gradientn(name = input$ecovar, colours = rev(topo.colors(33)), guide = F) +
      geom_point(data = vc_data(), 
                 aes(x = j_LON, y = j_LAT, shape = factor(POS_NEG)), 
                 size = 3, fill = "gray") + 
      scale_shape_discrete(solid = F, name = "Strain") +
      labs(x = "Longitude", y = "Latitude") + coord_fixed() + theme_gray() 
    
    # plot bioregions raster 
    bioreg_plot<-ggplot(data = hi_bioreg, aes(x = x, y = y)) +
      geom_raster(aes(fill = bioregions)) + 
      scale_fill_gradientn(name = input$ecovar, colours = rev(topo.colors(8)), guide = F) +
      geom_point(data = vc_data(), 
                 aes(x = j_LON, y = j_LAT, shape = factor(POS_NEG)), 
                 size = 3, fill = "gray") + 
      scale_shape_discrete(solid = F, name = "Strain") +
      labs(x = "Longitude", y = "Latitude") + coord_fixed() + theme_gray()    
    
    # barplot of vegetation class data
    ohia_v_plot<-ggplot(data = vc_table(), aes(x = landcover, y = all_ohia)) + 
      geom_bar(stat = "identity", fill = rev(topo.colors(33)), color = "grey") +
      labs(x = "", y = "count") + ggtitle("Entire Ohia Distribution") +
      coord_flip() + theme_gray()
    
    rod_v_plot<-ggplot(data = vc_table(), aes(x = landcover, y = strain)) + 
      geom_bar(stat = "identity", fill = rev(topo.colors(33)), color = "grey") +
      labs(x = "", y = "count") + ggtitle("ROD Only Distribution") +
      ylim(0, 30) + coord_flip() + theme_gray()
    
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
    
    #ggplot(t2, aes(x = Var1, y = Freq.y)) + 
    #geom_bar(na.rm = F, stat = "identity", fill = topo.colors(36), color = "grey") +
    #labs(x = "", y = "Count") + ggtitle("Entire Ohia Distribution") +
    #coord_flip() + theme_gray()
    
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