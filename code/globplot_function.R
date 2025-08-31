##############################################################
# Script to create spatial global plots  #####################
##############################################################

# Function takes in dataframe and creates a global plot who's breaks and 
# palettes depend on the variable given in input. Plots generates for a given 
# simulation and oxygen concentration.


globalplot <- function(mdf, diff = FALSE,sim="",plotlab="",o2="",var=var){
  
  # If a difference plot is defined as input:
  if(diff){
    pal = c("blue","white","red")
    lims = c(-1,1)
    
  }else{
  # list palettes that could be used 
  pal = c("white","gainsboro","darkslategray1","deepskyblue1","green2","yellow1","orange1","orangered2", "firebrick4", "black")}
  
  # set limits, breaks, palettes and units for different variables
  if(var == "mnfire"){
    breaks = TRUE
    mdf$breaks = cut(mdf$z.value,breaks= c(-1,0,1,2,4,6,8,10,Inf), 
                     labels=c("0","0-1", "1-2", "2-4", "4-6", "6-8", "8-10", "10+"))
    unit = "n fires"
    lims = c(-10,10)
    
  }else if (var == "burnedf"){
    breaks = TRUE
    mdf$breaks = cut(mdf$z.value,breaks= c(-1,0,0.0001,0.002,0.004,0.006,0.008,2), 
                     labels=c("0","0-10000000", "0.004-0.04", "0.01-0.07", "0.07-0.15", "0.15-0.3", "0.3-0.6"))
    leglab = c("0","0-0.01", "0.01-0.2", "0.2-0.4", "0.4-0.6", "0.6-0.8", "0.8+")
    unit = "burned area \n(% area)"
    lims = c(0,1)
    
  }else if (var == "tmp"){
      breaks = FALSE  # continuous scale
      unit = "Temperature (°C)"
      lims = c(-30, 40)
      pal = colorRampPalette(c("darkblue", "deepskyblue", "white", "orange", "red", "darkred"))(100)
      leglab = NULL  # no discrete legend labels for continuous scale
    
  }else if(var == "NPP"){
    breaks = TRUE
    mdf$breaks = mdf$z.value
    unit = "gm-2"
    if(diff){lims=c(-3000,3000)}else{lims = c(0,3000)
    pal = brewer.pal(6, "Blues")}
    mdf$breaks = cut(mdf$z.value,breaks= c(0,500,1000,1500,2000,2500,3000))
    
  }else if (var == "pftalbiomass"){
    breaks = TRUE
    # scale for aboveground biomass:
    brk=c(-1,0.01,0.2,0.5,1,1.5,2,2.5,3,6,9,12,15,18,21,24,27,Inf)
    if(diff){lims=c(-27,27)}else{lims = c(0,27)
    # colour scale for aboveground biomass:
    pal=c("white","#E6E9E2","#B2B6B2","#6E4A00","#A16706","#CB8200","#E0AB04","#FAD103","#FFF0A4",
          "#D6ED00","#A2DE03","#7FC108","#569A06","#308100","#006400","#014909","#080A04")}
    unit = "Aboveground biomass (kgm-2)"
    leglab = c("0-0.01","0.01-0.2","0.2-0.5","0.5-1","1-1.5","1.5-2","2-2.5","2.5-3","3-6","6-9","9-12",
               "12-15","15-18","18-21","21-24", "24-27", "27+")
    # make breaks
    mdf$z.value = mdf$z.value/1000
    mdf$breaks = cut(mdf$z.value,breaks= brk)
    
  }else if (var == "cover"){
    breaks = TRUE
    # Scale for fractional coverage
    bscale=c(-1,0.001,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,Inf)
    brk = c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1)
    lims = c(0,1)
    labs = c("0","10","20","30","40","50","60","70","80","90","100")
    pal=c("white","#B2B6B2","#A16706","#E0AB04","#FAD103","#FFF0A4","#D6ED00","#7FC108","#308100","#014909","#080A04")

    # make breaks
    mdf$breaks = cut(mdf$z.value,breaks= bscale)
    unit = "frac"
    
  }else if (var == "forestcov"){
    breaks = TRUE
    unit=""
    mdf$breaks = cut(mdf$z.value, breaks = c(0.6,1))
    leglab = c("forest cover")
    if(diff){lims=c(-1,1)}else{
    pal = c("forestgreen")}
    lims=c(0.6,1)
    
  }else if (var == "treecov"){
    breaks = TRUE
    unit= "Tree Cover (%)"
    mdf$breaks = cut(mdf$z.value, breaks = c(-1,0.1,0.2,0.3,0.4,0.5,0.6,1))
    leglab = c("0-10","10-20","20-30","30-40","40-50","50-60",">60 (forest)")
    pal = rev(terrain.colors(7))
    
  }else if (var == "height"){
    breaks = TRUE
    unit= "Tree height (m)"
    mdf$breaks = cut(mdf$z.value, breaks = c(-1,2,4,6,8,10,12,14,16,18,20,Inf))
    leglab = c("0-2","2-4","4-6","6-8","8-10","10-12","12-14","14-16","16-18","18-20","20+")
    pal = hsv(1, seq(0,1,length.out = 11) , 1)
  }
  
  # Create difference plot:
   if(diff){
    
     globplot = ggplot(mdf,aes_(mdf$lon.value, mdf$lat.value, fill=mdf$z.value))+geom_raster(interpolate = TRUE)+
       coord_equal(ylim = c())+
        scale_fill_gradientn(colours = c("blue","white","red"), name=paste(unit),limits = lims)+
       theme(  
         legend.key.size = unit(0.8, "cm"),
         legend.key = element_rect(colour = '#bdbdbd', size = 0.6),
         legend.title = element_text(size = 15, face="bold"),
         legend.box.margin = margin(0.1,0.1,0.1,0.1),
         legend.text = element_text(size=12),
         legend.box.background = element_rect(colour = "grey85", size = 1.2),
         panel.border = element_blank(),
         panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(),
         panel.background = element_blank(),
         axis.line = element_blank(),
         axis.text = element_blank(),
         axis.ticks = element_blank(),
         plot.margin = unit(c(0.1,0.1,0.1,0.1), "cm"))+
       xlab("")+ylab("")
  
  # Create global values plot:  
  }else{
  globplot = ggplot(mdf,aes_(mdf$lon.value, mdf$lat.value, fill=mdf$breaks))+geom_raster(interpolate = FALSE)+
    guides(fill = guide_legend(nrow = 1, label.position = "right"))+
    ggtitle(paste(plotlab))+
    labs(title = paste(o2, "% oxygen with " , sim," effects", sep=""),
         tag = paste(plotlab))+
    coord_equal(ylim = c())+
    scale_fill_manual(values=pal,labels = leglab,name=paste(unit),na.value = "white")+
    theme(  
      legend.key.size = unit(0.4, "cm"),
      legend.key = element_rect(colour = '#bdbdbd', size = 0.1),
      legend.title = element_text(size = 7),
      legend.position = "bottom",
      legend.box.margin = margin(0.1,0.1,0.05,0.1, unit="cm"),
      legend.text = element_text(size=7),
      legend.box.background = element_rect(colour = "grey33", size = 0.1),
      legend.margin = margin(0.1,0.1,0.05,0.1,unit = "cm"),
      panel.border = element_rect(fill=NA, colour="black", size=0.1),
      panel.background = element_rect(fill = "slategray1"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      plot.title = element_text(size=7,hjust = 0, vjust = -0.5,margin = margin(0,0,2,0)),
      plot.subtitle = element_text(size=6),
      plot.tag = element_text(size = 9, face="bold",hjust = -0.5, vjust = -1.),
      plot.tag.position = c(0, 1)
     )+
    xlab("")+ylab("")
  }
  
  # Return plot
  return(globplot)
  
}
