###########################################################################################
#                    
#
#    --- Last updated:  2020.05.28 By Daryl Yang <dediyang@bnl.gov>
###########################################################################################


#-------------------- Close all devices and delete all variables. ------------------------#
rm(list=ls(all=TRUE))   # clear workspace
graphics.off()          # close any open graphics
closeAllConnections()   # close any open connections to files
dlm <- .Platform$file.sep # <--- What is the platform specific delimiter?
#-----------------------------------------------------------------------------------------#

#---------------------------- Load required libraries ------------------------------------#
# Info: Loads required R libraries and warns if package is not availible.
ok <- require(raster) ; if (! ok) 
  #stop("*** Package pls is not available.  This is needed for model optimization ***")
  install.packages("raster")
ok <- require(readr) ; if (! ok) 
  install.packages("readr")
ok <- require(strucchange ) ; if (! ok) 
  install.packages("strucchange ")

library(ggplot2)
library(raster)
library(strucchange)
library(readr)

# Script options
options(digits.secs = 3)
options(digits = 15)
#-----------------------------------------------------------------------------------------#

#------------------------------ Set output directory -------------------------------------#
### Set ouput directory
out.dir <- '/Users/anncrumlish/Downloads/albedo analysis/step4_smoothed_tsstep4_smoothed_ts/'

### Create output folders
if (! file.exists(out.dir)) dir.create(out.dir,recursive=TRUE)
#-----------------------------------------------------------------------------------------#

#-------------------------------- Load in data--------------------------------------------#
data.dir <- '/Users/anncrumlish/Downloads/albedo analysis/step4_smoothed_ts/SG_Smoothed_council_lidar_area_albedo_ts_2013_2019.dat'
albedo.raster <- stack(data.dir)
#-----------------------------------------------------------------------------------------#

albedo.df <- as.data.frame(albedo.raster) #50 to 296 is the albedo on each individual day of the year for the summer
band.names <- names(albedo.df)
doy <- parse_number(band.names)

num.pixels <- nrow(albedo.df)

break.points.all.1 <- c() #do a different regression/model on different segments of the data, data is segmeted by finding the "breakpoints"
break.points.all.2 <- c()
for (i in 1:num.pixels)
{
  print(paste(i, 'of', num.pixels))
  albedo.pixel <- as.numeric(albedo.df[i,])
  data.combn <- data.frame(cbind(doy, albedo.pixel))
  names(data.combn) <- c('doy', 'albedo')
  data.combn <- na.omit(data.combn)
  
  bp <- breakpoints(albedo ~ doy, data = data.combn, breaks = 2)
  
  
  #---------------------
  # example plot for a single pixel
  data.combn <- data.combn[1:(nrow(data.combn)-1),]
  bp <- breakpoints(albedo ~ doy, data = data.combn, breaks = 2)
  
  ggplot(data = data.combn) +
    geom_point(aes(x = doy, y = albedo), shape = 16, colour = 'grey40', size = 2.5) +
    geom_line(aes(x = doy, y = fitted(bp, breaks = 2)), colour = 'red', size = 2) +
    xlim(c(50, 300)) + ylim(c(0, 1)) + labs(x = 'Day of Year', y = 'White Sky Albedo') +
    theme(axis.text = element_text(size=12, angle = 0), axis.title=element_text(size=12,face="bold")) +
    geom_point(aes(x = 185, y = 0.8), shape = 16, colour = 'grey40', size = 2.5) +
    annotate("text", x = 190 + 5,  y = 0.8, label = 'Savitzky-Golay Filterred', hjust=0) +
    geom_segment(aes(x = 180, xend = 190, y = 0.7, yend = 0.7), color = 'red', size = 1.2) +
    annotate("text", x = 190 + 5,  y = 0.7, label = 'Breakpoint regression', hjust=0)
  
  png.name <- paste0(out.dir, 'Breakpoint_Plot.png')
  ggsave(png.name, plot = last_plot(), width = 16, height = 10, units = 'cm')  
  #---------------------
  
  break.points.1 <- bp$breakpoints[1]
  break.points.2 <- bp$breakpoints[2]
  
  break.doy.1 <- data.combn$doy[break.points.1]
  break.doy.2 <- data.combn$doy[break.points.2]

  break.points.all.1 <- rbind(break.points.all.1, break.doy.1)
  break.points.all.2 <- rbind(break.points.all.2, break.doy.2)
}



break.point.arr <- matrix(break.points.all, nrow = albedo.raster@nrows, ncol = albedo.raster@ncols)

break.point.raster <- raster(break.points.all, template = albedo.raster)

writeRaster(break.point.raster, filename=file.path(out.dir, "break_point_30m.tif"), format="GTiff", overwrite=TRUE)

##Using previously created .tif to avoid running the long loop again
setwd("/Users/anncrumlish/Downloads/albedo analysis/step4_smoothed_ts")
break.point.raster <- raster('break_point_30m.tif')
plot(break.point.raster) #this shows the spatial distribution of doy break points, it is earliest near the bottomlands and later up on the hill slope (willow/alder)
break.point.raster #breakpoint ranges from May 1 to Sep. 16
#assuming that break point indicates some change in phenology, now you should plot breakpoint pixels against chm, or against pft fcover

stack <- rasterToPoints(break.point.raster)






