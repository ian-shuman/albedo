###########################################################################################
#
#        This script loads and cleans data and recreates spatial plots in Figure 1 of Shuman et al.
#
#    --- Last updated:  2025.02.06 By Ian Shuman <ins2109@columbia.edu>
###########################################################################################

#******************** close all devices and delete all variables *************************#
rm(list=ls(all=TRUE))   # clear workspace
graphics.off()          # close any open graphics
closeAllConnections()   # close any open connections to files
dlm <- .Platform$file.sep # <--- What is the platform specific delimiter?
options(warn = -1)
#*****************************************************************************************#

#****************************** load required libraries **********************************#
### install and load required R packages
list.of.packages <- c("raster", "caTools", "haven", 'terra', 'spatialEco', 'ggplo2', 'dplyr', 'cowplot', 'grid', 'tidyverse', 'ggfun', 'scatterpie', 'gtable')
# check for dependencies and install if needed
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, dependencies=c("Depends", "Imports",
                                                                       "LinkingTo"))
version_requirements <- c("3.3.2")
if (!packageVersion("ggplot2") >= version_requirements[1]) {
  remotes::install_version(package="ggplot2", version=paste0(">= ", version_requirements), 
                           dependencies=c("Depends", "Imports", "LinkingTo"), upgrade="ask",
                           quiet=TRUE)
}


# load libraries
invisible(lapply(list.of.packages, library, character.only = TRUE))

#Set wd to data directory
setwd("~/Downloads/albedo/Data")
#*****************************************************************************************#
#*
#*
#Load in data
summer_albedo <- rast('council_watershed_albedo_summer.tif') #1
winter_albedo <- rast('council_watershed_albedo_winter.tif') #1
canopy_height_model <- rast('council_watershed_chm_30m.tif') #4
canopy_height_model[canopy_height_model < 0] <- NA #there are about 6000 negative values which I think should be removed
aspect <- rast('council_watershed_aspect_30m.tif') #2
slope <- rast('council_watershed_slope_30m.tif') #2
topo <- rast('council_watershed_topo_30m.tif') #what is this #2
tpi <- rast('council_watershed_tpi_30m.tif')  #2 #topographic position index, positive values mean the cell is higher than its surrounding cells, negative mean it is lower
tri <- rast('council_watershed_tri_30m.tif') #terrain ruggedness index, higher values mean more heterogeneity in the surrounding cells #2
fcover <- rast('~/Downloads/albedo/Data/council_watershed_fcover_30m.dat') #3
pft <- rast('~/Downloads/albedo/Data/council_watershed_pft_30m.dat') #3
twi <- rast('~/Downloads/albedo/Data/twi.tif') #3
breakpoints1 <- rast("/Users/anncrumlish/Downloads/albedo analysis/step4_smoothed_ts/break_point1_30m.tif")
breakpoints <- rast("/Users/anncrumlish/Downloads/albedo analysis/step4_smoothed_ts/break_point_30m.tif")
breakpoints[breakpoints > 213] <- NA
breakpoints1[breakpoints1 < 92] <- NA


#All data have already been upscaled to a common 30m resolution in IDL. However, there are still
#very small differences in the extent and resolution (at the cm order of magnitide) of the products, 
#so we will resample again using bilinear interpolation. We will use fcover as the template
#because fcover values for all PFTs must necessarily sum to 1, which may not be true if resampled. 
chm.resamp <- resample(canopy_height_model, fcover)
topo.resamp <- resample(topo, fcover)
tpi.resamp <- resample(tpi, fcover)
tri.resamp <- resample(tri, fcover)
slope.resamp <- resample(slope, fcover)
aspect.resamp <- resample(aspect, fcover)
twi.resamp <- resample(twi, fcover)
summer_albedo.resamp <- resample(summer_albedo, fcover)
winter_albedo.resamp <- resample(winter_albedo, fcover)
breakpoints1.resamp <- resample(breakpoints1, fcover)
breakpoints.resamp <- resample(breakpoints, fcover)
stack <- c(chm.resamp, topo.resamp, tpi.resamp, tri.resamp, slope.resamp, aspect.resamp, twi.resamp, summer_albedo.resamp, winter_albedo.resamp, fcover, pft, breakpoints1.resamp, breakpoints.resamp)
stack_df <- as.data.frame(stack, xy = TRUE)
colnames(stack_df) <- c("x", "y" , "council_watershed_chm_30m", "council_watershed_topo_30m", "council_watershed_tpi_30m" , "council_watershed_tri_30m", "council_watershed_slope_30m", "council_watershed_aspect_30m", "twi", "council_watershed_albedo_summer", "council_watershed_albedo_winter", "Spruce_fcover", "Alder_fcover" , "Willow_fcover", "OtherTall_fcover" , "LowShrubs_fcover", "DwarfShrubs_fcover", "EvergreenShrubs_fcover", "Forb_fcover" , "DryGrass_fcover", "WetGrass_fcover", "Moss_fcover", "Lichen_fcover", "NPV_fcover", "category", "break_point1_30m", "break_point_30m")

#Check to ensure that physical laws are not broken after resampling
max(stack_df$council_watershed_albedo_winter, na.rm = T) <= 1 & 
  max(stack_df$council_watershed_albedo_summer, na.rm = T) <= 1 & 
  min(stack_df$council_watershed_albedo_winter, na.rm = T) >= 0 & 
  min(stack_df$council_watershed_albedo_summer, na.rm = T) >= 0 & 
  max(stack_df$council_watershed_slope_30m, na.rm = T) <= 90 &
  min(stack_df$council_watershed_slope_30m, na.rm = T) >= 0 &
  max(stack_df$council_watershed_aspect_30m, na.rm = T) <= 360 &
  min(stack_df$council_watershed_aspect_30m, na.rm = T) >= 0 &
  min(stack_df$council_watershed_chm_30m, na.rm = T) >= 0

save(stack_df, file = "stack_df.RData")


#Filter Categories for >75% fCover, use only for figures 2 and 3
stack_df2 <- stack_df[,12:24]
stack_df2[is.na(stack_df2)] <- 0
# Create a new column "HighCoverSpecies" which is a list of column names where the value is greater than 0.75
stack_df2$category2 <- apply(stack_df2, 1, function(row) {
  # Identify columns with values greater than 0.75
  category2 <- names(stack_df2)[row > 0.75]
  
  # Return the column names as a list
  if (length(category2) == 0) {
    return(NA)  # Return NA if no species have a cover greater than 0.75
  } else {
    return(list(category2))
  }
})

View(stack_df2)
stack_df2 <- cbind(stack_df[,1:11], stack_df2[,1:14], stack_df[,25:27])


stack_df2$category3 <- sapply(stack_df2$category2, function(x) {
  if (length(x) == 1) {
    return(unlist(x))
  } else {
    return(NA) # Placeholder for lists with more than one element
  }
})


stack_df2 <- stack_df2[!is.na(stack_df2$category3), ]
stack_df2 <- stack_df2[, !names(stack_df2) %in% c("category2", "category3")]
save(stack_df2, file = "stack_df2.RData")

#Make cursory plots to show distributions of PFTs and winter albedo
plot(pft)

plot(winter_albedo.resamp, col = gray.colors(255, start = 0, end = 1))

