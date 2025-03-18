###########################################################################################
#
#        This script loads and cleans data and recreates spatial plots in Figure 1 of Shuman et al.
#
#    --- Last updated:  2025.02.26 By Ian Shuman <ins2109@columbia.edu>
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
list.of.packages <- c("raster", "caTools", "haven", 'terra', 'spatialEco', 'ggplot2', 'dplyr', 'readr', 'cowplot', 'grid', 'tidyverse', 'ggfun', 'scatterpie', 'gtable', 'sf', 'rnaturalearth', 'rnaturalearthdata')
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


##All data have already been upscaled to a common 30m resolution in IDL. However, there are still
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
#save(stack_df, file = "stack_df.RData")

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

#View(stack_df2)
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
#save(stack_df2, file = "stack_df2.RData")

library(dplyr)
library(ggplot2)

# Sample data frame
pft_df <- as.data.frame(pft, xy = TRUE)

# Mutate the pft column and set the levels
pft_df <- pft_df %>%
  mutate(pft = case_when(
    category == "nonClass" ~ "Non-Vegetated Surface (NVS)",
    category == "Spruce" ~ "Evergreen Tree (ET)",
    category == "Alder" ~ "Decid. Tall Shrub- Alder (DTSA)",
    category == "Willow" ~ "Decid. Tall Shrub- Willow (DTSW)",
    category == "OtherTall" ~ "Decid. Tree (DT)",
    category == "LowS" ~ "Decid. Low Shrub (DLS)",
    category == "EvergS" ~ "Evergreen Shrub (ES)",
    category == "Forb" ~ "Forb (FO)",
    category == "DryG" ~ "Dry Graminoid (DG)",
    category == "WetG" ~ "Wet Graminoid (WG)",
    category == "Moss" ~ "Moss (MO)",
    category == "Lichen" ~ "Lichen (LI)",
    category == "NPV" ~ "Non-Photosynthetic Veg. (NPV)",
    TRUE ~ NA_character_
  )) %>%
  mutate(pft = factor(pft, levels = c("Evergreen Tree (ET)", "Decid. Tall Shrub- Alder (DTSA)", "Decid. Tall Shrub- Willow (DTSW)", 
                                      "Decid. Tree (DT)", "Decid. Low Shrub (DLS)", "Dry Graminoid (DG)", "Wet Graminoid (WG)", 
                                      "Moss (MO)", "Lichen (LI)", "Evergreen Shrub (ES)", "Forb (FO)", 
                                      "Non-Photosynthetic Veg. (NPV)", "Non-Vegetated Surface (NVS)")))

# Define colors for each level
colors <- c("#FF0000", "#008000", "#00CD00", "#669999", "#0066FF", "#00FF67", "#FF67FF", "#CC6600", "#FFFFFF", "#91FAFD", "#4FADEA", "#5A5A5A", "black")


# Create a sample plot
ggplot(pft_df, aes(x = x, y = y, fill = pft)) +
  geom_tile() +
  coord_equal() +
  scale_fill_manual(values = colors) +
  labs(fill = "Dominant PFT") +
  theme_minimal() +
  theme(
    legend.title = element_text(size = 28),  # Increase legend title text size
    legend.text = element_text(size = 20), 
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    plot.background = element_blank(), 
    legend.key.height = unit(0.8, "cm")
  )


wint_alb_df <- as.data.frame(winter_albedo.resamp, xy = TRUE)

ggplot(wint_alb_df, aes(x = x, y = y, fill = council_watershed_albedo_winter)) +
  geom_tile() +
  coord_equal() +
  scale_fill_gradient(low = "black", high = "white", 
                      breaks = seq(0, 1, by = 0.1)) +
  labs(fill = "Winter Albedo") +
  theme_minimal() +
  theme(
    legend.title = element_text(size = 28),  
    legend.text = element_text(size = 26), 
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    plot.background = element_blank(), 
    legend.key.height = unit(2, "cm"), 
    legend.key.width = unit(1.5, "cm")
  )

#Reproject pft for plotting 
pft_reprojected <- project(pft, 'EPSG:4326')
pft_reprojected_df <- as.data.frame(pft_reprojected, xy = T)

# Map of Alaska and major rivers
world <- ne_states(country = "United States of America", returnclass = "sf")
alaska <- world[world$name == "Alaska", ]
rivers <- ne_download(scale = 10, type = 'rivers_lake_centerlines', category = 'physical', returnclass = "sf")
rivers_c <- st_intersection(rivers, alaska)

council_extent <- ext(pft_reprojected)
center_x <- (council_extent[1] + council_extent[2]) / 2
center_y <- (council_extent[3] + council_extent[4]) / 2
center_point <- st_as_sf(data.frame(x = center_x, y = center_y), coords = c("x", "y"), crs = st_crs(pft_reprojected))


#Plot
ggplot() +
  ggplot2::geom_sf(data = alaska, color = 'black', linewidth = 1) +
  geom_tile(data = pft_reprojected_df, aes(x = x, y = y))+ 
  geom_sf(data = rivers_c, fill = "lightblue", color = "blue", alpha = 0.5) +
  geom_sf(data = center_point, shape = 24, color = "red", fill = "red", size = 5) +
  geom_sf(data = center_point, shape = 25, color = "red", fill = "red", size = 5) +
  theme_minimal() +
  ggplot2::coord_sf(crs = 'EPSG:4326', ylim = c(52, 71), xlim = c(-170, -131))+
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    plot.background = element_blank()
  )








