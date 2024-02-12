
#Create environment
rm(list=ls(all=TRUE))
library(raster)
library(caTools)
library(haven)
library(terra)
setwd("~/Downloads/albedo/Data")

#Load in data
summer_albedo <- rast('council_watershed_albedo_summer.tif') 
winter_albedo <- rast('council_watershed_albedo_winter.tif') 
canopy_height_model <- rast('council_watershed_chm_30m.tif') 
aspect <- rast('council_watershed_aspect_30m.tif')
slope <- rast('council_watershed_slope_30m.tif')
topo <- rast('council_watershed_topo_30m.tif') #what is this
tpi <- rast('council_watershed_tpi_30m.tif') #what is this
tri <- rast('council_watershed_tri_30m.tif') #what is this
fcover <- rast('~/Downloads/albedo/Data/council_watershed_fcover_30m.dat')
pft <- rast('~/Downloads/albedo/Data/council_watershed_pft_30m.dat')




#Define the minimum extent to which we should crop all rasters
summer_extent <- ext(summer_albedo)
winter_extent <- ext(winter_albedo)
chm_extent <- ext(canopy_height_model)
pft_extent <- ext(pft)
tri_extent <- ext(tri)
final_extent <- ext(chm_extent[1], winter_extent[2], chm_extent[3], winter_extent[4])

#Crop rasters
ext(summer_albedo) <- ext(final_extent)
ext(winter_albedo) <- ext(final_extent)
ext(fcover) <- ext(final_extent)
ext(pft) <- ext(final_extent)
ext(canopy_height_model) <- ext(final_extent)
ext(topo) <- ext(final_extent)
ext(tpi) <- ext(final_extent)
ext(tri) <- ext(final_extent)
ext(slope) <- ext(final_extent)
ext(aspect) <- ext(final_extent)

vegstack <- c(summer_albedo, winter_albedo, fcover, pft)
envistack <- c(canopy_height_model, topo, tpi, tri, slope, aspect)

vegstack2 <- vegstack
res(vegstack2) <- res(envistack)
vegstack2 <- resample(vegstack, vegstack2)
stack <- c(envistack, vegstack2)


#preliminary plots
plot(stack$council_watershed_fcover_30m_2, stack$council_watershed_albedo_summer) 
plot(stack$council_watershed_fcover_30m_2, stack$council_watershed_albedo_winter)

plot(stack$council_watershed_fcover_30m_3, stack$council_watershed_albedo_summer) 

plot(stack$council_watershed_fcover_30m_4, stack$council_watershed_albedo_summer)

barplot(pft, col = c("black", "red", "darkgreen", "lightgreen", "cyan", "blue", "yellow", "cyan","blue", "green", "pink", "orange", "white", "gray"))

