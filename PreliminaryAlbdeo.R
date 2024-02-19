
#Create environment
rm(list=ls(all=TRUE))
library(raster)
library(caTools)
library(haven)
library(terra)
library(spatialEco)
library(ggplot2)
setwd("~/Downloads/albedo/Data")

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
original_albedo <- read.table('~/Downloads/albedo/Data/council_watershed_albedo_ts_orignal.dat')
smoothed_albedo <- read.table('~/Downloads/albedo/Data/council_watershed_albedo_ts_smoothed.dat', skip = 5)


#Define the minimum extent to which we should crop all rasters
summer_extent <- ext(summer_albedo)
winter_extent <- ext(winter_albedo)
chm_extent <- ext(canopy_height_model)
pft_extent <- ext(pft)
tri_extent <- ext(tri)
final_extent <- ext(chm_extent[1], pft_extent[2], chm_extent[3], tri_extent[4])

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
names(fcover) <- c("Spruce_fcover" , "Alder_fcover" ,"Willow_fcover" , "OtherTall_fcover",  "LowShrubs_fcover" , "DwarfShrubs_fcover" ,"EvergreenShrubs_fcover" , "Forb_fcover" ,"DryGrass_fcover" , "WetGrass_fcover","Moss_fcover" ,"Lichen_fcover","NPV_fcover")
vegstack <- c(summer_albedo, winter_albedo, fcover, pft) #group together all layers with the same res/extent
envistack <- c(canopy_height_model, topo, tpi, tri, slope, aspect) #group together all layers with the same res/extent


envistack2 <- resample(envistack, vegstack) #maybe we need to put everything in the geometry of vegstack, so that fcovers still add up to 1, fcover and albdeo are our "most important" data
stack <- c(envistack2, vegstack)
stack_df <- as.data.frame(stack, xy = TRUE)

#preliminary plots
plot(stack$council_watershed_fcover_30m_2, stack$council_watershed_albedo_summer) 
plot(stack$council_watershed_fcover_30m_2, stack$council_watershed_albedo_winter)

hist(stack_df$council_watershed_chm_30m)
hist(stack_df$council_watershed_chm_30m)

stack_df$council_watershed_chm_30m[stack_df$council_watershed_chm_30m <0] <- NA
canopy_height_model[canopy_height_model < 0] <- NA
plot(canopy_height_model)


plot(stack$council_watershed_tri_30m, stack$council_watershed_tpi_30m)

barplot(pft, col = c("black", "red", "darkgreen", "lightgreen", "cyan", "blue", "yellow", "cyan","blue", "green", "pink", "orange", "white", "gray"))


#Pearson's correlation coefficient over space
r.cor <- rasterCorrelation(stack$council_watershed_albedo_summer, stack$Alder_fcover, type = "pearson")
plot(r.cor) #interesting to see where summer and winter albedo are strongly correlated

tri_tpi <- rasterCorrelation(stack$council_watershed_tri_30m, stack$council_watershed_tpi_30m, type = "pearson")
plot(tri_tpi)

#Linear regression - try with 
test <- lm(council_watershed_albedo_summer ~ council_watershed_chm_30m + council_watershed_topo_30m + council_watershed_tpi_30m + council_watershed_tri_30m + council_watershed_slope_30m + council_watershed_aspect_30m + Spruce_fcover + Alder_fcover + Willow_fcover + OtherTall_fcover + LowShrubs_fcover + EvergreenShrubs_fcover + DryGrass_fcover, data = as.data.frame(stack))#group together all things which might influence surface albedo
test2 <- lm(council_watershed_albedo_summer ~ council_watershed_chm_30m + council_watershed_topo_30m + council_watershed_tpi_30m + council_watershed_tri_30m + council_watershed_slope_30m + council_watershed_aspect_30m, data = as.data.frame(stack))#group together all things which might influence surface albedo besides fcover which has NAs



#DwarfShrubs_fcover, Forb_fcover, WetGrass_fcover, Moss_fcover, Lichen_fcover, NPV_fcoverall NA!
plot(test)
summary(test)

df <- as.data.frame(stack$council_watershed_albedo_summer)
drivers <- as.data.frame(driver_stack)




#try some stuff with ggplot




vegstack_df <- as.data.frame(vegstack, xy = TRUE)
ggplot(data = stack_df, aes(x = council_watershed_albedo_summer, y = Willow_fcover))+
  geom_point(aes(color = category, shape = "."))

ggplot(data = stack_df, aes(x = Willow_fcover, y = council_watershed_albedo_summer))+
  geom_point()+
  geom_smooth(method = "lm", se=FALSE)

ggplot(data = stack_df, aes(x = Alder_fcover, y = council_watershed_albedo_summer))+
  geom_point()+
  geom_smooth(method = "lm", se=FALSE)

ggplot(data = stack_df, aes(x = Willow_fcover, y = council_watershed_albedo_winter))+
  geom_point()+
  geom_smooth(method = "lm", se=FALSE)

ggplot(data = stack_df, aes(x = Alder_fcover, y = council_watershed_albedo_winter))+
  geom_point()+
  geom_smooth(method = "lm", se=FALSE)

stack_df$shrubcover <- stack_df$Alder_fcover + stack_df$Willow_fcover

ggplot(data = stack_df, aes(x = shrubcover, y = council_watershed_albedo_summer))+
  geom_point()+
  geom_smooth(method = "lm", se=FALSE)

ggplot(data = stack_df, aes(x = shrubcover, y = council_watershed_albedo_winter))+
  geom_point()+
  geom_smooth(method = "lm", se=FALSE)


ggplot(data = stack_df, aes(x = category, y = council_watershed_albedo_summer))+
  geom_boxplot()
ggplot(data = stack_df, aes(x = category, y = council_watershed_albedo_winter))+
  geom_boxplot()
ggplot(data = stack_df, aes(x = category, y = council_watershed_chm_30m))+
  geom_boxplot()

ggplot(data = stack_df, aes(x = council_watershed_chm_30m, y = council_watershed_albedo_summer))+
  geom_point()+
  geom_smooth(method = "lm", se=FALSE)

#Making simple models
winter_shrubmodel <- lm(stack_df$council_watershed_albedo_winter ~ stack_df$shrubcover)
summer_shrubmodel <- lm(stack_df$council_watershed_albedo_summer ~ stack_df$shrubcover)
summary(winter_shrubmodel)
summary(summer_shrubmodel)

spruce_model <- lm(stack_df$Spruce_fcover ~ stack_df$council_watershed_chm_30m)
summary(spruce_model) #Can chm predict spruce fcover? It should, because spruce should be the tallest, but there seems to be some issue visualizing chm in R

#these full models have pretty high R2 values
fullmodel_summer <- lm(council_watershed_albedo_summer ~ council_watershed_chm_30m + council_watershed_topo_30m + council_watershed_tpi_30m + council_watershed_tri_30m + council_watershed_slope_30m + council_watershed_aspect_30m + Spruce_fcover + Alder_fcover + Willow_fcover, data = stack_df)
summary(fullmodel_summer) #for some reason, including the less common PFTs causes the model to throw an error
fullmodel_winter <- lm(council_watershed_albedo_winter ~ council_watershed_chm_30m + council_watershed_topo_30m + council_watershed_tpi_30m + council_watershed_tri_30m + council_watershed_slope_30m + council_watershed_aspect_30m + Spruce_fcover + Alder_fcover + Willow_fcover, data = stack_df)
summary(fullmodel_winter) #for some reason, including the less common PFTs causes the model to throw an error










