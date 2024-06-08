
#Create environment
rm(list=ls(all=TRUE))
library(raster)
library(caTools)
library(haven)
library(terra)
library(spatialEco)
library(ggplot2)
library(dplyr)
library(cowplot)
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
twi <- rast('~/Downloads/albedo/Data/twi.tif') #3
breakpoints1 <- rast("/Users/anncrumlish/Downloads/albedo analysis/step4_smoothed_ts/break_point1_30m.tif")
breakpoints <- rast("/Users/anncrumlish/Downloads/albedo analysis/step4_smoothed_ts/break_point_30m.tif")
#original_albedo <- read.table('~/Downloads/albedo/Data/council_watershed_albedo_ts_orignal.dat')
#smoothed_albedo <- read.table('~/Downloads/albedo/Data/council_watershed_albedo_ts_smoothed.dat', skip = 5)


#Define the minimum extent to which we should crop all rasters
min_extent <- summer_albedo
min_extent <- crop(min_extent, ext(winter_albedo))
min_extent <- crop(min_extent, ext(canopy_height_model))
min_extent <- crop(min_extent, ext(pft))
min_extent <- crop(min_extent, ext(tri))
min_extent <- crop(min_extent, ext(breakpoints))

# Crop the rasters to the minimum extent
summer_albedo.c <- terra::crop(summer_albedo, min_extent)
winter_albedo.c <- terra::crop(winter_albedo, min_extent)
canopy_height_model.c <- terra::crop(canopy_height_model, min_extent)
aspect.c <- terra::crop(aspect, min_extent)
slope.c <- crop(slope, min_extent)
topo.c <- crop(topo, min_extent)
tpi.c <- crop(tpi, min_extent)
tri.c <- crop(tri, min_extent)
fcover.c <- crop(fcover, min_extent)
pft.c <- crop(pft, min_extent)
twi.c <- crop(twi, min_extent)
breakpoints1.c <- crop(breakpoints1, min_extent)
breakpoints.c <- crop(breakpoints, min_extent)

names(fcover.c) <- c("Spruce_fcover" , "Alder_fcover" ,"Willow_fcover" , "OtherTall_fcover",  "LowShrubs_fcover" , "DwarfShrubs_fcover" ,"EvergreenShrubs_fcover" , "Forb_fcover" ,"DryGrass_fcover" , "WetGrass_fcover","Moss_fcover" ,"Lichen_fcover","NPV_fcover")
vegstack <- c(summer_albedo.c, winter_albedo.c, fcover.c, pft.c) #group together all layers with the same res/extent
envistack <- c(canopy_height_model.c, topo.c, tpi.c, tri.c, slope.c, aspect.c, twi.c) #group together all layers with the same res/extent
breakpoints.c <- breakpoints.c #placeholder to keep track, we need to do the same procedure on breakpoints as we do on envistack

envistack2 <- resample(envistack, vegstack) #maybe we need to put everything in the geometry of vegstack, so that fcovers still add up to 1, fcover and albdeo are our "most important" data
breakpoints2 <- resample(mean(breakpoints), vegstack, method = "near")
breakpoints1 <- resample(mean(breakpoints1), vegstack)
stack <- c(envistack2, vegstack, breakpoints1, breakpoints2)
stack_df <- as.data.frame(stack, xy = TRUE)
colnames(stack_df) <- c("x", "y" , "council_watershed_chm_30m", "council_watershed_topo_30m", "council_watershed_tpi_30m" , "council_watershed_tri_30m", "council_watershed_slope_30m", "council_watershed_aspect_30m", "twi", "council_watershed_albedo_summer", "council_watershed_albedo_winter", "Spruce_fcover", "Alder_fcover" , "Willow_fcover", "OtherTall_fcover" , "LowShrubs_fcover", "DwarfShrubs_fcover", "EvergreenShrubs_fcover", "Forb_fcover" , "DryGrass_fcover", "WetGrass_fcover", "Moss_fcover", "Lichen_fcover", "NPV_fcover", "category", "break_point1_30m", "break_point_30m")



#preliminary plots
plot(stack$council_watershed_fcover_30m_2, stack$council_watershed_albedo_summer, color = ) 
plot(stack$council_watershed_fcover_30m_2, stack$council_watershed_albedo_winter)

hist(stack_df$council_watershed_chm_30m)
hist(stack_df$council_watershed_chm_30m)

stack_df$council_watershed_chm_30m[stack_df$council_watershed_chm_30m <0] <- NA
canopy_height_model[canopy_height_model < 0] <- NA
plot(canopy_height_model)


plot(stack$council_watershed_tri_30m, stack$council_watershed_tpi_30m)


barplot(pft, col = c("black", "red", "darkgreen", "lightgreen", "cyan", "blue", "yellow", "cyan","blue", "green", "pink", "orange", "white", "gray"))
barplot(stack_df$category, col = c("black", "red", "darkgreen", "lightgreen", "cyan", "blue", "yellow", "cyan","blue", "green", "pink", "orange", "white", "gray"))
barplot(table(stack_df$category),
        ylab = "Frequency",
        xlab = "Category",
        col = c("black", "red", "darkgreen", "lightgreen", "cyan", "blue", "yellow", "cyan","blue", "green", "pink", "orange", "white", "gray"))


#Pearson's correlation coefficient over space
r.cor <- rasterCorrelation(stack$council_watershed_albedo_summer, stack$Alder_fcover, type = "pearson")
plot(r.cor) #interesting to see where summer and winter albedo are strongly correlated

tri_tpi <- rasterCorrelation(stack$council_watershed_tri_30m, stack$council_watershed_tpi_30m, type = "pearson")
plot(tri_tpi)

#Linear regression - try with 
test <- lm(council_watershed_albedo_summer ~ council_watershed_chm_30m + council_watershed_topo_30m + council_watershed_tpi_30m + council_watershed_tri_30m + council_watershed_slope_30m + council_watershed_aspect_30m + Spruce_fcover + Alder_fcover + Willow_fcover + OtherTall_fcover + LowShrubs_fcover + EvergreenShrubs_fcover + DryGrass_fcover, data = stack_df, na.action = na.exclude)#group together all things which might influence surface albedo
test2 <- lm(council_watershed_albedo_summer ~ council_watershed_chm_30m + council_watershed_topo_30m + council_watershed_tpi_30m + council_watershed_tri_30m + council_watershed_slope_30m + council_watershed_aspect_30m, data = as.data.frame(stack))#group together all things which might influence surface albedo besides fcover which has NAs

#Using only the topography variables - R2 = 0.1521
noveg <- lm(council_watershed_albedo_summer ~ council_watershed_topo_30m + council_watershed_tpi_30m + council_watershed_tri_30m + council_watershed_slope_30m + council_watershed_aspect_30m, data = stack_df)#group together all things which might influence surface albedo besides fcover which has NAs

#Using only CHM and topography variables - R2 = 0.2192 
only.structure <- lm(council_watershed_albedo_summer ~ council_watershed_chm_30m + council_watershed_topo_30m + council_watershed_tpi_30m + council_watershed_tri_30m + council_watershed_slope_30m + council_watershed_aspect_30m, data = stack_df)#group together all things which might influence surface albedo besides fcover which has NAs

#Using only fcover of Spruce, Alder, and Willow and topography variables - R2 = 0.7227
only.SAWcomp <- lm(council_watershed_albedo_summer ~ council_watershed_topo_30m + council_watershed_tpi_30m + council_watershed_tri_30m + council_watershed_slope_30m + council_watershed_aspect_30m + Spruce_fcover + Alder_fcover + Willow_fcover, data = stack_df)#group together all things which might influence surface albedo besides fcover which has NAs

#Using fcover of Spruce, Alder, and Willow, CHM, and topography variables - R2 = 0.8996
structure.SAWcomp <- lm(council_watershed_albedo_summer ~ council_watershed_chm_30m + council_watershed_topo_30m + council_watershed_tpi_30m + council_watershed_tri_30m + council_watershed_slope_30m + council_watershed_aspect_30m + Spruce_fcover + Alder_fcover + Willow_fcover, data = stack_df)#group together all things which might influence surface albedo besides fcover which has NAs


#DwarfShrubs_fcover, Forb_fcover, WetGrass_fcover, Moss_fcover, Lichen_fcover, NPV_fcoverall NA!
plot(test)
summary(test)

df <- as.data.frame(stack$council_watershed_albedo_summer)
drivers <- as.data.frame(driver_stack)




#try some stuff with ggplot




vegstack_df <- as.data.frame(vegstack, xy = TRUE)
ggplot(data = stack_df, aes(x = council_watershed_albedo_summer, y = Willow_fcover))+
  geom_point(aes(color = category, shape = "."))

ggplot(data = stack_df, aes(x = Willow_fcover, y = council_watershed_albedo_summer, color = category))+
  geom_point()+
  #geom_smooth(method = "lm", se=FALSE)+ 
  scale_color_manual(values = c('#000000', '#FF0000', '#008000', '#00CD00', '#669999', '#0066FF', '#FFFF66', '#66FFFF', '#00B0F0', '#00FF67', '#FF67FF', '#CC6600', '#FFFFFF', '#5A5A5A'))+
  theme_classic()

ggplot(data = stack_df, aes(x = Alder_fcover, y = council_watershed_albedo_summer))+
  geom_point()+
  geom_smooth(method = "lm", se=FALSE)

ggplot(data = stack_df, aes(x = Willow_fcover, y = council_watershed_albedo_winter))+
  geom_point()+
  geom_smooth(method = "lm", se=FALSE)

ggplot(data = stack_df, aes(x = Alder_fcover, y = council_watershed_albedo_winter))+
  geom_point()+
  geom_smooth(method = "lm", se=FALSE)

ggplot(data = stack_df, aes(x = Alder_fcover, y = Spruce_fcover))+
  geom_point()+
  geom_smooth(method = "lm", se=FALSE)

stack_df$shrubcover <- stack_df$Alder_fcover + stack_df$Willow_fcover

ggplot(data = stack_df, aes(x = shrubcover, y = council_watershed_albedo_summer))+
  geom_point()+
  geom_smooth(method = "lm", se=FALSE)

ggplot(data = stack_df, aes(x = shrubcover, y = council_watershed_albedo_winter))+
  geom_point()+
  geom_smooth(method = "lm", se=FALSE)
ggplot(data = stack_df, aes(x = shrubcover, y = council_watershed_chm_30m))+
  geom_point()+
  geom_smooth(method = "lm", se=FALSE)
ggplot(data = stack_df, aes(x = council_watershed_chm_30m, y = break_point_30m))+
  geom_point()+
  geom_smooth(method = "lm", se=FALSE)

ggplot(data = stack_df, aes(x = category, y = council_watershed_albedo_summer))+
  geom_boxplot()
ggplot(data = stack_df, aes(x = category, y = council_watershed_albedo_winter))+
  geom_boxplot()
ggplot(data = stack_df, aes(x = category, y = council_watershed_chm_30m))+
  geom_boxplot()
ggplot(data = stack_df, aes(x = category, y = break_point_30m))+
  geom_violin()

ggplot(data = stack_df, aes(x = council_watershed_chm_30m, y = council_watershed_albedo_summer, color = category, shape = "."))+
  geom_point()+
  scale_color_manual(values = c("black", "red", "darkgreen", "lightgreen", "cyan", "blue", "yellow", "cyan","blue", "green", "pink", "orange", "white", "gray"))
ggplot(data = stack_df, aes(x = council_watershed_chm_30m, y = council_watershed_albedo_summer, color = category))+
  geom_point()+
  scale_color_manual(values = c("white", "red", "red", "red", "red", "blue", "blue", "blue", "blue", "blue", "blue", "blue", "blue"))
ggplot(data = stack_df, aes(x = council_watershed_chm_30m, y = council_watershed_albedo_summer, color = category))+
  geom_point()+
  scale_color_manual(values = colors)

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


### Try Random Forest, varpart in vegan
library(randomForest)
filtered_stack <- stack_df[,3:23]
filtered_stack <- filtered_stack[,c(1:7, 9:21)]
filtered_stack <- filtered_stack[-is.na(filtered_stack$council_watershed_albedo_summer),]
RF <- randomForest(council_watershed_albedo_summer ~ ., data = filtered_stack)

#Try varpart
library(vegan)
varmodel <- varpart(y = council_watershed_albedo_summer, X = council_watershed_chm_30m + council_watershed_slope_30m + Alder_fcover + Willow_fcover, data = stack_df)

#Try PCA
filtered_stack <- stack_df[275:37788,3:8]
pca <- prcomp(filtered_stack)




##breakpoint stuff


breakpoints2 <- crop(breakpoints, canopy_height_model)
extent(fcover) <- extent(canopy_height_model)
extent(breakpoints2) <- extent(canopy_height_model)
stack <- stack(breakpoints2, pft, canopy_height_model)
stack.df <- as.data.frame(stack)


ggplot(data = stack.df, aes(x = council_watershed_chm_30m, y = break_point_30m))+
  geom_point()+
  geom_smooth(method = "lm", se=FALSE)




#plot progression of breakpoints as chm increases

stack_alder <- stack_df |> subset(stack_df$Alder_fcover > 0.75)
stack_willow <- stack_df |> subset(stack_df$Willow_fcover > 0.75)

boxplot(break_point_30m ~ levels,data=stack_alder |> mutate(levels = cut(council_watershed_chm_30m,10)),horizontal=FALSE,las=2,cex.axis=0.6, xlab = "Canopy Height Bins")
title("Changes in breakpoint with changing canopy height for the alder PFT")

ggplot(data = stack_alder, aes(x=cut(council_watershed_chm_30m,seq(1, 12, 1)), y=break_point_30m)) +
  geom_boxplot(outlier.shape = NA) +
  geom_hline(yintercept = nonTSK, linetype="dashed") +
  geom_segment(aes(x = 7, xend = 8, y = 0.16, yend = 0.16),  size = 1, 
               linetype = "dashed") +
  annotate("text", x = 8.2,  y = 0.16, 
           label = stringr::str_wrap("Non-Woody Albedo", width = 10), hjust=0) +
  xlab("PFT") + ylab("Winter Albedo") + 
  scale_fill_manual(values = cols, name = "", labels = pfts[-c(6:8)]) +
  scale_x_discrete(labels = pfts[-c(6:8)]) +
  guides(fill=guide_legend(ncol = 2)) +
  theme(legend.position = c(0.8, 0.35), 
        legend.title = element_blank(), 
        legend.key.size = unit(0.75, 'cm'),
        legend.text = element_text(size = 12)) +
  theme(axis.text = element_text(size=13, color = 'black'),
        axis.title = element_text(size=13)) +
  theme(axis.line = element_line(colour = "black"),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())


boxplot(break_point_30m ~ levels,data=stack_df |> mutate(levels = cut(council_watershed_chm_30m,seq(0, 12, 0.5))),horizontal=FALSE,las=2,cex.axis=0.6, xlab = "Canopy Height Bins")
title("Changes in breakpoint with changing canopy height for all PFTs")

All_box <- stack_df %>%
  mutate(levels = cut(council_watershed_chm_30m, seq(0, 12, .5))) %>%
  drop_na(levels, break_point_30m) %>%
  ggplot(aes(x = levels, y = break_point_30m)) +
  geom_boxplot(fill = "gray") +
  scale_y_continuous(limits = c(120, 270), breaks = seq(120, 270, by = 10))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  labs(x = "Binned Canopy Height (m)", y = "Completion of Snowmelt (DOY)")+
  theme_classic() +
  theme(axis.text = element_text(size=21, color = 'black'),
        axis.title=element_text(size=25),
        axis.text.x = element_text(angle = 45, hjust = 1))
All_box_start <- stack_df %>%
  mutate(levels = cut(council_watershed_chm_30m, seq(0, 12, 1))) %>%
  drop_na(levels, break_point_30m) %>%
  ggplot(aes(x = levels, y = break_point1_30m)) +
  geom_boxplot(fill = "gray") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  labs(x = "Binned Canopy Height (m)", y = "Start of Snowmelt (DOY)")+
  theme_classic() +
  theme(axis.text = element_text(size=21, color = 'black'),
        axis.title=element_text(size=25),
        axis.text.x = element_text(angle = 45, hjust = 1))

stack_nonwoody <- stack_df |> filter((category != "Alder" & category != "Willow" & category != "Spruce" & category != "OtherTall" & category != "NPV" & category != "nonClass"))
nonwoody_box <- stack_nonwoody %>%
     mutate(levels = cut(council_watershed_chm_30m, seq(0, 12, 0.5))) %>%
     drop_na(levels, break_point_30m) %>%
     ggplot(aes(x = levels, y = break_point_30m)) +
     geom_boxplot(fill = 'gray') +
     theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
     labs(x = "Binned Canopy Height (m)", y = "Completion of Snowmelt (DOY)", title = "Completion of Snowmelt v.s. Canopy Height for Non-Woody PFTs")+
     theme_classic()
stack_alder <- stack_df |> subset(stack_df$Alder_fcover > 0.75)
stack_willow <- stack_df |> subset(stack_df$Willow_fcover > 0.75)
stack_ET <- stack_df |> subset(stack_df$Spruce_fcover > 0.75)
stack_DT <- stack_df |> subset(stack_df$OtherTall_fcover > 0.75)
DTSA_box <- stack_alder %>%
  mutate(levels = cut(council_watershed_chm_30m, seq(0, 12, 0.5))) %>%
  drop_na(levels, break_point_30m) %>%
  ggplot(aes(x = levels, y = break_point_30m)) +
  geom_boxplot(fill = '#008000') +
  scale_y_continuous(limits = c(120, 180))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  labs(x = "Binned Canopy Height (m)", y = "Completion of Snowmelt (DOY)", title = "Completion of Snowmelt v.s. Canopy Height for Alder Dominated Sites")+
  theme_classic()
DTSA_box_start <- stack_alder %>%
  mutate(levels = cut(council_watershed_chm_30m, seq(0, 12, 0.5))) %>%
  drop_na(levels, break_point_30m) %>%
  ggplot(aes(x = levels, y = break_point1_30m)) +
  geom_boxplot(fill = '#008000') +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  labs(x = "Binned Canopy Height (m)", y = "Start of Snowmelt (DOY)", title = "Completion of Snowmelt v.s. Canopy Height for Alder Dominated Sites")+
  theme_classic()
DTSW_box <- stack_willow %>%
  mutate(levels = cut(council_watershed_chm_30m, seq(0, 12, 0.5))) %>%
  drop_na(levels, break_point_30m) %>%
  ggplot(aes(x = levels, y = break_point_30m)) +
  geom_boxplot(fill = '#00CD00') +
  scale_y_continuous(limits = c(120, 180))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  labs(x = "Binned Canopy Height (m)", y = "Completion of Snowmelt (DOY)", title = "Completion of Snowmelt v.s. Canopy Height for Willow Dominated Sites")+
  theme_classic()
DTSW_box_start <- stack_willow %>%
  mutate(levels = cut(council_watershed_chm_30m, seq(0, 12, 0.5))) %>%
  drop_na(levels, break_point_30m) %>%
  ggplot(aes(x = levels, y = break_point1_30m)) +
  geom_boxplot(fill = '#00CD00') +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  labs(x = "Binned Canopy Height (m)", y = "Start of Snowmelt (DOY)", title = "Completion of Snowmelt v.s. Canopy Height for Willow Dominated Sites")+
  theme_classic()
ET_box <- stack_ET %>%
  mutate(levels = cut(council_watershed_chm_30m, seq(0, 12, 0.5))) %>%
  drop_na(levels, break_point_30m) %>%
  ggplot(aes(x = levels, y = break_point_30m)) +
  geom_boxplot(fill = '#FF0000') +
  scale_y_continuous(limits = c(120, 180))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  labs(x = "Binned Canopy Height (m)", y = "Completion of Snowmelt (DOY)", title = "Completion of Snowmelt v.s. Canopy Height for Spruce Dominated Sites")+
  theme_classic()
DT_box <- stack_DT %>%
  mutate(levels = cut(council_watershed_chm_30m, seq(0, 12, 0.5))) %>%
  drop_na(levels, break_point_30m) %>%
  ggplot(aes(x = levels, y = break_point_30m)) +
  geom_boxplot(fill = '#669999') +
  scale_y_continuous(limits = c(120, 180))+
  scale_x_binned()
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  labs(x = "Binned Canopy Height (m)", y = "Completion of Snowmelt (DOY)", title = "Completion of Snowmelt v.s. Canopy Height for Other Tall Shrub Dominated Sites")+
  theme_classic()
plot_grid(ET_box, DT_box, DTSA_box, DTSW_box, ncol = 2, align = "hv")


stack_DG <- stack_df |> subset(stack_df$DryGrass_fcover > 0.75)
DG_box <- stack_DG %>%
  mutate(levels = cut(council_watershed_chm_30m, seq(0, 12, 0.5))) %>%
  drop_na(levels, break_point_30m) %>%
  ggplot(aes(x = levels, y = break_point_30m)) +
  geom_boxplot(fill = '#008000') +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  labs(x = "Binned Canopy Height (m)", y = "Completion of Snowmelt (DOY)", title = "Completion of Snowmelt v.s. Canopy Height for Dry Graminoid Dominated Sites")+
  theme_classic()
stack_WG <- stack_df |> subset(stack_df$WetGrass_fcover > 0.75)
WG_box <- stack_WG %>%
  mutate(levels = cut(council_watershed_chm_30m, seq(0, 12, 0.5))) %>%
  drop_na(levels, break_point_30m) %>%
  ggplot(aes(x = levels, y = break_point_30m)) +
  geom_boxplot(fill = '#008000') +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  labs(x = "Binned Canopy Height (m)", y = "Completion of Snowmelt (DOY)", title = "Completion of Snowmelt v.s. Canopy Height for Wet Graminoid Dominated Sites")+
  theme_classic()

tall_df = filter(stack_df, council_watershed_chm_30m > 5)
mean(tall_df$break_point_30m)
short_df = filter(stack_df, council_watershed_chm_30m < 5)
mean(short_df$break_point_30m)
mean(short_df$break_point_30m) - mean(tall_df$break_point_30m)

stack_tall <- subset(stack_df, council_watershed_chm_30m > 5)
stack_short <- subset(stack_df, council_watershed_chm_30m < 5)
t.test(stack_tall$break_point_30m, stack_short$break_point_30m)

raw_cols <- c('#000000', '#FF0000', '#008000', '#00CD00', '#669999', '#0066FF', '#FFFF66', '#66FFFF', '#00B0F0', '#00FF67', '#FF67FF', '#CC6600', '#FFFFFF', '#5A5A5A')

stack_willow <- stack_df |> subset(stack_df$category == "Willow")
boxplot(break_point_30m ~ levels,data=stack_willow |> mutate(levels = cut(council_watershed_chm_30m,10)),horizontal=FALSE,las=2,cex.axis=0.6, xlab = "Canopy Height Bins")
title("Changes in breakpoint with changing canopy height for the willow PFT")

stack_spruce <- stack_df |> subset(stack_df$category == "Spruce")
boxplot(break_point_30m ~ levels,data=stack_spruce |> mutate(levels = cut(council_watershed_chm_30m,10)),horizontal=FALSE,las=2,cex.axis=0.6, xlab = "Canopy Height Bins")
title("Changes in breakpoint with changing canopy height for the sprcue PFT")

stack_OT <- stack_df |> subset(stack_df$category == "OtherTall")
boxplot(break_point_30m ~ levels,data=stack_OT |> mutate(levels = cut(council_watershed_chm_30m,10)),horizontal=FALSE,las=2,cex.axis=0.6, xlab = "Canopy Height Bins")
title("Changes in breakpoint with changing canopy height for the Other Tall PFT")

stack_lichen <- stack_df |> subset(stack_df$category == "Lichen")
boxplot(break_point_30m ~ levels,data=stack_lichen |> mutate(levels = cut(council_watershed_chm_30m,10)),horizontal=FALSE,las=2,cex.axis=0.6, xlab = "Canopy Height Bins")
title("Changes in breakpoint with changing canopy height for the lichen PFT")

stack_dryG <- stack_df |> subset(stack_df$category == "DryG")
boxplot(break_point_30m ~ levels,data=stack_dryG |> mutate(levels = cut(council_watershed_chm_30m,10)),horizontal=FALSE,las=2,cex.axis=0.6, xlab = "Canopy Height Bins")
title("Changes in breakpoint with changing canopy height for the dryG PFT")


ggplot(data = stack_alder, aes(x = council_watershed_chm_30m, y = break_point_30m))+
  geom_point()+
  geom_smooth(method = "lm", se=FALSE)

ggplot(data = stack_willow, aes(x = council_watershed_chm_30m, y = break_point_30m))+
  geom_point()+
  geom_smooth(method = "lm", se=FALSE)




Calculating twi
upslope <- function (dem, log = TRUE, atb = FALSE, deg = 0.12, fill.sinks = TRUE) 
{
  if (!all.equal(xres(dem), yres(dem))) {
    stop("Raster has differing x and y cell resolutions. Check that it is in a projected coordinate system (e.g. UTM) and use raster::projectRaster to reproject to one if not. Otherwise consider using raster::resample")
  }
  if (fill.sinks) {
    capture.output(dem <- invisible(raster::setValues(dem, topmodel::sinkfill(raster::as.matrix(dem), res = xres(dem), degree = deg))))
  }
  topidx <- topmodel::topidx(raster::as.matrix(dem), res = xres(dem))
  a <- raster::setValues(dem, topidx$area)
  if (log) {
    a <- log(a)
  }
  if (atb) {
    atb <- raster::setValues(dem, topidx$atb)
    a <- addLayer(a, atb)
    names(a) <- c("a", "atb")
  }
  return(a)
}

create_layers <- function (dem, fill.sinks = TRUE, deg = 0.1) 
{
  layers <- stack(dem)
  message("Building upslope areas...")
  a.atb <- upslope(dem, atb = TRUE, fill.sinks = fill.sinks, deg = deg)
  layers <- addLayer(layers, a.atb)
  names(layers) <- c("filled.elevations", "upslope.area", "twi")
  return(layers)
}

layers <- create_layers(topo)
twi.man <- log(layers$upslope.area / tan(slope / 180))
out.dir <- '/Users/anncrumlish/Downloads/albedo/data/'

writeRaster(twi.man, filename=file.path(out.dir, "twi.tif"), format="GTiff", overwrite=TRUE)


plot(twi.man, summer_albedo)


##Make legend for poster
colors = c("#000000", "#FF0000", "#008000", "#00CD00", "#669999", "#0066FF", "#FFFF66", "#66FFFF", "#00B0F0", "#00FF67", "#FF67FF", "#CC6600", "#FFFFFF", "#5A5A5A")
#except lichen needs to be white

PFT = c("Non-Vegetated", "Spruce", "Alder", "Willow", "Other Tall Vegetation", "Low Shrub", "Dwarf Shrub", "Evergreen Shrub", "Forb" ,"Dry Graminoid", "Wet Graminoid", "Moss", "Lichen", "Organic Debris")        
#what is NPV?

plot(0, 0, type = "n", xlim = c(0, 1), ylim = c(0, 1), axes = F, xlab = "", ylab = "")
legend("topright", legend = PFT, fill = colors, border = "black", bty = "n", cex = 1)


library(ggplot2)

# Define the colors and PFTs
colors <- c("#FF0000", "#008000", "#00CD00", "#669999", "#0066FF", "#FFFF66", "#66FFFF", "#00B0F0", "#00FF67", "#FF67FF", "#CC6600", "#FFFFFF", "#5A5A5A", "#000000")
PFT <- c("Evergreen Tree", "Alder", "Willow", "Other Decid. Tree", "Low Shrub", "Dwarf Shrub", "Evergreen Shrub", "Forb" ,"Dry Graminoid", "Wet Graminoid", "Moss", "Lichen", "NPV", "NVS")

library(ggplot2)

# Define the colors and PFTs
colors <- c("#000000", "#FF0000", "#008000", "#00CD00", "#669999", "#0066FF", "#FFFF66", "#66FFFF", "#00B0F0", "#00FF67", "#FF67FF", "#CC6600", "#FFFFFF", "#5A5A5A")
PFT <- c("Non-Vegetated Surface", "Spruce", "Alder", "Willow", "Other Tall Vegetation", "Low Shrub", "Dwarf Shrub", "Evergreen Shrub", "Forb" ,"Dry Graminoid", "Wet Graminoid", "Moss", "Lichen", "Non-Photosynthetic Vegetation")

# Create a data frame
df <- data.frame(PFT, colors)

# Convert PFT to a factor and specify the levels in the original order
df$PFT <- factor(df$PFT, levels = PFT)

# Create a named vector for the fill scale
fill_scale_values <- setNames(colors, PFT)

# Create the plot
ggplot(df, aes(x = 1, y = PFT, fill = PFT)) +
  geom_bar(width = 0.8, color = "black", stat = "identity") +
  scale_fill_manual(values = fill_scale_values) +
  theme_minimal() +
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        legend.key.size = unit(1.5, "cm"),
        legend.key.height = unit(0.7, "cm"),
        legend.key.width = unit(0.7, "cm"), 
        legend.text = element_text(size = 15), 
        legend.title  = element_text(size = 17))



correlation_matrix <- cor(sub_df[,2:7])

# Create the correlation plot
corrplot(correlation_matrix, method = "circle", type = "lower", diag = F)


#Breakpoint DOY and fcover

stack_alder_fcover <- as.data.frame(cbind(stack_df$Alder_fcover, stack_df$break_point_30m))
colnames(stack_alder_fcover) <- c("fcover", "DOY")
stack_alder_fcover <- na.omit(stack_alder_fcover)
DTSA_fbox <- stack_alder_fcover %>%
  mutate(levels = cut(fcover, seq(0, 1, 0.05))) %>%
  drop_na(levels, DOY) %>%
  ggplot(aes(x = levels, y = DOY)) +
  geom_boxplot(fill = '#008000') +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  labs(x = "Binned Alder Fractional Cover (%)", y = "Completion of Snowmelt (DOY)", title = "Completion of Snowmelt v.s. Fractional Cover of Alders")+
  theme_classic()

stack_willow_fcover <- as.data.frame(cbind(stack_df$Willow_fcover, stack_df$break_point_30m))
colnames(stack_willow_fcover) <- c("fcover", "DOY")
stack_willow_fcover <- na.omit(stack_willow_fcover)
DTSW_fbox <- stack_willow_fcover %>%
  mutate(levels = cut(fcover, seq(0, 1, 0.05))) %>%
  drop_na(levels, DOY) %>%
  ggplot(aes(x = levels, y = DOY)) +
  geom_boxplot(fill = "#00CD00") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  labs(x = "Binned Willow Fractional Cover (%)", y = "Completion of Snowmelt (DOY)", title = "Completion of Snowmelt v.s. Fractional Cover of Alders")+
  theme_classic()

stack_spruce_fcover <- as.data.frame(cbind(stack_df$Spruce_fcover, stack_df$break_point_30m))
colnames(stack_spruce_fcover) <- c("fcover", "DOY")
stack_spruce_fcover <- na.omit(stack_spruce_fcover)
ET_fbox <- stack_spruce_fcover %>%
  mutate(levels = cut(fcover, seq(0, 1, 0.05))) %>%
  drop_na(levels, DOY) %>%
  ggplot(aes(x = levels, y = DOY)) +
  geom_boxplot(fill = "#FF0000") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  labs(x = "Binned Willow Fractional Cover (%)", y = "Completion of Snowmelt (DOY)", title = "Completion of Snowmelt v.s. Fractional Cover of Alders")+
  theme_classic()
