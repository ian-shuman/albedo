###########################################################################################
#
#        This script recreates figure 3 in Shuman et al.
#
#    --- Last updated:  2025.03.08 By Ian Shuman <ins2109@columbia.edu>
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
list.of.packages <- c("dplyr", "ggplot2", "terra", "viridisLite", "cowplot")
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
#*****************************************************************************************#

#************************************** load data ****************************************#
#load in albedo timeseries and breakpoints
albedo.dir <- '/Users/anncrumlish/Downloads/albedo analysis/'
albedo.raster <- rast(paste0(albedo.dir, 'step4_smoothed_ts/SG_Smoothed_council_lidar_area_albedo_ts_2013_2019.dat'))
breakpoints <- rast(paste0(albedo.dir, 'step4_smoothed_ts/break_point_30m.tif'))
breakpoints[breakpoints > 213] <- NA

#load in fcover, chm, and breakpoints values over space
data.dir <- '~/Downloads/albedo/Data/'
fcover <- rast(paste0(data.dir, 'council_watershed_fcover_30m.dat')) 
pft <- rast(paste0(data.dir, 'council_watershed_pft_30m.dat'))
canopy_height_model <- rast(paste0(data.dir, 'council_watershed_chm_30m.tif'))
canopy_height_model[canopy_height_model < 0] <- NA #there are about 6000 negative values which should be removed

#make a dataframe for use when plotting the albedo of pfts over time
albedo.raster.pft = resample(albedo.raster, fcover) #resample albedo to preserve heterogeneity in fcover
breakpoints.raster.pft = resample(breakpoints, fcover)
timeseries.stack.pft <- c(albedo.raster.pft, fcover, pft, breakpoints.raster.pft)
timeseries.df.pft <- as.data.frame(timeseries.stack.pft)
#timeseries.df.pft <- na.omit(timeseries.df.pft)
colnames(timeseries.df.pft) <- c(colnames(timeseries.df.pft)[1:236], "Spruce_fcover", "Alder_fcover" , "Willow_fcover", "OtherTall_fcover" , "LowShrubs_fcover", "DwarfShrubs_fcover", "EvergreenShrubs_fcover", "Forb_fcover" , "DryGrass_fcover", "WetGrass_fcover", "Moss_fcover", "Lichen_fcover", "NPV_fcover", "category", "break_point_30m")

#make a dataframe for use when plotting the albedo of different chm bins over time
chm = canopy_height_model
#chm <- resample(canopy_height_model, albedo.raster)
albedo.raster.chm = resample(albedo.raster, chm) #resample albedo to preserve heterogeneity in chm
breakpoints.raster.chm = resample(breakpoints, chm)
timeseries.stack.chm <- c(albedo.raster.chm, chm, breakpoints.raster.chm)
timeseries.df.chm <- as.data.frame(timeseries.stack.chm)
timeseries.df.chm <- na.omit(timeseries.df.chm)
#***********************************************************************************************************************************************#




#********************************** plot albedo timeseries for PFTs with fCover >75% for figure 3a and 3c **************************************#
timeseries_df2 <- timeseries.df.pft[,237:249]
timeseries_df2[is.na(timeseries_df2)] <- 0
# Create a new column "HighCoverSpecies" which is a list of column names where the value is greater than 0.75
timeseries_df2$category2 <- apply(timeseries_df2, 1, function(row) {
  # Identify columns with values greater than 0.75
  category2 <- names(timeseries_df2)[row > 0.75]
  
  # Return the column names as a list
  if (length(category2) == 0) {
    return(NA)  # Return NA if no species have a cover greater than 0.75
  } else {
    return(list(category2))
  }
})

#View(timeseries_df2)
timeseries_df2 <- cbind(timeseries.df.pft[,1:236], timeseries_df2, timeseries.df.pft[,250:251])


timeseries_df2$category3 <- sapply(timeseries_df2$category2, function(x) {
  if (length(x) == 1) {
    return(unlist(x))
  } else {
    return(NA) # Placeholder for lists with more than one element
  }
})


timeseries_df2 <- timeseries_df2[!is.na(timeseries_df2$category3), ]
timeseries_df2 <- timeseries_df2[, !names(timeseries_df2) %in% c("category2", "category3")]

timeseries.ET <- subset(timeseries_df2, category == "Spruce")
timeseries.DTSA <- subset(timeseries_df2, category == "Alder")
timeseries.DTSW <- subset(timeseries_df2, category == "Willow")
timeseries.DT <- subset(timeseries_df2, category == "OtherTall")
timeseries.DLS <- subset(timeseries_df2, category == "LowS")
timeseries.DG <- subset(timeseries_df2, category == "DryG")
timeseries.WG <- subset(timeseries_df2, category == "WetG")
timeseries.MO <- subset(timeseries_df2, category == "Moss")
timeseries.LI <- subset(timeseries_df2, category == "Lichen")
timeseries.NPV <- subset(timeseries_df2, category == "NPV")

bp.ET <- mean(timeseries.ET$break_point_30m, na.rm = T)
bp.DTSA <- mean(timeseries.DTSA$break_point_30m, na.rm = T)
bp.DTSW <- mean(timeseries.DTSW$break_point_30m, na.rm = T)
bp.DT <- mean(timeseries.DT$break_point_30m, na.rm = T)
bp.DLS <- mean(timeseries.DLS$break_point_30m, na.rm = T)
bp.DG <- mean(timeseries.DG$break_point_30m, na.rm = T)
bp.WG <- mean(timeseries.WG$break_point_30m, na.rm = T)
bp.MO <- mean(timeseries.MO$break_point_30m, na.rm = T)
bp.LI <- mean(timeseries.LI$break_point_30m, na.rm = T)
bp.NPV <- mean(timeseries.NPV$break_point_30m, na.rm = T)

bp.pft.mean <- mean(timeseries.df.pft$break_point_30m, na.rm = T)
bp.pft.sd <- sd(timeseries.df.pft$break_point_30m, na.rm = T)
ribbon_data2 <- data.frame(
  x = c(bp.pft.mean - bp.pft.sd, bp.pft.mean + bp.pft.sd),
  ymin = -Inf,
  ymax = Inf
)


DOY  <- seq(50, 296, 1)
pft.series <- cbind(as.numeric(colnames(timeseries_df2[1:(ncol(timeseries_df2)-15)])), as.numeric(colMeans(timeseries.ET[1:236])), as.numeric(colMeans(timeseries.DTSA[1:236])), as.numeric(colMeans(timeseries.DTSW[1:236])), as.numeric(colMeans(timeseries.DT[1:236])), as.numeric(colMeans(timeseries.DLS[1:236])), as.numeric(colMeans(timeseries.DG[1:236])), as.numeric(colMeans(timeseries.WG[1:236])), as.numeric(colMeans(timeseries.MO[1:236])), as.numeric(colMeans(timeseries.LI[1:236])), as.numeric(colMeans(timeseries.NPV[1:236])))
colnames(pft.series) <- c('DOY', 'ET', 'DTSA', 'DTSW', 'DT', 'DLS', 'DG', 'WG', 'MO', 'LI', 'NPV')

doy.min <- 46
doy.max <- 286
bin.size <- 10 #changed from 15

bin.start <- doy.min
bin.end <- bin.start+bin.size
bin.stats.mean75 <- c()
bin.stats.sd75 <- c()
while (bin.end < doy.max) #calculate the 15 day average for albedo per pft over DOY
{
  bin.data <- pft.series[which(pft.series > bin.start & pft.series < bin.end), ]
  z.scores <- apply(bin.data, 2, scale)
  z.scores[is.na(z.scores)] <- 0
  bin.data[abs(z.scores) > 1] <- NaN
  bin.mean <- colMeans(bin.data, na.rm = T)
  bin.sd <- apply(bin.data, 2, sd, na.rm = TRUE)
  
  bin.name <- bin.start + 7.5
  
  bin.mean <- c(bin.name, bin.mean[-1])
  bin.stats.mean75 <- rbind(bin.stats.mean75, bin.mean)
  
  bin.sd <- c(bin.name, bin.sd[-1])
  bin.stats.sd75 <- rbind(bin.stats.sd75, bin.sd)
  
  bin.start <- bin.start + bin.size
  bin.end <- bin.start+bin.size
}

bin.stats.mean75 <- data.frame(bin.stats.mean75)
names(bin.stats.mean75) <- c('DOY', 'ET', 'DTSA', 'DTSW', 'DT', 'DLS', 'DG', 'WG', 'MO', 'LI', 'NPV')

bin.stats.sd75<- data.frame(bin.stats.sd75)
names(bin.stats.sd75) <- c('DOY', 'ET', 'DTSA', 'DTSW', 'DT', 'DLS', 'DG', 'WG', 'MO', 'LI', 'NPV')

xstart <- 200
xend <- 220

y1p <- 1
y2p <- 0.95
y3p <- 0.9
y4p <- 0.85
y5p <- 0.8
y6p <- 0.75
y7p <- 0.7
y8p <- 0.65
y9p <- 0.60
y10p <- 0.55
y11p <- 0.5

x_int = 10
y_r <- 0.015

cols <- c('#FF0000', '#008000', '#00CD00', '#669999', '#0066FF', '#00FF67', '#FF67FF', '#CC6600', '#FFFFFF', '#5A5A5A')
pft.series.plot = ggplot(data = bin.stats.mean75) +
  geom_point(aes(x = DOY, y = `ET`), colour = cols[1], size = 2.5) +
  geom_line(aes(x = DOY, y = `ET`), colour = cols[1], size = 1.2) +
  geom_ribbon(aes(x = DOY, ymin = bin.stats.mean75$`ET`-bin.stats.sd75$`ET`, ymax = bin.stats.mean75$`ET`+bin.stats.sd75$`ET`),
              fill = cols[1],alpha = 0.2) +
  geom_point(aes(x = DOY, y = `DTSA`), colour = cols[2], size = 2.5) +
  geom_line(aes(x = DOY, y = `DTSA`), colour = cols[2], size = 1.2,) +
  geom_ribbon(aes(x = DOY, ymin = bin.stats.mean75$`DTSA`-bin.stats.sd75$`DTSA`, ymax = bin.stats.mean75$`DTSA`+bin.stats.sd75$`DTSA`),
              fill = cols[2],alpha = 0.2)  +
  geom_point(aes(x = DOY, y = `DTSW`), colour = cols[3], size = 2.5) +
  geom_line(aes(x = DOY, y = `DTSW`), colour = cols[3], size = 1.2) +
  geom_ribbon(aes(x = DOY, ymin = bin.stats.mean75$`DTSW`-bin.stats.sd75$`DTSW`, ymax = bin.stats.mean75$`DTSW`+bin.stats.sd75$`DTSW`),
              fill = cols[3],alpha = 0.2) +
  geom_point(aes(x = DOY, y = `DT`), colour = cols[4], size = 2.5) +
  geom_line(aes(x = DOY, y = `DT`), colour = cols[4], size = 1.2,) +
  geom_ribbon(aes(x = DOY, ymin = bin.stats.mean75$`DT`-bin.stats.sd75$`DT`, ymax = bin.stats.mean75$`DT`+bin.stats.sd75$`DT`),
              fill = cols[4],alpha = 0.2)  +
  geom_point(aes(x = DOY, y = `DLS`), colour = cols[5], size = 2.5) +
  geom_line(aes(x = DOY, y = `DLS`), colour = cols[5], size = 1.2) +
  geom_ribbon(aes(x = DOY, ymin = bin.stats.mean75$`DLS`-bin.stats.sd75$`DLS`, ymax = bin.stats.mean75$`DLS`+bin.stats.sd75$`DLS`),
              fill = cols[5],alpha = 0.2) +
  geom_point(aes(x = DOY, y = `DG`), colour = cols[6], size = 2.5) +
  geom_line(aes(x = DOY, y = `DG`), colour = cols[6], size = 1.2,) +
  geom_ribbon(aes(x = DOY, ymin = bin.stats.mean75$`DG`-bin.stats.sd75$`DG`, ymax = bin.stats.mean75$`DG`+bin.stats.sd75$`DG`),
              fill = cols[6],alpha = 0.2)  +
  geom_point(aes(x = DOY, y = `WG`), colour = cols[7], size = 2.5) +
  geom_line(aes(x = DOY, y = `WG`), colour = cols[7], size = 1.2) +
  geom_ribbon(aes(x = DOY, ymin = bin.stats.mean75$`WG`-bin.stats.sd75$`WG`, ymax = bin.stats.mean75$`WG`+bin.stats.sd75$`WG`),
              fill = cols[7],alpha = 0.2) +
  geom_point(aes(x = DOY, y = `MO`), colour = cols[8], size = 2.5) +
  geom_line(aes(x = DOY, y = `MO`), colour = cols[8], size = 1.2) +
  geom_ribbon(aes(x = DOY, ymin = bin.stats.mean75$`MO`-bin.stats.sd75$`MO`, ymax = bin.stats.mean75$`MO`+bin.stats.sd75$`MO`),
              fill = cols[8],alpha = 0.2) +
  geom_line(aes(x = DOY, y = `LI`), colour = "black",  size = 1.5) +
  geom_point(aes(x = DOY, y = `LI`), colour = cols[9], size = 2.5) +
  geom_line(aes(x = DOY, y = `LI`), colour = cols[9], size = 1.2,) +
  geom_ribbon(aes(x = DOY, ymin = bin.stats.mean75$`LI`-bin.stats.sd75$`LI`, ymax = bin.stats.mean75$`LI`+bin.stats.sd75$`LI`),
              fill = "aliceblue",alpha = 0.2)  +
  geom_point(aes(x = DOY, y = `NPV`), colour = cols[10], size = 2.5) +
  geom_line(aes(x = DOY, y = `NPV`), colour = cols[10], size = 1.2,) +
  geom_ribbon(aes(x = DOY, ymin = bin.stats.mean75$`NPV`-bin.stats.sd75$`NPV`, ymax = bin.stats.mean75$`NPV`+bin.stats.sd75$`NPV`),
              fill = cols[10],alpha = 0.2)  +
  labs(x = NULL, y = NULL) + ylim(0.1, 1) +
  
  #geom_point(aes(x = bp.ET, y = 0.17), color = cols[1], pch = 8, size = 3)+
  #geom_point(aes(x = bp.DTSA, y = 0.19), color = cols[2], pch = 8, size = 3)+
  #geom_point(aes(x = bp.DTSW, y = 0.2388), color = cols[3], pch = 8, size = 3)+
  #geom_point(aes(x = bp.DT, y = 0.17), color = cols[4], pch = 8, size = 3)+
  #geom_point(aes(x = bp.DLS, y = 0.247), color = cols[5], pch = 8, size = 3)+
  #geom_point(aes(x = bp.DG, y = 0.2775), color = cols[6], pch = 8, size = 3)+
  #geom_point(aes(x = bp.WG, y = 0.2388), color = cols[7], pch = 8, size = 3)+
  #geom_point(aes(x = bp.MO, y = 0.255), color = cols[8], pch = 8, size = 3)+
  #geom_point(aes(x = bp.LI, y = 0.295), color = "black", pch = 8, size = 3)+
  #geom_point(aes(x = bp.NPV, y = 0.253), color = cols[10], pch = 8, size = 3)+
  geom_vline(xintercept = bp.pft.mean, color = "darkgrey")+
  geom_ribbon(data = ribbon_data2, aes(x = x, ymin = ymin, ymax = ymax), alpha = 0.2, fill = "darkgrey") +
  
  annotate("text", x = xstart + 2*x_int,  y = y1p, label = 'PFT', vjust=-0.2, size =6) +
  geom_segment(aes(x = xstart, xend = xend, y = y2p, yend = y2p), color = cols[1], size = 1.5) +
  annotate("rect", xmin = xstart, xmax = xend, ymin = y2p-y_r, ymax = y2p+y_r, alpha = .2, fill = cols[1]) +
  annotate("text", x = xend + x_int,  y = y2p, label = 'ET', hjust=0, size =6) +
  
  geom_segment(aes(x = xstart, xend = xend, y = y3p, yend = y3p), color = cols[2], size = 1.5) +
  annotate("rect", xmin = xstart, xmax = xend, ymin = y3p-y_r, ymax = y3p+y_r, alpha = .2, fill = cols[2]) +
  annotate("text", x = xend + x_int,  y = y3p, label = 'DTSA', hjust=0, size =6) +
  
  geom_segment(aes(x = xstart, xend = xend, y = y4p, yend = y4p), color = cols[3], size = 1.5) +
  annotate("rect", xmin = xstart, xmax = xend, ymin = y4p-y_r, ymax = y4p+y_r, alpha = .2, fill = cols[3]) +
  annotate("text", x = xend + x_int,  y = y4p, label = 'DTSW', hjust=0, size =6) +
  
  geom_segment(aes(x = xstart, xend = xend, y = y5p, yend = y5p), color = cols[4], size = 1.5) +
  annotate("rect", xmin = xstart, xmax = xend, ymin = y5p-y_r, ymax = y5p+y_r, alpha = .2, fill = cols[4]) +
  annotate("text", x = xend + x_int,  y = y5p, label = 'DT', hjust=0, size =6) +
  
  geom_segment(aes(x = xstart, xend = xend, y = y6p, yend = y6p), color = cols[5], size = 1.5) +
  annotate("rect", xmin = xstart, xmax = xend, ymin = y6p-y_r, ymax = y6p+y_r, alpha = .2, fill = cols[5]) +
  annotate("text", x = xend + x_int,  y = y6p, label = 'DLS', hjust=0, size =6) +
  
  geom_segment(aes(x = xstart, xend = xend, y = y7p, yend = y7p), color = cols[6], size = 1.5) +
  annotate("rect", xmin = xstart, xmax = xend, ymin = y7p-y_r, ymax = y7p+y_r, alpha = .2, fill = cols[6]) +
  annotate("text", x = xend + x_int,  y = y7p, label = 'DG', hjust=0, size =6) +
  
  geom_segment(aes(x = xstart, xend = xend, y = y8p, yend = y8p), color = cols[7], size = 1.5) +
  annotate("rect", xmin = xstart, xmax = xend, ymin = y8p-y_r, ymax = y8p+y_r, alpha = .2, fill = cols[7]) +
  annotate("text", x = xend + x_int,  y = y8p, label = 'WG', hjust=0, size =6) +
  
  geom_segment(aes(x = xstart, xend = xend, y = y9p, yend = y9p), color = cols[8], size = 1.5) +
  annotate("rect", xmin = xstart, xmax = xend, ymin = y9p-y_r, ymax = y9p+y_r, alpha = .2, fill = cols[8]) +
  annotate("text", x = xend + x_int,  y = y9p, label = 'MO', hjust=0, size =6) +
  
  geom_segment(aes(x = xstart, xend = xend, y = y10p, yend = y10p), color = "black", size = 1.8) +
  geom_segment(aes(x = xstart, xend = xend, y = y10p, yend = y10p), color = cols[9], size = 1.5) +
  annotate("rect", xmin = xstart, xmax = xend, ymin = y10p-y_r, ymax = y10p-0.005, alpha = .7, fill = "aliceblue") +
  annotate("rect", xmin = xstart, xmax = xend, ymin = y10p+0.005, ymax = y10p+y_r, alpha = .7, fill = "aliceblue") +
  annotate("text", x = xend + x_int,  y = y10p, label = 'LI', hjust=0, size = 6) +
  
  geom_segment(aes(x = xstart, xend = xend, y = y11p, yend = y11p), color = cols[10], size = 1.5) +
  annotate("rect", xmin = xstart, xmax = xend, ymin = y11p-y_r, ymax = y11p+y_r, alpha = .2, fill = cols[10]) +
  annotate("text", x = xend + x_int,  y = y11p, label = 'NPV', hjust=0, size =6) +
  
  theme_classic()+
  theme(axis.text = element_text(size=27), axis.title=element_text(size=29))+
  labs(x = "Day of Year", y = "Albedo") + ylim(0.1, 1)
pft.series.plot

#Zoom into DOY 150-250 to show summer dynamics with more detail (Fig. 3c)
bin.stats.mean752 = filter(bin.stats.mean75, DOY > 150)
bin.stats.sd752 = filter(bin.stats.sd75, DOY > 150)
pft.bottom = ggplot(data = bin.stats.mean752) +
  geom_point(aes(x = DOY, y = `ET`), colour = cols[1], size = 2.5) +
  geom_line(aes(x = DOY, y = `ET`), colour = cols[1], size = 1.2) +
  geom_ribbon(aes(x = DOY, ymin = bin.stats.mean752$`ET`-bin.stats.sd752$`ET`, ymax = bin.stats.mean752$`ET`+bin.stats.sd752$`ET`),
              fill = cols[1],alpha = 0.2) +
  geom_point(aes(x = DOY, y = `DTSA`), colour = cols[2], size = 2.5) +
  geom_line(aes(x = DOY, y = `DTSA`), colour = cols[2], size = 1.2,) +
  geom_ribbon(aes(x = DOY, ymin = bin.stats.mean752$`DTSA`-bin.stats.sd752$`DTSA`, ymax = bin.stats.mean752$`DTSA`+bin.stats.sd752$`DTSA`),
              fill = cols[2],alpha = 0.2)  +
  geom_point(aes(x = DOY, y = `DTSW`), colour = cols[3], size = 2.5) +
  geom_line(aes(x = DOY, y = `DTSW`), colour = cols[3], size = 1.2) +
  geom_ribbon(aes(x = DOY, ymin = bin.stats.mean752$`DTSW`-bin.stats.sd752$`DTSW`, ymax = bin.stats.mean752$`DTSW`+bin.stats.sd752$`DTSW`),
              fill = cols[3],alpha = 0.2) +
  geom_point(aes(x = DOY, y = `DT`), colour = cols[4], size = 2.5) +
  geom_line(aes(x = DOY, y = `DT`), colour = cols[4], size = 1.2,) +
  geom_ribbon(aes(x = DOY, ymin = bin.stats.mean752$`DT`-bin.stats.sd752$`DT`, ymax = bin.stats.mean752$`DT`+bin.stats.sd752$`DT`),
              fill = cols[4],alpha = 0.2)  +
  geom_point(aes(x = DOY, y = `DLS`), colour = cols[5], size = 2.5) +
  geom_line(aes(x = DOY, y = `DLS`), colour = cols[5], size = 1.2) +
  geom_ribbon(aes(x = DOY, ymin = bin.stats.mean752$`DLS`-bin.stats.sd752$`DLS`, ymax = bin.stats.mean752$`DLS`+bin.stats.sd752$`DLS`),
              fill = cols[5],alpha = 0.2) +
  geom_point(aes(x = DOY, y = `DG`), colour = cols[6], size = 2.5) +
  geom_line(aes(x = DOY, y = `DG`), colour = cols[6], size = 1.2,) +
  geom_ribbon(aes(x = DOY, ymin = bin.stats.mean752$`DG`-bin.stats.sd752$`DG`, ymax = bin.stats.mean752$`DG`+bin.stats.sd752$`DG`),
              fill = cols[6],alpha = 0.2)  +
  geom_point(aes(x = DOY, y = `WG`), colour = cols[7], size = 2.5) +
  geom_line(aes(x = DOY, y = `WG`), colour = cols[7], size = 1.2) +
  geom_ribbon(aes(x = DOY, ymin = bin.stats.mean752$`WG`-bin.stats.sd752$`WG`, ymax = bin.stats.mean752$`WG`+bin.stats.sd752$`WG`),
              fill = cols[7],alpha = 0.2) +
  geom_point(aes(x = DOY, y = `MO`), colour = cols[8], size = 2.5) +
  geom_line(aes(x = DOY, y = `MO`), colour = cols[8], size = 1.2) +
  geom_ribbon(aes(x = DOY, ymin = bin.stats.mean752$`MO`-bin.stats.sd752$`MO`, ymax = bin.stats.mean752$`MO`+bin.stats.sd752$`MO`),
              fill = cols[8],alpha = 0.2) +
  geom_line(aes(x = DOY, y = `LI`), colour = "black",  size = 1.5) +
  geom_point(aes(x = DOY, y = `LI`), colour = cols[9], size = 2.5) +
  geom_line(aes(x = DOY, y = `LI`), colour = cols[9], size = 1.2,) +
  geom_ribbon(aes(x = DOY, ymin = bin.stats.mean752$`LI`-bin.stats.sd752$`LI`, ymax = bin.stats.mean752$`LI`+bin.stats.sd752$`LI`),
              fill = "aliceblue",alpha = 0.2)  +
  geom_point(aes(x = DOY, y = `NPV`), colour = cols[10], size = 2.5) +
  geom_line(aes(x = DOY, y = `NPV`), colour = cols[10], size = 1.2,) +
  geom_ribbon(aes(x = DOY, ymin = bin.stats.mean752$`NPV`-bin.stats.sd752$`NPV`, ymax = bin.stats.mean752$`NPV`+bin.stats.sd752$`NPV`),
              fill = cols[10],alpha = 0.2)  +
  labs(x = 'Day of Year', y = 'Albedo') + ylim(0.1, 0.25) +
  
  theme_classic()+
  theme(axis.text = element_text(size=27), axis.title=element_text(size=29)) 
pft.bottom
#***********************************************************************************************************************************************#



#********************************** plot albedo timeseries for 1m canopy height bins for figure 3b and 3d **************************************#
timeseries.s <- subset(timeseries.df.chm, council_watershed_chm_30m <5)
timeseries.l <- subset(timeseries.df.chm, council_watershed_chm_30m >5)

timeseries.0 <- subset(timeseries.df.chm, council_watershed_chm_30m <1)
timeseries.1 <- subset(timeseries.df.chm, council_watershed_chm_30m <2 & council_watershed_chm_30m >1)
timeseries.2 <- subset(timeseries.df.chm, council_watershed_chm_30m <3 & council_watershed_chm_30m >2)
timeseries.3 <- subset(timeseries.df.chm, council_watershed_chm_30m <4 & council_watershed_chm_30m >3)
timeseries.4 <- subset(timeseries.df.chm, council_watershed_chm_30m <5 & council_watershed_chm_30m >4)
timeseries.5 <- subset(timeseries.df.chm, council_watershed_chm_30m <6 & council_watershed_chm_30m >5)
timeseries.6 <- subset(timeseries.df.chm, council_watershed_chm_30m <7 & council_watershed_chm_30m >6)
timeseries.7 <- subset(timeseries.df.chm, council_watershed_chm_30m <8 & council_watershed_chm_30m >7)
timeseries.8 <- subset(timeseries.df.chm, council_watershed_chm_30m >8)

bp.0 <- mean(timeseries.0$break_point_30m)
bp.1 <- mean(timeseries.1$break_point_30m)
bp.2 <- mean(timeseries.2$break_point_30m)
bp.3 <- mean(timeseries.3$break_point_30m)
bp.4 <- mean(timeseries.4$break_point_30m)
bp.5 <- mean(timeseries.5$break_point_30m)
bp.6 <- mean(timeseries.6$break_point_30m)
bp.7 <- mean(timeseries.7$break_point_30m)
bp.8 <- mean(timeseries.8$break_point_30m)

bp.chm.mean <- mean(timeseries.df.chm$break_point_30m, na.rm = T)
bp.chm.sd <- sd(timeseries.df.chm$break_point_30m, na.rm = T)
ribbon_data <- data.frame(
  x = c(bp.chm.mean - bp.chm.sd, bp.chm.mean + bp.chm.sd),
  ymin = -Inf,
  ymax = Inf
)


DOY <- seq(50, 296, 1)
chm.series <- cbind(as.numeric(colnames(timeseries.df.chm[1:(ncol(timeseries.df.chm)-2)])), as.numeric(colMeans(timeseries.0)[1:236]), as.numeric(colMeans(timeseries.1)[1:236]), as.numeric(colMeans(timeseries.2)[1:236]), as.numeric(colMeans(timeseries.3)[1:236]), as.numeric(colMeans(timeseries.4)[1:236]), as.numeric(colMeans(timeseries.5)[1:236]), as.numeric(colMeans(timeseries.6)[1:236]), as.numeric(colMeans(timeseries.7)[1:236]), as.numeric(colMeans(timeseries.8)[1:236]))
colnames(chm.series) <- c('DOY', '<1 m', '1-2 m', '2-3 m', '3-4 m', '4-5 m', '5-6 m', '6-7 m', '7-8 m', '8+ m')

doy.min <- 46
doy.max <- 286
bin.size <- 10 #changed from 15

bin.start <- doy.min
bin.end <- bin.start+bin.size
bin.stats.mean2 <- c()
bin.stats.sd2 <- c()
while (bin.end < doy.max) #calculate the 15 day average for albedo per pft over DOY
{
  bin.data <- chm.series[which(chm.series > bin.start & chm.series < bin.end), ]
  z.scores <- apply(bin.data, 2, scale)
  z.scores[is.na(z.scores)] <- 0
  bin.data[abs(z.scores) > 1] <- NaN
  bin.mean <- colMeans(bin.data, na.rm = T)
  bin.sd <- apply(bin.data, 2, sd, na.rm = TRUE)
  
  bin.name <- bin.start + 7.5
  
  bin.mean <- c(bin.name, bin.mean[-1])
  bin.stats.mean2 <- rbind(bin.stats.mean2, bin.mean)
  
  bin.sd <- c(bin.name, bin.sd[-1])
  bin.stats.sd2 <- rbind(bin.stats.sd2, bin.sd)
  
  bin.start <- bin.start + bin.size
  bin.end <- bin.start+bin.size
}

bin.stats.mean2 <- data.frame(bin.stats.mean2)
names(bin.stats.mean2) <- c('DOY', '<1 m', '1-2 m', '2-3 m', '3-4 m', '4-5 m', '5-6 m', '6-7 m', '7-8 m', '8+ m')

bin.stats.sd2 <- data.frame(bin.stats.sd2)
names(bin.stats.sd2) <- c('DOY', '<1 m', '1-2 m', '2-3 m', '3-4 m', '4-5 m', '5-6 m', '6-7 m', '7-8 m', '8+ m')

xstart <- 200
xend <- 220

y1b <- 1
y2b <- 0.95
y3b <- 0.9
y4b <- 0.85
y5b <- 0.8
y6b <- 0.75
y7b <- 0.7
y8b <- 0.65
y9b <- 0.60
y10b <- 0.55

x_int = 10
y_r <- 0.015

cols <- viridis(n = 9)

#cols <- RColorBrewer::brewer.pal(10, 'Spectral')
#cols <- rev(cols)
chm.series.plot = ggplot(data = bin.stats.mean2) +
  geom_point(aes(x = DOY, y = `<1 m`), colour = cols[1], size = 2.5) +
  geom_line(aes(x = DOY, y = `<1 m`), colour = cols[1], size = 1.2) +
  geom_ribbon(aes(x = DOY, ymin = bin.stats.mean2$`<1 m`-bin.stats.sd2$`<1 m`, ymax = bin.stats.mean2$`<1 m`+bin.stats.sd2$`<1 m`),
              fill = cols[1],alpha = 0.2) +
  geom_point(aes(x = DOY, y = `1-2 m`), colour = cols[2], size = 2.5) +
  geom_line(aes(x = DOY, y = `1-2 m`), colour = cols[2], size = 1.2,) +
  geom_ribbon(aes(x = DOY, ymin = bin.stats.mean2$`1-2 m`-bin.stats.sd2$`1-2 m`, ymax = bin.stats.mean2$`1-2 m`+bin.stats.sd2$`1-2 m`),
              fill = cols[2],alpha = 0.2)  +
  geom_point(aes(x = DOY, y = `2-3 m`), colour = cols[3], size = 2.5) +
  geom_line(aes(x = DOY, y = `2-3 m`), colour = cols[3], size = 1.2) +
  geom_ribbon(aes(x = DOY, ymin = bin.stats.mean2$`2-3 m`-bin.stats.sd2$`2-3 m`, ymax = bin.stats.mean2$`2-3 m`+bin.stats.sd2$`2-3 m`),
              fill = cols[3],alpha = 0.2) +
  geom_point(aes(x = DOY, y = `3-4 m`), colour = cols[4], size = 2.5) +
  geom_line(aes(x = DOY, y = `3-4 m`), colour = cols[4], size = 1.2,) +
  geom_ribbon(aes(x = DOY, ymin = bin.stats.mean2$`3-4 m`-bin.stats.sd2$`3-4 m`, ymax = bin.stats.mean2$`3-4 m`+bin.stats.sd2$`3-4 m`),
              fill = cols[4],alpha = 0.2)  +
  geom_point(aes(x = DOY, y = `4-5 m`), colour = cols[5], size = 2.5) +
  geom_line(aes(x = DOY, y = `4-5 m`), colour = cols[5], size = 1.2) +
  geom_ribbon(aes(x = DOY, ymin = bin.stats.mean2$`4-5 m`-bin.stats.sd2$`4-5 m`, ymax = bin.stats.mean2$`4-5 m`+bin.stats.sd2$`4-5 m`),
              fill = cols[5],alpha = 0.2) +
  geom_point(aes(x = DOY, y = `5-6 m`), colour = cols[6], size = 2.5) +
  geom_line(aes(x = DOY, y = `5-6 m`), colour = cols[6], size = 1.2,) +
  geom_ribbon(aes(x = DOY, ymin = bin.stats.mean2$`5-6 m`-bin.stats.sd2$`5-6 m`, ymax = bin.stats.mean2$`5-6 m`+bin.stats.sd2$`5-6 m`),
              fill = cols[6],alpha = 0.2)  +
  geom_point(aes(x = DOY, y = `6-7 m`), colour = cols[7], size = 2.5) +
  geom_line(aes(x = DOY, y = `6-7 m`), colour = cols[7], size = 1.2) +
  geom_ribbon(aes(x = DOY, ymin = bin.stats.mean2$`6-7 m`-bin.stats.sd2$`6-7 m`, ymax = bin.stats.mean2$`6-7 m`+bin.stats.sd2$`6-7 m`),
              fill = cols[7],alpha = 0.2) +
  geom_point(aes(x = DOY, y = `7-8 m`), colour = cols[8], size = 2.5) +
  geom_line(aes(x = DOY, y = `7-8 m`), colour = cols[8], size = 1.2) +
  geom_ribbon(aes(x = DOY, ymin = bin.stats.mean2$`7-8 m`-bin.stats.sd2$`7-8 m`, ymax = bin.stats.mean2$`7-8 m`+bin.stats.sd2$`7-8 m`),
              fill = cols[8],alpha = 0.2) +
  geom_point(aes(x = DOY, y = `8+ m`), colour = cols[9], size = 2.5) +
  geom_line(aes(x = DOY, y = `8+ m`), colour = cols[9], size = 1.2,) +
  geom_ribbon(aes(x = DOY, ymin = bin.stats.mean2$`8+ m`-bin.stats.sd2$`8+ m`, ymax = bin.stats.mean2$`8+ m`+bin.stats.sd2$`8+ m`),
              fill = cols[9],alpha = 0.2)  +
  labs(x = NULL, y = NULL) + ylim(0.1, 1) +
  annotate("text", x = xstart + 2*x_int,  y = y1b, label = 'Canopy Height', vjust=-0.2, size =6) +
  geom_segment(aes(x = xstart, xend = xend, y = y2b, yend = y2b), color = cols[1], size = 1.5) +
  #geom_point(aes(x = bp.0, y = 0.26), color = cols[1], pch = 8, size = 3)+
  #geom_point(aes(x = bp.1, y = 0.2125), color = cols[2], pch = 8, size = 3)+
  #geom_point(aes(x = bp.2, y = 0.199), color = cols[3], pch = 8, size = 3)+
  #geom_point(aes(x = bp.3, y = 0.185), color = cols[4], pch = 8, size = 3)+
  #geom_point(aes(x = bp.4, y = 0.1765), color = cols[5], pch = 8, size = 3)+
  #geom_point(aes(x = bp.5, y = 0.1575), color = cols[6], pch = 8, size = 3)+
  #geom_point(aes(x = bp.6, y = 0.1485), color = cols[7], pch = 8, size = 3)+
  #geom_point(aes(x = bp.7, y = 0.1475), color = cols[8], pch = 8, size = 3)+
  #geom_point(aes(x = bp.8, y = 0.1435), color = cols[9], pch = 8, size = 3)+
  geom_vline(xintercept = bp.chm.mean, color = "darkgrey")+
  geom_ribbon(data = ribbon_data, aes(x = x, ymin = ymin, ymax = ymax), alpha = 0.2, fill = "darkgrey") +
  annotate("rect", xmin = xstart, xmax = xend, ymin = y2b-y_r, ymax = y2b+y_r, alpha = .2, fill = cols[1]) +
  annotate("text", x = xend + x_int,  y = y2b, label = '<1 m', hjust=0, size =6) +
  
  geom_segment(aes(x = xstart, xend = xend, y = y3b, yend = y3b), color = cols[2], size = 1.5) +
  annotate("rect", xmin = xstart, xmax = xend, ymin = y3b-y_r, ymax = y3b+y_r, alpha = .2, fill = cols[2]) +
  annotate("text", x = xend + x_int,  y = y3b, label = '1-2 m', hjust=0, size =6) +
  
  geom_segment(aes(x = xstart, xend = xend, y = y4b, yend = y4b), color = cols[3], size = 1.5) +
  annotate("rect", xmin = xstart, xmax = xend, ymin = y4b-y_r, ymax = y4b+y_r, alpha = .2, fill = cols[3]) +
  annotate("text", x = xend + x_int,  y = y4b, label = '2-3 m', hjust=0, size =6) +
  
  geom_segment(aes(x = xstart, xend = xend, y = y5b, yend = y5b), color = cols[4], size = 1.5) +
  annotate("rect", xmin = xstart, xmax = xend, ymin = y5b-y_r, ymax = y5b+y_r, alpha = .2, fill = cols[4]) +
  annotate("text", x = xend + x_int,  y = y5b, label = '3-4 m', hjust=0, size =6) +
  
  geom_segment(aes(x = xstart, xend = xend, y = y6b, yend = y6b), color = cols[5], size = 1.5) +
  annotate("rect", xmin = xstart, xmax = xend, ymin = y6b-y_r, ymax = y6b+y_r, alpha = .2, fill = cols[5]) +
  annotate("text", x = xend + x_int,  y = y6b, label = '4-5 m', hjust=0, size =6) +
  
  geom_segment(aes(x = xstart, xend = xend, y = y7b, yend = y7b), color = cols[6], size = 1.5) +
  annotate("rect", xmin = xstart, xmax = xend, ymin = y7b-y_r, ymax = y7b+y_r, alpha = .2, fill = cols[6]) +
  annotate("text", x = xend + x_int,  y = y7b, label = '5-6 m', hjust=0, size =6) +
  
  geom_segment(aes(x = xstart, xend = xend, y = y8b, yend = y8b), color = cols[7], size = 1.5) +
  annotate("rect", xmin = xstart, xmax = xend, ymin = y8b-y_r, ymax = y8b+y_r, alpha = .2, fill = cols[7]) +
  annotate("text", x = xend + x_int,  y = y8b, label = '6-7 m', hjust=0, size =6) +
  
  geom_segment(aes(x = xstart, xend = xend, y = y9b, yend = y9b), color = cols[8], size = 1.5) +
  annotate("rect", xmin = xstart, xmax = xend, ymin = y9b-y_r, ymax = y9b+y_r, alpha = .2, fill = cols[8]) +
  annotate("text", x = xend + x_int,  y = y9b, label = '7-8 m', hjust=0, size =6) +
  geom_segment(aes(x = xstart, xend = xend, y = y10b, yend = y10b), color = cols[9], size = 1.5) +
  annotate("rect", xmin = xstart, xmax = xend, ymin = y10b-y_r, ymax = y10b+y_r, alpha = .2, fill = cols[9]) +
  annotate("text", x = xend + x_int,  y = y10b, label = '8+ m', hjust=0, size =6) +
  theme_classic()+
  theme(axis.text = element_text(size=27), 
        axis.title=element_text(size=29), 
        legend.position = c(1,0.5))+
  labs(x = "Day of Year", y = NULL) + ylim(0.1, 1) 
chm.series.plot

#Zoom into DOY 150-250 to show summer dynamics with more detail (Fig. 3d)
bin.stats.mean3 = filter(bin.stats.mean2, DOY > 150)
bin.stats.sd3 = filter(bin.stats.sd2, DOY > 150)
chm.bottom = ggplot(data = bin.stats.mean3) +
  geom_point(aes(x = DOY, y = `<1 m`), colour = cols[1], size = 2.5) +
  geom_line(aes(x = DOY, y = `<1 m`), colour = cols[1], size = 1.2) +
  geom_ribbon(aes(x = DOY, ymin = bin.stats.mean3$`<1 m`-bin.stats.sd3$`<1 m`, ymax = bin.stats.mean3$`<1 m`+bin.stats.sd3$`<1 m`),
              fill = cols[1],alpha = 0.2) +
  geom_point(aes(x = DOY, y = `1-2 m`), colour = cols[2], size = 2.5) +
  geom_line(aes(x = DOY, y = `1-2 m`), colour = cols[2], size = 1.2,) +
  geom_ribbon(aes(x = DOY, ymin = bin.stats.mean3$`1-2 m`-bin.stats.sd3$`1-2 m`, ymax = bin.stats.mean3$`1-2 m`+bin.stats.sd3$`1-2 m`),
              fill = cols[2],alpha = 0.2)  +
  geom_point(aes(x = DOY, y = `2-3 m`), colour = cols[3], size = 2.5) +
  geom_line(aes(x = DOY, y = `2-3 m`), colour = cols[3], size = 1.2) +
  geom_ribbon(aes(x = DOY, ymin = bin.stats.mean3$`2-3 m`-bin.stats.sd3$`2-3 m`, ymax = bin.stats.mean3$`2-3 m`+bin.stats.sd3$`2-3 m`),
              fill = cols[3],alpha = 0.2) +
  geom_point(aes(x = DOY, y = `3-4 m`), colour = cols[4], size = 2.5) +
  geom_line(aes(x = DOY, y = `3-4 m`), colour = cols[4], size = 1.2,) +
  geom_ribbon(aes(x = DOY, ymin = bin.stats.mean3$`3-4 m`-bin.stats.sd3$`3-4 m`, ymax = bin.stats.mean3$`3-4 m`+bin.stats.sd3$`3-4 m`),
              fill = cols[4],alpha = 0.2)  +
  geom_point(aes(x = DOY, y = `4-5 m`), colour = cols[5], size = 2.5) +
  geom_line(aes(x = DOY, y = `4-5 m`), colour = cols[5], size = 1.2) +
  geom_ribbon(aes(x = DOY, ymin = bin.stats.mean3$`4-5 m`-bin.stats.sd3$`4-5 m`, ymax = bin.stats.mean3$`4-5 m`+bin.stats.sd3$`4-5 m`),
              fill = cols[5],alpha = 0.2) +
  geom_point(aes(x = DOY, y = `5-6 m`), colour = cols[6], size = 2.5) +
  geom_line(aes(x = DOY, y = `5-6 m`), colour = cols[6], size = 1.2,) +
  geom_ribbon(aes(x = DOY, ymin = bin.stats.mean3$`5-6 m`-bin.stats.sd3$`5-6 m`, ymax = bin.stats.mean3$`5-6 m`+bin.stats.sd3$`5-6 m`),
              fill = cols[6],alpha = 0.2)  +
  geom_point(aes(x = DOY, y = `6-7 m`), colour = cols[7], size = 2.5) +
  geom_line(aes(x = DOY, y = `6-7 m`), colour = cols[7], size = 1.2) +
  geom_ribbon(aes(x = DOY, ymin = bin.stats.mean3$`6-7 m`-bin.stats.sd3$`6-7 m`, ymax = bin.stats.mean3$`6-7 m`+bin.stats.sd3$`6-7 m`),
              fill = cols[7],alpha = 0.2) +
  geom_point(aes(x = DOY, y = `7-8 m`), colour = cols[8], size = 2.5) +
  geom_line(aes(x = DOY, y = `7-8 m`), colour = cols[8], size = 1.2) +
  geom_ribbon(aes(x = DOY, ymin = bin.stats.mean3$`7-8 m`-bin.stats.sd3$`7-8 m`, ymax = bin.stats.mean3$`7-8 m`+bin.stats.sd3$`7-8 m`),
              fill = cols[8],alpha = 0.2) +
  geom_point(aes(x = DOY, y = `8+ m`), colour = cols[9], size = 2.5) +
  geom_line(aes(x = DOY, y = `8+ m`), colour = cols[9], size = 1.2,) +
  geom_ribbon(aes(x = DOY, ymin = bin.stats.mean3$`8+ m`-bin.stats.sd3$`8+ m`, ymax = bin.stats.mean3$`8+ m`+bin.stats.sd3$`8+ m`),
              fill = cols[9],alpha = 0.2)  +
  labs(x = 'Day of Year', y = NULL) + ylim(0.1, 0.25) +
  theme_classic()+
  theme(axis.text = element_text(size=27), axis.title=element_text(size=29)) 
chm.bottom
#***********************************************************************************************************************************************#


#********************************** plot all albedo timeseries panels together *****************************************************************#

plot_grid(pft.series.plot, chm.series.plot, pft.bottom, chm.bottom, ncol = 2, rel_heights = c(2, 1))
#***********************************************************************************************************************************************#
