library(ggplot2)
library(matrixStats)

out.dir <- "/Users/anncrumlish/Downloads/albedo analysis/pft_stats"

data.dir <- '/Users/anncrumlish/Downloads/albedo analysis/pft_stats/PFT_Albedo_Mean_and_SD.csv'
data.original <- read.csv(data.dir, header = T)
data.original[data.original > 10000] <- 10000

EVT.mean <- data.original$Spruce/10000
DTSA.mean <- data.original$Alder/10000
DTSW.mean <- data.original$Willow/10000
DTSO.mean <- data.original$OtherTall/10000
DLS.mean <- data.original$LowS/10000
EVS.mean <- data.original$EvergS/10000
DGR.mean <- data.original$DryG/10000
WGT.mean <- data.original$WetG/10000
MOS.mean <- data.original$Moss/10000
LIC.mean <- data.original$Lichen/10000

data.combn <- cbind(data.original$DOY, EVT.mean, DTSA.mean, DTSW.mean, DTSO.mean, DLS.mean, EVS.mean,
                    DGR.mean, WGT.mean, MOS.mean, LIC.mean)
data.combn <- data.frame(data.combn)
names(data.combn) <- c('DOY', 'EVT', 'DTSA', 'DTSW', 'DTSO', 'DLS', 'EVS', 'DGR', 'WGR', 'MOS', 'LIC')

doy.min <- 46
doy.max <- 286
bin.size <- 15 #changed from 15

bin.start <- doy.min
bin.end <- bin.start+bin.size
bin.stats.mean <- c()
bin.stats.sd <- c()
while (bin.end < doy.max) #calculate the 15 day average for albedo per pft over DOY
{
  bin.data <- data.combn[which(data.combn > bin.start & data.combn < bin.end), ]
  z.scores <- apply(bin.data, 2, scale)
  z.scores[is.na(z.scores)] <- 0
  bin.data[abs(z.scores) > 1] <- NaN
  bin.mean <- colMeans(bin.data, na.rm = T)
  bin.sd <- apply(bin.data, 2, sd, na.rm = TRUE)
  
  bin.name <- bin.start + 7.5
  
  bin.mean <- c(bin.name, bin.mean[-1])
  bin.stats.mean <- rbind(bin.stats.mean, bin.mean)
  
  bin.sd <- c(bin.name, bin.sd[-1])
  bin.stats.sd <- rbind(bin.stats.sd, bin.sd)
  
  bin.start <- bin.start + bin.size
  bin.end <- bin.start+bin.size
}

bin.stats.mean <- data.frame(bin.stats.mean)
names(bin.stats.mean) <- c('DOY', 'EVT', 'DTSA', 'DTSW', 'DTSO', 'DLS', 'EVS', 'DGR', 'WGR', 'MOS', 'LIC')

bin.stats.sd <- data.frame(bin.stats.sd)
names(bin.stats.sd) <- c('DOY', 'EVT', 'DTSA', 'DTSW', 'DTSO', 'DLS', 'EVS', 'DGR', 'WGR', 'MOS', 'LIC')

xstart <- 200
xend <- 220

y0 <- 1
y1 <- 1 - 0.05
y2 <- 0.95 - 0.05
y3 <- 0.9 - 0.05
y4 <- 0.85 - 0.05
y5 <- 0.8 - 0.05
y6 <- 0.75 - 0.05
y7 <- 0.7 - 0.05
y8 <- 0.65 - 0.05
y9 <- 0.60 - 0.05
y10 <- 0.55 - 0.05

x_int = 10
y_r <- 0.015

cols <- c('#FF0000', '#008000', '#669999', '#00CD00', '#FF67FF', '#0066FF', '#CC6600', '#00FF67', '#FFFFFF', '#5A5A5A')
  
#cols <- RColorBrewer::brewer.pal(10, 'Spectral')
#cols <- rev(cols)
pft.series = ggplot(data = bin.stats.mean) +
  geom_point(aes(x = DOY, y = EVT), colour = cols[1], size = 2.5) +
  geom_line(aes(x = DOY, y = EVT), colour = cols[1], size = 1.2) +
  geom_ribbon(aes(x = DOY, ymin = bin.stats.mean$EVT-bin.stats.sd$EVT, ymax = bin.stats.mean$EVT+bin.stats.sd$EVT),
    fill = cols[1],alpha = 0.1) +
  geom_point(aes(x = DOY, y = DTSA), colour = cols[2], size = 2.5) +
  geom_line(aes(x = DOY, y = DTSA), colour = cols[2], size = 1.2,) +
  geom_ribbon(aes(x = DOY, ymin = bin.stats.mean$DTSA-bin.stats.sd$DTSA, ymax = bin.stats.mean$DTSA+bin.stats.sd$DTSA),
              fill = cols[2],alpha = 0.2) +
  geom_point(aes(x = DOY, y = DTSO), colour = cols[3], size = 2.5) +
  geom_line(aes(x = DOY, y = DTSO), colour = cols[3], size = 1.2) +
  geom_ribbon(aes(x = DOY, ymin = bin.stats.mean$DTSO-bin.stats.sd$DTSO, ymax = bin.stats.mean$DTSO+bin.stats.sd$DTSO),
              fill = cols[3],alpha = 0.2) +
  geom_point(aes(x = DOY, y = DTSW), colour = cols[4], size = 2.5) +
  geom_line(aes(x = DOY, y = DTSW), colour = cols[4], size = 1.2)  +
  geom_ribbon(aes(x = DOY, ymin = bin.stats.mean$DTSW-bin.stats.sd$DTSW, ymax = bin.stats.mean$DTSW+bin.stats.sd$DTSW),
              fill = cols[4],alpha = 0.2) +
  geom_point(aes(x = DOY, y = WGR), colour = cols[5], size = 2.5) +
  geom_line(aes(x = DOY, y = WGR), colour = cols[5],  size = 1.2) +
  geom_ribbon(aes(x = DOY, ymin = bin.stats.mean$WGR-bin.stats.sd$WGR, ymax = bin.stats.mean$WGR+bin.stats.sd$WGR),
              fill = cols[5],alpha = 0.2) +
  geom_point(aes(x = DOY, y = DLS), colour = cols[6], size = 2.5) +
  geom_line(aes(x = DOY, y = DLS), colour = cols[6], size = 1.2)  +
  geom_ribbon(aes(x = DOY, ymin = bin.stats.mean$DLS-bin.stats.sd$DLS, ymax = bin.stats.mean$DLS+bin.stats.sd$DLS),
              fill = cols[6],alpha = 0.2) +
  geom_point(aes(x = DOY, y = MOS), colour = cols[7], size = 2.5) +
  geom_line(aes(x = DOY, y = MOS), colour = cols[7],  size = 1.2) +
  geom_ribbon(aes(x = DOY, ymin = bin.stats.mean$MOS-bin.stats.sd$MOS, ymax = bin.stats.mean$MOS+bin.stats.sd$MOS),
              fill = cols[7],alpha = 0.2) +
  geom_point(aes(x = DOY, y = DGR), colour = cols[8], size = 2.5) +
  geom_line(aes(x = DOY, y = DGR), colour = cols[8],  size = 1.2) +
  geom_ribbon(aes(x = DOY, ymin = bin.stats.mean$DGR-bin.stats.sd$DGR, ymax = bin.stats.mean$DGR+bin.stats.sd$DGR),
              fill = cols[8],alpha = 0.2) +
  geom_ribbon(aes(x = DOY, ymin = bin.stats.mean$LIC-bin.stats.sd$LIC, ymax = bin.stats.mean$LIC+bin.stats.sd$LIC),
              fill = "aliceblue",alpha = 0.3) +
  geom_line(aes(x = DOY, y = LIC), colour = "black",  size = 1.5) +
  geom_line(aes(x = DOY, y = LIC), colour = cols[9],  size = 1.2) +
  geom_point(aes(x = DOY, y = LIC), colour = cols[9], size = 2.5) +
  geom_point(aes(x = DOY, y = EVS), colour = cols[10], size = 2.5) +
  geom_line(aes(x = DOY, y = EVS), colour = cols[10],  size = 1.2) +
  geom_ribbon(aes(x = DOY, ymin = bin.stats.mean$EVS-bin.stats.sd$EVS, ymax = bin.stats.mean$EVS+bin.stats.sd$EVS),
              fill = cols[10],alpha = 0.2) +
  labs(x = 'Day of Year', y = 'White Sky Albedo') + ylim(0.1, 1) +
  geom_segment(aes(x = xstart, xend = xend, y = y1, yend = y1), color = cols[1], size = 1.5) +
  annotate("rect", xmin = xstart, xmax = xend, ymin = y1-y_r, ymax = y1+y_r, alpha = .2, fill = cols[1]) +
  annotate("text", x = xend + x_int,  y = y1, label = 'ET', hjust=0, size =6) +
  geom_segment(aes(x = xstart, xend = xend, y = y2, yend = y2), color = cols[2], size = 1.5) +
  annotate("rect", xmin = xstart, xmax = xend, ymin = y2-y_r, ymax = y2+y_r, alpha = .2, fill = cols[2]) +
  annotate("text", x = xend + x_int,  y = y2, label = 'DTSA', hjust=0, size =6) +
  geom_segment(aes(x = xstart, xend = xend, y = y3, yend = y3), color = cols[3], size = 1.5) +
  annotate("rect", xmin = xstart, xmax = xend, ymin = y3-y_r, ymax = y3+y_r, alpha = .2, fill = cols[3]) +
  annotate("text", x = xend + x_int,  y = y3, label = 'DT', hjust=0, size = 6) +
  geom_segment(aes(x = xstart, xend = xend, y = y4, yend = y4), color = cols[4], size = 1.5) +
  annotate("rect", xmin = xstart, xmax = xend, ymin = y4-y_r, ymax = y4+y_r, alpha = .2, fill = cols[4]) +
  annotate("text", x = xend + x_int,  y = y4, label = 'DTSW', hjust=0, size =6) +
  geom_segment(aes(x = xstart, xend = xend, y = y5, yend = y5), color = cols[5], size = 1.5) +
  annotate("rect", xmin = xstart, xmax = xend, ymin = y5-y_r, ymax = y5+y_r, alpha = .2, fill = cols[5]) +
  annotate("text", x = xend + x_int,  y = y5, label = 'WG', hjust=0, size = 6) +
  geom_segment(aes(x = xstart, xend = xend, y = y6, yend = y6), color = cols[6], size = 1.5) +
  annotate("rect", xmin = xstart, xmax = xend, ymin = y6-y_r, ymax = y6+y_r, alpha = .2, fill = cols[6]) +
  annotate("text", x = xend + x_int,  y = y6, label = 'DLS', hjust=0, size = 6) +
  geom_segment(aes(x = xstart, xend = xend, y = y7, yend = y7), color = cols[7], size = 1.5) +
  annotate("rect", xmin = xstart, xmax = xend, ymin = y7-y_r, ymax = y7+y_r, alpha = .2, fill = cols[7]) +
  annotate("text", x = xend + x_int,  y = y7, label = 'MO', hjust=0, size = 6) +
  geom_segment(aes(x = xstart, xend = xend, y = y8, yend = y8), color = cols[8], size = 1.5) +
  annotate("rect", xmin = xstart, xmax = xend, ymin = y8-y_r, ymax = y8+y_r, alpha = .2, fill = cols[8]) +
  annotate("text", x = xend + x_int,  y = y8, label = 'DG', hjust=0, size = 6) +
  geom_segment(aes(x = xstart, xend = xend, y = y9, yend = y9), color = "black", size = 1.8) +
  geom_segment(aes(x = xstart, xend = xend, y = y9, yend = y9), color = cols[9], size = 1.5) +
  annotate("rect", xmin = xstart, xmax = xend, ymin = y9-y_r, ymax = y9-0.005, alpha = .7, fill = "aliceblue") +
  annotate("rect", xmin = xstart, xmax = xend, ymin = y9+0.005, ymax = y9+y_r, alpha = .7, fill = "aliceblue") +
  annotate("text", x = xend + x_int,  y = y9, label = 'LI', hjust=0, size = 6) +
  geom_segment(aes(x = xstart, xend = xend, y = y10, yend = y10), color = cols[10], size = 1.5) +
  annotate("rect", xmin = xstart, xmax = xend, ymin = y10-y_r, ymax = y10+y_r, alpha = .2, fill = cols[10]) +
  annotate("text", x = xend + x_int,  y = y10, label = 'ES', hjust=0, size = 6)+
  annotate("text", x = xstart + 2*x_int,  y = y0, label = 'PFT', vjust=-0.2, size =6) +
  theme_classic()+
  theme(axis.text = element_text(size=27), axis.title=element_text(size=29)) 


#png.name <- paste0(out.dir, '/', 'PFT_Albedo_Time_Series.png')
#ggsave(png.name, plot = last_plot(), width = 20, height = 15, units = 'cm')


#### Timeseries by CHM######



library(terra)
library(dplyr)
library(ggplot2)
library(viridisLite)

#load in albedo timeseries
data.dir <- '/Users/anncrumlish/Downloads/albedo analysis/step4_smoothed_ts/SG_Smoothed_council_lidar_area_albedo_ts_2013_2019.dat'
albedo.raster <- rast(data.dir)

#load in CHM values over space
setwd("~/Downloads/albedo/Data")
canopy_height_model <- rast('council_watershed_chm_30m.tif') #4
canopy_height_model[canopy_height_model < 0] <- NA #there are about 6000 negative values which I think should be removed

#chm <- resample(canopy_height_model, albedo.raster)
chm = canopy_height_model
albedo.raster = resample(albedo.raster, chm) #resample albedo to preserve heterogeneity in chm
timeseries.stack <- c(albedo.raster, chm)
timeseries.df <- as.data.frame(timeseries.stack)
timeseries.df <- na.omit(timeseries.df)


timeseries.s <- subset(timeseries.df, council_watershed_chm_30m <5)
timeseries.l <- subset(timeseries.df, council_watershed_chm_30m >5)

timeseries.0 <- subset(timeseries.df, council_watershed_chm_30m <1)
timeseries.1 <- subset(timeseries.df, council_watershed_chm_30m <2 & council_watershed_chm_30m >1)
timeseries.2 <- subset(timeseries.df, council_watershed_chm_30m <3 & council_watershed_chm_30m >2)
timeseries.3 <- subset(timeseries.df, council_watershed_chm_30m <4 & council_watershed_chm_30m >3)
timeseries.4 <- subset(timeseries.df, council_watershed_chm_30m <5 & council_watershed_chm_30m >4)
timeseries.5 <- subset(timeseries.df, council_watershed_chm_30m <6 & council_watershed_chm_30m >5)
timeseries.6 <- subset(timeseries.df, council_watershed_chm_30m <7 & council_watershed_chm_30m >6)
timeseries.7 <- subset(timeseries.df, council_watershed_chm_30m <8 & council_watershed_chm_30m >7)
timeseries.8 <- subset(timeseries.df, council_watershed_chm_30m >8)



DOY <- seq(50, 296, 1)
chm.series <- cbind(as.numeric(colnames(timeseries.df[1:ncol(timeseries.df)-1])), as.numeric(colMeans(timeseries.0)[1:236]), as.numeric(colMeans(timeseries.1)[1:236]), as.numeric(colMeans(timeseries.2)[1:236]), as.numeric(colMeans(timeseries.3)[1:236]), as.numeric(colMeans(timeseries.4)[1:236]), as.numeric(colMeans(timeseries.5)[1:236]), as.numeric(colMeans(timeseries.6)[1:236]), as.numeric(colMeans(timeseries.7)[1:236]), as.numeric(colMeans(timeseries.8)[1:236]))
colnames(chm.series) <- c('DOY', '<1 m', '1-2 m', '2-3 m', '3-4 m', '4-5 m', '5-6 m', '6-7 m', '7-8 m', '8+ m')

doy.min <- 46
doy.max <- 286
bin.size <- 15 #changed from 15

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

cols <- viridis(n = 8)

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
  geom_point(aes(x = DOY, y = `8+ m`), colour = cols[8], size = 2.5) +
  geom_line(aes(x = DOY, y = `8+ m`), colour = cols[8], size = 1.2,) +
  geom_ribbon(aes(x = DOY, ymin = bin.stats.mean2$`8+ m`-bin.stats.sd2$`8+ m`, ymax = bin.stats.mean2$`8+ m`+bin.stats.sd2$`8+ m`),
              fill = cols[8],alpha = 0.2)  +
  labs(x = 'Day of Year', y = NULL) + ylim(0.1, 1) +
  annotate("text", x = xstart + 2*x_int,  y = y1b, label = 'Canopy Height', vjust=-0.2, size =6) +
  geom_segment(aes(x = xstart, xend = xend, y = y2b, yend = y2b), color = cols[1], size = 1.5) +
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
  annotate("text", x = xend + x_int,  y = y9b, label = '8+ m', hjust=0, size =6) +
  theme_classic()+
  theme(axis.text = element_text(size=27), axis.title=element_text(size=29)) 

library(cowplot)
plot_grid(pft.series, chm.series.plot, labels = c("A", "B"))






