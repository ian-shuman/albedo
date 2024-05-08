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
bin.size <- 7 #changed from 15

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

y1 <- 1
y2 <- 0.95
y3 <- 0.9
y4 <- 0.85
y5 <- 0.8
y6 <- 0.75
y7 <- 0.7
y8 <- 0.65
y9 <- 0.60
y10 <- 0.55

x_int = 10
y_r <- 0.015

cols <- c('#FF0000', '#008000', '#669999', '#00CD00', '#FF67FF', '#0066FF', '#CC6600', '#00FF67', '#FFFFFF', '#5A5A5A')
  
#cols <- RColorBrewer::brewer.pal(10, 'Spectral')
#cols <- rev(cols)
ggplot(data = bin.stats.mean) +
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
  geom_point(aes(x = DOY, y = LIC), colour = cols[9], size = 2.5) +
  geom_line(aes(x = DOY, y = LIC), colour = cols[9],  size = 1.2) +
  geom_ribbon(aes(x = DOY, ymin = bin.stats.mean$LIC-bin.stats.sd$LIC, ymax = bin.stats.mean$LIC+bin.stats.sd$LIC),
              fill = cols[9],alpha = 0.2) +
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
  geom_segment(aes(x = xstart, xend = xend, y = y9, yend = y9), color = cols[9], size = 1.5) +
  annotate("rect", xmin = xstart, xmax = xend, ymin = y9-y_r, ymax = y9+y_r, alpha = .2, fill = cols[9]) +
  annotate("text", x = xend + x_int,  y = y9, label = 'LI', hjust=0, size = 6) +
  geom_segment(aes(x = xstart, xend = xend, y = y10, yend = y10), color = cols[10], size = 1.5) +
  annotate("rect", xmin = xstart, xmax = xend, ymin = y10-y_r, ymax = y10+y_r, alpha = .2, fill = cols[10]) +
  annotate("text", x = xend + x_int,  y = y10, label = 'ES', hjust=0, size = 6)+
  theme_classic()+
  theme(axis.text = element_text(size=27), axis.title=element_text(size=29)) 


png.name <- paste0(out.dir, '/', 'PFT_Albedo_Time_Series.png')
ggsave(png.name, plot = last_plot(), width = 20, height = 15, units = 'cm')

















