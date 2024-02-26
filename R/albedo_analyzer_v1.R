###########################################################################################
#
#        this script does analyze albedo data and make plot
#
#    --- Last updated:  2022.04.04 By Daryl Yang <dediyang@bnl.gov>
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
list.of.packages <- c("ggplot2", "readr", "GMCM", "reshape2", "raster", "ggpmisc", 
                      "randomForest", "caTools", "pls", "spectratrait", "terra")
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

#************************************ user parameters ************************************#
outDIR <- "/Users/anncrumlish/Downloads/albedo analysis/map_stats"
### Create output folders
if (! file.exists(outDIR)) dir.create(outDIR,recursive=TRUE)

# define pft names
pfts <- c("ET","DT","DTSA","DTSW","DLS","DDS","ES","FO","DG","WG","MO","LI","NVS" )

# define outlier removal function
VIP <- function(object) {
  SS <- c(object$Yloadings)^2 * colSums(object$scores^2)
  Wnorm2 <- colSums(object$loading.weights^2)
  SSW <- sweep(object$loading.weights^2, 2, SS / Wnorm2, "*") # Replace with matrix mult.
  sqrt(nrow(SSW) * apply(SSW, 1, cumsum) / cumsum(SS))
}
#*****************************************************************************************#

#************************************** load data ****************************************#
dataDIR <- "/Users/anncrumlish/Downloads/albedo analysis/map_analysis/image_stats.csv"
dataORIG <- read.csv(dataDIR, header = TRUE)
dataORIG$PFT[dataORIG$PFT == 2] <- 100
dataORIG$PFT[dataORIG$PFT == 3] <- 101
dataORIG$PFT[dataORIG$PFT == 4] <- 2
dataORIG$PFT[dataORIG$PFT == 100] <- 3
dataORIG$PFT[dataORIG$PFT == 101] <- 4
#*****************************************************************************************#
# load in fcover map
mapDIR <- "/Users/anncrumlish/Downloads/albedo analysis/veg/council_watershed_pft_30m.dat"
mapORIG <- raster(mapDIR)
raw_cols <- c('#000000', '#FF0000', '#008000', '#00CD00', '#669999', '#0066FF', '#FFFF66', '#66FFFF', '#00B0F0', '#00FF67', '#FF67FF', '#CC6600', '#FFFFFF', '#5A5A5A')
#When I load in mapORIG, the colortable is blank. I can get the proper colors if I load mapORIG as a spatRaster, but do not know how to using raster. I have manually created a vector of the colors by converting the spatRaster coltab into hex here
#mapORIG2 <- terra::rast(mapDIR)
#coltab(mapORIG2)
cols <- raw_cols[-c(1, 6:8)]
cols <- c(cols[1], cols[4], cols[2:3], cols[-c(1:4)])
#********************************** plot pft albedo **************************************#
# extract pft winter albedo
albedoWT <- cbind(dataORIG[4:17], dataORIG[c(22)])
albedoWT[albedoWT > 10000] <- 10000
albedoWT[albedoWT < 0] <- 0
albedoWT[which(albedoWT$PFT == 0),] <- NA
albedoWT[which(albedoWT$PFT == 6),] <- NA
albedoWT[which(albedoWT$PFT == 7),] <- NA
albedoWT[which(albedoWT$PFT == 8),] <- NA
albedoWT <- na.omit(albedoWT)
# calculate mean of non tall vegetation
nonTSK <- mean(albedoWT[which(albedoWT$PFT>5),15])
# make a plot for winter albedo
ggplot(data = albedoWT, aes(x=as.factor(PFT), y=winter_albedo, fill = as.factor(PFT))) +
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

pdfNAME = paste0(outDIR, "/", "pft_albedo_winter_v1.pdf")
ggsave(pdfNAME, plot = last_plot(), width = 16, height = 13, units = 'cm')

# 
# 
# ader.data <- albedoWT[which(albedoWT$PFT == 3),15]
# will.data <- albedoWT[which(albedoWT$PFT == 4),15]
# 
# alder.ratio <- 1- mean(ader.data, na.rm = TRUE)/nonTSK
# will.ratio <- 1- mean(will.data, na.rm = TRUE)/nonTSK
# 
# data.comn <- rbind(ader.data, will.data)
# 
# low.range <- quantile(data.comn, 0.25)
# high.range <- quantile(data.comn, 0.75)
# 
# low.perc<- 1- low.range/nonTSK
# hihg.perc <- 1 - high.range/nonTSK


# extract pft summer albedo
albedoSM <- cbind(dataORIG[4:17], dataORIG[c(23)])
albedoSM[albedoSM > 10000] <- 10000
albedoSM[albedoSM < 0] <- 0
albedoSM[which(albedoSM$PFT == 0),] <- NA
albedoSM[which(albedoSM$PFT == 6),] <- NA
albedoSM[which(albedoSM$PFT == 7),] <- NA
albedoSM[which(albedoSM$PFT == 8),] <- NA
albedoSM <- na.omit(albedoSM)
# calculate mean of non tall vegetation
nonTSK <- mean(albedoSM[which(albedoSM$PFT>5),15])
# make a plot for winter albedo
ggplot(data = albedoSM, aes(x=as.factor(PFT), y=summer_albedo, fill = as.factor(PFT))) +
  geom_boxplot(outlier.shape = NA) +
  geom_hline(yintercept = nonTSK, linetype="dashed") +
  xlab("PFT") + ylab("Summer Albedo") + ylim(0, 0.22) +
  scale_fill_manual(values = cols, name = "", labels = pfts[-c(6:8)]) +
  scale_x_discrete(labels = pfts[-c(6:8)]) +
  geom_segment(aes(x = 7, xend = 8, y = 0.01, yend = 0.01),  size = 1, 
               linetype = "dashed") +
  annotate("text", x = 8.2,  y = 0.01, 
           label = stringr::str_wrap("Non-Woody Albedo", width = 10), hjust=0) +
  guides(fill=guide_legend(ncol = 2)) +
  theme(legend.position = c(0.8, 0.35), 
        legend.title = element_blank(), 
        legend.key.size = unit(0.75, 'cm'),
        legend.text = element_text(size = 12)) +
  theme(axis.text = element_text(size=13, color = 'black'),
        axis.title=element_text(size=13)) +
  theme(axis.line = element_line(colour = "black"),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

pdfNAME = paste0(outDIR, "/", "pft_albedo_summer.pdf")
ggsave(pdfNAME, plot = last_plot(), width = 16, height = 13, units = 'cm')
#*****************************************************************************************#


#******************************* plot  structure-albedo **********************************#
# extract pft winter albedo
albedoSTR <- cbind(dataORIG[3], dataORIG[c(18)], dataORIG[c(22, 23)])
albedoSTR[albedoSTR > 10] <- NA
albedoSTR[albedoSTR < 0] <- NA
albedoSTR[which(albedoSTR$PFT == 0),] <- NA
albedoSTR[which(albedoSTR$PFT == 6),] <- NA
albedoSTR[which(albedoSTR$PFT == 7),] <- NA
albedoSTR[which(albedoSTR$PFT == 8),] <- NA
albedoSTR <- na.omit(albedoSTR)
# make a structure-albedo plot for winter
my.formula <- y ~ x
ggplot(data = albedoSTR, aes(x=CHM, y=winter_albedo)) +
  #geom_point(size = 1) +
  geom_hex(aes(fill = stat(log(count))), bins = 70, 
           breaks = log(c(0, 1, 2, 4, 6))) +
  #scale_fill_gradient(low = 'purple4', high = 'yellow') +
  scale_fill_viridis_c() +
  geom_smooth(color = 'red', method = "lm", formula = my.formula) +
  #scale_color_manual(values = cols, labels = pfts[-c(6:8)]) +
  #guides(color=guide_legend(ncol = 2)) +
  stat_poly_eq(aes(label = paste(..eq.label.., sep = "~~~")), 
               label.x.npc = 6, label.y.npc = 0.95,
               eq.with.lhs = "italic(hat(y))~`=`~",
               eq.x.rhs = "~italic(x)",
               formula = my.formula, parse = TRUE, size = 4) + 
  stat_poly_eq(aes(label = paste(..rr.label.., sep = "~~~")), 
               label.x.npc = 6, label.y.npc = 0.85,
               formula = my.formula, parse = TRUE, size = 4) +
  xlab("Canopy Height (m)") + ylab("Winter Albedo") + xlim(0, 8) +
  theme(legend.position = "", 
        legend.title = element_blank(), 
        legend.key.size = unit(0.4, 'cm'),
        legend.text = element_text(size = 10),
        legend.spacing.x = unit(0.2, 'cm'),
        legend.spacing.y = unit(0.2, 'cm')) +
  theme(axis.text = element_text(size=12, color = 'black'),
        axis.title=element_text(size=12)) +
  theme(axis.line = element_line(colour = "black"),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

pdfNAME = paste0(outDIR, "/", "chm_albedo_winter.pdf")
ggsave(pdfNAME, plot = last_plot(), width = 12, height = 10, units = 'cm')

# make a structure-albedo plot for summer
ggplot(data = albedoSTR, aes(x=CHM, y=summer_albedo)) +
  #geom_point(size = 1) +
  geom_hex(aes(fill = stat(log(count))), bins = 80, 
           breaks = log(c(0, 1, 2, 4, 6))) +
  #scale_fill_gradient(low = 'purple4', high = 'yellow') +
  scale_fill_viridis_c() +
  #geom_smooth() +
  geom_smooth(color = 'red', method = "lm", formula = my.formula) +
  #scale_color_manual(values = cols, labels = pfts[-c(6:8)]) +
  #guides(color=guide_legend(ncol = 2)) +
  stat_poly_eq(aes(label = paste(..eq.label.., sep = "~~~")), 
               label.x.npc = 6, label.y.npc = 0.95,
               eq.with.lhs = "italic(hat(y))~`=`~",
               eq.x.rhs = "~italic(x)",
               formula = my.formula, parse = TRUE, size = 4) + 
  stat_poly_eq(aes(label = paste(..rr.label.., sep = "~~~")), 
               label.x.npc = 6, label.y.npc = 0.85,
               formula = my.formula, parse = TRUE, size = 4) +
  xlab("Canopy Height (m)") + ylab("Summer Albedo") + xlim(0, 8) +
  theme(legend.position = "", 
        legend.title = element_blank(), 
        legend.key.size = unit(0.4, 'cm'),
        legend.text = element_text(size = 10),
        legend.spacing.x = unit(0.2, 'cm'),
        legend.spacing.y = unit(0.2, 'cm')) +
  theme(axis.text = element_text(size=12, color = 'black'),
        axis.title=element_text(size=12)) +
  theme(axis.line = element_line(colour = "black"),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

pdfNAME = paste0(outDIR, "/", "chm_albedo_summer.pdf")
ggsave(pdfNAME, plot = last_plot(), width = 12, height = 10, units = 'cm')
#*****************************************************************************************#


#******************************* plot  pft structure *************************************#
pftSTR <- cbind(dataORIG[4:17], dataORIG[c(18)])
pftSTR[pftSTR > 20] <- NA
pftSTR[pftSTR < 0] <- NA
pftSTR[which(pftSTR$PFT == 0),] <- NA
pftSTR[which(pftSTR$PFT == 6),] <- NA
pftSTR[which(pftSTR$PFT == 7),] <- NA
pftSTR[which(pftSTR$PFT == 8),] <- NA
pftSTR <- na.omit(pftSTR)
# make a density plot for pft height
ggplot(data = pftSTR, aes(x = factor(PFT), y = CHM, fill = factor(PFT)))+
  geom_boxplot(outlier.shape = NA) +
  scale_fill_manual(values = cols, name = "", labels = pfts[-c(6:8)]) +
  scale_x_discrete(labels = pfts[-c(6:8)]) +
  guides(fill=guide_legend(ncol = 2)) +
  xlab("CHM") + ylab("Density") + ylim(0, 8) +
  theme(legend.position = c(0.7, 0.7), 
        legend.title = element_blank(), 
        legend.key.size = unit(0.75, 'cm'),
        legend.text = element_text(size = 12)) +
  theme(axis.text = element_text(size=13, color = 'black'),
        axis.title=element_text(size=13)) +
  theme(axis.line = element_line(colour = "black"),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

pdfNAME = paste0(outDIR, "/", "pft_height.pdf")
ggsave(pdfNAME, plot = last_plot(), width = 16, height = 13, units = 'cm')
#*****************************************************************************************#


#******************************* plot  fcover albedo *************************************#
fcoverALBEDO <- cbind(dataORIG[4:17], dataORIG[c(18)], dataORIG[c(22, 23)])
fcoverALBEDO[fcoverALBEDO > 15] <- NA
fcoverALBEDO[fcoverALBEDO < 0] <- NA
fcoverALBEDO[which(fcoverALBEDO$PFT == 0),] <- NA
fcoverALBEDO[which(fcoverALBEDO$PFT == 6),] <- NA
fcoverALBEDO[which(fcoverALBEDO$PFT == 7),] <- NA
fcoverALBEDO[which(fcoverALBEDO$PFT == 8),] <- NA
fcoverALBEDO <- na.omit(fcoverALBEDO)
for (pft in pfts[-c(6:8)])
{
  # detemine pdf id
  pftID <- which(pfts == pft)
  # extract pft data frame
  pftDF <- fcoverALBEDO[which(fcoverALBEDO$PFT == pftID),]
  # keep data where pft fcover > 0.5
  pftDF <- data.frame(pftDF[which(names(pftDF) == pft)], pftDF[-c(1:14)])
  names(pftDF) <- c('fcover', 'chm', 'albedoWT', 'albedoSM')
  pftDF <- na.omit(pftDF)
  # make a pft struccture scatter plot
  ggplot(data=pftDF, aes(x=fcover, y=albedoSM)) +
    geom_hex(aes(fill=stat(log(count))), bins = 50, 
             breaks = log(c(0, 1, 2, 4, 6))) +
    scale_fill_viridis_c() +
    geom_smooth(color = 'red', method = "lm", formula = my.formula) +
    stat_poly_eq(formula = my.formula, 
                 aes(label = paste(..eq.label.., ..rr.label..,..p.value.label..,
                                   sep = "*`,`~")), 
                 label.x.npc = 'right', label.y.npc = 'bottom',
                 parse = T) +
    #stat_poly_eq(aes(label = paste(..eq.label.., sep = "~~~")), 
    #             label.x.npc = 6, label.y.npc = 0.95,
    #             eq.with.lhs = "italic(hat(y))~`=`~",
    #             eq.x.rhs = "~italic(x)",
    #             formula = my.formula, parse = TRUE, size = 4) + 
    #stat_poly_eq(aes(label = paste(..rr.label.., sep = "~~~")), 
    #             label.x.npc = 6, label.y.npc = 0.85,
    #             formula = my.formula, parse = TRUE, size = 4) +
    xlab("fCover") + ylab("Summer Albedo") +
    theme(legend.position = "", 
          legend.title = element_blank(), 
          legend.key.size = unit(0.4, 'cm'),
          legend.text = element_text(size = 10),
          legend.spacing.x = unit(0.2, 'cm'),
          legend.spacing.y = unit(0.2, 'cm')) +
    theme(axis.text = element_text(size=12, color = 'black'),
          axis.title=element_text(size=12)) +
    theme(axis.line = element_line(colour = "black"),
          panel.background = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
  
  pdfNAME = paste0(outDIR, "/", pft, "_fcover_summer_albedo_regression.pdf")
  ggsave(pdfNAME, plot = last_plot(), width = 10, height = 9, units = 'cm')

  
  ggplot(data=pftDF, aes(x=chm, y=albedoSM)) +
    geom_hex(aes(fill=stat(log(count))), bins = 50, 
             breaks = log(c(0, 1, 2, 4, 6))) +
    scale_fill_viridis_c() +
    geom_smooth(color = 'red', method = "lm", formula = my.formula) +
    stat_poly_eq(formula = my.formula, 
                 aes(label = paste(..eq.label.., ..rr.label..,..p.value.label..,
                                   sep = "*`,`~")), 
                 label.x.npc = 'right', label.y.npc = 'bottom',
                 parse = T) +
    #stat_poly_eq(aes(label = paste(..eq.label.., sep = "~~~")), 
    #             label.x.npc = 6, label.y.npc = 0.95,
    #             eq.with.lhs = "italic(hat(y))~`=`~",
    #             eq.x.rhs = "~italic(x)",
    #             formula = my.formula, parse = TRUE, size = 4) + 
    #stat_poly_eq(aes(label = paste(..rr.label.., sep = "~~~")), 
    #             label.x.npc = 6, label.y.npc = 0.85,
    #             formula = my.formula, parse = TRUE, size = 4) +
    xlab("CHM") + ylab("Summer Albedo") + 
    theme(legend.position = "", 
          legend.title = element_blank(), 
          legend.key.size = unit(0.4, 'cm'),
          legend.text = element_text(size = 10),
          legend.spacing.x = unit(0.2, 'cm'),
          legend.spacing.y = unit(0.2, 'cm')) +
    theme(axis.text = element_text(size=12, color = 'black'),
          axis.title=element_text(size=12)) +
    theme(axis.line = element_line(colour = "black"),
          panel.background = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
  
  pdfNAME = paste0(outDIR, "/", pft, "_chm_summer_albedo_regression.pdf")
  ggsave(pdfNAME, plot = last_plot(), width = 10, height = 9, units = 'cm')
}
#*****************************************************************************************#

#******************************** run a PLSR analysis ************************************#
# run RF on 
#dataRF <- dataORIG[, -c(1,2, 22)]
dataRF <- dataORIG[, -c(1,2, 3, 21, 23, 24:ncol(dataORIG))]
dataRF$CHM[dataRF$CHM > 15] <- NA
dataRF$CHM[dataRF$CHM < 0] <- NA
dataRF[which(dataRF$PFT == 0),] <- NA
dataRF[which(dataRF$PFT == 6),] <- NA
dataRF[which(dataRF$PFT == 7),] <- NA
dataRF[which(dataRF$PFT == 8),] <- NA
dataRF <- dataRF[-c(7:9)]
dataRF <- na.omit(dataRF)

dataIN <- data.frame()
for (pft in pfts[-c(6:8)])
{
  # detemine pdf id
  pftID <- which(pfts == pft)
  # extract pft data
  pftDF <- dataRF[which(dataRF$PFT == pftID), ]
  if(nrow(pftDF) > 0)
  {
    sample = sample(1:nrow(pftDF), size = 300, replace = TRUE)
    pftSMPL <- pftDF[sample, ]
    dataIN = rbind(dataIN, pftSMPL)
  }
}

### plsr model on winter albedo
dataIN_WT <- na.omit(dataIN[-c(1,17)])
dataIN_WT$winter_albedo <- log(dataIN_WT$winter_albedo/(1-dataIN_WT$winter_albedo))
# create training and test samples
sample = sample.split(dataIN_WT$winter_albedo, SplitRatio = .7, group = dataIN_WT$PFT)
dataTRAIN = subset(dataIN_WT, sample == TRUE)
dataTEST  = subset(dataIN_WT, sample == FALSE)
# create training and testing samples
vips <- c()
coeffs <- c()
preds <- c()
for (i in 1:100)
{
  print(i)
  sample = sample.split(dataTRAIN$winter_albedo, SplitRatio = .7, group = dataTRAIN$PFT)
  train = subset(dataTRAIN, sample == TRUE)
  test  = subset(dataTRAIN, sample == FALSE)
  plsr.out <- plsr(winter_albedo ~ ., data = train, ncomp = 10, validation = "LOO")
  # extract vip and coefficient
  coeff <- coef(plsr.out, ncomp = 10)[,1,1]
  vip <- VIP(plsr.out)[10,]
  # store vip and coefficient values
  vips <- cbind(vips, vip)
  coeffs <- cbind(coeffs, coeff)
  
  #predict on new data
  pred <- predict(plsr.out, ncomp=10, newdata = dataTEST)
  preds <- cbind(preds, pred)
}
# calculate mean and standard deviation of vip and coefficient
vipMEAN <- apply(vips, 1, mean)
vipSD <- apply(vips, 1, sd)
coeffMEAN <- apply(coeffs, 1, mean)
#names(coeffMEAN) <- names(vipMEAN)
coeffSD <- apply(coeffs, 1, sd)
#names(coeffSD) <- names(vipMEAN)
# make a vip plot
vipPLOT <- data.frame(names(vipMEAN), as.vector(vipMEAN), 
                            as.vector(vipSD))
names(vipPLOT) <- c('pft', 'mean', 'sd')

vipPLOT$pft <- factor(vipPLOT$pft,levels = names(dataIN_WT[-c(14,15)]))
ggplot(data = vipPLOT, aes(x=pft, y=mean, fill=pft)) +
  geom_bar(stat = 'identity') +
  xlab("PFT") + ylab("VIP") + 
  theme(legend.position = "") +
  theme(axis.text = element_text(size=11, color = 'black'),
        axis.title=element_text(size=11),
        axis.text.x = element_text(angle = 45, vjust = 0.5)) +
  theme(axis.line = element_line(colour = "black"),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

pdfNAME = paste0(outDIR, "/", "winter_vip.pdf")
ggsave(pdfNAME, plot = last_plot(), width = 14, height = 7, units = 'cm')

# make a coefficient plot
coeffPLOT <- data.frame(names(coeffMEAN), as.vector(coeffMEAN), 
                      as.vector(coeffSD))
names(coeffPLOT) <- c('pft', 'mean', 'sd')
coeffPLOT$pft <- factor(coeffPLOT$pft,levels = names(dataIN_WT[-c(14,15)]))
ggplot(data = coeffPLOT, aes(x=pft, y=mean, fill=pft)) +
  geom_bar(stat = 'identity') +
  xlab("PFT") + ylab("Coefficient") +
  theme(legend.position = "") +
  theme(axis.text = element_text(size=11, color = 'black'),
        axis.title=element_text(size=11),
        axis.text.x = element_text(angle = 45, vjust = 0.5)) +
  theme(axis.line = element_line(colour = "black"),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

pdfNAME = paste0(outDIR, "/", "winter_coefficient.pdf")
ggsave(pdfNAME, plot = last_plot(), width = 14, height = 7, units = 'cm')

## plot model prediction result
predMEAN <- apply(preds, 1, mean)
predSD <- apply(preds, 1, sd)
predPLOT <- data.frame(dataTEST$winter_albedo, predMEAN, predSD)
names(predPLOT) <- c('ref', 'pred', 'unc')
ggplot(data = predPLOT, aes(x=ref, y=pred)) +
  geom_point(size = 2, fill='grey20', color='grey20') +
  geom_errorbar(aes(ymin = pred-unc, ymax = pred+unc), color='grey20') +
  geom_smooth(method = lm, formula = my.formula, color='black') +
  stat_poly_eq(formula = my.formula, 
               aes(label = paste(..eq.label.., ..rr.label..,..p.value.label..,
                                 sep = "*`,`~")), 
               label.x.npc = 'right', label.y.npc = 'bottom',
               parse = T) +
  geom_abline(intercept = 0, slope = 1, linetype='dashed', color='red') +
  xlab("Observed (logodd)") + ylab("Predicted (logodd)") + 
  xlim(-1, 4) + ylim(-1,4) +
  theme(legend.position = "") +
  theme(axis.text = element_text(size=11, color = 'black'),
        axis.title=element_text(size=11)) +
  theme(axis.line = element_line(colour = "black"),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

pdfNAME = paste0(outDIR, "/", "winter_plsr_prediction.pdf")
ggsave(pdfNAME, plot = last_plot(), width = 12, height = 12, units = 'cm')




### plsr model on summer albedo
dataIN_SM <- na.omit(dataIN[-c(1,15)])
dataIN_SM$summer_albedo <- log(dataIN_SM$summer_albedo/(1-dataIN_SM$summer_albedo))
# create training and test samples
sample = sample.split(dataIN_SM$summer_albedo, SplitRatio = .7, group = dataIN_SM$PFT)
dataTRAIN = subset(dataIN_SM, sample == TRUE)
dataTEST  = subset(dataIN_SM, sample == FALSE)
# create training and testing samples
vips <- c()
coeffs <- c()
preds <- c()
for (i in 1:100)
{
  print(i)
  sample = sample.split(dataTRAIN$summer_albedo, SplitRatio = .7, group = dataTRAIN$PFT)
  train = subset(dataTRAIN, sample == TRUE)
  test  = subset(dataTRAIN, sample == FALSE)
  plsr.out <- plsr(summer_albedo ~ ., data = train, ncomp = 10, validation = "LOO")
  # extract vip and coefficient
  coeff <- coef(plsr.out, ncomp = 10)[,1,1]
  vip <- VIP(plsr.out)[10,]
  # store vip and coefficient values
  vips <- cbind(vips, vip)
  coeffs <- cbind(coeffs, coeff)
  
  #predict on new data
  pred <- predict(plsr.out, ncomp=10, newdata = dataTEST)
  preds <- cbind(preds, pred)
}
# calculate mean and standard deviation of vip and coefficient
vipMEAN <- apply(vips, 1, mean)
vipSD <- apply(vips, 1, sd)
coeffMEAN <- apply(coeffs, 1, mean)
#names(coeffMEAN) <- names(vipMEAN)
coeffSD <- apply(coeffs, 1, sd)
#names(coeffSD) <- names(vipMEAN)
# make a vip plot
vipPLOT <- data.frame(names(vipMEAN), as.vector(vipMEAN), 
                      as.vector(vipSD))
names(vipPLOT) <- c('pft', 'mean', 'sd')

vipPLOT$pft <- factor(vipPLOT$pft,levels = names(dataIN_SM[-c(14)]))
ggplot(data = vipPLOT, aes(x=pft, y=mean, fill=pft)) +
  geom_bar(stat = 'identity') +
  xlab("PFT") + ylab("VIP") + 
  theme(legend.position = "") +
  theme(axis.text = element_text(size=11, color = 'black'),
        axis.title=element_text(size=11),
        axis.text.x = element_text(angle = 45, vjust = 0.5)) +
  theme(axis.line = element_line(colour = "black"),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

pdfNAME = paste0(outDIR, "/", "summer_vip.pdf")
ggsave(pdfNAME, plot = last_plot(), width = 14, height = 7, units = 'cm')

# make a coefficient plot
coeffPLOT <- data.frame(names(coeffMEAN), as.vector(coeffMEAN), 
                        as.vector(coeffSD))
names(coeffPLOT) <- c('pft', 'mean', 'sd')
coeffPLOT$pft <- factor(coeffPLOT$pft,levels = names(dataIN_SM[-c(14)]))
ggplot(data = coeffPLOT, aes(x=pft, y=mean, fill=pft)) +
  geom_bar(stat = 'identity') +
  xlab("PFT") + ylab("Coefficient") +
  theme(legend.position = "") +
  theme(axis.text = element_text(size=11, color = 'black'),
        axis.title=element_text(size=11),
        axis.text.x = element_text(angle = 45, vjust = 0.5)) +
  theme(axis.line = element_line(colour = "black"),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

pdfNAME = paste0(outDIR, "/", "summer_coefficient.pdf")
ggsave(pdfNAME, plot = last_plot(), width = 14, height = 7, units = 'cm')

## plot model prediction result
predMEAN <- apply(preds, 1, mean)
predSD <- apply(preds, 1, sd)
predPLOT <- data.frame(dataTEST$summer_albedo, predMEAN, predSD)
names(predPLOT) <- c('ref', 'pred', 'unc')
ggplot(data = predPLOT, aes(x=ref, y=pred)) +
  geom_point(size = 2, fill='grey20', color='grey20') +
  geom_errorbar(aes(ymin = pred-unc, ymax = pred+unc), color='grey20') +
  geom_smooth(method = lm, formula = my.formula, color='black') +
  stat_poly_eq(formula = my.formula, 
               aes(label = paste(..eq.label.., ..rr.label..,..p.value.label..,
                                 sep = "*`,`~")), 
               label.x.npc = 'right', label.y.npc = 'bottom',
               parse = T) +
  geom_abline(intercept = 0, slope = 1, linetype='dashed', color='red') +
  xlab("Observed (logodd)") + ylab("Predicted (logodd)") + 
  xlim(-2.5, -1.2) + ylim(-2.5,-1.2) +
  theme(legend.position = "") +
  theme(axis.text = element_text(size=11, color = 'black'),
        axis.title=element_text(size=11)) +
  theme(axis.line = element_line(colour = "black"),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

pdfNAME = paste0(outDIR, "/", "summer_plsr_prediction.pdf")
ggsave(pdfNAME, plot = last_plot(), width = 12, height = 12, units = 'cm')














