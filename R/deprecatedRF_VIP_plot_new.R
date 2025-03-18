#Follow script from albedo_analyzer_v1.R to set up, this is a temporary document
rm(list=ls(all=TRUE))   # clear workspace
graphics.off()          # close any open graphics
closeAllConnections()   # close any open connections to files
dlm <- .Platform$file.sep # <--- What is the platform specific delimiter?
options(warn = -1)
#*****************************************************************************************#

#****************************** load required libraries **********************************#
### install and load required R packages
list.of.packages <- c("ggplot2", "readr", "GMCM", "reshape2", "raster", "ggpmisc", 
                      "randomForest", "caTools", "pls", "spectratrait", "terra", "pdp")
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
my.formula <- y ~ x
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
colors <- c("#FF0000", "#008000", "#00CD00", "#669999", "#0066FF", "#00FF67", "#FF67FF", "#CC6600", "#FFFFFF", "#5A5A5A", "lightgrey", "lightgrey", "lightgrey", "lightgrey")
#When I load in mapORIG, the colortable is blank. I can get the proper colors if I load mapORIG as a spatRaster, but do not know how to using raster. I have manually created a vector of the colors by converting the spatRaster coltab into hex here
#mapORIG2 <- terra::rast(mapDIR)
#coltab(mapORIG2)
cols <- raw_cols[-c(1, 6:8)]
cols <- c(cols[1], cols[4], cols[2:3], cols[-c(1:4)])

# run RF on 
#dataRF <- dataORIG[, -c(1,2, 22)]
dataRF <- dataORIG[, -c(1,2, 3, 21, 24:ncol(dataORIG))]
dataRF <- df2_subset[, -c(1,2, 3, 26:ncol(df2_subset)-2)]
dataRF$CHM[dataRF$CHM > 15] <- NA
dataRF$CHM[dataRF$CHM < 0] <- NA
dataRF[which(dataRF$PFT == 0),] <- NA
dataRF[which(dataRF$PFT == 6),] <- NA
dataRF[which(dataRF$PFT == 7),] <- NA
dataRF[which(dataRF$PFT == 8),] <- NA
dataRF <- dataRF[-c(7:9)]
dataRF <- na.omit(dataRF)

dataIN <- data.frame()
for (pft in pfts[-c(6:8)]) #randomly sample from the council site to create a dataframe which contains all pfts (except DDS, ES, and FO)
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

### RF model on winter albedo
dataIN_WT <- na.omit(dataIN[-c(1,16)])
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
  RF.out <- randomForest(winter_albedo ~ ., data = train, importance = T)
  # extract vip and coefficient
  imp <- importance(RF.out, type = 1)
  # store vip and coefficient values
  vips <- cbind(vips, imp)
  
  #predict on new data
  pred <- predict(RF.out, ncomp=10, newdata = dataTEST)
  preds <- cbind(preds, pred)
}

vip.winter <- vips
vipMEAN.winter <- apply(vip.winter, 1, mean)
vipSD.winter <- apply(vip.winter, 1, sd)
#coeffMEAN <- apply(coeffs, 1, mean)
#names(coeffMEAN) <- names(vipMEAN)
#coeffSD <- apply(coeffs, 1, sd)
#names(coeffSD) <- names(vipMEAN)
# make a vip plot
vipPLOT.winter <- data.frame(names(vipMEAN.winter), as.vector(vipMEAN.winter), 
                             as.vector(vipSD.winter))
names(vipPLOT.winter) <- c('pft', 'mean', 'sd')

vipPLOT.winter$pft <- factor(vipPLOT.winter$pft,levels = names(dataIN_WT[-c(14,15)]))
xfactor.w <- rownames(vip.winter)
xfactor.w <- factor(xfactor.w, levels = names(dataIN_WT[-c(14,15)]))
fill_scale_values <- setNames(colors, levels(vipPLOT.winter$pft))
bar.winter <- ggplot(data = vipPLOT.winter, aes(x=xfactor.w, y=mean, fill=pft)) +
  geom_bar(stat = 'identity', color = "black") +
  scale_fill_manual(values = fill_scale_values) +
  xlab("Predictor Variable") + ylab("VIP") + ggtitle("Winter")+
  theme(legend.position = "") +
  theme(axis.text = element_text(size=19, color = 'black'),
        axis.title=element_text(size=21),
        axis.text.x = element_text(angle = 45, vjust = 0.5), 
        plot.title = element_text(face = "bold", size = 25)) +
  theme(axis.line = element_line(colour = "black"),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  scale_y_continuous(breaks = seq(0, 50, by = 10))+
  geom_errorbar(aes(ymin=mean-2*sd, ymax=mean+2*sd), width=.2,
                position=position_dodge(.9))


pdfNAME = paste0(outDIR, "/", "RF_winter_vip_2sd_poster.pdf")
ggsave(pdfNAME, plot = last_plot(), width = 14, height = 7, units = 'cm')
#why was aspect taken out? when you include it, it's very high in importance!

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

pdfNAME = paste0(outDIR, "/", "winter_RF_prediction.pdf")
ggsave(pdfNAME, plot = last_plot(), width = 12, height = 12, units = 'cm')




### RF model on summer albedo
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
  RF.out <- randomForest(summer_albedo ~ ., data = train, importance = T)
  # extract vip and coefficient
  imp <- importance(RF.out, type = 1)
  # store vip and coefficient values
  vips <- cbind(vips, imp)
  
  #predict on new data
  pred <- predict(RF.out, ncomp=10, newdata = dataTEST)
  preds <- cbind(preds, pred)
}

vips.summer <- vips
vipMEAN.summer <- apply(vips.summer, 1, mean)
vipSD.summer <- apply(vips.summer, 1, sd)
#coeffMEAN <- apply(coeffs, 1, mean)
#names(coeffMEAN) <- names(vipMEAN)
#coeffSD <- apply(coeffs, 1, sd)
#names(coeffSD) <- names(vipMEAN)
# make a vip plot
vipPLOT.summer <- data.frame(names(vipMEAN.summer), as.vector(vipMEAN.summer), 
                             as.vector(vipSD.summer))
names(vipPLOT.summer) <- c('pft', 'mean', 'sd')


# Create a named vector for the fill scale
fill_scale_values <- setNames(colors, levels(vipPLOT.summer$pft))


vipPLOT.summer$pft <- factor(vipPLOT.summer$pft,levels = names(dataIN_WT[-c(14,15)]))
xfactor <- rownames(vips.summer)
xfactor <- factor(xfactor, levels = names(dataIN_WT[-c(14,15)]))
bar.summer <- ggplot(data = vipPLOT.summer, aes(x=xfactor, y=mean, fill=pft)) +
  geom_bar(stat = 'identity', color = "black") +
  scale_fill_manual(values = fill_scale_values) +
  xlab("Predictor Variable") + ylab("VIP") + ggtitle("Summer")+
  theme(legend.position = "") +
  theme(axis.text = element_text(size=19, color = 'black'),
        axis.title=element_text(size=21),
        axis.text.x = element_text(angle = 45, vjust = 0.5), 
        plot.title = element_text(face = "bold", size = 25)) +
  theme(axis.line = element_line(colour = "black"),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  geom_errorbar(aes(ymin=mean-2*sd, ymax=mean+2*sd), width=.2,
                position=position_dodge(.9))


pdfNAME = paste0(outDIR, "/", "RF_summer_vip_2SD_poster.pdf")
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
  #xlim(-1, 4) + ylim(-1,4) +
  theme(legend.position = "") +
  theme(axis.text = element_text(size=11, color = 'black'),
        axis.title=element_text(size=11)) +
  theme(axis.line = element_line(colour = "black"),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

pdfNAME = paste0(outDIR, "/", "summer_RF_prediction.pdf")
ggsave(pdfNAME, plot = last_plot(), width = 12, height = 12, units = 'cm')

plot_grid(bar.summer, bar.winter, nrow = 2)
plot_grid(
  bar.summer,
  ggdraw() + theme_void(),  # Empty plot for the line
  bar.winter,
  ncol = 1,
  rel_heights = c(1, 0.05, 1)  # Adjust this to change the thickness of the line
)

# Add the line
grid::grid.lines(x = c(0, 1), y = c(0.5, 0.5), gp = grid::gpar(col = "black"))

#look at the differences of winter - summer, positive means more important in the winter
vip.differences = vipPLOT.winter - vipPLOT.summer
vip.differences$pft <- vipPLOT.winter$pft
vip.differences$sd <- sqrt((((vipPLOT.winter$sd)^2)/100)+(((vipPLOT.summer$sd)^2)/100))
View(vip.differences)

ggplot(data = vip.differences, aes(x=pft, y=mean, fill=pft)) +
  geom_bar(stat = 'identity') +
  xlab("PFT") + ylab("VIP") + 
  theme(legend.position = "") +
  theme(axis.text = element_text(size=11, color = 'black'),
        axis.title=element_text(size=11),
        axis.text.x = element_text(angle = 45, vjust = 0.5)) +
  theme(axis.line = element_line(colour = "black"),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  geom_errorbar(aes(ymin=mean-2*sd, ymax=mean+2*sd), width=.2,
                position=position_dodge(.9))



#look at partial dependence plots too?

RF_SM <- randomForest(summer_albedo ~ ET + DTSA + DTSW + DT + DLS + DG + WG + MO + LI + NVS + CHM + DEM + slope, data = dataRF, importance = T)
RF_WT <- randomForest(winter_albedo ~ ET + DTSA + DTSW + DT + DLS + DG + WG + MO + LI + NVS + CHM + DEM + slope, data = dataRF, importance = T)

imp_SM <- as.data.frame(importance(RF_SM, type = 1))
imp_WT <- as.data.frame(importance(RF_WT, type = 1))
x <- rownames(imp_SM)
imp_SM_df <- cbind(x, imp_SM)
ax <- rownames(imp_WT)
imp_WT_df <- cbind(x, imp_WT)


library(forcats)
ggplot(data = imp_SM_df, aes(x = x, y = `%IncMSE`, fill = x, color = x)) +
  geom_bar(stat = 'identity') +
  xlab("PFT") + ylab("VIP") + 
  theme(legend.position = "") +
  theme(axis.text = element_text(size=11, color = 'black'),
        axis.title=element_text(size=11),
        axis.text.x = element_text(angle = 45, vjust = 0.5)) +
  theme(axis.line = element_line(colour = "black"),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  scale_fill_manual(values = c("#FF0000", "#669999", "#008000" ,"#00CD00", "#00B0F0" ,"#00FF67", "#FF67FF", "#CC6600",
                               "#FFFFFF", "#5A5A5A", "black", "black", "black"))


##Code from biometry which demonstrates how to use factor to properly order the columns
treatment <- c("Control", "Low-level of infection", "Intermediate-level of infection", "Heavy-level of infection")
abnormalities <- c(0, 22, 17, 18)
examined <- c(31, 31, 18, 18)
prop <- abnormalities/examined
lower <- prop - qnorm(0.975)* sqrt( (prop * (1 - prop))/ examined)
upper <- prop + qnorm(0.975)* sqrt( (prop * (1 - prop))/ examined)
frogs.df <- as.data.frame(cbind(treatment, abnormalities, examined, prop, lower, upper))

library(ggplot2)
library(forcats)
ggplot(frogs.df, aes(x=fct_reorder(treatment, prop), y=prop)) + 
  geom_bar(stat = "identity")+
  geom_errorbar(data=frogs.df,aes(x=treatment,ymin=lower,ymax=upper))+
  labs(x = "Treatment", y = "Proportion of Abnormalities")+
  theme(axis.text.x = element_text(angle = 30, vjust = 0.5))

####DO IT WITH MY DATASET
rm(list=ls(all=TRUE))
setwd("~/Downloads/albedo/Data")
load("stack_df.RData")
library(pdp)

colnames(stack_df) <- c("x", "y", "CHM", "DEM", "TPI", "TRI", "slope", "aspect", "TWI", "summer_albedo", "winter_albedo", "ET", "DTSA", "DTSW", "DT", "DLS", "DDS", "ES", "FO", "DG", "WG", "MO", "LI", "NPV", "category", "bp1", "bp2")
my.formula <- y ~ x
stack_df[which(stack_df$category == 0),] <- NA
stack_df[which(stack_df$category == "DDS"),] <- NA
stack_df[which(stack_df$category == "EvergS"),] <- NA
stack_df[which(stack_df$category == "Forb"),] <- NA
cols_to_fill <- c("ET", "DTSA", "DTSW", "DT", "DLS", "DDS", "ES", "FO", "DG", "WG", "MO", "LI", "NPV")
colors <- c("#FF0000", "#008000", "#00CD00", "#669999", "#0066FF", "#00FF67", "#FF67FF", "#CC6600", "#FFFFFF", "#5A5A5A", "lightgrey", "lightgrey", "lightgrey", "lightgrey", "lightgrey", "lightgrey", "lightgrey", "lightgrey")



for (col in cols_to_fill) {
  stack_df[is.na(stack_df[[col]]), col] <- 0
}
sub_df <- subset(stack_df, select = -c(x, y, DDS, ES, FO, bp1, bp2))
sub_df$category2 <- sapply(as.character(sub_df$category), function(x) if (x %in% c('nonClass', 'NPV')) 'NVS' else x)
sub_df$CHM[sub_df$CHM > 15] <- NA
sub_df$CHM[sub_df$CHM < 0] <- NA
sub_df <- na.omit(sub_df)
sub_df <- subset(sub_df, select = -c(category))

pfts = unique(sub_df$category2)
dataIN <- data.frame()
for (pft in pfts) #randomly sample from the council site to create a dataframe which contains all pfts (except DDS, ES, and FO)
{
  # detemine pdf id
  pftID <- which(pfts == pft)
  # extract pft data
  pftDF <- sub_df[which(sub_df$category2 == pft), ]
  if(nrow(pftDF) > 0)
  {
    sample = sample(1:nrow(pftDF), size = 300, replace = TRUE)
    pftSMPL <- pftDF[sample, ]
    dataIN = rbind(dataIN, pftSMPL)
  }
}

dataIN_WT <- subset(dataIN, select = -c(`summer_albedo`, `category2`))

dataIN_WT$winter_albedo <- log(dataIN_WT$winter_albedo/(1-dataIN_WT$winter_albedo))
# create training and test samples
sample = sample.split(dataIN_WT$winter_albedo, SplitRatio = .7, group = dataIN_WT$PFT)
dataTRAIN = subset(dataIN_WT, sample == TRUE)
dataTEST  = subset(dataIN_WT, sample == FALSE)
# create training and testing samples
vips <- c()
coeffs <- c()
preds <- c()
outputs.w <- vector("list", length = 100)
for (i in 1:100)
{
  print(i)
  sample = sample.split(dataTRAIN$winter_albedo, SplitRatio = .7, group = dataTRAIN$PFT)
  train = subset(dataTRAIN, sample == TRUE)
  test  = subset(dataTRAIN, sample == FALSE)
  RF.out <- randomForest(winter_albedo ~ ., data = train, importance = T)
  
  # extract vip and coefficient
  imp <- importance(RF.out, type = 1)
  # store vip and coefficient values
  vips <- cbind(vips, imp)
  outputs.w[[i]] <- RF.out
  #predict on new data
  pred <- predict(RF.out, ncomp=10, newdata = dataTEST)
  preds <- cbind(preds, pred)
}

vip.winter <- vips
vipMEAN.winter <- apply(vip.winter, 1, mean)
vipSD.winter <- apply(vip.winter, 1, sd)
# make a vip plot
vipPLOT.winter <- data.frame(names(vipMEAN.winter), as.vector(vipMEAN.winter), 
                             as.vector(vipSD.winter))
names(vipPLOT.winter) <- c('pft', 'mean', 'sd')
rownames(vipPLOT.winter) <- unique(vipPLOT.winter$pft)
new_order <- c("ET", "DTSA", "DTSW", "DT", "DLS", "DG", "WG", "MO", "LI", "NPV", "CHM", "DEM", "slope", "aspect", "TPI", "TRI", "TWI")
row_order <- match(new_order, vipPLOT.winter$pft)
vipPLOT.winter <- vipPLOT.winter[row_order,]
vipPLOT.winter$pattern_group <- ifelse(vipPLOT.winter$pft %in% c("ET", "DT", "DTSA", "DTSW"), "Woody", "Non-woody")
xfactor.w <- new_order
xfactor.w <- factor(xfactor.w, levels = unique(vipPLOT.winter$pft))
fill_scale_values <- setNames(colors, levels(xfactor.w))
bar.winter <- ggplot(data = vipPLOT.winter, aes(x=xfactor.w, y=mean, fill=pft, pattern = pattern_group)) +
  geom_bar(stat = 'identity', color = "black") +
  geom_col_pattern(pattern_density = 0.1, pattern_fill = "black", pattern_angle = 45, fill = NA, color = "black") +
  scale_pattern_manual(values = c("Woody" = "stripe", "Non-woody" = NA), name = "Pattern") +
  xlab("PFT") + ylab("VIP") + ggtitle("Winter")+
  scale_fill_manual(values = fill_scale_values)+
  theme(legend.position = "") +
  theme(axis.text = element_text(size=19, color = 'black'),
        axis.title=element_text(size=21),
        axis.text.x = element_text(angle = 45, vjust = 0.5), 
        plot.title = element_text(face = "bold", size = 25)) +
  theme(axis.line = element_line(colour = "black"),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  scale_y_continuous(breaks = seq(0, 70, by = 10))+
  geom_errorbar(aes(ymin=mean-2*sd, ymax=mean+2*sd), width=.2,
                position=position_dodge(.9))
bar.winter 
predMEAN.w <- apply(preds, 1, mean)
predSD.w <- apply(preds, 1, sd)
predPLOT.w <- data.frame(dataTEST$winter_albedo, predMEAN.w, predSD.w)
names(predPLOT.w) <- c('ref', 'pred', 'unc')
pred.winter = ggplot(data = predPLOT.w, aes(x=ref, y=pred)) +
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
pred.winter
#Make the partial dependence plots
variables_of_interest <- c("CHM", "DEM", "TPI", "TRI", "slope", "aspect", "TWI", "ET", "DTSA", "DTSW", "DT", "DLS", "DG", "WG", "MO", "LI", "NPV")
for (v in variables_of_interest){
  variable_of_interest <- v
  yhats.w <- vector("list", length = length(outputs.w))
  for (i in 1:length(outputs.w)){
    pd <- pdp::partial(outputs.w[[i]], pred.var = v)
    a = pd$yhat
    b = exp(a) / (1 + exp(a)) #need to reverse the log odds that was taken of albedo for RF prediction
    yhats.w[[i]] <- b
  }
  yhat.all.w <- pd[,1]
  for (i in 1:length(yhats.w)){
    object <- yhats.w[i]
    yhat.all.w <- cbind(yhat.all.w, unlist(object))  
  }
  pd.all.w <- as.data.frame(yhat.all.w)
  pd.mean.w <- cbind(pd.all.w$yhat.all.w, rowMeans(yhat.all.w[,-1]))
  pd.final.w <- cbind(pd.mean.w, apply(yhat.all.w[,-1], 1, sd, na.rm = TRUE))
  colnames(pd.final.w) <- c(variable_of_interest, "Mean", "SD")
  pd.final <- as.data.frame(pd.final.w)
  plot.name <- paste0("pdp.winter.", v)
  plot <- ggplot(pd.final.w, aes(x = pd.final.w[,1], y = Mean)) +
    geom_line() +
    labs(x = v)+
    geom_ribbon(aes(ymin = Mean - SD, ymax = Mean + SD), alpha = 0.2)
  assign(plot.name, plot)
}








######### My data summer ################


dataIN_SM <- subset(dataIN, select = -c(`winter_albedo`, `category2`))

dataIN_SM$summer_albedo <- log(dataIN_SM$summer_albedo/(1-dataIN_SM$summer_albedo))
# create training and test samples
sample = sample.split(dataIN_SM$summer_albedo, SplitRatio = .7, group = dataIN_SM$PFT)
dataTRAIN = subset(dataIN_SM, sample == TRUE)
dataTEST  = subset(dataIN_SM, sample == FALSE)
# create training and testing samples
vips <- c()
coeffs <- c()
preds <- c()
outputs <- vector("list", length = 100)
for (i in 1:100)
{
  print(i)
  sample = sample.split(dataTRAIN$summer_albedo, SplitRatio = .7, group = dataTRAIN$PFT)
  train = subset(dataTRAIN, sample == TRUE)
  test  = subset(dataTRAIN, sample == FALSE)
  RF.out <- randomForest(summer_albedo ~ ., data = train, importance = T)
  outputs[[i]] <- RF.out
  # extract vip and coefficient
  imp <- importance(RF.out, type = 1)
  # store vip and coefficient values
  vips <- cbind(vips, imp)
  
  
  
  #predict on new data
  pred <- predict(RF.out, ncomp=10, newdata = dataTEST)
  preds <- cbind(preds, pred)
}

vip.summer <- vips
vipMEAN.summer <- apply(vip.summer, 1, mean)
vipSD.summer <- apply(vip.summer, 1, sd)
# make a vip plot
vipPLOT.summer <- data.frame(names(vipMEAN.summer), as.vector(vipMEAN.summer), 
                             as.vector(vipSD.summer))
names(vipPLOT.summer) <- c('pft', 'mean', 'sd')
rownames(vipPLOT.summer) <- unique(vipPLOT.summer$pft)
new_order <- c("ET", "DTSA", "DTSW", "DT", "DLS", "DG", "WG", "MO", "LI", "NPV", "CHM", "DEM", "slope", "aspect", "TPI", "TRI", "TWI")
row_order <- match(new_order, vipPLOT.summer$pft)
vipPLOT.summer <- vipPLOT.summer[row_order,]
vipPLOT.summer$pattern_group <- ifelse(vipPLOT.summer$pft %in% c("ET", "DT", "DTSA", "DTSW"), "Woody", "Non-woody")
xfactor.s <- new_order
xfactor.s <- factor(xfactor.s, levels = unique(vipPLOT.summer$pft))
fill_scale_values <- setNames(colors, levels(xfactor.s))
bar.summer <- ggplot(data = vipPLOT.summer, aes(x=xfactor.s, y=mean, fill=pft, pattern = pattern_group)) +
  geom_bar(stat = 'identity', color = "black") +
  geom_col_pattern(pattern_density = 0.1, pattern_fill = "black", pattern_angle = 45, fill = NA, color = "black") +
  scale_pattern_manual(values = c("Woody" = "stripe", "Non-woody" = NA), name = "Pattern") +
  xlab("PFT") + ylab("VIP") + ggtitle("Summer")+
  scale_fill_manual(values = fill_scale_values)+
  theme(legend.position = "") +
  theme(axis.text = element_text(size=19, color = 'black'),
        axis.title=element_text(size=21),
        axis.text.x = element_text(angle = 45, vjust = 0.5), 
        plot.title = element_text(face = "bold", size = 25)) +
  theme(axis.line = element_line(colour = "black"),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  scale_y_continuous(breaks = seq(0, 70, by = 10))+
  geom_errorbar(aes(ymin=mean-2*sd, ymax=mean+2*sd), width=.2,
                position=position_dodge(.9))
bar.summer
predMEAN.s <- apply(preds, 1, mean)
predSD.s <- apply(preds, 1, sd)
predPLOT.s <- data.frame(dataTEST$summer_albedo, predMEAN.s, predSD.s)
names(predPLOT.s) <- c('ref', 'pred', 'unc')
pred.summer = ggplot(data = predPLOT.s, aes(x=ref, y=pred)) +
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
  #xlim(-1, 4) + ylim(-1,4) +
  theme(legend.position = "") +
  theme(axis.text = element_text(size=11, color = 'black'),
        axis.title=element_text(size=11)) +
  theme(axis.line = element_line(colour = "black"),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
pred.summer
#Make the partial dependence plots
variables_of_interest <- c("CHM", "DEM", "TPI", "TRI", "slope", "aspect", "TWI", "ET", "DTSA", "DTSW", "DT", "DLS", "DG", "WG", "MO", "LI", "NPV") #c("CHM", "ET", "DG", "aspect", "TRI", "TWI")
for (v in variables_of_interest){
  variable_of_interest <- v
  yhats <- vector("list", length = length(outputs))
  for (i in 1:length(outputs.w)){
    pd <- pdp::partial(outputs[[i]], pred.var = v)
    a = pd$yhat
    b = exp(a) / (1 + exp(a)) #need to reverse the log odds that was taken of albedo for RF prediction
    yhats[[i]] <- b
  }
  yhat.all <- pd[,1]
  for (i in 1:length(yhats)){
    object <- yhats[i]
    yhat.all <- cbind(yhat.all, unlist(object))  
  }
  pd.all <- as.data.frame(yhat.all)
  pd.mean <- cbind(pd.all$yhat.all, rowMeans(yhat.all[,-1]))
  pd.final <- cbind(pd.mean, apply(yhat.all[,-1], 1, sd, na.rm = TRUE))
  colnames(pd.final) <- c(variable_of_interest, "Mean", "SD")
  pd.final <- as.data.frame(pd.final)
  plot.name <- paste0("pdp.summer.", v)
  plot <- ggplot(pd.final, aes(x = pd.final[,1], y = Mean)) +
    geom_line() +
    labs(x = v)+
    geom_ribbon(aes(ymin = Mean - SD, ymax = Mean + SD), alpha = 0.2)
  assign(plot.name, plot)
}



library(ggbreak)
library(ggplot2)
library(tidyr)

#Define a function to split and sort a dataframe so that it can be plotted with geom_path
splitsort <- function(df, var_col, mean_summer_col, mean_winter_col, sd_summer_col, sd_winter_col) {
  result_df <- data.frame()
  for (i in 1:nrow(df)) {
    row <- df[i, ]
    if (!is.na(row[[mean_summer_col]]) & !is.na(row[[mean_winter_col]])) {
      # Create two new rows, one for summer and one for winter
      summer_row <- row
      winter_row <- row
      summer_row[[mean_winter_col]] <- NA
      summer_row[[sd_winter_col]] <- NA
      winter_row[[mean_summer_col]] <- NA
      winter_row[[sd_summer_col]] <- NA
      result_df <- rbind(result_df, summer_row, winter_row)
    } else {
      result_df <- rbind(result_df, row)
    }
  }
  # Sort the resulting data frame
  result_df <- result_df %>%
    arrange(
      desc(!is.na(result_df[[mean_summer_col]])), 
      result_df[[var_col]]
    )
  return(result_df)
}

#Plot the Summer and Winter PDPs both using a break and using separate Y axes in the same graph

pdp.all.aspect.df <- rbind(pdp.summer.aspect$data, pdp.winter.aspect$data)
pdp.all.aspect.df$Model <- c(rep("Summer", nrow(pdp.summer.aspect$data)), rep("Winter", nrow(pdp.winter.aspect$data)))
pdp.all.aspect <- ggplot(pdp.all.aspect.df, aes(x = aspect, y = Mean, fill = Model))+
  geom_line(aes(color = Model)) +
  labs(x = "Aspect (°)", y = "Mean Seasonal Albedo")+ 
  geom_ribbon(aes(ymin = Mean - SD, ymax = Mean + SD, fill = Model), alpha = 0.2)+
  scale_y_break(breaks = c(0.175, 0.84), scale = "fixed")+
  theme_classic()+
  theme(legend.position = "none")

coeff <- 5.5
pdp.all.aspect.w <- pivot_wider(data = pdp.all.aspect.df, names_from = Model, values_from = c(Mean, SD), names_sep = "_")
pdp.all.aspect.w <- splitsort(pdp.all.aspect.w, "aspect", "Mean_Summer", "Mean_Winter", "SD_Summer", "SD_Winter")
pdpw.aspect <- ggplot(pdp.all.aspect.w, aes(x = aspect)) +
  geom_path(aes(y=Mean_Winter, color = "Winter"), na.rm = T) + 
  geom_path( aes(y=Mean_Summer*coeff, color = "Summer"), na.rm = T) + 
  scale_y_continuous(name = "Winter Albedo", sec.axis = sec_axis(~./coeff, name="Summer Albedo"))+
  scale_color_manual(values = c("red", "blue"))+
  scale_fill_manual(values = c("red", "blue"))+
  labs(x = "Aspect (°)")+
  geom_ribbon(aes(ymin = Mean_Winter - SD_Winter, ymax = Mean_Winter + SD_Winter, fill = "Winter"), alpha = 0.2, na.rm = T)+
  geom_ribbon(aes(ymin = Mean_Summer*coeff - SD_Summer*coeff, ymax = Mean_Summer*coeff + SD_Summer*coeff, fill = "Summer"), alpha = 0.2, na.rm = T)+
  labs(fill = "Model", color = "Model")+
  theme_classic()+
  theme(legend.position = "none", 
        axis.text = element_text(size=15, color = 'black'),
        axis.title=element_text(size=18))


pdp.all.CHM.df <- rbind(pdp.summer.CHM$data, pdp.winter.CHM$data)
pdp.all.CHM.df$Model <- c(rep("Summer", nrow(pdp.summer.CHM$data)), rep("Winter", nrow(pdp.winter.CHM$data)))
pdp.all.CHM <- ggplot(pdp.all.CHM.df, aes(x = CHM, y = Mean, fill = Model))+
  geom_line(aes(color = Model)) +
  labs(x = "Canopy Height (m)", y = "Mean Seasonal Albedo")+ 
  geom_ribbon(aes(ymin = Mean - SD, ymax = Mean + SD, fill = Model), alpha = 0.2)+
  scale_y_break(breaks = c(0.2, 0.75), scale = "fixed")+
  theme_classic()+
  theme(legend.position = "none")

coeff <- 5.5
pdp.all.CHM.w <- pivot_wider(data = pdp.all.CHM.df, names_from = Model, values_from = c(Mean, SD), names_sep = "_")
pdp.all.CHM.w <- splitsort(pdp.all.CHM.w, "CHM", "Mean_Summer", "Mean_Winter", "SD_Summer", "SD_Winter")
pdpw.CHM <- ggplot(pdp.all.CHM.w, aes(x = CHM)) +
  geom_path(aes(y=Mean_Winter, color = "Winter"), na.rm = T) + 
  geom_path( aes(y=Mean_Summer*coeff, color = "Summer"), na.rm = T) + 
  scale_y_continuous(name = "Winter Albedo", sec.axis = sec_axis(~./coeff, name="Summer Albedo"))+
  scale_color_manual(values = c("red", "blue"))+
  scale_fill_manual(values = c("red", "blue"))+
  labs(x = "Canopy Height (m)")+
  geom_ribbon(aes(ymin = Mean_Winter - SD_Winter, ymax = Mean_Winter + SD_Winter, fill = "Winter"), alpha = 0.2, na.rm = T)+
  geom_ribbon(aes(ymin = Mean_Summer*coeff - SD_Summer*coeff, ymax = Mean_Summer*coeff + SD_Summer*coeff, fill = "Summer"), alpha = 0.2, na.rm = T)+
  labs(fill = "Model", color = "Model")+
  theme_classic()+
  theme(legend.position = "none", 
        axis.text = element_text(size=15, color = 'black'),
        axis.title=element_text(size=18))


pdp.all.DEM.df <- rbind(pdp.summer.DEM$data, pdp.winter.DEM$data)
pdp.all.DEM.df$Model <- c(rep("Summer", nrow(pdp.summer.DEM$data)), rep("Winter", nrow(pdp.winter.DEM$data)))
pdp.all.DEM <- ggplot(pdp.all.DEM.df, aes(x = DEM, y = Mean, fill = Model))+
  geom_line(aes(color = Model)) +
  labs(x = "Elevation (m)", y = NULL)+ 
  geom_ribbon(aes(ymin = Mean - SD, ymax = Mean + SD, fill = Model), alpha = 0.2)+
  scale_y_break(breaks = c(0.2, 0.75), scale = "fixed")+
  theme_classic()+
  theme(legend.position = "none")

coeff <- 5.5
pdp.all.DEM.w <- pivot_wider(data = pdp.all.DEM.df, names_from = Model, values_from = c(Mean, SD), names_sep = "_")
pdp.all.DEM.w <- splitsort(pdp.all.DEM.w, "DEM", "Mean_Summer", "Mean_Winter", "SD_Summer", "SD_Winter")
pdpw.DEM <- ggplot(pdp.all.DEM.w, aes(x = DEM)) +
  geom_path(aes(y=Mean_Winter, color = "Winter"), na.rm = T) + 
  geom_path( aes(y=Mean_Summer*coeff, color = "Summer"), na.rm = T) + 
  scale_y_continuous(name = "Winter Albedo", sec.axis = sec_axis(~./coeff, name="Summer Albedo"))+
  scale_color_manual(values = c("red", "blue"))+
  scale_fill_manual(values = c("red", "blue"))+
  labs(x = "Elevation (m)")+
  geom_ribbon(aes(ymin = Mean_Winter - SD_Winter, ymax = Mean_Winter + SD_Winter, fill = "Winter"), alpha = 0.2, na.rm = T)+
  geom_ribbon(aes(ymin = Mean_Summer*coeff - SD_Summer*coeff, ymax = Mean_Summer*coeff + SD_Summer*coeff, fill = "Summer"), alpha = 0.2, na.rm = T)+
  labs(fill = "Model", color = "Model")+
  theme_classic()+
  theme(legend.position = "none", 
        axis.text = element_text(size=15, color = 'black'),
        axis.title=element_text(size=18))


pdp.all.DG.df <- rbind(pdp.summer.DG$data, pdp.winter.DG$data)
pdp.all.DG.df$Model <- c(rep("Summer", nrow(pdp.summer.DG$data)), rep("Winter", nrow(pdp.winter.DG$data)))
pdp.all.DG <- ggplot(pdp.all.DG.df, aes(x = DG, y = Mean))+
  geom_line(aes(color = Model)) +
  labs(x = "DG Fractional Cover (%)", y = NULL)+ 
  geom_ribbon(aes(ymin = Mean - SD, ymax = Mean + SD, fill = Model), alpha = 0.2)+
  scale_y_break(breaks = c(0.2, 0.8), scale = "fixed")+
  theme_classic()+
  theme(legend.position = "none")

coeff <- 5.5
pdp.all.DG.w <- pivot_wider(data = pdp.all.DG.df, names_from = Model, values_from = c(Mean, SD), names_sep = "_")
pdp.all.DG.w <- splitsort(pdp.all.DG.w, "DG", "Mean_Summer", "Mean_Winter", "SD_Summer", "SD_Winter")
pdpw.DG <- ggplot(pdp.all.DG.w, aes(x = DG)) +
  geom_path(aes(y=Mean_Winter, color = "Winter"), na.rm = T) + 
  geom_path( aes(y=Mean_Summer*coeff, color = "Summer"), na.rm = T) + 
  scale_y_continuous(name = "Winter Albedo", sec.axis = sec_axis(~./coeff, name="Summer Albedo"))+
  scale_color_manual(values = c("red", "blue"))+
  scale_fill_manual(values = c("red", "blue"))+
  labs(x = "DG Fractional Cover (%)")+
  geom_ribbon(aes(ymin = Mean_Winter - SD_Winter, ymax = Mean_Winter + SD_Winter, fill = "Winter"), alpha = 0.2, na.rm = T)+
  geom_ribbon(aes(ymin = Mean_Summer*coeff - SD_Summer*coeff, ymax = Mean_Summer*coeff + SD_Summer*coeff, fill = "Summer"), alpha = 0.2, na.rm = T)+
  labs(fill = "Model", color = "Model")+
  theme_classic()+
  theme(legend.position = "inside", 
        legend.position.inside = c(.75, .25), 
        axis.text = element_text(size=15, color = 'black'),
        axis.title= element_text(size=18), 
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 15))


pdp.all.DLS.df <- rbind(pdp.summer.DLS$data, pdp.winter.DLS$data)
pdp.all.DLS.df$Model <- c(rep("Summer", nrow(pdp.summer.DLS$data)), rep("Winter", nrow(pdp.winter.DLS$data)))
pdp.all.DLS <- ggplot(pdp.all.DLS.df, aes(x = DLS, y = Mean, fill = Model))+
  geom_line(aes(color = Model)) +
  labs(x = "DLS Fractional Cover (%)", y = NULL)+ 
  geom_ribbon(aes(ymin = Mean - SD, ymax = Mean + SD, fill = Model), alpha = 0.2)+
  scale_y_break(breaks = c(0.2, 0.75), scale = "fixed")+
  theme_classic()+
  theme(legend.position = "none")

coeff <- 5.5
pdp.all.DLS.w <- pivot_wider(data = pdp.all.DLS.df, names_from = Model, values_from = c(Mean, SD), names_sep = "_")
pdp.all.DLS.w <- splitsort(pdp.all.DLS.w, "DLS", "Mean_Summer", "Mean_Winter", "SD_Summer", "SD_Winter")
pdpw.DLS <- ggplot(pdp.all.DLS.w, aes(x = DLS)) +
  geom_path(aes(y=Mean_Winter, color = "Winter"), na.rm = T) + 
  geom_path( aes(y=Mean_Summer*coeff, color = "Summer"), na.rm = T) + 
  scale_y_continuous(name = "Winter Albedo", sec.axis = sec_axis(~./coeff, name="Summer Albedo"))+
  scale_color_manual(values = c("red", "blue"))+
  scale_fill_manual(values = c("red", "blue"))+
  labs(x = "DLS Fractional Cover (%)")+
  geom_ribbon(aes(ymin = Mean_Winter - SD_Winter, ymax = Mean_Winter + SD_Winter, fill = "Winter"), alpha = 0.2, na.rm = T)+
  geom_ribbon(aes(ymin = Mean_Summer*coeff - SD_Summer*coeff, ymax = Mean_Summer*coeff + SD_Summer*coeff, fill = "Summer"), alpha = 0.2, na.rm = T)+
  labs(fill = "Model", color = "Model")+
  theme_classic()+
  theme(legend.position = "none", 
        axis.text = element_text(size=15, color = 'black'),
        axis.title=element_text(size=18))


pdp.all.DT.df <- rbind(pdp.summer.DT$data, pdp.winter.DT$data)
pdp.all.DT.df$Model <- c(rep("Summer", nrow(pdp.summer.DT$data)), rep("Winter", nrow(pdp.winter.DT$data)))
pdp.all.DT <- ggplot(pdp.all.DT.df, aes(x = DT, y = Mean, fill = Model))+
  geom_line(aes(color = Model)) +
  labs(x = "DT Fractional Cover (%)", y = NULL)+ 
  geom_ribbon(aes(ymin = Mean - SD, ymax = Mean + SD, fill = Model), alpha = 0.2)+
  scale_y_break(breaks = c(0.2, 0.75), scale = "fixed")+
  theme_classic()+
  theme(legend.position = "none")

coeff <- 5.5
pdp.all.DT.w <- pivot_wider(data = pdp.all.DT.df, names_from = Model, values_from = c(Mean, SD), names_sep = "_")
pdp.all.DT.w <- splitsort(pdp.all.DT.w, "DT", "Mean_Summer", "Mean_Winter", "SD_Summer", "SD_Winter")
pdpw.DT <- ggplot(pdp.all.DT.w, aes(x = DT)) +
  geom_path(aes(y=Mean_Winter, color = "Winter"), na.rm = T) + 
  geom_path( aes(y=Mean_Summer*coeff, color = "Summer"), na.rm = T) + 
  scale_y_continuous(name = "Winter Albedo", sec.axis = sec_axis(~./coeff, name="Summer Albedo"))+
  scale_color_manual(values = c("red", "blue"))+
  scale_fill_manual(values = c("red", "blue"))+
  labs(x = "DT Fractional Cover (%)")+
  geom_ribbon(aes(ymin = Mean_Winter - SD_Winter, ymax = Mean_Winter + SD_Winter, fill = "Winter"), alpha = 0.2, na.rm = T)+
  geom_ribbon(aes(ymin = Mean_Summer*coeff - SD_Summer*coeff, ymax = Mean_Summer*coeff + SD_Summer*coeff, fill = "Summer"), alpha = 0.2, na.rm = T)+
  labs(fill = "Model", color = "Model")+
  theme_classic()+
  theme(legend.position = "none", 
        axis.text = element_text(size=15, color = 'black'),
        axis.title=element_text(size=18))


pdp.all.DTSA.df <- rbind(pdp.summer.DTSA$data, pdp.winter.DTSA$data)
pdp.all.DTSA.df$Model <- c(rep("Summer", nrow(pdp.summer.DTSA$data)), rep("Winter", nrow(pdp.winter.DTSA$data)))
pdp.all.DTSA <- ggplot(pdp.all.DTSA.df, aes(x = DTSA, y = Mean, fill = Model))+
  geom_line(aes(color = Model)) +
  labs(x = "DTSA Fractional Cover (%)", y = NULL)+ 
  geom_ribbon(aes(ymin = Mean - SD, ymax = Mean + SD, fill = Model), alpha = 0.2)+
  scale_y_break(breaks = c(0.2, 0.75), scale = "fixed")+
  theme_classic()+
  theme(legend.position = "none")

coeff <- 5.5
pdp.all.DTSA.w <- pivot_wider(data = pdp.all.DTSA.df, names_from = Model, values_from = c(Mean, SD), names_sep = "_")
pdp.all.DTSA.w <- splitsort(pdp.all.DTSA.w, "DTSA", "Mean_Summer", "Mean_Winter", "SD_Summer", "SD_Winter")
pdpw.DTSA <- ggplot(pdp.all.DTSA.w, aes(x = DTSA)) +
  geom_path(aes(y=Mean_Winter, color = "Winter"), na.rm = T) + 
  geom_path( aes(y=Mean_Summer*coeff, color = "Summer"), na.rm = T) + 
  scale_y_continuous(name = "Winter Albedo", sec.axis = sec_axis(~./coeff, name="Summer Albedo"))+
  scale_color_manual(values = c("red", "blue"))+
  scale_fill_manual(values = c("red", "blue"))+
  labs(x = "DTSA Fractional Cover (%)")+
  geom_ribbon(aes(ymin = Mean_Winter - SD_Winter, ymax = Mean_Winter + SD_Winter, fill = "Winter"), alpha = 0.2, na.rm = T)+
  geom_ribbon(aes(ymin = Mean_Summer*coeff - SD_Summer*coeff, ymax = Mean_Summer*coeff + SD_Summer*coeff, fill = "Summer"), alpha = 0.2, na.rm = T)+
  labs(fill = "Model", color = "Model")+
  theme_classic()+
  theme(legend.position = "none", 
        axis.text = element_text(size=15, color = 'black'),
        axis.title=element_text(size=18))

pdp.all.DTSW.df <- rbind(pdp.summer.DTSW$data, pdp.winter.DTSW$data)
pdp.all.DTSW.df$Model <- c(rep("Summer", nrow(pdp.summer.DTSW$data)), rep("Winter", nrow(pdp.winter.DTSW$data)))
pdp.all.DTSW <- ggplot(pdp.all.DTSW.df, aes(x = DTSW, y = Mean, fill = Model))+
  geom_line(aes(color = Model)) +
  labs(x = "DTSW Fractional Cover (%)", y = NULL)+ 
  geom_ribbon(aes(ymin = Mean - SD, ymax = Mean + SD, fill = Model), alpha = 0.2)+
  scale_y_break(breaks = c(0.2, 0.75), scale = "fixed")+
  theme_classic()+
  theme(legend.position = "none")

coeff <- 5.5
pdp.all.DTSW.w <- pivot_wider(data = pdp.all.DTSW.df, names_from = Model, values_from = c(Mean, SD), names_sep = "_")
pdp.all.DTSW.w <- splitsort(pdp.all.DTSW.w, "DTSW", "Mean_Summer", "Mean_Winter", "SD_Summer", "SD_Winter")
pdpw.DTSW <- ggplot(pdp.all.DTSW.w, aes(x = DTSW)) +
  geom_path(aes(y=Mean_Winter, color = "Winter"), na.rm = T) + 
  geom_path( aes(y=Mean_Summer*coeff, color = "Summer"), na.rm = T) + 
  scale_y_continuous(name = "Winter Albedo", sec.axis = sec_axis(~./coeff, name="Summer Albedo"))+
  scale_color_manual(values = c("red", "blue"))+
  scale_fill_manual(values = c("red", "blue"))+
  labs(x = "DTSW Fractional Cover (%)")+
  geom_ribbon(aes(ymin = Mean_Winter - SD_Winter, ymax = Mean_Winter + SD_Winter, fill = "Winter"), alpha = 0.2, na.rm = T)+
  geom_ribbon(aes(ymin = Mean_Summer*coeff - SD_Summer*coeff, ymax = Mean_Summer*coeff + SD_Summer*coeff, fill = "Summer"), alpha = 0.2, na.rm = T)+
  labs(fill = "Model", color = "Model")+
  theme_classic()+
  theme(legend.position = "none", 
        axis.text = element_text(size=15, color = 'black'),
        axis.title=element_text(size=18))


pdp.all.ET.df <- rbind(pdp.summer.ET$data, pdp.winter.ET$data)
pdp.all.ET.df$Model <- c(rep("Summer", nrow(pdp.summer.ET$data)), rep("Winter", nrow(pdp.winter.ET$data)))
pdp.all.ET <- ggplot(pdp.all.ET.df, aes(x = ET, y = Mean, fill = Model))+
  geom_line(aes(color = Model)) +
  labs(x = "ET Fractional Cover (%)", y = NULL)+ 
  geom_ribbon(aes(ymin = Mean - SD, ymax = Mean + SD, fill = Model), alpha = 0.2)+
  scale_y_break(breaks = c(0.2, 0.75), scale = "fixed")+
  theme_classic()+
  theme(legend.position = "none")

coeff <- 5.5
pdp.all.ET.w <- pivot_wider(data = pdp.all.ET.df, names_from = Model, values_from = c(Mean, SD), names_sep = "_")
pdp.all.ET.w <- splitsort(pdp.all.ET.w, "ET", "Mean_Summer", "Mean_Winter", "SD_Summer", "SD_Winter")
pdpw.ET <- ggplot(pdp.all.ET.w, aes(x = ET)) +
  geom_path(aes(y=Mean_Winter, color = "Winter"), na.rm = T) + 
  geom_path( aes(y=Mean_Summer*coeff, color = "Summer"), na.rm = T) + 
  scale_y_continuous(name = "Winter Albedo", sec.axis = sec_axis(~./coeff, name="Summer Albedo"))+
  scale_color_manual(values = c("red", "blue"))+
  scale_fill_manual(values = c("red", "blue"))+
  labs(x = "ET Fractional Cover (%)")+
  geom_ribbon(aes(ymin = Mean_Winter - SD_Winter, ymax = Mean_Winter + SD_Winter, fill = "Winter"), alpha = 0.2, na.rm = T)+
  geom_ribbon(aes(ymin = Mean_Summer*coeff - SD_Summer*coeff, ymax = Mean_Summer*coeff + SD_Summer*coeff, fill = "Summer"), alpha = 0.2, na.rm = T)+
  labs(fill = "Model", color = "Model")+
  theme_classic()+
  theme(legend.position = "none", 
        axis.text = element_text(size=15, color = 'black'),
        axis.title=element_text(size=18))


pdp.all.LI.df <- rbind(pdp.summer.LI$data, pdp.winter.LI$data)
pdp.all.LI.df$Model <- c(rep("Summer", nrow(pdp.summer.LI$data)), rep("Winter", nrow(pdp.winter.LI$data)))
pdp.all.LI <- ggplot(pdp.all.LI.df, aes(x = LI, y = Mean, fill = Model))+
  geom_line(aes(color = Model)) +
  labs(x = "LI Fractional Cover (%)", y = NULL)+ 
  geom_ribbon(aes(ymin = Mean - SD, ymax = Mean + SD, fill = Model), alpha = 0.2)+
  scale_y_break(breaks = c(0.2, 0.75), scale = "fixed")+
  theme_classic()+
  theme(legend.position = "none")

coeff <- 5.5
pdp.all.LI.w <- pivot_wider(data = pdp.all.LI.df, names_from = Model, values_from = c(Mean, SD), names_sep = "_")
pdp.all.LI.w <- splitsort(pdp.all.LI.w, "LI", "Mean_Summer", "Mean_Winter", "SD_Summer", "SD_Winter")
pdpw.LI <- ggplot(pdp.all.LI.w, aes(x = LI)) +
  geom_path(aes(y=Mean_Winter, color = "Winter"), na.rm = T) + 
  geom_path( aes(y=Mean_Summer*coeff, color = "Summer"), na.rm = T) + 
  scale_y_continuous(name = "Winter Albedo", sec.axis = sec_axis(~./coeff, name="Summer Albedo"))+
  scale_color_manual(values = c("red", "blue"))+
  scale_fill_manual(values = c("red", "blue"))+
  labs(x = "LI Fractional Cover (%)")+
  geom_ribbon(aes(ymin = Mean_Winter - SD_Winter, ymax = Mean_Winter + SD_Winter, fill = "Winter"), alpha = 0.2, na.rm = T)+
  geom_ribbon(aes(ymin = Mean_Summer*coeff - SD_Summer*coeff, ymax = Mean_Summer*coeff + SD_Summer*coeff, fill = "Summer"), alpha = 0.2, na.rm = T)+
  labs(fill = "Model", color = "Model")+
  theme_classic()+
  theme(legend.position = "none", 
        axis.text = element_text(size=15, color = 'black'),
        axis.title=element_text(size=18))


pdp.all.MO.df <- rbind(pdp.summer.MO$data, pdp.winter.MO$data)
pdp.all.MO.df$Model <- c(rep("Summer", nrow(pdp.summer.MO$data)), rep("Winter", nrow(pdp.winter.MO$data)))
pdp.all.MO <- ggplot(pdp.all.MO.df, aes(x = MO, y = Mean, fill = Model))+
  geom_line(aes(color = Model)) +
  labs(x = "MO Fractional Cover (%)", y = NULL)+ 
  geom_ribbon(aes(ymin = Mean - SD, ymax = Mean + SD, fill = Model), alpha = 0.2)+
  scale_y_break(breaks = c(0.2, 0.75), scale = "fixed")+
  theme_classic()+
  theme(legend.position = "none")

coeff <- 5.5
pdp.all.MO.w <- pivot_wider(data = pdp.all.MO.df, names_from = Model, values_from = c(Mean, SD), names_sep = "_")
pdp.all.MO.w <- splitsort(pdp.all.MO.w, "MO", "Mean_Summer", "Mean_Winter", "SD_Summer", "SD_Winter")
pdpw.MO <- ggplot(pdp.all.MO.w, aes(x = MO)) +
  geom_path(aes(y=Mean_Winter, color = "Winter"), na.rm = T) + 
  geom_path( aes(y=Mean_Summer*coeff, color = "Summer"), na.rm = T) + 
  scale_y_continuous(name = "Winter Albedo", sec.axis = sec_axis(~./coeff, name="Summer Albedo"))+
  scale_color_manual(values = c("red", "blue"))+
  scale_fill_manual(values = c("red", "blue"))+
  labs(x = "MO Fractional Cover (%)")+
  geom_ribbon(aes(ymin = Mean_Winter - SD_Winter, ymax = Mean_Winter + SD_Winter, fill = "Winter"), alpha = 0.2, na.rm = T)+
  geom_ribbon(aes(ymin = Mean_Summer*coeff - SD_Summer*coeff, ymax = Mean_Summer*coeff + SD_Summer*coeff, fill = "Summer"), alpha = 0.2, na.rm = T)+
  labs(fill = "Model", color = "Model")+
  theme_classic()+
  theme(legend.position = "none", 
        axis.text = element_text(size=15, color = 'black'),
        axis.title=element_text(size=18))


pdp.all.NPV.df <- rbind(pdp.summer.NPV$data, pdp.winter.NPV$data)
pdp.all.NPV.df$Model <- c(rep("Summer", nrow(pdp.summer.NPV$data)), rep("Winter", nrow(pdp.winter.NPV$data)))
pdp.all.NPV <- ggplot(pdp.all.NPV.df, aes(x = NPV, y = Mean, fill = Model))+
  geom_line(aes(color = Model)) +
  labs(x = "NPV Fractional Cover (%)", y = NULL)+ 
  geom_ribbon(aes(ymin = Mean - SD, ymax = Mean + SD, fill = Model), alpha = 0.2)+
  scale_y_break(breaks = c(0.2, 0.75), scale = "fixed")+
  theme_classic()+
  theme(legend.position = "none")

coeff <- 5.5
pdp.all.NPV.w <- pivot_wider(data = pdp.all.NPV.df, names_from = Model, values_from = c(Mean, SD), names_sep = "_")
pdp.all.NPV.w <- splitsort(pdp.all.NPV.w, "NPV", "Mean_Summer", "Mean_Winter", "SD_Summer", "SD_Winter")
pdpw.NPV <- ggplot(pdp.all.NPV.w, aes(x = NPV)) +
  geom_path(aes(y=Mean_Winter, color = "Winter"), na.rm = T) + 
  geom_path( aes(y=Mean_Summer*coeff, color = "Summer"), na.rm = T) + 
  scale_y_continuous(name = "Winter Albedo", sec.axis = sec_axis(~./coeff, name="Summer Albedo"))+
  scale_color_manual(values = c("red", "blue"))+
  scale_fill_manual(values = c("red", "blue"))+
  labs(x = "NPV Fractional Cover (%)")+
  geom_ribbon(aes(ymin = Mean_Winter - SD_Winter, ymax = Mean_Winter + SD_Winter, fill = "Winter"), alpha = 0.2, na.rm = T)+
  geom_ribbon(aes(ymin = Mean_Summer*coeff - SD_Summer*coeff, ymax = Mean_Summer*coeff + SD_Summer*coeff, fill = "Summer"), alpha = 0.2, na.rm = T)+
  labs(fill = "Model", color = "Model")+
  theme_classic()+
  theme(legend.position = "none", 
        axis.text = element_text(size=15, color = 'black'),
        axis.title=element_text(size=18))


pdp.all.slope.df <- rbind(pdp.summer.slope$data, pdp.winter.slope$data)
pdp.all.slope.df$Model <- c(rep("Summer", nrow(pdp.summer.slope$data)), rep("Winter", nrow(pdp.winter.slope$data)))
pdp.all.slope <- ggplot(pdp.all.slope.df, aes(x = slope, y = Mean, fill = Model))+
  geom_line(aes(color = Model)) +
  labs(x = "Slope (°)", y = NULL)+ 
  geom_ribbon(aes(ymin = Mean - SD, ymax = Mean + SD, fill = Model), alpha = 0.2)+
  scale_y_break(breaks = c(0.2, 0.8), scale = "fixed")+
  theme_classic()+
  theme(legend.position = "none")

coeff <- 5.5
pdp.all.slope.w <- pivot_wider(data = pdp.all.slope.df, names_from = Model, values_from = c(Mean, SD), names_sep = "_")
pdp.all.slope.w <- splitsort(pdp.all.slope.w, "slope", "Mean_Summer", "Mean_Winter", "SD_Summer", "SD_Winter")
pdpw.slope <- ggplot(pdp.all.slope.w, aes(x = slope)) +
  geom_path(aes(y=Mean_Winter, color = "Winter"), na.rm = T) + 
  geom_path( aes(y=Mean_Summer*coeff, color = "Summer"), na.rm = T) + 
  scale_y_continuous(name = "Winter Albedo", sec.axis = sec_axis(~./coeff, name="Summer Albedo"))+
  scale_color_manual(values = c("red", "blue"))+
  scale_fill_manual(values = c("red", "blue"))+
  labs(x = "Slope (°)")+
  geom_ribbon(aes(ymin = Mean_Winter - SD_Winter, ymax = Mean_Winter + SD_Winter, fill = "Winter"), alpha = 0.2, na.rm = T)+
  geom_ribbon(aes(ymin = Mean_Summer*coeff - SD_Summer*coeff, ymax = Mean_Summer*coeff + SD_Summer*coeff, fill = "Summer"), alpha = 0.2, na.rm = T)+
  labs(fill = "Model", color = "Model")+
  theme_classic()+
  theme(legend.position = "none", 
        axis.text = element_text(size=15, color = 'black'),
        axis.title=element_text(size=18))


pdp.all.TPI.df <- rbind(pdp.summer.TPI$data, pdp.winter.TPI$data)
pdp.all.TPI.df$Model <- c(rep("Summer", nrow(pdp.summer.TPI$data)), rep("Winter", nrow(pdp.winter.TPI$data)))
pdp.all.TPI <- ggplot(pdp.all.TPI.df, aes(x = TPI, y = Mean, fill = Model))+
  geom_line(aes(color = Model)) +
  labs(x = "TPI", y = NULL)+ 
  geom_ribbon(aes(ymin = Mean - SD, ymax = Mean + SD, fill = Model), alpha = 0.2)+
  scale_y_break(breaks = c(0.2, 0.8), scale = "fixed")+
  theme_classic()+
  theme(legend.position = "none")

coeff <- 5.5
pdp.all.TPI.w <- pivot_wider(data = pdp.all.TPI.df, names_from = Model, values_from = c(Mean, SD), names_sep = "_")
pdp.all.TPI.w <- splitsort(pdp.all.TPI.w, "TPI", "Mean_Summer", "Mean_Winter", "SD_Summer", "SD_Winter")
pdpw.TPI <- ggplot(pdp.all.TPI.w, aes(x = TPI)) +
  geom_path(aes(y=Mean_Winter, color = "Winter"), na.rm = T) + 
  geom_path( aes(y=Mean_Summer*coeff, color = "Summer"), na.rm = T) + 
  scale_y_continuous(name = "Winter Albedo", sec.axis = sec_axis(~./coeff, name="Summer Albedo"))+
  scale_color_manual(values = c("red", "blue"))+
  scale_fill_manual(values = c("red", "blue"))+
  labs(x = "TPI")+
  geom_ribbon(aes(ymin = Mean_Winter - SD_Winter, ymax = Mean_Winter + SD_Winter, fill = "Winter"), alpha = 0.2, na.rm = T)+
  geom_ribbon(aes(ymin = Mean_Summer*coeff - SD_Summer*coeff, ymax = Mean_Summer*coeff + SD_Summer*coeff, fill = "Summer"), alpha = 0.2, na.rm = T)+
  labs(fill = "Model", color = "Model")+
  theme_classic()+
  theme(legend.position = "none", 
        axis.text = element_text(size=15, color = 'black'),
        axis.title=element_text(size=18))

pdp.all.TRI.df <- rbind(pdp.summer.TRI$data, pdp.winter.TRI$data)
pdp.all.TRI.df$Model <- c(rep("Summer", nrow(pdp.summer.TRI$data)), rep("Winter", nrow(pdp.winter.TRI$data)))
pdp.all.TRI <- ggplot(pdp.all.TRI.df, aes(x = TRI, y = Mean, fill = Model))+
  geom_line(aes(color = Model)) +
  labs(x = "TRI", y = NULL)+ 
  geom_ribbon(aes(ymin = Mean - SD, ymax = Mean + SD, fill = Model), alpha = 0.2)+
  scale_y_break(breaks = c(0.2, 0.8), scale = "fixed")+
  theme_classic()+
  theme(legend.position = "none")

coeff <- 5.5
pdp.all.TRI.w <- pivot_wider(data = pdp.all.TRI.df, names_from = Model, values_from = c(Mean, SD), names_sep = "_")
pdp.all.TRI.w <- splitsort(pdp.all.TRI.w, "TRI", "Mean_Summer", "Mean_Winter", "SD_Summer", "SD_Winter")
pdpw.TRI <- ggplot(pdp.all.TRI.w, aes(x = TRI)) +
  geom_path(aes(y=Mean_Winter, color = "Winter"), na.rm = T) + 
  geom_path(aes(y=Mean_Summer*coeff, color = "Summer"), na.rm = T) + 
  scale_y_continuous(name = "Winter Albedo", sec.axis = sec_axis(~./coeff, name="Summer Albedo"))+
  scale_color_manual(values = c("red", "blue"))+
  scale_fill_manual(values = c("red", "blue"))+
  labs(x = "TRI")+
  geom_ribbon(aes(ymin = Mean_Winter - SD_Winter, ymax = Mean_Winter + SD_Winter, fill = "Winter"), alpha = 0.2, na.rm = T)+
  geom_ribbon(aes(ymin = Mean_Summer*coeff - SD_Summer*coeff, ymax = Mean_Summer*coeff + SD_Summer*coeff, fill = "Summer"), alpha = 0.2, na.rm = T)+
  labs(fill = "Model", color = "Model")+
  theme_classic()+
  theme(legend.position = "none", 
        axis.text = element_text(size=15, color = 'black'),
        axis.title=element_text(size=18))

pdp.all.TWI.df <- rbind(pdp.summer.TWI$data, pdp.winter.TWI$data)
pdp.all.TWI.df$Model <- c(rep("Summer", nrow(pdp.summer.TWI$data)), rep("Winter", nrow(pdp.winter.TWI$data)))
pdp.all.TWI <- ggplot(pdp.all.TWI.df, aes(x = TWI, y = Mean, fill = Model))+
  geom_line(aes(color = Model)) +
  labs(x = "TWI", y = NULL)+ 
  geom_ribbon(aes(ymin = Mean - SD, ymax = Mean + SD, fill = Model), alpha = 0.2)+
  scale_y_break(breaks = c(0.2, 0.8), scale = "fixed")+
  theme_classic()+
  theme(legend.position = "none")

coeff <- 5.5
pdp.all.TWI.w <- pivot_wider(data = pdp.all.TWI.df, names_from = Model, values_from = c(Mean, SD), names_sep = "_")
pdp.all.TWI.w <- splitsort(pdp.all.TWI.w, "TWI", "Mean_Summer", "Mean_Winter", "SD_Summer", "SD_Winter")
pdpw.TWI <- ggplot(pdp.all.TWI.split, aes(x = TWI)) +
  geom_path(aes(y=Mean_Winter, color = "Winter"), na.rm = T) + 
  geom_path(aes(y=Mean_Summer*coeff, color = "Summer"), na.rm = T) + 
  scale_y_continuous(name = "Winter Albedo", sec.axis = sec_axis(~./coeff, name="Summer Albedo"))+
  scale_color_manual(values = c("red", "blue"))+
  scale_fill_manual(values = c("red", "blue"))+
  labs(x = "TWI")+
  geom_ribbon(aes(ymin = Mean_Winter - SD_Winter, ymax = Mean_Winter + SD_Winter, fill = "Winter"), alpha = 0.2, na.rm = T)+
  geom_ribbon(aes(ymin = Mean_Summer*coeff - SD_Summer*coeff, ymax = Mean_Summer*coeff + SD_Summer*coeff, fill = "Summer"), alpha = 0.2, na.rm = T)+
  labs(fill = "Model", color = "Model")+
  theme_classic()+
  theme(legend.position = "none", 
        axis.text = element_text(size=15, color = 'black'),
        axis.title=element_text(size=18))


pdp.all.WG.df <- rbind(pdp.summer.WG$data, pdp.winter.WG$data)
pdp.all.WG.df$Model <- c(rep("Summer", nrow(pdp.summer.WG$data)), rep("Winter", nrow(pdp.winter.WG$data)))
pdp.all.WG <- ggplot(pdp.all.WG.df, aes(x = WG, y = Mean, fill = Model))+
  geom_line(aes(color = Model)) +
  labs(x = "WG Fractional Cover (%)", y = NULL)+ 
  geom_ribbon(aes(ymin = Mean - SD, ymax = Mean + SD, fill = Model), alpha = 0.2)+
  scale_y_break(breaks = c(0.2, 0.75), scale = "fixed")+
  theme_classic()+
  theme(legend.position = "none")

coeff <- 5.5
pdp.all.WG.w <- pivot_wider(data = pdp.all.WG.df, names_from = Model, values_from = c(Mean, SD), names_sep = "_")
pdp.all.WG.w <- splitsort(pdp.all.WG.w, "WG", "Mean_Summer", "Mean_Winter", "SD_Summer", "SD_Winter")
pdpw.WG <- ggplot(pdp.all.WG.w, aes(x = WG)) +
  geom_path(aes(y=Mean_Winter, color = "Winter"), na.rm = T) + 
  geom_path( aes(y=Mean_Summer*coeff, color = "Summer"), na.rm = T) + 
  scale_y_continuous(name = "Winter Albedo", sec.axis = sec_axis(~./coeff, name="Summer Albedo"))+
  scale_color_manual(values = c("red", "blue"))+
  scale_fill_manual(values = c("red", "blue"))+
  labs(x = "WG Fractional Cover (%)")+
  geom_ribbon(aes(ymin = Mean_Winter - SD_Winter, ymax = Mean_Winter + SD_Winter, fill = "Winter"), alpha = 0.2, na.rm = T)+
  geom_ribbon(aes(ymin = Mean_Summer*coeff - SD_Summer*coeff, ymax = Mean_Summer*coeff + SD_Summer*coeff, fill = "Summer"), alpha = 0.2, na.rm = T)+
  labs(fill = "Model", color = "Model")+
  theme_classic()+
  theme(legend.position = "none", 
        axis.text = element_text(size=15, color = 'black'),
        axis.title=element_text(size=18))


#Plot the partial dependence plots

#Plot top 6 VIPs (breaks and overlaid) for Figure 6
#plot_grid(pdp.all.CHM, pdp.all.ET, pdp.all.DG, pdp.all.aspect, pdp.all.TRI, pdp.all.TWI)
#aplot::plot_list(pdp.all.CHM, pdp.all.ET, pdp.all.DG, pdp.all.aspect, pdp.all.TRI, pdp.all.TWI, labels = c("A", "B", "C", "D", "E", "F"))
#This next one uses the old top 6 from a previous run
#aplot::plot_list(pdpw.CHM, pdpw.ET, pdpw.DG, pdpw.aspect, pdpw.TRI, pdpw.TWI, labels = c("A", "B", "C", "D", "E", "F"))
#This uses the new top 6 from the current run, DEM instead of TRI
aplot::plot_list(pdpw.CHM, pdpw.ET, pdpw.DG, pdpw.DEM, pdpw.aspect, pdpw.TWI, labels = c("A", "B", "C", "D", "E", "F"))
outDIR <- "/Users/anncrumlish/Downloads/albedo analysis/map_stats"
pdfNAME = paste0(outDIR, "/", "Figure6.pdf")
ggsave(pdfNAME, plot = last_plot(), width = 35, height = 13, units = 'cm')

#Plot all variables for Supplemental Figure
pdpw.DG2 <- pdpw.DG + theme(legend.position = "none") #need to remove the legend here, as it was there for inclusion in Figure 6
pdpw.TWI2  <- pdpw.TWI + theme(legend.position = "inside", legend.position.inside = c(0.25, 0.25))
aplot::plot_list(pdpw.ET, pdpw.DTSA, pdpw.DTSW, pdpw.DT, pdpw.DLS, pdpw.DG2, pdpw.WG, pdpw.MO, pdpw.LI, pdpw.NPV, pdpw.CHM, pdpw.DEM, pdpw.slope, pdpw.aspect, pdpw.TPI, pdpw.TRI, pdpw.TWI2, ncol = 3, labels = c("A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N", "O", "P", "Q"))
pdfNAME = paste0(outDIR, "/", "FigureS7.pdf")
ggsave(pdfNAME, plot = last_plot(), width = 40, height = 50, units = 'cm')


load("/Users/anncrumlish/Downloads/albedo analysis/map_stats/FinalRFrun.RData")
###########
#DO NOT TOUCH
#NAME = paste0(outDIR, "/", "FinalRFrun.RData")
#save.image(NAME)





#plot_grid(pdp.all.CHM, pdp.all.ET, pdp.all.DG, pdp.all.aspect, pdp.all.TRI, pdp.all.TWI)
aplot::plot_list(pdp.all.CHM, pdp.all.ET, pdp.all.DG, pdp.all.aspect, pdp.all.TRI, pdp.all.TWI, labels = c("A", "B", "C", "D", "E", "F"))



pdp.summer.DG
plot_grid(pdp.summer.DG, pdp.summer.ET, pdp.summer.chm, pdp.summer.DEM)

#Plot both VIP plots together

plot_grid(bar.summer, bar.winter, nrow = 2, labels = c("A", "B"))
plot_grid(
  bar.summer,
  ggdraw() + theme_void(),  # Empty plot for the line
  bar.winter,
  ncol = 1,
  rel_heights = c(1, 0.05, 1)  # Adjust this to change the thickness of the line
)

# Add the line
grid::grid.lines(x = c(0, 1), y = c(0.5, 0.5), gp = grid::gpar(col = "black"))

plot_grid(bar.summer, bar.winter, labels = c("A", "B"), nrow = 2)
outDIR <- "/Users/anncrumlish/Downloads/albedo analysis/map_stats"
pdfNAME = paste0(outDIR, "/", "combined_RF_prediction_2SD_ms.pdf")
ggsave(pdfNAME, plot = last_plot(), width = 12, height = 12, units = 'cm')
plot_grid(pred.summer, pred.winter, labels = c("A", "B"), nrow = 2)
plot_grid(pdp.summer.ET, pdp.summer.DG, pdp.summer.TRI, pdp.winter.CHM, pdp.winter.aspect, pdp.winter.ET, nrow = 2, labels = c("A", "B", "C", "D", "E", "F"))
