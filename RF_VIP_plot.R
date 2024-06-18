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
  xlab("PFT") + ylab("VIP") + ggtitle("Winter")+
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
  xlab("PFT") + ylab("VIP") + ggtitle("Summer")+
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
xfactor.w <- new_order
xfactor.w <- factor(xfactor.w, levels = unique(vipPLOT.winter$pft))
fill_scale_values <- setNames(colors, levels(xfactor.w))
bar.winter <- ggplot(data = vipPLOT.winter, aes(x=xfactor.w, y=mean, fill=pft)) +
  geom_bar(stat = 'identity', color = "black") +
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
variables_of_interest <- c("CHM", "ET", "DG", "aspect", "TRI", "TWI")
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
xfactor.s <- new_order
xfactor.s <- factor(xfactor.s, levels = unique(vipPLOT.summer$pft))
fill_scale_values <- setNames(colors, levels(xfactor.s))
bar.summer <- ggplot(data = vipPLOT.summer, aes(x=xfactor.s, y=mean, fill=pft)) +
  geom_bar(stat = 'identity', color = "black") +
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
variables_of_interest <- c("CHM", "ET", "DG", "aspect", "TRI", "TWI")
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
pdp.all.ET.df <- rbind(pdp.summer.ET$data, pdp.winter.ET$data)
pdp.all.ET.df$Model <- c(rep("Summer", nrow(pdp.summer.ET$data)), rep("Winter", nrow(pdp.winter.ET$data)))
pdp.all.ET <- ggplot(pdp.all.ET.df, aes(x = ET, y = Mean, fill = Model))+
  geom_line(aes(color = Model)) +
  labs(x = "ET Fractional Cover (%)", y = NULL)+ 
  geom_ribbon(aes(ymin = Mean - SD, ymax = Mean + SD, fill = Model), alpha = 0.2)+
  scale_y_break(breaks = c(0.2, 0.75), scale = "fixed")+
  theme_classic()+
  theme(legend.position = "none")

pdp.all.DG.df <- rbind(pdp.summer.DG$data, pdp.winter.DG$data)
pdp.all.DG.df$Model <- c(rep("Summer", nrow(pdp.summer.DG$data)), rep("Winter", nrow(pdp.winter.DG$data)))
pdp.all.DG <- ggplot(pdp.all.DG.df, aes(x = DG, y = Mean))+
  geom_line(aes(color = Model)) +
  labs(x = "DG Fractional Cover (%)", y = NULL)+ 
  geom_ribbon(aes(ymin = Mean - SD, ymax = Mean + SD, fill = Model), alpha = 0.2)+
  scale_y_break(breaks = c(0.18, 0.825), scale = "fixed")+
  theme_classic()+
  theme(legend.position = "none")
  

pdp.all.CHM.df <- rbind(pdp.summer.CHM$data, pdp.winter.CHM$data)
pdp.all.CHM.df$Model <- c(rep("Summer", nrow(pdp.summer.CHM$data)), rep("Winter", nrow(pdp.winter.CHM$data)))
pdp.all.CHM <- ggplot(pdp.all.CHM.df, aes(x = CHM, y = Mean, fill = Model))+
  geom_line(aes(color = Model)) +
  labs(x = "CHM (m)", y = "Mean Seasonal Albedo")+ 
  geom_ribbon(aes(ymin = Mean - SD, ymax = Mean + SD, fill = Model), alpha = 0.2)+
  scale_y_break(breaks = c(0.2, 0.75), scale = "fixed")+
  theme_classic()+
  theme(legend.position = "none")

pdp.all.TWI.df <- rbind(pdp.summer.TWI$data, pdp.winter.TWI$data)
pdp.all.TWI.df$Model <- c(rep("Summer", nrow(pdp.summer.TWI$data)), rep("Winter", nrow(pdp.winter.TWI$data)))
pdp.all.TWI <- ggplot(pdp.all.TWI.df, aes(x = TWI, y = Mean, fill = Model))+
  geom_line(aes(color = Model)) +
  labs(x = "TWI", y = NULL)+ 
  geom_ribbon(aes(ymin = Mean - SD, ymax = Mean + SD, fill = Model), alpha = 0.2)+
  scale_y_break(breaks = c(0.2, 0.84), scale = "fixed")+
  theme_classic()+
  theme(legend.position = "none")

pdp.all.aspect.df <- rbind(pdp.summer.aspect$data, pdp.winter.aspect$data)
pdp.all.aspect.df$Model <- c(rep("Summer", nrow(pdp.summer.aspect$data)), rep("Winter", nrow(pdp.winter.aspect$data)))
pdp.all.aspect <- ggplot(pdp.all.aspect.df, aes(x = aspect, y = Mean, fill = Model))+
  geom_line(aes(color = Model)) +
  labs(x = "Aspect (°)", y = "Mean Seasonal Albedo")+ 
  geom_ribbon(aes(ymin = Mean - SD, ymax = Mean + SD, fill = Model), alpha = 0.2)+
  scale_y_break(breaks = c(0.175, 0.84), scale = "fixed")+
  theme_classic()+
  theme(legend.position = "none")

pdp.all.TRI.df <- rbind(pdp.summer.TRI$data, pdp.winter.TRI$data)
pdp.all.TRI.df$Model <- c(rep("Summer", nrow(pdp.summer.TRI$data)), rep("Winter", nrow(pdp.winter.TRI$data)))
pdp.all.TRI <- ggplot(pdp.all.TRI.df, aes(x = TRI, y = Mean, fill = Model))+
  geom_line(aes(color = Model)) +
  labs(x = "TRI", y = NULL)+ 
  geom_ribbon(aes(ymin = Mean - SD, ymax = Mean + SD, fill = Model), alpha = 0.2)+
  scale_y_break(breaks = c(0.2, 0.84), scale = "fixed")+
  theme_classic()+
  theme(legend.position = "none")


pdp.all.ET.df <- rbind(pdp.summer.ET$data, pdp.winter.ET$data)
pdp.all.ET.df$Model <- c(rep("Summer", nrow(pdp.summer.ET$data)), rep("Winter", nrow(pdp.winter.ET$data)))
pdp.all.ET <- ggplot(pdp.all.ET.df, aes(x = ET, y = Mean, fill = Model))+
  geom_line(aes(color = Model)) +
  labs(x = "ET Fractional Cover (%)", y = NULL)+ 
  geom_ribbon(aes(ymin = Mean - SD, ymax = Mean + SD, fill = Model), alpha = 0.2)+
  scale_y_break(breaks = c(0.2, 0.75), scale = "fixed")+
  theme_classic()+
  theme(legend.position = "none")

pdp.all.DG.df <- rbind(pdp.summer.DG$data, pdp.winter.DG$data)
pdp.all.DG.df$Model <- c(rep("Summer", nrow(pdp.summer.DG$data)), rep("Winter", nrow(pdp.winter.DG$data)))
pdp.all.DG <- ggplot(pdp.all.DG.df, aes(x = DG, y = Mean))+
  geom_line(aes(color = Model)) +
  labs(x = "DG Fractional Cover (%)", y = NULL)+ 
  geom_ribbon(aes(ymin = Mean - SD, ymax = Mean + SD, fill = Model), alpha = 0.2)+
  scale_y_break(breaks = c(0.2, 0.8), scale = "fixed")+
  theme_classic()+
  theme(legend.position = "none")


pdp.all.CHM.df <- rbind(pdp.summer.CHM$data, pdp.winter.CHM$data)
pdp.all.CHM.df$Model <- c(rep("Summer", nrow(pdp.summer.CHM$data)), rep("Winter", nrow(pdp.winter.CHM$data)))
pdp.all.CHM <- ggplot(pdp.all.CHM.df, aes(x = CHM, y = Mean, fill = Model))+
  geom_line(aes(color = Model)) +
  labs(x = "CHM (m)", y = "Mean Seasonal Albedo")+ 
  geom_ribbon(aes(ymin = Mean - SD, ymax = Mean + SD, fill = Model), alpha = 0.2)+
  scale_y_break(breaks = c(0.2, 0.75), scale = "fixed")+
  theme_classic()+
  theme(legend.position = "none")

pdp.all.TWI.df <- rbind(pdp.summer.TWI$data, pdp.winter.TWI$data)
pdp.all.TWI.df$Model <- c(rep("Summer", nrow(pdp.summer.TWI$data)), rep("Winter", nrow(pdp.winter.TWI$data)))
pdp.all.TWI <- ggplot(pdp.all.TWI.df, aes(x = TWI, y = Mean, fill = Model))+
  geom_line(aes(color = Model)) +
  labs(x = "TWI", y = NULL)+ 
  geom_ribbon(aes(ymin = Mean - SD, ymax = Mean + SD, fill = Model), alpha = 0.2)+
  scale_y_break(breaks = c(0.2, 0.8), scale = "fixed")+
  theme_classic()+
  theme(legend.position = "none")

pdp.all.aspect.df <- rbind(pdp.summer.aspect$data, pdp.winter.aspect$data)
pdp.all.aspect.df$Model <- c(rep("Summer", nrow(pdp.summer.aspect$data)), rep("Winter", nrow(pdp.winter.aspect$data)))
pdp.all.aspect <- ggplot(pdp.all.aspect.df, aes(x = aspect, y = Mean, fill = Model))+
  geom_line(aes(color = Model)) +
  labs(x = "Aspect (°)", y = "Mean Seasonal Albedo")+ 
  geom_ribbon(aes(ymin = Mean - SD, ymax = Mean + SD, fill = Model), alpha = 0.2)+
  scale_y_break(breaks = c(0.175, 0.84), scale = "fixed")+
  theme_classic()+
  theme(legend.position = "none")

pdp.all.TRI.df <- rbind(pdp.summer.TRI$data, pdp.winter.TRI$data)
pdp.all.TRI.df$Model <- c(rep("Summer", nrow(pdp.summer.TRI$data)), rep("Winter", nrow(pdp.winter.TRI$data)))
pdp.all.TRI <- ggplot(pdp.all.TRI.df, aes(x = TRI, y = Mean, fill = Model))+
  geom_line(aes(color = Model)) +
  labs(x = "TRI", y = NULL)+ 
  geom_ribbon(aes(ymin = Mean - SD, ymax = Mean + SD, fill = Model), alpha = 0.2)+
  scale_y_break(breaks = c(0.2, 0.8), scale = "fixed")+
  theme_classic()+
  theme(legend.position = "none")
#plot_grid(pdp.all.CHM, pdp.all.ET, pdp.all.DG, pdp.all.aspect, pdp.all.TRI, pdp.all.TWI)
aplot::plot_list(pdp.all.CHM, pdp.all.ET, pdp.all.DG, pdp.all.aspect, pdp.all.TRI, pdp.all.TWI, labels = c("A", "B", "C", "D", "E", "F"))



pdp.summer.DG
plot_grid(pdp.summer.DG, pdp.summer.ET, pdp.summer.chm, pdp.summer.DEM)

#Plot both VIP plots together

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

plot_grid(bar.summer, bar.winter, labels = c("A", "B"), nrow = 2)
plot_grid(pred.summer, pred.winter, labels = c("A", "B"), nrow = 2)
plot_grid(pdp.summer.ET, pdp.summer.DG, pdp.summer.TRI, pdp.winter.CHM, pdp.winter.aspect, pdp.winter.ET, nrow = 2, labels = c("A", "B", "C", "D", "E", "F"))
