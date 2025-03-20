###########################################################################################
#
#        This script runs the random forest model to create the VIP and PDP figures (Figures 5 and 6 in Shuman et al.)
#
#    --- Last updated:  2025.03.20 By Ian Shuman <ins2109@columbia.edu>
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
                      "randomForest", "caTools", "pls", "spectratrait", "terra", "pdp", "ggbreak", "tidyr", "cowplot", "Metrics")
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

#Set output directory to save graphs to
outDIR <- "/Users/anncrumlish/Downloads/albedo analysis/map_stats"
#*****************************************************************************************#
#*
load("stack_df.RData") #Load from the Figure 1 script


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

#Plot Figure 5
bar.summer.plot <- bar.summer + xlab("Predictor Variable") 
bar.winter.plot <- bar.winter + xlab("Predictor Variable")
plot_grid(bar.summer.plot, bar.winter.plot, nrow = 2)
pdfNAME = paste0(outDIR, "/", "Figure5.pdf")
ggsave(pdfNAME, plot = last_plot(), width = 40, height = 20, units = 'cm')

#Plot Figure S6
pred.summer.plot <- pred.summer + 
  ggtitle("Summer RF Model")+ 
  theme(title = element_text(face = "bold", size = 22), 
        axis.text = element_text(size = 18), 
        axis.title = element_text(size = 20)) +
  stat_poly_eq(formula = my.formula, 
               aes(label = paste(..eq.label.., ..rr.label..,..p.value.label..,
                                 sep = "*`,`~")), 
               label.x = 1, label.y = 0.1,
               parse = T, size = 6) 
pred.winter.plot <- pred.winter +
  ggtitle("Winter RF Model")+ 
  theme(title = element_text(face = "bold", size = 22), 
        axis.text = element_text(size = 18), 
        axis.title = element_text(size = 20))+
  stat_poly_eq(formula = my.formula, 
               aes(label = paste(..eq.label.., ..rr.label..,..p.value.label..,
                                 sep = "*`,`~")), 
               label.x = 1, label.y = 0.1,
               parse = T, size = 6) 

plot_grid(pred.summer.plot, pred.winter.plot, nrow = 1)
pdfNAME = paste0(outDIR, "/", "FigureS6.pdf")
ggsave(pdfNAME, plot = last_plot(), width = 40, height = 20, units = 'cm')

#Calculate RMSE
Metrics::rmse(predPLOT.s$ref, predPLOT.s$pred) #summer RF
Metrics::rmse(predPLOT.w$ref, predPLOT.w$pred) #winter RF


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

coeff <- 5.5
pdp.all.ET.w <- pivot_wider(data = pdp.all.ET.df, names_from = Model, values_from = c(Mean, SD), names_sep = "_")
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


pdp.all.CHM.df <- rbind(pdp.summer.CHM$data, pdp.winter.CHM$data)
pdp.all.CHM.df$Model <- c(rep("Summer", nrow(pdp.summer.CHM$data)), rep("Winter", nrow(pdp.winter.CHM$data)))
pdp.all.CHM <- ggplot(pdp.all.CHM.df, aes(x = CHM, y = Mean, fill = Model))+
  geom_line(aes(color = Model)) +
  labs(x = "CHM (m)", y = "Mean Seasonal Albedo")+ 
  geom_ribbon(aes(ymin = Mean - SD, ymax = Mean + SD, fill = Model), alpha = 0.2)+
  scale_y_break(breaks = c(0.2, 0.75), scale = "fixed")+
  theme_classic()+
  theme(legend.position = "none")

coeff <- 5.5
pdp.all.CHM.w <- pivot_wider(data = pdp.all.CHM.df, names_from = Model, values_from = c(Mean, SD), names_sep = "_")
pdpw.CHM <- ggplot(pdp.all.CHM.w, aes(x = CHM)) +
  geom_path(aes(y=Mean_Winter, color = "Winter"), na.rm = T) + 
  geom_path( aes(y=Mean_Summer*coeff, color = "Summer"), na.rm = T) + 
  scale_y_continuous(name = "Winter Albedo", sec.axis = sec_axis(~./coeff, name="Summer Albedo"))+
  scale_color_manual(values = c("red", "blue"))+
  scale_fill_manual(values = c("red", "blue"))+
  labs(x = "CHM (m)")+
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
pdpw.TWI <- ggplot(pdp.all.TWI.w, aes(x = TWI)) +
  geom_path(aes(y=Mean_Winter, color = "Winter"), na.rm = T) + 
  geom_path( aes(y=Mean_Summer*coeff, color = "Summer"), na.rm = T) + 
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
pdpw.TRI <- ggplot(pdp.all.TRI.w, aes(x = TRI)) +
  geom_path(aes(y=Mean_Winter, color = "Winter"), na.rm = T) + 
  geom_path( aes(y=Mean_Summer*coeff, color = "Summer"), na.rm = T) + 
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

#plot_grid(pdp.all.CHM, pdp.all.ET, pdp.all.DG, pdp.all.aspect, pdp.all.TRI, pdp.all.TWI)

#Plot Figure 6
aplot::plot_list(pdpw.CHM, pdpw.ET, pdpw.DG, pdpw.DEM , pdpw.aspect, pdpw.TWI)
pdfNAME = paste0(outDIR, "/", "Figure6.pdf")
ggsave(pdfNAME, plot = last_plot(), width = 40, height = 20, units = 'cm')

#Plot Figure S7
aplot::plot_list(pdpw.ET, pdpw.DTSA, pdpw.DTSW, pdpw.DT, pdpw.DLS, pdpw.DG, pdpw.WG, pdpw.MO, pdpw.LI, pdpw.NPV, pdpw.CHM, pdpw.DEM, pdpw.slope, pdpw.aspect, pdpw.TPI, pdpw.TRI, pdpw.TWI, ncol = 4)
pdfNAME = paste0(outDIR, "/", "FigureS7.png")
ggsave(pdfNAME, plot = last_plot(), width = 60, height = 60, units = 'cm')
