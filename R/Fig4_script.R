###########################################################################################
#
#        This script recreates figure 4 in Shuman et al.
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
list.of.packages <- c('raster', 'caTools', 'haven', 'terra', 'ggplot2', 'dplyr', 'cowplot', 
                      'grid', 'gtable', 'ggfun', 'scatterpie', 'viridis', 'multcompView', 'tidyverse', 'ggpattern')
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
#set directories
data.dir <- "~/Downloads/albedo/Data/"
breakpoints.dir <- "~/Downloads/albedo analysis/step4_smoothed_ts/"

#Load in data
summer_albedo <- rast(paste0(data.dir, 'council_watershed_albedo_summer.tif')) 
winter_albedo <- rast(paste0(data.dir, 'council_watershed_albedo_winter.tif')) 
canopy_height_model <- rast(paste0(data.dir, 'council_watershed_chm_30m.tif')) 
canopy_height_model[canopy_height_model < 0] <- NA #there are about 6000 negative values which I think should be removed
aspect <- rast(paste0(data.dir, 'council_watershed_aspect_30m.tif')) 
slope <- rast(paste0(data.dir, 'council_watershed_slope_30m.tif')) 
topo <- rast(paste0(data.dir, 'council_watershed_topo_30m.tif')) 
tpi <- rast(paste0(data.dir, 'council_watershed_tpi_30m.tif'))  #topographic position index, positive values mean the cell is higher than its surrounding cells, negative mean it is lower
tri <- rast(paste0(data.dir, 'council_watershed_tri_30m.tif')) #terrain ruggedness index, higher values mean more heterogeneity in the surrounding cells 
fcover <- rast(paste0(data.dir, 'council_watershed_fcover_30m.dat')) 
pft <- rast(paste0(data.dir, 'council_watershed_pft_30m.dat')) 
twi <- rast(paste0(data.dir, 'twi.tif'))
breakpoints1 <- rast(paste0(breakpoints.dir, "break_point1_30m.tif"))
breakpoints <- rast(paste0(breakpoints.dir, "break_point_30m.tif"))
breakpoints[breakpoints > 213] <- NA
breakpoints1[breakpoints1 < 92] <- NA

#resample data to match fcover
chm.resamp <- resample(canopy_height_model, fcover)
topo.resamp <- resample(topo, fcover)
tpi.resamp <- resample(tpi, fcover)
tri.resamp <- resample(tri, fcover)
slope.resamp <- resample(slope, fcover)
aspect.resamp <- resample(aspect, fcover)
twi.resamp <- resample(twi, fcover)
summer_albedo.resamp <- resample(summer_albedo, fcover)
winter_albedo.resamp <- resample(winter_albedo, fcover)
breakpoints1.resamp <- resample(breakpoints1, fcover)
breakpoints.resamp <- resample(breakpoints, fcover)
stack <- c(chm.resamp, topo.resamp, tpi.resamp, tri.resamp, slope.resamp, aspect.resamp, twi.resamp, summer_albedo.resamp, winter_albedo.resamp, fcover, pft, breakpoints1.resamp, breakpoints.resamp)
stack_df <- as.data.frame(stack, xy = TRUE)
colnames(stack_df) <- c("x", "y" , "council_watershed_chm_30m", "council_watershed_topo_30m", "council_watershed_tpi_30m" , "council_watershed_tri_30m", "council_watershed_slope_30m", "council_watershed_aspect_30m", "twi", "council_watershed_albedo_summer", "council_watershed_albedo_winter", "Spruce_fcover", "Alder_fcover" , "Willow_fcover", "OtherTall_fcover" , "LowShrubs_fcover", "DwarfShrubs_fcover", "EvergreenShrubs_fcover", "Forb_fcover" , "DryGrass_fcover", "WetGrass_fcover", "Moss_fcover", "Lichen_fcover", "NPV_fcover", "category", "break_point1_30m", "break_point_30m")
#***********************************************************************************************************************************************#

#********************************** filter stack_df for PFT fCover >75% for figure 4 ***********************************************************#
stack_df2 <- stack_df[,12:24]
stack_df2[is.na(stack_df2)] <- 0
# Create a new column "HighCoverSpecies" which is a list of column names where the value is greater than 0.75
stack_df2$category2 <- apply(stack_df2, 1, function(row) {
  # Identify columns with values greater than 0.75
  category2 <- names(stack_df2)[row > 0.75]
  
  # Return the column names as a list
  if (length(category2) == 0) {
    return(NA)  # Return NA if no species have a cover greater than 0.75
  } else {
    return(list(category2))
  }
})
#View(stack_df2)
stack_df2 <- cbind(stack_df[,1:11], stack_df2[,1:14], stack_df[,25:27])
stack_df2$category3 <- sapply(stack_df2$category2, function(x) {
  if (length(x) == 1) {
    return(unlist(x))
  } else {
    return(NA) # Placeholder for lists with more than one element
  }
})
stack_df2 <- stack_df2[!is.na(stack_df2$category3), ]
stack_df2 <- stack_df2[, !names(stack_df2) %in% c("category2", "category3")]
#***********************************************************************************************************************************************#





#********************************** create figure 4 *********************************************************************************************#
##### Figure 4a- Completion of albedo transition (breakpoint) v.s. PFT with significance ######
#Use stack_df2 so that categories represent only pixels >75% fCover of a certain PFT

pfts <- c("ET","DT","DTSA","DTSW","DLS","DDS","ES","FO","DG","WG","MO","LI","NPV" )
pfts <- pfts[-c(6:8)]
raw_cols <- c('#000000', '#FF0000', '#008000', '#00CD00', '#669999', '#0066FF', '#FFFF66', '#66FFFF', '#00B0F0', '#00FF67', '#FF67FF', '#CC6600', '#FFFFFF', '#5A5A5A')
cols <- raw_cols[-c(1, 6:8)]
cols <- c(cols[1], cols[4], cols[2:3], cols[-c(1:4)])
stack_df_4b <- as.data.frame(filter(stack_df2, !category %in% c("EvergS", "Forb", "nonClass")))
stack_df_4b$category <- as.character(stack_df_4b$category)
stack_df_4b$category <- replace(stack_df_4b$category, stack_df_4b$category == "Spruce", "ET")
stack_df_4b$category <- replace(stack_df_4b$category, stack_df_4b$category == "OtherTall", "DT")
stack_df_4b$category <- replace(stack_df_4b$category, stack_df_4b$category == "Alder", "DTSA")
stack_df_4b$category <- replace(stack_df_4b$category, stack_df_4b$category == "Willow", "DTSW")
stack_df_4b$category <- replace(stack_df_4b$category, stack_df_4b$category == "LowS", "DLS")
stack_df_4b$category <- replace(stack_df_4b$category, stack_df_4b$category == "DryG", "DG")
stack_df_4b$category <- replace(stack_df_4b$category, stack_df_4b$category == "WetG", "WG")
stack_df_4b$category <- replace(stack_df_4b$category, stack_df_4b$category == "Moss", "MO")
stack_df_4b$category <- replace(stack_df_4b$category, stack_df_4b$category == "Lichen", "LI")
stack_df_4b$category <- factor(stack_df_4b$category,levels = c("ET", "DTSA", "DTSW", "DT", "DLS", "DG", "WG", "MO", "LI", "NPV"))
stack_df_4b$pattern_group <- ifelse(stack_df_4b$category %in% c("ET", "DT", "DTSA", "DTSW"), "Woody", "Non-woody")

#First just make the boxplots
Fig4a = ggplot(data = stack_df_4b, aes(pattern = pattern_group))+
  geom_boxplot(aes(x = category, y = break_point_30m, fill = category), outlier.shape = NA)+
  geom_boxplot_pattern(
    aes(x = category, y = break_point_30m),
    pattern_density = 0.1, 
    pattern_fill = "black", 
    pattern_angle = 45, 
    fill = NA, 
    color = "black", 
    outlier.shape = NA
  ) +
  scale_pattern_manual(values = c("Woody" = "stripe", "Non-woody" = NA), name = "Pattern") +
  xlab("PFT") + ylab("Completion of Albedo Transition (DOY)")+
  scale_fill_manual(values = cols)+
  theme_classic() +
  theme(axis.text = element_text(size=21, color = 'black'),
        axis.title=element_text(size=25),
        axis.text.x = element_text(angle = 60, hjust = 1), 
        legend.position = "none", 
        aspect.ratio = 1) 

#Calculate the mean break point by PFT category for reporting in the text
mean(stack_df$break_point_30m, na.rm = T)
2*sd(stack_df$break_point_30m, na.rm = T)
for (i in unique(stack_df2$category)){
  print(i)
  a <- filter(stack_df2, stack_df2$category == i)
  print(paste0("Mean breakpoint: ", mean(a$break_point_30m, na.rm = T)))
  print(paste0("2 SD of breakpoint: ", 2*sd(a$break_point_30m, na.rm = T)))
  print(paste0("Difference between ", i, " breakpoint and all pixels: ", mean(stack_df$break_point_30m, na.rm = T) - mean(a$break_point_30m, na.rm = T))) #subtract dominant PFT mean from the mean of all data/PFTs (including non-dominant pixels)
  print("_______________________")
}

##Perform ANOVA and Tukey's HSD post-hoc to test for significant differences in albedo transition between PFTs #####
pfts <- c("ET","DT","DTSA","DTSW","DLS","DDS","ES","FO","DG","WG","MO","LI","NVS" )
pfts <- pfts[-c(6:8)]
cat <- as.vector(stack_df_4b$category)
anova.df <- as.data.frame(na.omit(cbind(stack_df_4b$break_point_30m, cat)))
colnames(anova.df) <- c("break_point_30m", "category")
anova.df$break_point_30m <- as.numeric(anova.df$break_point_30m)
anova.df$category <- as.factor(anova.df$category)

model=lm(anova.df$break_point_30m ~ anova.df$category)
ANOVA=aov(model)

# Tukey test to study each pair of treatment :
TUKEY <- TukeyHSD(x=ANOVA, 'anova.df$category', conf.level=0.95)

# Tukey test representation :
plot(TUKEY , las=1 , col="brown")

#Generate labels to depict **Non-significant** pairs on the boxplot
generate_label_df <- function(TUKEY, variable){
  # Extract labels and factor levels from Tukey post-hoc 
  Tukey.levels <- TUKEY[[variable]][,4]
  Tukey.labels <- data.frame(multcompLetters(Tukey.levels)['Letters'])
  #Order them
  Tukey.labels$treatment=rownames(Tukey.labels)
  Tukey.labels=Tukey.labels[order(Tukey.labels$treatment) , ]
  return(Tukey.labels)
}

LABELS <- generate_label_df(TUKEY , "anova.df$category")
LABELS$Letters <- sapply(LABELS$Letters, function(str) {
  paste(strsplit(str, NULL)[[1]], collapse = " ")
})
#Plot Figure 4a with labels from TUKEY post hoc test
stack_df_4b$category <- factor(stack_df_4b$category,levels = c("ET", "DTSA", "DTSW", "DT", "DLS", "DG", "WG", "MO", "LI", "NPV"))

over <- 0.1 * max(anova.df$break_point_30m)
boxplot_stats <- boxplot(anova.df$break_point_30m ~ anova.df$category, plot = FALSE)
label_positions <- c(160, 170, 170, 160, 160, 150, 165, 160, 165, 160) #boxplot_stats$stats[nrow(boxplot_stats$stats), ] + over
label_df <- data.frame(
  category = levels(anova.df$category),
  y = label_positions,
  label = LABELS$Letters, 
  pattern_group = c("Woody", "Woody", "Woody", "Woody", "Non-woody", "Non-woody", "Non-woody", "Non-woody", "Non-woody", "Non-woody")
)

anova.df$category <- factor(anova.df$category,levels = c("ET", "DTSA", "DTSW", "DT", "DLS", "DG", "WG", "MO", "LI", "NPV"))
label_df$category <- factor(label_df$category,levels = c("ET", "DTSA", "DTSW", "DT", "DLS", "DG", "WG", "MO", "LI", "NPV"))
Fig4a_Tukey = Fig4a + geom_text(data = label_df, aes(x = category, y = y, label = label), vjust = -0.5,  size = 7)
Fig4a_Tukey

##### Figure 4b- Completion of albedo transition (breakpoint) v.s. CHM with PFT fractional cover ######
#Create tall and short dfs for regression and testing the significance of the visually observed 5.5m threshold
tall_df = filter(stack_df, council_watershed_chm_30m > 5.5 & council_watershed_chm_30m < 12.5)
tall_df = as.data.frame(tall_df)
tall_df = filter(tall_df, is.na(break_point_30m) == F)
mean(tall_df$break_point_30m)
short_df = filter(stack_df, council_watershed_chm_30m < 5.5)
short_df = as.data.frame(short_df)
short_df = filter(short_df, is.na(break_point_30m) == F)
mean(short_df$break_point_30m)
mean(short_df$break_point_30m) - mean(tall_df$break_point_30m)
t.test(short_df$break_point_30m, tall_df$break_point_30m) #Student's t-test to show that <5.5m and >5.5m breakpoint DOY are significantly different

#Calculate mean breakpoint DOY for the 1m CHM bins
bin_means = data.frame(matrix(vector("logical"), ncol = 14))
colnames(bin_means) = colnames (stack_df[,c(3, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24)])
for (i in 0:11){
  start = i
  finish = i + 1
  bin_init = filter(stack_df, council_watershed_chm_30m > start & council_watershed_chm_30m < finish)
  bin = bin_init[,c(3, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24)]
  colnames(bin) = colnames(stack_df[,c(3, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24)])
  bin[] <- lapply(bin, function(x) replace(x, is.na(x), 0))
  bin_means = rbind(bin_means, as.vector(colMeans(bin)))
}
colnames(bin_means) = colnames(stack_df[,c(3, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24)])
#bin_means$council_watershed_chm_30m <- c("0-1", "1-2", "2-3", "3-4", "4-5", "5-6", "6-7 m", "7-8 m", "8-9 m", "9-10 m", "10-11 m", "11-12 m", "12-13 m", "13-14 m", "14-15 m", "15-16 m", "16+ m")
colnames(bin_means) = c("CHM (m)", "ET", "DTSA", "DTSW", "DT", "LS", "DS", "ES", "FO", "DG", "WG", "MO", "LI", "NPV")

bin_means$Y <- rep(as.numeric(175), nrow(bin_means))
bin_means$group <- c(0:(nrow(bin_means)-1))
#bin_means <- bin_means[1:16, ]


#Create df to use for pie chart inset
long <- seq(1, 12, 1)
lat <- rep(130, 12)
#seq((130 + 15*2), 130, -2)
d <- data.frame(long=long, lat=lat)
n <- nrow(d)
d$region <- factor(1:n)
d$ET <- bin_means$ET * 100
d$DTSA <- bin_means$DTSA * 100
d$DTSW <- bin_means$DTSW * 100
d$DT <- bin_means$DT * 100
d$DLS <- bin_means$LS * 100
d$DG <- bin_means$DG * 100
d$WG <- bin_means$WG* 100
d$MO <- bin_means$MO * 100
d$LI <- bin_means$LI * 100
d$NPV <- bin_means$NPV * 100

inset <-  ggplot() + #make the pft pie charts to show fCover for each CHM bin
  coord_equal() + 
  geom_scatterpie(aes(x=long, y=lat, group=region, r = .45), data=d,cols=c("ET", "DTSA", "DTSW", "DT", "DLS", "DG", "WG", "MO", "LI", "NPV"))+
  scale_fill_manual(values = c('#FF0000', '#008000', '#00CD00', '#669999', '#0066FF', '#00FF67', '#FF67FF', '#CC6600', '#FFFFFF', '#5A5A5A'))+
  theme_classic()+
  theme(legend.position = "none",
        axis.text=element_blank(), 
        axis.title=element_blank(), 
        axis.ticks=element_blank(), 
        axis.line=element_blank())
ldf <- data.frame(matrix(vector("logical"), ncol = 3, nrow = 10))
ldf$PFT <- c("ET", "DTSA", "DTSW", "DT", "DLS", "DG", "WG", "MO", "LI", "NPV")
ldf$X <- seq(1, nrow(ldf), 1)
ldf$Y <- seq(1, nrow(ldf),1)
ldf$PFT <- factor(ldf$PFT,levels = c("ET", "DTSA", "DTSW", "DT", "DLS", "DG", "WG", "MO", "LI", "NPV"))

#Create separate legend for the pie chart insert
legend.plot <-  ggplot(data = ldf, aes(x=X, y=Y, fill = PFT)) + #legend for PFT pie charts
  geom_point(pch = 21, color = "black", size = 8)+
  scale_fill_manual(values = c('#FF0000', '#008000', '#00CD00', '#669999', '#0066FF', '#00FF67', '#FF67FF', '#CC6600', '#FFFFFF', '#5A5A5A'), 
                    guide = guide_legend(nrow = 5))+
  theme_classic()+
  theme(axis.text=element_blank(), 
        axis.title=element_blank(), 
        axis.ticks=element_blank(), 
        axis.line=element_blank(), 
        legend.title = element_text(size = 18), 
        legend.text = element_text(size = 14))
legend <- get_legend(legend.plot)

#Make the actual DOY v.s CHM boxplot with optional regression lines (not included in final manuscript due to non-linear relationship)
cols = viridis(12)
All_box <- stack_df %>%
  mutate(levels = cut(council_watershed_chm_30m, seq(0, 12, 1))) %>%
  drop_na(levels, break_point_30m) %>%
  ggplot(aes(x = levels, y = break_point_30m)) +
  geom_boxplot(aes(fill = levels), outlier.shape = NA, width = .9) +
  #geom_smooth(data = short_df, method = "lm", aes(x = council_watershed_chm_30m, y = break_point_30m), color = "blue")+
  #geom_smooth(data = tall_df, method = "lm", aes(x = council_watershed_chm_30m, y = break_point_30m), color = "blue")+
  #geom_smooth(data = filter(stack_df, council_watershed_chm_30m < 12.5), method = "gam", aes(x = council_watershed_chm_30m, y = break_point_30m), color = "red")+ #GAM does a good job of showing the overall trend within the boxplots
  scale_y_continuous(limits = c(125, 170), breaks = seq(120, 170, by = 10))+
  scale_fill_manual(values = cols)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  labs(x = "Binned Canopy Height (m)", y = "Completion of Albedo Transition (DOY)")+
  theme_classic() +
  theme(axis.text = element_text(size=21, color = 'black'),
        axis.title=element_text(size=25),
        axis.text.x = element_text(angle = 60, hjust = 1), 
        legend.position = "none", 
        aspect.ratio = 1)

summary(lm(data = short_df, break_point_30m ~ council_watershed_chm_30m))
summary(lm(data = tall_df, break_point_30m ~ council_watershed_chm_30m))

Fig4b<- All_box + 
  annotation_custom(
    grob = ggplotGrob(inset), 
    xmin = -0.25, xmax = 13.5, ymin = 120, ymax = 130)+
  annotation_custom(
    grob = legend, 
    xmin = 10, xmax = 12, ymin = 150, ymax = 170)
Fig4b


##### Combine and save Figures 4a and 4b ######

plot_grid(Fig4a_Tukey, Fig4b)

outDIR <- "/Users/anncrumlish/Downloads/albedo analysis/map_stats"
pdfNAME = paste0(outDIR, "/", "Figure4.pdf")
ggsave(pdfNAME, plot = last_plot(), width = 40, height = 20, units = 'cm')
#***********************************************************************************************************************************************#

























