###########################################################################################
#
#        This script recreates figures 2 and S1-S5 in Shuman et al.
#
#    --- Last updated:  2025.02.06 By Ian Shuman <ins2109@columbia.edu>
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
                      "randomForest", "caTools", "pls", "spectratrait", "terra", "cowplot", "viridis", "multcompView", "ggpattern", "tidyr", "quantreg")
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

# define outlier removal function and Tukey labeling function
VIP <- function(object) {
  SS <- c(object$Yloadings)^2 * colSums(object$scores^2)
  Wnorm2 <- colSums(object$loading.weights^2)
  SSW <- sweep(object$loading.weights^2, 2, SS / Wnorm2, "*") # Replace with matrix mult.
  sqrt(nrow(SSW) * apply(SSW, 1, cumsum) / cumsum(SS))
}

generate_label_df <- function(TUKEY, variable){
  # Extract labels and factor levels from Tukey post-hoc 
  Tukey.levels <- TUKEY[[variable]][,4]
  Tukey.labels <- data.frame(multcompLetters(Tukey.levels)['Letters'])
  #Order them
  Tukey.labels$treatment=rownames(Tukey.labels)
  Tukey.labels=Tukey.labels[order(Tukey.labels$treatment) , ]
  return(Tukey.labels)
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

#make breakpoints dataset
breakpoints <- raster("/Users/anncrumlish/Downloads/albedo analysis/step4_smoothed_ts/break_point_30m.tif")
test <- raster(nrow = 154, ncol = 242)
breakpoints.df <- as.data.frame(breakpoints, xy = T)
breakpoints2 <- resample(breakpoints, mapORIG)

#********************************** filter dataORIG for PFT fCover >75% for figures 2 + 3 **************************************#

dataORIG2 <- dataORIG[,5:17]
dataORIG2[is.na(dataORIG2)] <- 0
# Create a new column "HighCoverSpecies" which is a list of column names where the value is greater than 0.75
dataORIG2$PFT2 <- apply(dataORIG2, 1, function(row) {
  # Identify columns with values greater than 0.75
  PFT2 <- names(dataORIG2)[row > 0.75]
  
  # Return the column names as a list
  if (length(PFT2) == 0) {
    return(NA)  # Return NA if no species have a cover greater than 0.75
  } else {
    return(list(PFT2))
  }
})
dataORIG2 <- cbind(dataORIG[,1:4], dataORIG2[,1:14], dataORIG[,18:ncol(dataORIG)])


dataORIG2$PFT3 <- sapply(dataORIG2$PFT2, function(x) {
  if (length(x) == 1) {
    return(unlist(x))
  } else {
    return(NA) # Placeholder for lists with more than one element
  }
})
dataORIG2 <- dataORIG2[!is.na(dataORIG2$PFT3), ]
dataORIG2 <- dataORIG2[, !names(dataORIG2) %in% c("PFT2", "PFT3")]



#********************************** plot pft albedo **************************************#
# extract pft winter albedo
albedoWT <- cbind(dataORIG2[4:17], dataORIG2[c(22)])
albedoWT[albedoWT > 10000] <- 10000
albedoWT[albedoWT < 0] <- 0
albedoWT[which(albedoWT$PFT == 0),] <- NA
albedoWT[which(albedoWT$PFT == 6),] <- NA
albedoWT[which(albedoWT$PFT == 7),] <- NA
albedoWT[which(albedoWT$PFT == 8),] <- NA
albedoWT <- na.omit(albedoWT)
albedoWT_plot <- albedoWT
albedoWT_plot$pattern_group <- ifelse(albedoWT$PFT %in% c(1, 2, 3, 4), "Woody", "Non-woody")
pft_mapping <- c("1" = "ET","2" = "DT", "3" = "DTSA", "4" = "DTSW", "5" = "DLS", "9" = "DG", "10" = "WG", "11" = "MO", "12" = "LI", "13" = "NVS")
pft_levels <- c("ET", "DT", "DTSA", "DTSW", "DLS", "DG", "WG", "MO", "LI", "NVS")
albedoWT_plot <- albedoWT_plot %>%
  mutate(pft = case_when(
    PFT == 1 ~ "ET",
    PFT == 2 ~ "DT",
    PFT == 3 ~ "DTSA",
    PFT == 4 ~ "DTSW",
    PFT == 5 ~ "DLS",
    PFT == 9 ~ "DG",
    PFT == 10 ~ "WG",
    PFT == 11 ~ "MO",
    PFT == 12 ~ "LI",
    PFT == 13 ~ "NVS",
    TRUE ~ NA_character_))
albedoWT_plot$pft <- factor(albedoWT_plot$pft, levels = pft_levels)


# calculate mean of non tall vegetation
nonTSK <- mean(albedoWT[which(albedoWT$PFT>5),15])
# make a plot for winter albedo
winter.pft = ggplot(data = albedoWT_plot, aes(x=as.factor(PFT), y=winter_albedo, fill = as.factor(PFT), pattern = as.factor(pft))) +
  geom_boxplot(outlier.shape = NA) +
  geom_hline(yintercept = nonTSK, linetype="dashed") +
  #geom_segment(aes(x = 7.9, xend = 8.9, y = 0.03, yend = 0.03),  size = 1, 
  #linetype = "dashed") +
  #annotate("text", x = 9.1,  y = 0.03, 
  #label = stringr::str_wrap("Non-Woody Albedo", width = 10), hjust=0) +
  xlab("PFT") + ylab("Winter Albedo") + ylim(0, 1) + 
  scale_fill_manual(values = cols, name = "", labels = pfts[-c(6:8)]) +
  scale_x_discrete(labels = pfts[-c(6:8)]) +
  geom_boxplot_pattern(
    pattern_density = 0.1, 
    pattern_fill = "black", 
    pattern_angle = 45, 
    fill = NA, 
    color = "black", 
    outliers = F
  ) +
  scale_pattern_manual(values = c("stripe", "stripe", "stripe", "stripe", NA, NA, NA, NA, NA, NA), name = "Pattern") +
  guides(fill = "none", pattern = guide_legend(override.aes = list(fill = cols))) +
  theme(legend.position = "none",
        #c(0.87, 0.37), 
        legend.title = element_blank(), 
        legend.key.size = unit(0.75, 'cm'),
        legend.text = element_text(size = 12)) +
  theme(axis.text = element_text(size=15, color = 'black'),
        axis.title = element_text(size=18)) +
  theme(axis.line = element_line(colour = "black"),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
winter.pft
pdfNAME = paste0(outDIR, "/", "pft_albedo_winter_v1.pdf")
ggsave(pdfNAME, plot = last_plot(), width = 16, height = 13, units = 'cm')


#ANOVA and Tukey's test on winter albedo between PFTs
albedoWT2 <- albedoWT
albedoWT2$PFT <- as.character(albedoWT2$PFT)
albedoWT2$PFT <- replace(albedoWT2$PFT, albedoWT2$PFT == "1", "ET")
albedoWT2$PFT <- replace(albedoWT2$PFT, albedoWT2$PFT == "2", "DT")
albedoWT2$PFT <- replace(albedoWT2$PFT, albedoWT2$PFT == "3", "DTSA")
albedoWT2$PFT <- replace(albedoWT2$PFT, albedoWT2$PFT == "4", "DTSW")
albedoWT2$PFT <- replace(albedoWT2$PFT, albedoWT2$PFT == "5", "DLS")
albedoWT2$PFT <- replace(albedoWT2$PFT, albedoWT2$PFT == "9", "DG")
albedoWT2$PFT <- replace(albedoWT2$PFT, albedoWT2$PFT == "10", "WG")
albedoWT2$PFT <- replace(albedoWT2$PFT, albedoWT2$PFT == "11", "MO")
albedoWT2$PFT <- replace(albedoWT2$PFT, albedoWT2$PFT == "12", "LI")
albedoWT2$PFT <- replace(albedoWT2$PFT, albedoWT2$PFT == "13", "NVS")
cat <- albedoWT2$PFT
anova.df <- as.data.frame(na.omit(cbind(albedoWT$winter_albedo, cat)))
colnames(anova.df) <- c("albedo", "category")
anova.df$albedo <- as.numeric(anova.df$albedo)
anova.df$category <- factor(anova.df$category, levels = c("ET", "DT", "DTSA", "DTSW", "DLS", "DG", "WG", "MO", "LI", "NVS"))
model=lm(anova.df$albedo ~ anova.df$category)
ANOVA=aov(model)
TUKEY <- TukeyHSD(x=ANOVA, 'anova.df$category', conf.level=0.95)
plot(TUKEY , las=1 , col="brown")

#Add results of Tukey test to the boxplot
LABELS <- generate_label_df(TUKEY , "anova.df$category")
order <- c("ET", "DT", "DTSA", "DTSW", "DLS", "DG", "WG", "MO", "LI", "NVS")
LABELS <- LABELS %>%
  slice(match(order, treatment))
over <- 0.01 * max(anova.df$albedo)
boxplot_stats <- boxplot(anova.df$albedo ~ anova.df$category, plot = FALSE)
label_positions <- boxplot_stats$stats[nrow(boxplot_stats$stats), ] + over
label_df <- data.frame(
  category = levels(anova.df$category),
  y = label_positions,
  label = LABELS$Letters,
  pattern_group = c("stripe", "stripe", "stripe", "stripe", NA, NA, NA, NA, NA, NA),
  pft = order
)
label_df$category <- as.character(label_df$category)
label_df <- label_df %>%
  slice(match(order, category))
label_df$category <- factor(label_df$category,levels = c("ET", "DT", "DTSA", "DTSW", "DLS", "DG", "WG", "MO", "LI", "NVS"))
label_df$PFT <- c("1", "2", "3", "4", "5", "9", "10", "11", "12", "13")
winter.pft_Tukey = winter.pft + geom_text(data = label_df, aes(x = PFT, y = y, label = label), vjust = -0.5) + ylim(0, 1.099)  
winter.pft_Tukey

# Perform pairwise t-tests comparing each PFT to nonTSK
pft_levels <- unique(albedoWT2$PFT)
p_values <- sapply(pft_levels, function(pft) {
  t.test(albedoWT2$winter_albedo[albedoWT2$PFT == pft], mu = nonTSK)$p.value
})

# Adjust p-values for multiple comparisons
p_adjusted <- p.adjust(p_values, method = "bonferroni")

# Generate labels based on significance
significance_labels <- ifelse(p_adjusted < 0.05, "*", "ns")

# Create a data frame for labels
label_df <- data.frame(
  PFT = pft_levels,
  y = rep(1.05, length(pft_levels)),  # Position labels above the plot
  label = significance_labels
)

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
albedoSM <- cbind(dataORIG2[4:17], dataORIG2[c(23)])
albedoSM[albedoSM > 10000] <- 10000
albedoSM[albedoSM < 0] <- 0
albedoSM[which(albedoSM$PFT == 0),] <- NA
albedoSM[which(albedoSM$PFT == 6),] <- NA
albedoSM[which(albedoSM$PFT == 7),] <- NA
albedoSM[which(albedoSM$PFT == 8),] <- NA
albedoSM <- na.omit(albedoSM)
albedoSM_plot <- albedoSM
albedoSM_plot$pattern_group <- ifelse(albedoSM_plot$PFT %in% c(1, 2, 3, 4), "Woody", "Non-woody")
pft_mapping <- c("1" = "ET","2" = "DT", "3" = "DTSA", "4" = "DTSW", "5" = "DLS", "9" = "DG", "10" = "WG", "11" = "MO", "12" = "LI", "13" = "NVS")
pft_levels <- c("ET", "DT", "DTSA", "DTSW", "DLS", "DG", "WG", "MO", "LI", "NVS")
albedoSM_plot <- albedoSM_plot %>%
  mutate(pft = case_when(
    PFT == 1 ~ "ET",
    PFT == 2 ~ "DT",
    PFT == 3 ~ "DTSA",
    PFT == 4 ~ "DTSW",
    PFT == 5 ~ "DLS",
    PFT == 9 ~ "DG",
    PFT == 10 ~ "WG",
    PFT == 11 ~ "MO",
    PFT == 12 ~ "LI",
    PFT == 13 ~ "NVS",
    TRUE ~ NA_character_))
albedoSM_plot$pft <- factor(albedoSM_plot$pft, levels = pft_levels)


# calculate mean of non tall vegetation
nonTSK <- mean(albedoSM[which(albedoSM$PFT>5),15])
# make a plot for winter albedo
summer.pft = ggplot(data = albedoSM_plot, aes(x=as.factor(PFT), y=summer_albedo, fill = as.factor(PFT), pattern = as.factor(pft))) +
  geom_boxplot(outlier.shape = NA) +
  geom_hline(yintercept = nonTSK, linetype="dashed") +
  geom_boxplot_pattern(
    pattern_density = 0.1, 
    pattern_fill = "black", 
    pattern_angle = 45, 
    fill = NA, 
    color = "black", 
    outliers = F
  ) +
  xlab("PFT") + ylab("Summer Albedo") + ylim(0, 0.22) +
  scale_fill_manual(values = cols, name = "", labels = pfts[-c(6:8)]) +
  scale_pattern_manual(values = c("stripe", "stripe", "stripe", "stripe", NA, NA, NA, NA, NA, NA), name = "Pattern") +
  guides(fill = "none", pattern = guide_legend(ncol = 5, override.aes = list(fill = cols)))+
  scale_x_discrete(labels = pfts[-c(6:8)]) +
  geom_segment(aes(x = 1.3, xend = 2.3, y = 0.03, yend = 0.03),  size = 1, 
               linetype = "dashed") +
  annotate("text", x = 2.5,  y = 0.03, 
           label = stringr::str_wrap("Non-Woody Albedo", width = 10), hjust=0) +
  theme(legend.position = c(0.65, 0.15), 
        legend.title = element_blank(), 
        legend.key.size = unit(0.75, 'cm'),
        legend.text = element_text(size = 12)) +
  theme(axis.text = element_text(size=15, color = 'black'),
        axis.title=element_text(size=18)) +
  theme(axis.line = element_line(colour = "black"),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
summer.pft
pdfNAME = paste0(outDIR, "/", "pft_albedo_summer.pdf")
ggsave(pdfNAME, plot = last_plot(), width = 16, height = 13, units = 'cm')

#ANOVA and Tukey's test on summer albedo between PFTs
albedoSM2 <- albedoSM_plot
albedoSM2$PFT <- as.character(albedoSM2$PFT)
albedoSM2$PFT <- replace(albedoSM2$PFT, albedoSM2$PFT == "1", "ET")
albedoSM2$PFT <- replace(albedoSM2$PFT, albedoSM2$PFT == "2", "DT")
albedoSM2$PFT <- replace(albedoSM2$PFT, albedoSM2$PFT == "3", "DTSA")
albedoSM2$PFT <- replace(albedoSM2$PFT, albedoSM2$PFT == "4", "DTSW")
albedoSM2$PFT <- replace(albedoSM2$PFT, albedoSM2$PFT == "5", "DLS")
albedoSM2$PFT <- replace(albedoSM2$PFT, albedoSM2$PFT == "9", "DG")
albedoSM2$PFT <- replace(albedoSM2$PFT, albedoSM2$PFT == "10", "WG")
albedoSM2$PFT <- replace(albedoSM2$PFT, albedoSM2$PFT == "11", "MO")
albedoSM2$PFT <- replace(albedoSM2$PFT, albedoSM2$PFT == "12", "LI")
albedoSM2$PFT <- replace(albedoSM2$PFT, albedoSM2$PFT == "13", "NVS")
cat <- albedoSM2$PFT
anova.df <- as.data.frame(na.omit(cbind(albedoSM$summer_albedo, cat)))
colnames(anova.df) <- c("albedo", "category")
anova.df$albedo <- as.numeric(anova.df$albedo)
anova.df$category <- anova.df$category <- factor(anova.df$category, levels = c("ET", "DT", "DTSA", "DTSW", "DLS", "DG", "WG", "MO", "LI", "NVS"))
model=lm(anova.df$albedo ~ anova.df$category)
ANOVA=aov(model)
TUKEY <- TukeyHSD(x=ANOVA, 'anova.df$category', conf.level=0.95)
plot(TUKEY , las=1 , col="brown")

#Add results of Tukey test to the boxplot
LABELS <- generate_label_df(TUKEY , "anova.df$category")
order <- c("ET", "DT", "DTSA", "DTSW", "DLS", "DG", "WG", "MO", "LI", "NVS")
LABELS <- LABELS %>%
  slice(match(order, treatment))
over <- 0.01 * max(anova.df$albedo)
boxplot_stats <- boxplot(anova.df$albedo ~ anova.df$category, plot = FALSE)
label_positions <- boxplot_stats$stats[nrow(boxplot_stats$stats), ] + over
label_df <- data.frame(
  category = levels(anova.df$category),
  y = label_positions,
  label = LABELS$Letters, 
  pattern_group = c("stripe", "stripe", "stripe", "stripe", NA, NA, NA, NA, NA, NA),
  pft = order
)
label_df$category <- as.character(label_df$category)
label_df$category <- factor(label_df$category,levels = c("ET", "DT", "DTSA", "DTSW", "DLS", "DG", "WG", "MO", "LI", "NVS"))
label_df$PFT <- c("1", "2", "3", "4", "5", "9", "10", "11", "12", "13")
summer.pft_Tukey = summer.pft + geom_text(data = label_df, aes(x = PFT, y = y, label = label), vjust = -0.5)
summer.pft_Tukey

# Perform pairwise t-tests comparing each PFT to nonTSK
pft_levels <- unique(albedoSM2$PFT)
p_values <- sapply(pft_levels, function(pft) {
  t.test(albedoSM2$winter_albedo[albedoSM2$PFT == pft], mu = nonTSK)$p.value
})

# Adjust p-values for multiple comparisons
p_adjusted <- p.adjust(p_values, method = "bonferroni")

# Generate labels based on significance
significance_labels <- ifelse(p_adjusted < 0.05, "*", "ns")

# Create a data frame for labels
label_df <- data.frame(
  PFT = pft_levels,
  y = rep(1.05, length(pft_levels)),  # Position labels above the plot
  label = significance_labels
)

# extract chm summer albedo
albedoCHM <- cbind(dataORIG[18], dataORIG[c(22, 23)])
albedoCHM[albedoCHM > 10000] <- 10000
albedoCHM[albedoCHM < 0] <- 0
albedoCHM <- na.omit(albedoCHM)

for (i in 1:nrow(albedoCHM)){
  if(albedoCHM$CHM[i] > 0 && albedoCHM$CHM[i] < 1){
    albedoCHM$Size[i] <- "<1 m"
  }
  if(albedoCHM$CHM[i] > 1 && albedoCHM$CHM[i] < 2){
    albedoCHM$Size[i] <- "1-2 m"
  }
  if(albedoCHM$CHM[i] > 2 && albedoCHM$CHM[i] < 3){
    albedoCHM$Size[i] <- "2-3 m"
  }
  if(albedoCHM$CHM[i] > 3 && albedoCHM$CHM[i] < 4){
    albedoCHM$Size[i] <- "3-4 m"
  }
  if(albedoCHM$CHM[i] > 4 && albedoCHM$CHM[i] < 5){
    albedoCHM$Size[i] <- "4-5 m"
  }
  if(albedoCHM$CHM[i] > 5 && albedoCHM$CHM[i] < 6){
    albedoCHM$Size[i] <- "5-6 m"
  }
  if(albedoCHM$CHM[i] > 6 && albedoCHM$CHM[i] < 7){
    albedoCHM$Size[i] <- "6-7 m"
  }
  if(albedoCHM$CHM[i] > 7 && albedoCHM$CHM[i] < 8){
    albedoCHM$Size[i] <- "7-8 m"
  }
  if(albedoCHM$CHM[i] > 8){
    albedoCHM$Size[i] <- "8+ m"
  }
  
}

cols2 = viridis(n = 12)
# make a plot for CHM and winter albedo
nonTSK <- mean(albedoWT[which(albedoWT$PFT>5),15])
albedoCHM$Size <- as.factor(albedoCHM$Size)
albedoCHM$Size <- factor(albedoCHM$Size, levels = (levels(albedoCHM$Size)))
winter.chm = albedoCHM %>% 
  mutate(levels = cut(CHM, seq(0, 12, 1))) %>%
  drop_na(levels, winter_albedo) %>%
  ggplot(aes(x=levels, y=winter_albedo)) +
  geom_boxplot(outlier.shape = NA, aes(fill = levels)) +
  #geom_smooth(data = albedoCHM, method = "lm", aes(x = CHM, y = winter_albedo), color = "blue")+
  #geom_hline(yintercept = nonTSK, linetype="dashed") +
  xlab("Canopy Height") + ylab("Winter Albedo") + ylim(0, 1) + coord_cartesian(xlim = c(1, 12)) +
  #scale_fill_manual(values = c("red", "blue"), name = "", labels = c("Short", "Tall")) +
  #scale_x_discrete(labels = pfts[-c(6:8)]) +
  #geom_segment(aes(x = 2.05, xend = 2.25, y = 0.0, yend = 0.0),  size = 1, 
  #linetype = "dashed") +
  #annotate("text", x = 2.3,  y = 0.015, 
  #label = stringr::str_wrap("Non-Woody Albedo", width = 10), hjust=0) +
  scale_fill_manual(values = cols2)+
  guides(fill=guide_legend(ncol = 7)) +
  theme(legend.position = "none",
        #c(0.63, 0.09), 
        legend.title = element_blank(), 
        legend.key.size = unit(0.75, 'cm'),
        legend.text = element_text(size = 12)) +
  theme(axis.text = element_text(size=15, color = 'black'),
        axis.text.x = element_text(, angle = 60, hjust = 1),
        axis.title=element_text(size=18)) +
  theme(axis.line = element_line(colour = "black"),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

#summer albedo and CHM
nonTSK <- mean(albedoSM[which(albedoSM$PFT>5),15])
summer.chm = albedoCHM %>% 
  mutate(levels = cut(CHM, seq(0, 12, 1))) %>%
  drop_na(levels, summer_albedo) %>%
  ggplot(aes(x=levels, y=summer_albedo)) +
  geom_boxplot(outlier.shape = NA, aes(fill = levels)) +
  #geom_smooth(data = albedoCHM, method = "lm", aes(x = CHM, y = summer_albedo), color = "blue")+
  #geom_hline(yintercept = nonTSK, linetype="dashed") +
  xlab("Canopy Height") + ylab("Summer Albedo") + ylim(0, 0.22) + coord_cartesian(xlim = c(1, 12)) +
  #scale_fill_manual(values = c("red", "blue"), name = "", labels = c("Short", "Tall")) +
  #scale_x_discrete(labels = pfts[-c(6:8)]) +
  #geom_segment(aes(x = 2.5, xend = 2.95, y = 0.06, yend = 0.06),  size = 1, 
  #linetype = "dashed") +
  #annotate("text", x = 2.5,  y = 0.045, 
  #label = stringr::str_wrap("Non-Woody Albedo", width = 10), hjust=0) +
  scale_fill_manual(values = cols2)+
  guides(fill=guide_legend(ncol = 5)) +
  theme(legend.position = c(0.64, 0.2), 
        legend.title = element_blank(), 
        legend.key.size = unit(0.75, 'cm'),
        legend.text = element_text(size = 12)) +
  theme(axis.text = element_text(size=15, color = 'black'),
        axis.text.x = element_text(, angle = 60, hjust = 1),
        axis.title=element_text(size=18)) +
  theme(axis.line = element_line(colour = "black"),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

summary(lm(data = albedoCHM, formula = summer_albedo ~ CHM))
summary(lm(data = albedoCHM, formula = winter_albedo ~ CHM))
summary(rq(summer_albedo ~ CHM, data = albedoCHM, tau = 0.5))
summary(rq(winter_albedo ~ CHM, data = albedoCHM, tau = 0.5))


plot_grid(summer.pft, winter.pft, summer.chm, winter.chm, labels = c("A", "B", "C", "D"))
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
chm.winter = ggplot(data = albedoSTR, aes(x=CHM, y=winter_albedo)) +
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
chm.winter
pdfNAME = paste0(outDIR, "/", "chm_albedo_winter.pdf")
ggsave(pdfNAME, plot = last_plot(), width = 12, height = 10, units = 'cm')

# make a structure-albedo plot for summer
chm.summer = ggplot(data = albedoSTR, aes(x=CHM, y=summer_albedo)) +
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
chm.summer
pdfNAME = paste0(outDIR, "/", "chm_albedo_summer.pdf")
ggsave(pdfNAME, plot = last_plot(), width = 12, height = 10, units = 'cm')

plot_grid(winter.pft, summerr.pft, chm.winter, chm.summer, labels = c("A", "B", "C", "D"))

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
  xlab("Dominant PFT") + ylab("CHM") + ylim(0, 8) +
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
  ggplot(data=pftDF, aes(x=fcover, y=albedoWT)) +
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
    xlab("fCover") + ylab("Winter Albedo") +
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
  
  pdfNAME = paste0(outDIR, "/", pft, "_fcover_winter_albedo_regression.pdf")
  ggsave(pdfNAME, plot = last_plot(), width = 10, height = 9, units = 'cm')
  
  
  chm.regression = ggplot(data=pftDF, aes(x=chm, y=albedoWT)) +
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
    xlab("CHM") + ylab("Winter Albedo") + 
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
          panel.grid.minor = element_blank(), 
          strip.background = element_rect(fill = "white", color = "black"),
          strip.text = element_text(color = "black", size = 12))
  chm.regression
  pdfNAME = paste0(outDIR, "/", pft, "_chm_winter_albedo_regression.pdf")
  ggsave(pdfNAME, plot = last_plot(), width = 10, height = 9, units = 'cm')
}

#plot all on one plot
combined_df <- data.frame()

for (pft in pfts[-c(6:8)]) {
  # Determine pdf id
  pftID <- which(pfts == pft)
  # Extract pft data frame
  pftDF <- fcoverALBEDO[which(fcoverALBEDO$PFT == pftID),]
  # Keep data where pft fcover > 0.5
  pftDF <- data.frame(pftDF[which(names(pftDF) == pft)], pftDF[-c(1:14)])
  names(pftDF) <- c('fcover', 'chm', 'albedoWT', 'albedoSM')
  pftDF <- na.omit(pftDF)
  # Add a column to indicate the pft
  pftDF$pft <- pft
  # Combine with the main data frame
  combined_df <- rbind(combined_df, pftDF)
}

desired_order <-c("ET", "DT", "DTSA", "DTSW", "DLS", "DG", "WG", "MO", "LI", "NVS")
combined_df$pft <- factor(combined_df$pft, levels = desired_order)


ggplot(data = combined_df, aes(x = fcover, y = albedoWT)) +
  geom_hex(aes(fill = stat(log(count))), bins = 50, breaks = log(c(0, 1, 2, 4, 6))) +
  scale_fill_viridis_c() +
  geom_smooth(color = 'red', method = "lm", formula = my.formula) +
  stat_poly_eq(formula = my.formula, 
               aes(label = paste(..eq.label.., ..rr.label.., ..p.value.label.., sep = "*`,`~")), 
               label.x.npc = 'right', label.y = 0, parse = TRUE) +
  xlab("Fractional Cover (fCover)") + ylab("Winter Albedo") +
  theme(legend.position = "", 
        legend.title = element_blank(), 
        legend.key.size = unit(0.4, 'cm'),
        legend.text = element_text(size = 10),
        legend.spacing.x = unit(0.2, 'cm'),
        legend.spacing.y = unit(0.2, 'cm')) +
  theme(axis.text = element_text(size = 12, color = 'black'),
        axis.title = element_text(size = 12)) +
  theme(axis.line = element_line(colour = "black"),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        strip.background = element_blank(),
        strip.text = element_text(color = "black", size = 12)) +
  facet_wrap(~ pft)

# Plot the data using facet_wrap for CHM
ggplot(data = combined_df, aes(x = chm, y = albedoWT)) +
  geom_hex(aes(fill = stat(log(count))), bins = 50, breaks = log(c(0, 1, 2, 4, 6))) +
  scale_fill_viridis_c() +
  geom_smooth(color = 'red', method = "lm", formula = my.formula) +
  stat_poly_eq(formula = my.formula, 
               aes(label = paste(..eq.label.., ..rr.label.., ..p.value.label.., sep = "*`,`~")), 
               label.x.npc = 'right', label.y= 0, parse = TRUE) +
  xlab("Canopy Height (m)") + ylab("Winter Albedo") + 
  theme(legend.position = "", 
        legend.title = element_blank(), 
        legend.key.size = unit(0.4, 'cm'),
        legend.text = element_text(size = 10),
        legend.spacing.x = unit(0.2, 'cm'),
        legend.spacing.y = unit(0.2, 'cm')) +
  theme(axis.text = element_text(size = 12, color = 'black'),
        axis.title = element_text(size = 12)) +
  theme(axis.line = element_line(colour = "black"),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        strip.background = element_blank(),
        strip.text = element_text(color = "black", size = 12)) +
  facet_wrap(~ pft)
