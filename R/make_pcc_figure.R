###########################################################################################
#
#                  quick process and plot weather station data
#
#    --- Last updated:  2020.01.16 By Daryl Yang <dediyang@bnl.gov>
###########################################################################################

#-------------------- close all devices and delete all variables. ------------------------#
rm(list=ls(all=TRUE))   # clear workspace
graphics.off()          # close any open graphics
closeAllConnections()   # close any open connections to files
dlm <- .Platform$file.sep # <--- What is the platform specific delimiter?
#-----------------------------------------------------------------------------------------#

#---------------------------- load required libraries ------------------------------------#
# Info: Loads required R libraries and warns if package is not availible.
ok <- require(ggplot2) ; if (! ok) 
  #stop("*** Package pls is not available.  This is needed for model optimization ***")
  install.packages("ggplot2")
ok <- require(stats) ; if (! ok) 
  install.packages("stats")
ok <- require(ggcorrplot) ; if (! ok) 
  install.packages("ggcorrplot")

# Script options
options(digits.secs = 3)
options(digits = 5)
#-----------------------------------------------------------------------------------------#

#------------------------------ set output directory -------------------------------------#
### Set ouput directory
out.dir <- file.path('G:\\MyWork\\Albedo_Scaling\\Figures')

### Create output folders
if (! file.exists(out.dir)) dir.create(out.dir,recursive=TRUE)
#-----------------------------------------------------------------------------------------#

#------------------------------ define user parameters -----------------------------------#
# number of lines to skip in the csv file
skip <- 0
# column name for time
tvar <- 'X'
# environmental variable to display
var <- 'Avg'
#-----------------------------------------------------------------------------------------#

#-------------------------------- load in data--------------------------------------------#
data.dir <- 'G:\\MyWork\\Albedo_Scaling\\Data_Analysis\\map_analysis\\map_data_v2.csv'
data.original <- read.csv(data.dir, header = TRUE, skip = skip)

data <- data.original[, -1]
data <- na.omit(data)

corr = round(cor(data), 1)
p.mat = cor_pmat(data)


ggcor = ggcorrplot(corr, type = "lower", lab = TRUE, lab_size = 2, method="square", 
                   colors = c("blue", "white", "red"), p.mat = p.mat, sig.level = 0.05/22, outline.col = "white", ggtheme = ggplot2::theme_gray) +
  theme(axis.text.x = element_text(color="black", size=9, angle=45, hjust = 1)) +
  theme(axis.text.y = element_text(color="black", size=9, angle=0)) +
  theme(legend.position = c(0.1, 0.7), legend.title = element_blank(), legend.background = element_blank(), legend.key = element_blank())



png.name = paste0(out.dir, "\\", 'pcc_corr_plot.png')
ggsave(png.name, plot = last_plot(), width = 20, height = 20, units = 'cm')






























