library(ggplot2)

out.dir <- 'G:\\MyWork\\Albedo_Scaling\\Figures'

data.dir <- 'G:\\MyWork\\Albedo_Scaling\\Data_Processing\\Step5_PFT_Albedo_Analysis\\PFT_Summer_Albedo'
file.list <- list.files(data.dir, pattern = '.csv', recursive = TRUE, full.names = TRUE)
#this script seems like an old version which creates a boxplot of albedo for pfts, which is also created in albedo_analyzer_v1.R
data.combn <- c()
for (file in file.list)
{
  data.original <- read.csv(file, header = T)
  data.size <- length(data.original)
  row.name <- rep(names(data.original), data.size)
  data.update <- cbind(row.name, data.original)
  names(data.update) <- c('pft', 'albedo')
  
  data.combn <- rbind(data.combn, data.update)
}

data.combn$albedo[data.combn$albedo > 2500] <- 0
data.combn <- na.omit(data.combn)
data.combn$albedo <- data.combn$albedo/10000.0

ggplot(data = data.combn, aes(x = reorder(pft, albedo, median), y = albedo)) + 
  geom_boxplot(aes(fill = reorder(pft, albedo, median)), outlier.size = -1, varwidth = FALSE, position =  position_dodge2(width = 1, preserve = "single")) + 
  scale_fill_brewer(palette="Spectral", direction = -1) +
  ylim(c(0, 0.3)) + labs(x = 'Plant Functional Type', y = 'White Sky Albedo') + 
  theme(legend.position = c(0.7, 0.2), legend.direction = "horizontal", legend.title = element_blank(), legend.text = element_text(colour="black", size=12), legend.background = element_blank(), legend.key = element_blank()) +
  theme(axis.text = element_text(size=12, angle = 0), axis.title=element_text(size=12,face="bold"))
  

png.name <- paste0(out.dir, '\\', 'Boxplot_PFT_Summer_Albedo.png')
ggsave(png.name, plot = last_plot(), width = 20, height = 15, units = 'cm')





















