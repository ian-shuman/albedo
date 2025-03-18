###########################################################################################
#
#        This script calculates topographic wetness index (TWI) from a DEM for Shumen et al. (in prep)
#
#    --- Last updated:  2025.03.18 By Ian Shuman <ins2109@columbia.edu>
###########################################################################################


#******************** close all devices and delete all variables *************************#
rm(list=ls(all=TRUE))   # clear workspace
graphics.off()          # close any open graphics
closeAllConnections()   # close any open connections to files
dlm <- .Platform$file.sep # <--- What is the platform specific delimiter?
options(warn = -1)
#*****************************************************************************************#

#************************* define user parameters and load data **************************#
data.dir <- '~/Downloads/albedo/Data'
out.dir <- '/Users/anncrumlish/Downloads/albedo/data/'
topo <- rast(paste0(data.dir, 'council_watershed_topo_30m.tif')) 
#*****************************************************************************************#


#Define twi functions by modifying the R package "whitebox" (Lindsay, 2016; Q. Wu & Brown, 2022)
upslope <- function (dem, log = TRUE, atb = FALSE, deg = 0.12, fill.sinks = TRUE) 
{
  if (!all.equal(xres(dem), yres(dem))) {
    stop("Raster has differing x and y cell resolutions. Check that it is in a projected coordinate system (e.g. UTM) and use raster::projectRaster to reproject to one if not. Otherwise consider using raster::resample")
  }
  if (fill.sinks) {
    capture.output(dem <- invisible(raster::setValues(dem, topmodel::sinkfill(raster::as.matrix(dem), res = xres(dem), degree = deg))))
  }
  topidx <- topmodel::topidx(raster::as.matrix(dem), res = xres(dem))
  a <- raster::setValues(dem, topidx$area)
  if (log) {
    a <- log(a)
  }
  if (atb) {
    atb <- raster::setValues(dem, topidx$atb)
    a <- addLayer(a, atb)
    names(a) <- c("a", "atb")
  }
  return(a)
}

create_layers <- function (dem, fill.sinks = TRUE, deg = 0.1) 
{
  layers <- stack(dem)
  message("Building upslope areas...")
  a.atb <- upslope(dem, atb = TRUE, fill.sinks = fill.sinks, deg = deg)
  layers <- addLayer(layers, a.atb)
  names(layers) <- c("filled.elevations", "upslope.area", "twi")
  return(layers)
}

#Calculate TWI
layers <- create_layers(topo)
twi.man <- log(layers$upslope.area / tan(slope / 180))

#Save raster
writeRaster(twi.man, filename=file.path(out.dir, "twi.tif"), format="GTiff", overwrite=TRUE)
