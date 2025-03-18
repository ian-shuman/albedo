# Overview

This repository contains code written by I Shuman and D Yang for Shuman et al. (in prep). The code processes a Landsat-derived albedo product developed by Erb et al. (2022) to create a high resolution time series of broadband white sky albedo across a 33.7 square kilometer site near Council, AK on the Seward Peninsula. Additional code analyzes the timeseries to determine how fine scale plant functional type fractional cover and canopy height affects spatiotemporal variation in surface albedo and implements a random forest model to determine the importance of vegetation characteristics relative to topography when predicting surface albedo. 

# License

This repository holds an MIT License, as described in the LICENSE file.

# Software versions

This repository is primarily built in the R environment using R version 4.4.0.

# Package versions
* `caTools` v. 1.18.2
* `cowplot` v. 1.1.3
* `dplyr` v. 1.1.4
* `ggbreak` v. 0.1.2
* `ggfun` v. 0.1.5
* `ggpattern` v. 1.1.1
* `ggplot2` v. 3.5.1
* `ggpmisc` v. 0.5.5
* `GMCM` v. 1.4
* `grid` v. 4.4.0
* `gtable` v. 0.3.5
* `haven` v. 2.5.4
* `multcompView` v. 0.1.10
* `pdp` v. 0.8.1
* `pls` v. 2.8.3
* `quantreg` v. 5.98
* `randomForest` v. 4.7.1.1
* `raster` v. 3.6.26
* `readr` v. 2.1.5
* `reshape2` v. 1.4.4
* `scatterpie` v. 0.2.3
* `spatialEco` v. 2.0.2
* `spectratrait` v. 1.2.1
* `strucchange` v. 1.5.3
* `terra` v. 1.7.78
* `tidyr` v. 1.3.1
* `tidyverse` v. 2.0.0 
* `viridis` v. 0.6.5
* `viridisLite` v. 0.4.2

# Directory structure

* IDL: The IDL subdirectory contains IDL code developed by D Yang for processing Erb et al.'s (2022) surface albedo product for use in the analyses of Shuman et al. (in prep). Description of each file is provided in the following subsection, scripts should be run in the order presented there. 
* R: The bulk of the repository is housed in the R subdirectory. These are the scripts used for analyzing the surface albedo product with respect to vegetation characteristics and topography. The scripts within the R subdirectory are ordered according to the figures presented in Shuman et al. (in prep). Descriptions of each figure's script are presented in the following subsection, and scripts should be run in squential order (i.e. 1 > 2 > 3 > 4 > 56).
  * deprecated: Contains old versions of code presented in the figure scripts with preliminary or additional analyses that are not presented in Shuman et al. (in prep). These scripts are for the authors' reference only.
* ROI: Contains the region of interest (ROI) used to define the 33.7 square kilometer site near Council, AK. Three file formats (.tif, .enp, and .hdr) are provided.

## Code organization

* **IDL/step1_landsat_clip.pro**: This script clips Landsat surface reflectance images to the Council, AK site using shapefiles. 

  * Inputs: Unzipped landsate surface reflectance images (GeoTIFF, .tif) from Erb et al. (2022) and shapefiles for the region of interest (ROI, .shp). The ROI used in Shuman et al. (in prep) can be found in the ROI subdirectory.  
  
  * Outputs: Clipped Landsat surface reflectance files (GeoTIFF, .tif).

 * **IDL/step2_layerstack_images.pro**: This script builds a layer-stack of time series data from the albedo product. All Landsat images of the Council site acquired between 2015 and 2019 are stacked in chronological order by day of year (DOY). 

   * Inputs: Clipped Landsat surface reflectance files (GeoTIFF, .tif) from IDL/step1_landsat_clip.pro.  
  
   * Outputs: Layer-stack of Landsat-derived surface albedo rasters (GeoTIFF, .tif). This time series is a "representative year," in which all images collected between 2015 and 2019 have been combined sequentially by DOY as if in a single year because there were not enough images collected in year individual year to analyze temporal variation in albedo. 

* **IDL/step3_time_series_smooth.pro**: This script smooths the albedo time series using a Savitzky-Golay filter and removes outlier values (e.g. cloudy pixels) based on NDVI. 

   * Inputs: Layer-stack of Landsat-derived surface albedo rasters (GeoTIFF, .tif) from IDL/step2_layerstack_images.pro.
  
   * Outputs: Smoothed and filtered layer-stack of Landsat-derived surface albedo rasters (GeoTIFF, .tif).
 
* **IDL/generate_albedo_datasets_for_analysis.pro**: This script summarizes the albedo timeseries over the winter and summer seasons by calculating the mean winter (DOY 1 - 121) albedo and mean summer (DOY 182 – 243) albedo of images collected during those respective date spans. 

   * Inputs: Smoothed and filtered layer-stack of Landsat-derived surface albedo rasters (GeoTIFF, .tif) from IDL/step3_time_series_smooth.pro.
  
   * Outputs: Two raster images summarizing the mean winter surface albedo and mean summer surface albedo across the Council, AK site (GeoTIFF, .tif).

 * **IDL/albedo_image_statistics.pro**: This script combines the seasonal albedo data with third-party plant functional type and topographic data into a single tabular data file for further analyses in R. 

  * Inputs:
   * Mean summer albedo raster from IDL/generate_albedo_datasets_for_analysis.pro (GeoTIFF, .tif)
   * Mean winter albedo raster from IDL/generate_albedo_datasets_for_analysis.pro (GeoTIFF, .tif)
   * 30 meter resolution plant functional type (PFT) fractional cover raster clipped to the Council, AK site; from Yang et al. (2024) (GeoTIFF, .tif)
   * 30 meter resolution topography rasters clipped to the Council, AK site; from Singhania et al. (2022) (GeoTIFF, .tif); includes the LiDAR collected digital elevation model (DEM), slope, aspect, and hillslope
   * **Optional** raster file describing the transition date from winter to summer albedo based on breakpoint regression; see R/breakpoint_regression for calculation
  
  * Outputs: Single tabular file (.csv) summarizing mean summer albedo, mean winter albedo, DEM, slope, aspect, hillslope, and transition date at 30 meter spatial resolution across the Council, AK site
  

  
# Third-Party Data Availability

## Landsat Surface Albedo Timeseries

* Main article:
  * Erb, A. M., Li, Z., Sun, Q., Paynter, I., Wang, Z., & Schaaf, C. (2022). Evaluation of the Landsat-8 Albedo Product across the Circumpolar Domain. Remote Sensing, 14(21), 5320. https://doi.org/10.3390/rs14215320

* More details:
  * Li, Z., Erb, A., Sun, Q., Liu, Y., Shuai, Y., Wang, Z., Boucher, P., & Schaaf, C. (2018). Preliminary assessment of 20-m surface albedo retrievals from sentinel-2A surface reflectance and MODIS/VIIRS surface anisotropy measures. Remote Sensing of Environment, 217, 352–365. https://doi.org/10.1016/j.rse.2018.08.025
  * Wang, Z., Erb, A. M., Schaaf, C. B., Sun, Q., Liu, Y., Yang, Y., Shuai, Y., Casey, K. A., & Román, M. O. (2016). Early spring post-fire snow albedo dynamics in high latitude boreal forests using Landsat-8 OLI data. Remote Sensing of Environment, 185, 71–83. https://doi.org/10.1016/j.rse.2016.02.059
  * Shuai, Y., Masek, J. G., Gao, F., & Schaaf, C. B. (2011). An algorithm for the retrieval of 30-m snow-free albedo from Landsat surface reflectance and MODIS BRDF. Remote Sensing of Environment, 115(9), 2204–2216. https://doi.org/10.1016/j.rse.2011.04.019
    
## Plant Functional Type Fractional Cover

* Raw data:
  * Yang D; Serbin S (2024): Maps of plant functional type (PFT), PFT fractional cover, and uncertainty derived from AVIRIS-NG data, 2019, Seward Peninsula. Next-Generation Ecosystem Experiments (NGEE) Arctic, ESS-DIVE repository. Dataset. doi:10.15485/2441506 accessed via https://data.ess-dive.lbl.gov/datasets/doi:10.15485/2441506 on 2025-03-18
  
* More details:
  * Yang, D., Morrison, B. D., Hanston, W., McMahon, A., Baskaran, L., Hayes, D. J., Miller, C. E., & Serbin, S. P. (2023). Integrating very-high-resolution UAS data and airborne imaging spectroscopy to map the fractional composition of Arctic plant functional types in Western Alaska. Remote Sensing of Environment, 286, 113430. https://doi.org/10.1016/j.rse.2022.113430
 
## Topography and Canopy Height

* Raw data:
  * Singhania A; Glennie C; Fernandez-Diaz J; Hauser D (2023): National Center for Airborne Laser Mapping (NCALM) LiDAR, Imagery, and DEM data from five NGEE Arctic Sites, Seward Peninsula, Alaska, August 2021. Next-Generation Ecosystem Experiments (NGEE) Arctic, ESS-DIVE repository. Dataset. doi:10.5440/1832016 accessed via https://data.ess-dive.lbl.gov/datasets/doi:10.5440/1832016 on 2025-03-18
    
