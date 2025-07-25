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
* `Metrics` v. 0.1.4
* `multcompView` v. 0.1.10
* `ncdf4` v. 1.24
* `pdp` v. 0.8.1
* `pls` v. 2.8.3
* `quantreg` v. 5.98
* `randomForest` v. 4.7.1.1
* `raster` v. 3.6.26
* `readr` v. 2.1.5
* `reshape2` v. 1.4.4
* `rnaturalearth` v.1.0.1
* `rnaturalearthdata` v. 1.0.0
* `scatterpie` v. 0.2.3
* `sf` v. 1.0.16
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
* R: The bulk of the repository is housed in the R subdirectory. These are the scripts used for analyzing the surface albedo product with respect to vegetation characteristics and topography. The scripts within the R subdirectory are ordered according to the figures presented in Shuman et al. (in prep). Descriptions of each figure's script are presented in the following subsection, and scripts should be run in the order presented there.
  * deprecated: Contains old versions of code presented in the figure scripts with preliminary or additional analyses that are not presented in Shuman et al. (in prep). These scripts are for the authors' reference only.
* ROI: Contains the region of interest (ROI) used to define the 33.7 square kilometer site near Council, AK. Three file formats (.tif, .enp, and .hdr) are provided.

## Code organization

To fully replicate the analyses used in Shuman et al. (in prep), code should be run in the order shown here.

IDL Subdirectory:

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
 
R Subdirectory:
 
* **R/make_pcc_figure.R**: This script generates a quick correlation plot to visualize the correlations between vegetation, topography, and albedo variables. 

   * Inputs: Tabular file (.csv) of variables similar to that produced by IDL/albedo_image_statistics.pro.
  
   * Outputs: No data outputs saved. A correlation plot figure is saved that was used to check correlations within the data, but this figure is not presented in the final version of the Shuman et al. (in prep) manuscript.

* **R/breakpoint_regression.R**: This script uses breakpoint regression to calculate the DOY on which each pixel in the Council, AK site begins and ends the transition from a "winter albedo regime" (i.e. dominated by high albedo snow) to a "summer albedo regime" (i.e. dominated by low albedo leafy vegetation). 

   * Inputs: Time series layer-stack of white-sky albedo raster images from IDL/step3_time_series_smooth.pro.
  
   * Outputs: A raster image representing the completion DOY for the spring albedo transition (second breakpoint) for each pixel across the Council, AK site (GeoTIFF, .tif). Low-light conditions and frequent cloud cover during early spring created significant uncertainty in the identification of the first breakpoint, so while it is calculated in this script, it is not used in the subsequent analyses of Shuman et al. (in prep).
 
* **R/calculate_topography.R**: This script calculates the topographic wetness index (TPI), terrain ruggedness index (TRI), and topographic position index (TPI). TWI is calculated by modifying the "whitebox" R package and is used as a proxy for moisture when predicting surface albeodo in Shuman et a. (in prep). TRI and TPI are calculated using the "terrain" and "focal" functions in the raster or terra R packages.

   * Inputs:
     * For TWI: 30 meter DEM and slope raster images clipped to the Council, AK site; from Singhania et al. (2022) (GeoTIFF, .tif)
     * For TRI and TPI: Tabular file (.csv) of variables similar to that produced by IDL/albedo_image_statistics.pro, must contain at least the DEM and slope values
  
   * Outputs: 30 meter TWI, TRI, and TPI raster images (GeoTIFF, .tif) 

* **R/Fig1_script.R**: This script loads and re-processes the vegetation, topography, and albedo data to optimize it for random forest implementation in the R environment, and also generates the plots used in Figures 1 and S8 of Shuman et al. (in prep). Re-processing includes bi-linear resampling to make sure all data sources are exactly the same extent and resoltion, as well as identifying all 30 meter pixels which contain >75% fractional cover for a single PFT (i.e. are dominated by one PFT).   

   * Inputs: Tabular file (.csv) of variables similar to that produced by IDL/albedo_image_statistics.pro.
  
   * Outputs:
     * An .RData object containing a dataframe with all variables analyzed in Shuman et al. (in prep), including for each 30m pixel:
       * X coordinate
       * Y coordinate
       * Canopy Height (m)
       * Elevation (m)
       * Topographic Position Index
       * Terrain Ruggedness Index
       * Slope (degrees)
       * Aspect (degrees)
       * Topographic Wetness Index
       * Mean summer albedo
       * Mean winter albedo
       * Evergreen tree fractional cover 
       * Alder fractional cover 
       * Willow fractional cover 
       * Other deciduous tall tree fractional cover 
       * Low shrub fractional cover 
       * Dwarf shrub fractional cover 
       * Evergreen shrub fractional cover 
       * Forb fractional cover 
       * Dry graminoid fractional cover 
       * Wet graminoid fractional cover 
       * Moss fractional cover 
       * Lichen fractional cover 
       * Non-photosynthetic vegetation fractional cover
       * PFT with highest fractional cover in pixel
       * Start date of spring albedo transition (DOY)
       * End date of spring albedo transition (DOY)
     * Plots representing panels (A), (B), and (C) of Figure 1 in Shuman et al. (in prep) depicting the location of the Council, AK site (A), the PFT with the highest fractional cover for each pixel across the site (B), and the variation in winter albedo across the site (C).  
  
* **R/Fig2_script.R**: This script analyzes variation in summer and winter albedo versus multiple formulations of plant functional type and canopy height and uses that information to create Figures 2, S4-S7, and S11 from Shuman et al. (in prep). Significance tests comparing seasonal albedo between pixels dominated by different PFTs and by different 1 meter bins of canopy height are conducted, but are not presented in Shuman et al. (in prep). 

   * Inputs: Tabular file (.csv) of variables similar to that produced by IDL/albedo_image_statistics.pro. The file is reprocessed in this script, the .RData file developed in R/Fig1_script.R is not used here. 
  
   * Outputs: No data outputs saved. All plots representing Figures 2, S4-S7, and S11 in Shuman et al. (in prep) are saved to a user-specified directory.
 
* **R/Fig3_script.R**: This script visualizes the surface albedo time series by dominant PFT (>75% fractional cover) and canopy height (binned in 1 meter increments) using a 15 day moving-window average. These visualizations are used to create Figure 3 in Shuman et al., (in prep). 

   * Inputs:
     * Smoothed and filtered layer-stack of Landsat-derived surface albedo rasters (30 m resolution); from IDL/step3_time_series_smooth.pro (.dat) 
     * Plant functional type (PFT) fractional cover map clipped to the Council, AK site (30 m resolution); from Yang et al. (2024) (.dat)
     * Raster image of canopy height clipped to the Council, AK site (30 m resolution); from processing Singhania et al. (2022) (GeoTIFF, .tif)
     * 30 meter resolution raster image of the breakpoint regression identified date of completion for the spring albedo transition across the Council, AK site (30 m resolution); from R/breakpoint_regression.R (GeoTIFF, .tif)
  
   * Outputs: No data outputs saved. All panels represented in Figure 3 in Shuman et al. (in prep) are plotted but not explicitly saved.
 
 * **R/Fig4_script.R**: This script plots variation in the completion date of the spring albedo transition (calculated in R/breakpoint_regression.R) by dominant PFT (>75% fractional cover) and by canopy height (in 1 m bins). Significance tests are conduced on the differences in transition completion DOY between dominant PFTs, and the fractional composition of dominant PFTs is calculated for each canopy height bin. These analyses are then used to create Figure 4 in Shuman et al. (in prep). This script also ingests Daymet snow water equivalent data (SWE) to produce figure S12.   

   * Inputs: Data frame which is the output of R/Fig1_script.R. Alternatively, 30 meter raster images (GeoTIFF, .tif) of mean summer albedo, mean winter albedo, canopy height, aspect, slope, DEM, TPI, TRI, fractional cover of each PFT, dominant PFT at a pixel, TWI, and the first and second (optional) DOY breakpoints of the spring albedo transition which are generated as outputs in the scripts above can be used as inputs to generate a nearly identical dataframe in this script. Daymet analyses require Daymet data (see third party data) and a 30 meter raster image (GeoTIFF, .tif) for cropping Daymet to the extent of the Council site. 

   * Outputs: No data outputs saved. All plots representing Figures 4 and S12 in Shuman et al. (in prep) are saved to a user-specified directory.

 * **R/Fig56_script.R**: This script uses a random forest model to predict winter and summer surface albedo from the vegetation composition, vegetation structure, topography, and moisture variables. The output of the random forest model is analyzed for variable importance, partial least squares regression (PLSR) is used to determine the goodness of fit for the random forest models, and partial dependence plots are created for each predictor variable. The variable importance analyses are used to create Figure 5, the PLSR check are used to create Figure S6, and the partial dependence plots are used to creat Figures 6 and S7 in Shuman et al. (in prep).

   * Inputs: Data frame which is the output of R/Fig1_script.R. 

   * Outputs: No data outputs saved. All panels represented in Figures 5-6 and S6-S7 in Shuman et al. (in prep) are plotted and saved to a user-specified directory. We recommend saving the workspace image of the R environment to avoid re-running this computationally intensive script.

* **R/FigS2_script.R**: This script describes the summarizing, filtering, and smoothing process done in IDL to create supplemental figure S2. 

   * Inputs: A raster file (.dat and .hdr) representing the raw albedo observations produced by IDL/step2_layerstack_images.pro and a raster file (.dat and .hdr) representing the filtered and smoothed albedo observatiosn produced by IDL/step3_time_series_smooth.pro.
  
   * Outputs: No data outputs are saved. All plots representing Figure S2 are saved to a user-specified directory.
  
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

## Snow Water Equivalent (SWE)

* Raw data:
  * Thornton, M. M., Shrestha, R., Wei, Y., Thornton, P. E., & Kao, S.-C. (2022). JDaymet: Daily Surface Weather Data on a 1-km Grid for North America, Version 4 R1 (Version 4.1, p. 0 MB) [netCDF]. ORNL Distributed Active Archive Center. https://doi.org/10.3334/ORNLDAAC/2129
    
