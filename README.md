# Overview

This repository contains code written by I Shuman and D Yang for Shuman et al. (in prep). The code processes a Landsat-derived albedo product developed by Erb et al. (2022) to create a high resolution time series of broadband white sky albedo across a 33 square kilometer site near Council, AK on the Seward Peninsula. Additional code analyzes the timeseries to determine how fine scale plant functional type fractional cover and canopy height affects spatiotemporal variation in surface albedo and implements a random forest model to determine the importance of vegetation characteristics relative to topography when predicting surface albedo. 

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
* `terra` v. 1.7.78
* `tidyr` v. 1.3.1
* `tidyverse` v. 2.0.0 
* `viridis` v. 0.6.5
* `viridisLite` v. 0.4.2

# Directory structure

* IDL: The IDL subdirectory contains IDL code developed by D Yang for processing Erb et al.'s (2022) surface albedo product for use in the analyses of Shuman et al. (in prep). Description of each file is provided in the following subsection. 
* R: The bulk of the repository is housed in the R subdirectory. These are the scripts used for analyzing the surface albedo product with respect to vegetation characteristics and topography. The scripts within the R subdirectory are ordered according to the figures presented in Shuman et al. (in prep). Descriptions of each figure's script are presented in the following subsection, and scripts should be run in squential order (i.e. 1 > 2 > 3 > 4 > 56).
  * deprecated: Contains old versions of code presented in the figure scripts with preliminary or additional analyses that are not presented in Shuman et al. (in prep). These scripts are for the authors' reference only. 

## Code organization

* **0.SummaryFigures.R**: This script produces the figures given in the supplement of Shuman et al. (in prep) that summarize the response and drive variables in space

  * Inputs: Environmental covariates (xdata), and community- and biome-level response data (ydata) from the GJAMDATA/ directory. **The inputs that are loaded at the top of the script are saved from step 4. Step 0 is the only out-of-sequence step because it is not part of the analysis and only serves as data (not analysis) visualization**
  
  * Outputs: None saved. Figures of reconstructed vegetation and environmental covariates produced
  
* **0.SummaryFigures_OOS.R**: This script produces the figures given in the supplement of Shuman et al. (in prep) that summarize the response and driver variables of withheld data for performing out-of-sample validation. The script is an exact mirror of 0.SummaryFigures.R but for out-of-sample validation data
  
  * Inputs: Community- and biome-level out-of-sample data from the GJAMDATA/Withheld for Validation/ subdirectory. **The inputs that are loaded at the top of the script are saved from step 5. Step 0 is the only out-of-sequence step because it is not part of the analysis and only serves as data (not analysis) visualization**
  
  * Outputs: None saved. Figures of reconstructed vegetation and environmental covariates produced

* **1.Process.R**: Takes tables of response (ydata) and driver (xdata) variables and formats them according to the format required by the `gjam()` function. The script conducts simple data checks before saving the output in an RData file.

  * Inputs: xdata (drivers) and ydata (response variables) from CSVs specific to each management area as described in the methods of Shuman et al. (in prep) and in the Input Data section of the README
  * Outputs: GJAMDATA/processed_xydata.RData: a RData file with the xdata and ydata databases with all management areas combined
    * xdata: 78224 observations of 17 variables
      * long: longitude (decimal degrees, EPSG4326)
      * lat: latitude (decimal degrees, EPSG4326)
      * Slope: topographic slope, as defined in Input data
      * Aspect: topographic aspect, as defined in Input data
      * CAC: soil CaCO3 concentration, as defined in Input data
      * CEC: cation exchange capacity, as defined in Input data
      * CLA: soil clay content, as defined in Input data
      * SAN: soil sand content, as defined in Input data
      * SIL: soil silt content, as defined in Input data
      * WAT: soil available water content, as defined in Input data
      * mean.SWI: SAGA Wetness Index, as defined in Input data
      * Hydric: presence of hydric soils, as defined in Input data
      * Floodplain: presence of floodplain, as defined in Input data
      * totalPPT: mean total annual precipitation, as defined in Input data
      * MeanTEMP: mean annual temperature, as defined in Input data
      * direction: cardinal direction of aspect, as defined in Input data
      * marea: nickname of the management area derived from the file names of the inputs
    * ydata: 78225 observations of 32 variables
      * each column contains presence or absence of a given taxon at a given corner. No.tree as defined above is the only non-taxon column
      
* **2.Process_OOS.R**: This script is identical to 1.Process.R but formats the out-of-sample data that will be used to validate the model in step 8.

  * Inputs: xdata (drivers) and ydata (response variables) from CSVs in the GJAMDATA/Withheld for Validation/ subdirectory specific to each management area as described in the methods of Shuman et al. (in prep). Structure is explained in the section Input data of the README.
  * Outputs: GJAMDATA/Withheld for Validation/validation_processed_xydata.RData: a RData file with the xdata_oos and ydata_oos (oos denoting out-of-sample data) dataframes with all management areas combined
    * xdata_oos: 24699 out-of-sample observations of 17 variables. Variables as described in 1.Process.R
    * ydata_oos: 24699 out-of-sample observations of 28 variables. Variables as described in 1.Process.R. The fewer number of variables is because some rare taxa that are present in the in-sample data and later combined into "other hardwood" and "other conifer" categories are absent
    
* **3.Combine_marea.R**: This step takes the output of 2.Process_OOS.R and matches each management area from the out-of-sample data with the nearest in-sample management area using Euclidean distance. This is necessary because the out-of-sample prediction must have only random effect groups that were present in the in-sample dataset.
  
  * Inputs:
    * GJAMDATA/Withheld for Validation/validation_processed_xydata.RData
    * GJAMDATA/processed_xydata.RData
  * Outputs:
    * GJAMDATA/Withheld for Validation/validation_processed_xydata_fixmarea.RData
      * Identical to the xdata_oos and ydata_oos in GJAM/Withheld for Validation/validation_processed_xydata.RData except that the contents of the marea (management area) column were replaced with management areas in the in-sample data
      * type: new column denoting that this is out-of-sample data
      
* **4.Reduce.R**: This script takes the output of 1.Process.R and reduces the number of individual taxa by grouping black gum, sweet gum, and black gum/sweet gum; poplar, tulip poplar, and poplar/tulip poplar; and any uncommon taxon with either other hardwood or other conifer. This reduces the dimensionality of the response variables to reduce the probability of overfitting and improve parameter inference.
  * Inputs:
    * GJAMDATA/processed_xydata.RData
  * Outputs:
    * GJAMDATA/processed_xydata_2.RData
      * xdata: 78224 observations of 17 variables. Identical to the xdata output saved in GJAMDATA/processed_xydata.RData. Saved again here for convenience down the pipeline
      * ydata: 78224 observations of 15 variables. Presence or absence (1/0) of each taxonomic group
        * No.tree: corners where no trees were observed, as defined in Input data
        * Poplar.tulip.poplar: presence/absence of taxonomic groups poplar, tulip poplar, and poplar/tulip poplar in the original data
        * Black.gum.sweet.gum: presence/absence of taxonomic groups black gum, sweet gum, and black gum/sweet gum in the original data
        * Other.conifer: presence/absence of taxonomic groups bald cypress, pine, tamarack, cedar/juniper in the original data
        * Other.hardwood: presence/absence of taxonomic groups birch, locust, willow, cherry, sycamore, buckeye, hackberry, mulberry, other.hardwood, alder, and chestnut in the original data
        * All other columns as described in ydata in Input data

* **4.Reduce_ecosystem.R**: This script takes the output of 4.Reduce.R and further reduces the dimensions of the response variable to include only three ecosystem states: prairie, savanna, and forest. This allows an analysis of both biome-level and community-level drivers of vegetation distributions.
  * Inputs: GJAMDATA/processed_xydata_2.RData
  * Outputs: GJAMDATA/processed_xydata_2_ecosystem.RData
    * xdata: 78224 observations of 17 variables. Identical to the xdata output saved in GJAMDATA/processed_xydata.RData. Saved again here for convenience down the pipeline
    * ydata: 78224 observations of 3 variables. Presence or absence (1/0) of each biome
      * Prairie: present at any corner only including presence of No.tree in the original data
      * Savanna: present at any corner including presence of oak and/or hickory and no other tree taxa in the original data
      * Forest: present at any corner including any tree taxon other than exclusively oak and/or hickory
      
* **5.Reduce_OOS.R**: This script is identical to 4.Reduce.R for the out-of-sample data. The script takes the output of 3.Combine_marea.R as the input
  * Inputs: GJAMDATA/Withheld for Validation/validation_processed_xydata_fixmarea.RData
  * Outputs: GJAMDATA/Withheld for Validation/validation_processed_xydata_fixmarea_reduced.RData
    * xdata_oos: 24699 observations of 18 variables. Identical to the xdata saved in GJAMDATA/Withheld for Validation/validation_processed_xydata_fixmarea.RData. Saved again for convenience down the pipeline
    * ydata_oos: 24699 observations of 15 variables. Presence/absence (1/0) of each taxonomic group
      * No.tree: corners where no trees were observed, as defined in Input data
      * Poplar.tulip.poplar: presence/absence of poplar, tulip poplar, and poplar/tulip poplar as defined above for in-sample data
      * Black.gum.sweet.gum: presence/absence of black gum, sweet gum, black gum/sweet gum as defined above for in-sample data
      * Other.conifer: presence/absence of taxonomic groups bald cypress, pine, tamarack, cedar/juniper as defined above for in-sample data
      * Other.hardwood: presence/absence of taxonomic groups birch, locust, willow, cherry, sycamore, buckeye, hackberry, mulberry, other.hardwood, alder, and chestnut as defined above for in-sample data

* **5.Reduce_OOS_ecosystem.R**: This script is identical to 4.Reduce_ecosystem.R for the out-of-sample data. It takes the output of 5.Reduce_OOS.R as the input
  * Inputs: GJAMDATA/Withheld for Validation/validation_processed_xydata_fixmarea_reduced.RData
  * Outputs: GJAMDATA/Withheld for Validation/validation_processed_xydata_fixmarea_reduced_ecosystem.RData
    * xdata_oos: 24699 observations of 18 variables. Identical to the xdata saved in GJAMDATA/Withheld for Validation/validation_processed_xydata_fixmarea.RData. Saved again for conveniences down the pipeline
    * ydata_oos: 24699 observations of 3 variables. Presence/absence (1/0) of each biome
      * Prairie: present at any corner only including presence of No.tree in the original dataset
      *  Savanna: present at any corner including presence of oak and/or hickory but not any other tree taxon in the original dataset
      * Forest: present at any corner including any tree taxon other than exclusively oak and/or hickory

* **6.Run**: This step is a folder with the run scripts for each of the four simulation types. There are four subdirectories of 6.Run, each of which contains both the R script and a job submission script for running the simulations on the University of Notre Dame cluster computing system. The R scripts are set up in the same way and differ in two ways: (1) the data files that are called as inputs and (2) the covariates that are named as independent variables in the gjam formula. Each job should be submitted three to four times for four independent MCMC chains. The jobs differ in two ways: using the community- or biome-level response variables (All_taxa vs Reduced_taxa) and including or not including aspect as an environmental covariate (ASPECT vs NOASPECT)
  * Sub directories: 
    * All_taxa~all_cov_ASPECT
      * Input: GJAMDATA/processed_xydata_2.RData
      * Output:
        * out/All_taxa~all_cov_ASPECT/all_taxa-all_cov_ASPECT_1.RData
        * out/All_taxa~all_cov_ASPECT/all_taxa-all_cov_ASPECT_2.RData
        * out/All_taxa~all_cov_ASPECT/all_taxa-all_cov_ASPECT_3.RData
        * out/All_taxa~all_cov_ASPECT/all_taxa-all_cov_ASPECT_4.RData
    * All_taxa~all_cov_NOASPECT
      * Input: GJAMDATA/processed_xydata_2.RData
      * Output:
        * out/All_taxa~all_cov_NOASPECT/all_taxa-all_cov_NOASPECT_1.RData
        * out/All_taxa~all_cov_NOASPECT/all_taxa-all_cov_NOAPSECT_2.RData
        * out/All_taxa~all_cov_NOASPECT/all_taxa-all_cov_NOAPSECT_3.RData
    * Reduced_taxa~all_cov_ASPECT
      * Input: GJAMDATA/processed_xydata_2_ecosystem.RData
      * Output:
        * out/Reduced_taxa~all_cov_ASPECT/reduced_taxa-all_cov_ASPECT_1.RData
        * out/Reduced_taxa~all_cov_ASPECT/reduced_taxa-all_cov_ASPECT_2.RData
        * out/Reduced_taxa~all_cov_ASPECT/reduced_taxa-all_cov_ASPECT_3.RData
        * out/Reduced_taxa~all_cov_ASPECT/reduced_taxa-all_cov_ASPECT_4.RData
    * Reduced_taxa~all_cov_NOASPECT
      * Input: GJAMDATA/processed_xydata_2_ecosystem.RData
      * Output:
        * out/Reduced_taxa~all_cov_NOASPECT/reduced_taxa-all_cov_NOASPECT_1.RData
        * out/Reduced_taxa~all_cov_NOASPECT/reduced_taxa-all_cov_NOASPECT_2.Rdata
        * out/Reduced_taxa~all_cov_NOASPECT/reduced_taxa-all_cov_NOASPECT_3.RData
        * out/Reduced_taxa~all_cov_NOASPECT/reduced_taxa-all_cov_NOASPECT_4.RData
  * The formula includes all covariates: mean precipitation, mean temperature, topographic slope, SAGA Wetness Index, presence of hydric soils, presence of a floodplain, soil CaCO3 concentration, cation exchange capacity, soil sand content, soil clay content, soil water content, and topographic direction (for ASPECT runs). The response variable is the joint presence or absence of each taxon (All_taxa) or biome (Reduced_taxa)
  * All the submit.sh files have identical specifications as follows:
    * -M denotes which email information about the status of the job should go to
    * -m abe specifies to send emails about the job being aborted, beginning, and ending
    * -pe smp 10 denotes submitting the job to 10 parallel cores on a single node
    * -q long specifies the job should be submitted to the long queue
    * -N specifies the name of the job
    
* **7.Combine.R**: GJAM only allows for one chain to be run at a time. Given that the model is a Bayesian model relying on MCMC, multiple chains are useful for estimating the full uncertainty of the model and to reduce the impact of initial chain values on inference. We therefore ran the model 3-4 times using identical specifications in each of the simulations given in 6.Run and we then combine the chains into a more usable format here.
  * Inputs: individual files within one of the subfolders of out/ specific to the type of model run. The type being loaded can be changed by modifying the `type` parameter using the four available options specified at the top of the script. The input is always the entire global environment saved for each run of GJAM, but the specific input depends on the type of run you want to manipulate
  * Outputs: RData object named combined.RData in the specified subdirectory of the out/ folder. The subdirectory is specified in the `type` parameter at the beginning of the script. Each combined.RData object contains the following:
    * bFacGibbs: dataframe with 3200 samples. Number of columns differs according to the number of parameters in the specified model. Columns are beta coefficient estimates standardized for X with correlation scale for W. See Clark et al. (2017) for more information. Includes intercepts and coefficient estimates for "yes" and "no" of binary variables (hydric and floodplain). Also includes chain (1-4) and iter (201-1000) columns for the MCMC chain and MCMC iteration. iter starts at 201 because burn-in has already been removed. 
    * bgibbs: dataframe with 3200 samples. Number of columns differs according to the number of parameters in the specified model. Columns are beta coefficient estimates standardized for X. See Clark et al. (2017) for more information. Includes intercepts but only "yes" of binary variables (hydric and floodplain). Also includes chain (1-4) and iter (201-1000) columns for the MCMC chain and MCMC iteration. iter starts at 201 because burn-in has already been removed.
    * bgibbsUn: dataframe with 3200 samples. Number of columns differs according to the number of parameters in the specified model. Columns are unstandardized beta coefficient estimates. Includes intercepts but only "yes" of binary variables (hydric and floodplain). Also includes chain (1-4) and iter (201-1000) columns for the MCMC chain and MCMC iteration. iter starts at 201 because burn-in has already been removed. 
    * fSensGibbs: dataframe with 3200 samples. Number of columns differs according to the number of parameters in the specified model. Columns are the covariance between X and the responses they elicit from Y. See Clark et al. (2017) for more information. Also includes chain (1-4) and iter (201-1000) columns for the MCMC chain and MCMC iteration. iter starts at 201 because burn-in has already been removed.
    * sgibbs: dataframe with 3200 samples. Number of columns differs according to the number of parameters in the specified model. Columns are covariances between response variables. See Clark et al. (2017) for more information. Also includes chain (1-4) and iter (201-1000) columns for the MCMC chain and MCMC iteration. iter starts at 201 because burn-in has already been removed.
    
* **8.Visualize.R**: This script produces figures used in Shuman et al. (in prep) for both community-level and biome-level runs, depending on the input called into the script and the designation of "all" (community-level) or "reduced" (biome-level) at the beginning of the script. The script uses the output of 7.Combine.R, but the specific input depends on the simulation of interest.
  * Inputs: combined.RData file from the specified subdirectory of the out/ folder. The file is specified by modifying line 11.
  * Outputs: none. Figures produced

* **9.Predict_OOS.R**: This script uses the model fit of one of the simulations from 6.Run to predict the out-of-sample data. The type of prediction here predicts the response variable only from the environmental covariate sand does not take into account the covariance between taxa or biomes. Visualization is included in the same script.
  * Inputs:
    * The global environment RData file for the first chain of any of the four model runs. The first chain is specified as follows, using the All_taxa\~all_cov_ASPECT model run type as an example: out/All_taxa\~all_cov_ASPECT/all_taxa-all_cov_ASPECT_1.RD ta
    * The out-of-sample data that fits the model run type. One of the following:
      * GJAMDATA/Withheld for Validation/validation_processed_xydata_fixmarea_reduced.R ata: for All_taxa model runs 
      * GJAMDATA/Withheld for Validation/validation_processed_xydata_fixmarea_reduced_e ecosystem.RData for Reduced_taxa model runs
  * Outputs: none saved. This is the last step, so all analyses of the validation are done in the same step

* **9.Predict_OOS_conditional.R**: This script is similar to 9.Predict_OOS.R and uses the same inputs. The difference is that the type of prediction implemented in this script accounts for the covariance between taxa or ecosystem types, which we hypothesize will improve prediction.
  * Inputs:
    * The global environment RData file for the first chain of any of the four model runs. The first chain is specified as follows, using the All_taxa\~all_cov_ASPECT model run type as an example: out/All_taxa\~all_cov_ASPECT/all_taxa-all_cov_ASPECT_1.RData
    * The out-of-sample data that fits the model run type. One of the following
      * GJAMDATA/Withheld for Validation/validation_processed_xydata_fixmarea_reduced.R ata: for All_taxa model runs
      * GJAMDATA/Withheld for Validation/validation_processed_xydata_fixmarea_reduced_ecosystem.Rdata for Reduced_taxa model runs
  * Outputs: the validation is computationally intensive, so the out-of-sample prediction for the All_taxa~all_cov_NOASPECT and Reduced_taxa~all_cov_NOASPECT model run type are saved as intermediate outputs as follows:
    * out/cond_pred_all_taxa.RData: conditional prediction with All_taxa model run type
    * out/cond_pred_reduced_taxa.RData: conditional prediction with Reduced_taxa model run type
    
* **utils.R**: This script contains utility functions for the code. Specifically, there is a function for manually calculating the Gelman Rubin diagnostic for assessing chain convergence because our output is not in the proper format to use the default functions available in R. I additionally removed any identical chains (usually 2/4) from the output prior to calculating the diagnostic statistic. The identical chains are an artifact of the gjam function and cannot be avoided to the authors' knowledge. Removing the identical chains offers a more conservative view of chain convergence.
  
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
    
