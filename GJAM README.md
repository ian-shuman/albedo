# Overview

This repository contains code written by AM Willson and I Shuman for Shuman et al. (in prep). The code uses GJAM to investigate the vegetation-environment relationship during the period of European settlement in the U.S. Midwest, as documented via the Public Land Surveys of the 19th century.

# License

This repository holds an MIT License, as described in the LICENSE file.

# Software versions

This repository is entirely built in the R environment using R version 4.3.1.

# Package versions

* `corrplot` v. 0.92
* `cowplot` v. 1.1.1
* `dplyr` v. 1.1.3
* `fields` v. 14.1
* `ggplot2` v. 3.4.4
* `gjam` v. 2.6.2
* `lme4` v. 1.1.34
* `maps` v. 3.4.1
* `piecewiseSEM` v. 2.3.0
* `sf` v. 1.0.14
* `sfheaders` v. 0.4.3
* `stats` v. 4.3.1
* `stringr` v. 1.5.0
* `tibble` v. 3.2.1
* `tidyr` v. 1.3.0
* `utils` v 4.3.1
* `virids` v. 0.6.3

# Directory structure

* GJAMDATA: **Data** and **intermediate outputs** from processing steps required to run the analysis code are available here. Data are separated from intermediate outputs because data are always contained within subdirectories named "X" and "Y."
  * In-sample data: 
    * GJAMDATA/X
    * GJAMDATA/Y
  * Out-of-sample data: 
    * GJAMDATA/Withheld for Validation/X 
    * GJAMDATA/Withheld for Validation/Y
  * Intermediate outputs are available here and described after the steps producing them below:
    * GJAMDATA/processed_xydata.RData
    * GJAMDATA/processed_xydata_2.RData
    * GJAMDATA/processed_xydata_2_ecosystem.RData
    * GJAMDATA/Withheld for Validation/validation_processed_xydata.RData
    * GJAMDATA/Withheld for Validation/validation_processed_xydata_fixmarea.RData
    * GJAMDATA/Withheld for Validation/validation_processed_xydata_fixmarea_reduced.RData
    * GJAMDATA/Withheld for Validation/validation_processed_xydata_fixmarea_reduced_ecosystem.RData
* out: **Outputs** of the GJAM model are available in this directory. Any data or outputs not included in the repository are available upon request at this time. Processing for the data inputs to this repository are available from the following repository:
  * <https://github.com/ian-shuman/GJAMCovariates.git>: produces point-level reconstructions of environmental covariates
* R: The bulk of the repository is housed in the R subdirectory. The scripts within the R subdirectory are ordered according to the sequence of steps that should be taken to recreate the analysis. Descriptions of each step are presented in the following subsection, and should be run in the order presented, except step 0 must follow step 1.
* Review: Documents related to an code review conducted prior to submitting the associated manuscript for publication. The documents are as follows:
  * Code Review Checklist - GJAM - AW-CK - Notre Dame Checklist.pdf: A filled out rubric pertaining to reviewing the code, workflow, and statistical methodology of the manuscript associated with this paper. The checklist was jointly written by current members of the McLachlan Lab at the University of Notre Dame as of December 2023 (Jason McLachlan, Alyssa Willson, Jody Peters, Cazimir Kowalski, Nate Kroeze, Hannah O'Grady). The "author checklist" and "author notes" sections were filled out by Alyssa Willson and the "reviewer checklist" and "reviewer notes" were filled out by Cazimir Kowalski.
  * ck_review.Rmd: additional notes related to the code review conducted by C. Kowalski

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
  
# Input data

## xdata

* 15 in-sample and 7 out-of-sample CSVs. Each CSV contains point-level (corner) data for one management area
  * GJAMDATA/X/IL_Mixed1_X_Ian10mDataStructure_RandomEffects_Fixed.csv
  * GJAMDATA/X/IL_Prairie1_X_Ian10mDataStructure_RandomEffects_Fixed.csv
  * GJAMDATA/X/IL_Prairie2_X_Ian10mDataStructure_RandomEffects_Fixed.csv
  * GJAMDATA/X/IL_Prairie3_X_Ian10mDataStructure_RandomEffects_Fixed.csv
  * GJAMDATA/X/IL_River2_X_Ian10mDataStructure_RandomEffects_Fixed.csv
  * GJAMDATA/X/IL_Shawnee_X_Ian10mDataStructure_RandomEffects_Fixed.csv
  * GJAMDATA/X/IL_Small1_X_Ian10mDataStructure_RandomEffects_Fixed.csv
  * GJAMDATA/X/IL_Small3_X_Ian10mDataStructure_RandomEffects_Fixed.csv
  * GJAMDATA/X/IL_StudyArea_X_Ian10mDataStructure_RandomEffects_Fixed.csv
  * GJAMDATA/X/IN_Dunes_X_Ian10mDataStructure_RandomEffects_Fixed.csv
  * GJAMDATA/X/IN_Forest1_X_Ian10mDataStructure_RandomEffects_Fixed.csv
  * GJAMDATA/X/IN_Forest2_X_Ian10mDataStructure_RandomEffects_Fixed.csv
  * GJAMDATA/X/IN_HoosierNorth_X_Ian10mDataStructure_RandomEffects_Fixed.csv
  * GJAMDATA/X/IN_HoosierSouth_X_Ian10mDataStructure_RandomEffects_Fixed.csv
  * GJAMDATA/X/IN_NE_X_Ian10mDataStructure_RandomEffects_Fixed.csv
  * GJAMDATA/Withheld for Validation/X/IL_Forest1_X_Ian10mDataStructure_RandomEffects_Fixed.csv
  * GJAMDATA/Withheld for Validation/X/IL_River1_X_Ian10mDataStructure_RandomEffects_Fixed.csv
  * GJAMDATA/Withheld for Validation/X/IL_Small2_X_Ian10mDataStructure_RandomEffects_Fixed.csv
  * GJAMDATA/Withheld for Validation/X/IN_Forest3_X_Ian10mDataStructure_RandomEffects_Fixed.csv
  * GJAMDATA/Withheld for Validation/X/IN_Forest4_X_Ian10mDataStructure_RandomEffects_Fixed.csv
  * GJAMDATA/Withheld for Validation/X/IN_Indianapolis_X_Ian10mDataStructure_RandomEffects_Fixed.csv
  * GJAMDATA/Withheld for Validation/X/IN_Prairie1_X_Ian10mDataStructure_RandomEffects_Fixed.csv

* Dimensions:
  * Rows correspond to individual corners in the PLS dataset. Each CSV has a different number of rows corresponding to the size of the management area
  * Columns are as follows:
    * uniqueID: an identifier for each corner within a given dataframe. When binded together, the uniqueID must be combined with the management area to get a unique corner identifier across management areas
    * x: easting (m)
    * y: northing (m)
    * long: longitude (decimal degrees, EPSG:4326)
    * lat: latitude (decimal degrees, EPSG:4326)
    * Slope: topographic slope (m)
    * Aspect: topographic aspect (degrees). NA indicates slope of 0, when aspect cannot be calculated
    * CAC: soil CaCO3 concentration (%)
    * CEC: cation exchange capacity (meq/100g)
    * CLA: soil clay content (%)
    * SAN: soil sand content (%)
    * SIL: soil silt content (%)
    * WAT: soil available water content (cm)
    * mean.SWI: SAGA Wetness Index (unit area/radian)
    * Hydric: presence or absence of hydric soils (1/0 binary)
    * Floodplain: presence or absence of floodplain (1/0 binary)
    * totalPPT: mean total annual precipitation (mm)
    * MeanTEMP: mean annual temperature (deg F)
    * direction: cardinal direction of aspect. Simplification of aspect to improve use of aspect with the high number of points with aspect NA because slope = 0 (N = north, S = south, E = east, W = west, NS = no slope)
    
## ydata

* 15 in-sample and 7 out-of-sample CSVs corresponding to the xdata files. Each file contains point-level vegetation data
  * GJAMDATA/Y/IL_Mixed1_Y_DataStructure.csv
  * GJAMDATA/Y/IL_Prairie2_Y_DataStructure.csv
  * GJAMDATA/Y/IL_River2_Y_DataStructure.csv
  * GJAMDATA/Y/IL_Small1_Y_DataStructure.csv
  * GJAMDATA/Y/IL_StudyArea_Y_DataStructure.csv
  * GJAMDATA/Y/IN_Forest1_Y_DataStructure.csv
  * GJAMDATA/Y/IN_HoosierNorth_Y_DataStructure.csv
  * GJAMDATA/Y/IN_NE_Y_DataStructure.csv
  * GJAMDATA/Y/IL_Prairie1_Y_DataStructure.csv
  * GJAMDATA/Y/IL_Prairie3_Y_DataStructure.csv
  * GJAMDATA/Y/IL_Shawnee_Y_DataStructure.csv
  * GJAMDATA/Y/IL_Small3_Y_DataStructure.csv
  * GJAMDATA/Y/IN_Dunes_Y_DataStructure.csv
  * GJAMDATA/Y/IN_Forest2_Y_DataStructure.csv
  * GJAMDATA/Y/IN_HoosierSouth_Y_DataStructure.csv
  * GJAMDATA/Withheld for Validation/Y/IL_Forest1_Y_DataStructure.csv
  * GJAMDATA/Withheld for Validation/Y/IL_Small2_Y_DataStructure.csv
  * GJAMDATA/Withheld for Validation/Y/IN_Forest4_Y_DataStructure.csv
  * GJAMDATA/Withheld for Validation/Y/IN_Prairie1_Y_DataStructure.csv
  * GJAMDATA/Withheld for Validation/Y/IL_River1_Y_DataStructure.csv
  * GJAMDATA/Withheld for Validation/Y/IN_Forest3_Y_DataStructure.csv
  * GJAMDATA/Withheld for Validation/Y/IN_Indianapolis_Y_DataStructure.csv
  
* Dimensions:
  * Rows correspond to individual corners in the PLS dataset. Each CSV has a different number of rows corresponding to the size of the management area
  * Columns are as follows:
    * uniqueID: an identifier for each corner within a given dataframe. When binded together, the uniqueID must be combined with the management area to get a unique corner identifier across management areas
    * chainstree, chainstree2, chainstree3, chainstree4: distance of each tree at a given corner to the corner. There can be between 1 and 4 trees at a given corner, so chainstree2-chainstree4 can be NA. Units are chains. Currently not used in the analysis
    * No.tree, Oak, Elm, Hickory, Ash, Unknown.tree, Poplar, Maple, Sycamore, Other.hardwood, Mulberry, Basswood, Walnut, Cherry, Locust, No.data, Hackberry, Willow, Buckeye, Birch, Black.gum.sweet.gum, Water, Sweet.gum, Black.gum, Ironwood, Poplar.tulip.poplar, Beech, Dogwood, Wet, NA.: Columns indicating the presence or absence of each taxon/group at each corner (1/0 presence/absence). Columns that are not specific to a taxonomic group are as follows:
      * No.tree: indicates no tree was present within a given distance of the corner
      * Unknown.tree: the surveyor's notes were impossible to read to determine the tree taxon present at the site. Removed during data processing
      * No.data: No data available at a given corner. Removed during data processing
      * Water: The corner is within a body of water. Removed during data processing
      * Wet: The corner is in a wetland area. Removed during data processing
      * NA: Miscellaneous column removed during data processing
    * No.tree_dist, Oak_dist, Elm_dist, Hickory_dist, Ash_dist, Unknown.tree_dist, Poplar_dist, Maple_dist, Sycamore_dist, Other.hardwood_dist, Mulberry_dist, Basswood_dist, Walnut_dist, Cherry_dist, Locust_dist: distances (in chains) of the recorded tree of each taxon to the corner. Currently not used in the analysis
    
