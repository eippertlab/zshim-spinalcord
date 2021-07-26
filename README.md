# zshim-spinalcord
# under construction
The code in this repository is associated with the manuscript “Kaptan, M., Vannesjo, S.J., Mildner, T., Horn, U., Hartley-Davies, R., Oliva, V., Brooks, J.C.W., Weiskopf, N., Finsterbusch, J., Eippert, F. (2021). Automated slice-specific z-shimming for fMRI of the human spinal cord.” (biorxiv link of our paper) and the corresponding dataset (https://openneuro.org/datasets/ds003743).

Please note that this repository contains code for two different purposes: while the code contained in directory ZShim_OnlineCalculation can be used for automated calculation of slice-specific z-shims during an experiment, the code contained in directory ZShim_Results can be used to reproduce the results of the above-mentioned manuscript.

## Required software
MATLAB (https://www.mathworks.com/products/matlab.html), version 2018a or higher

Spinal Cord Toolbox (SCT; https://spinalcordtoolbox.com/en/latest/), version 3.2.7 or higher

FSL (https://fsl.fmrib.ox.ac.uk/fsl/fslwiki), version 5.0 or higher

dcm2niix (https://github.com/rordenlab/dcm2niix), version 1.0.20180622 or higher

## Automated selection of slice-specific z-shims (directory "ZShim_OnlineCalculation")
The code in directory ZShim_OnlineCalculation and the respective subdirectories (EPI_based for EPI-based selection and FM_based for field map based selection) was used to determine the slice-specific z-shims during data acquisition.

Please see section "2.4.2 Automated selection" in the aforementioned manuscript for a detailed explanation of the automated z-shim methods.

SCT version 3.2.7, FSL version 5.0, dcm2niix version 1.0.20180622, and MATLAB version 2019a were used for automated selection.

For the EPI-based calculation, there is a single Matlab script (*ZShim_OnlineCalculation_EPIbased.m*) that needs to be run. For the FM-based calculation there is also a main Matlab script (*ZShim_OnlineCalculation_FMbased.m*), but this requires the functions *compute_shims.m, compute_map_coords.m, compute_SH_basis.m, and raised_cosine.m* to be run.

If you use one of these methods in your research, please cite the paper “Kaptan, M., Vannesjo, S.J., Mildner, T., Horn, U., Hartley-Davies, R., Oliva, V., Brooks, J.C.W., Weiskopf, N., Finsterbusch, J., Eippert, F. (2021). Automated slice-specific z-shimming for fMRI of the human spinal cord” and add a link to this repository https://github.com/eippertlab/zshim-spinalcord

If you have any questions, please contact Merve Kaptan (mkaptan@cbs.mpg.de) or Falk Eippert (eippert@cbs.mpg.de).

## ZShim_Results
For preprocessing and results, SCT version 4.2.2, FSL version 6.0.3, and MATLAB version 2021a were used.

### Preprocessing, Directory: "Step1_Preprocessing"

Please see the section "2.5 Preprocessing" in the aforementioned paper for a detailed explanation of the preprocessing pipeline. 

ZShim_Preprocessing_I.m: This script contains the preprocessing steps that were employed to preprocess the following dataset: https://openneuro.org/datasets/ds003743
Please note that files that were created or modified manually (e.g. “disc_labels.nii.gz”) and files that were created based on 250 EPI volumes under different sequence variants (e.g., “moco2_volumes_mean.nii.gz”) are shared under participant-specific directories of the "derivatives" directory. 

After preprocessing, the extracted signals (see section "2.5.4 EPI signal extraction" in the aforementioned manuscript) were saved in the "extracted_signal" subdirectory of the "derivatives" parent directory. For the extracted signal (e.g., tSNR) in native space, see the subdirectory "signal_nativespace" and for the extracted signal in template space, see the subdirectory "signal_templatespace". 

Please note that a copy of the PAM50 template that was used for normalization is also shared for reproducibility of the results ("derivatives", "template", "SCT", "PAM50").


### Results, Directory: "Step2_CalculateResults"

Please see section "2.6 Statistical analysis" in the aforementioned paper for a detailed explanation of the statistical tests and results. 

The results of the statistical tests reported in the manuscript can be calculated based on the matrices (.mat files) contained in the "extracted_signal" directory (under the "derivatives" parent directory) and relevant subdirectories without running the preprocessing pipeline. Alternatively, the script *ZShim_Preprocessing_I.m* can be used to preprocess the data and recreate the “extracted_signal” directory to run the following scripts.

*ZShimResults_III_I.m*: calculate all the results & supplementary results that were reported under the section "3.1 Replication and extension of previous findings" of the aforementioned manuscript.

*ZShimResults_III_II.m*: calculate all the results & supplementary results that were reported under the section "3.2 Automation of z-shimming" of the aforementioned manuscript.

Please note that we are not sharing code for section “3.3 Validation of EPI-based automation approach”, because this is based on an (externally acquired) independent data-set, which we cannot publicly share.

#### Helper_FunctionsScripts
This directory contains helper functions and scripts to produce results & figures, called via *ZShimResults_III_I.m* or *ZShimResults_III_II.m*

##### Figure_Code
This directory contains functions & scripts to create the figures & supplementary figures in the aforementioned manuscript. Please note that there are no scripts for Figure 1 (as it is mostly based on images created in fsleyes), and the code to produce Figure 4 cannot be shared as the underlying data is an externally acquired (independent) data-set.

*ZShimFigures_FigureII.m*: create the plots in Figure 2 (Figure 2. Replication and extension of previous results) of the aforementioned manuscript, called at the end the script *ZShimResults_III_I.m*

*ZShimFigures_FiguresIII.m*: create most of the plots in Figure 3 (Figure 3. Performance of both automated methods) of the aforementioned manuscript, called at the end of the script *ZShimResults_III_II.m*

*ZShimFigures_SupplementaryFigures.m*: run this script to create the supplementary figures (Supplementary Figures 1,2 and 5). Note that Supplementary Figure 4 is created by the function named *ZShim_SuppResults_III_II_II_IV.m* (called from the script *ZShimResults_III_II.m*). Please note that there are no scripts for Supplementary Figure 3 (as it does not require any code).

###### distributionPlot
Reference for the “distributionPlot” :
Jonas (2021). Violin Plots for plotting multiple distributions (distributionPlot.m) (https://www.mathworks.com/matlabcentral/fileexchange/23661-violin-plots-for-plotting-multiple-distributions-distributionplot-m), MATLAB Central File Exchange. 


