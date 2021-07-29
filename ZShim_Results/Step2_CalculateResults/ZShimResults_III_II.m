% In this script, the calculation of the results for the section "3.2
% Automation of z-shimming" are listed for the
% following manuscript:Kaptan, M., Vannesjo, S. J., Mildner, T., Horn, U., Hartley-Davies, R., 
% Oliva, V., Brooks, J. C. W., Weiskopf, N., Finsterbusch, J., & Eippert, F. (2021). 
% Automated slice-specific z-shimming for fMRI of the human spinal cord. BioRxiv, 
% 2021.07.27.454049. https://doi.org/10.1101/2021.07.27.454049

% The organization of the results were kept consistent with the
% manuscript.

% Merve Kaptan, mkaptan@cbs.mpg.de
% 04.07.2021


% 3.2 Automation of z-shimming
clc;clear all;close all
% path to the provided folder with extracted signal data
datapath = '/data/pt_02098/ZShim_BIDS_0807/derivatives/extracted_signal/';
% path to analysis code
codepath = '/data/pt_02098/ZShim_BIDS_0807/derivatives/code/';
% add code folder & subfolders to path 
addpath(genpath(codepath))
% 1 = save the results as a matrix, 0 = do not save
saveResults = 0;
% 1 = run all the & reconstrcution of zref EPI image 
% supplementary functions again (can take a long time!)
% the readily extracted signal can be found under 
% "extracted_signal/signal_templatespace/ReconstructedSignal"
recalculateResults = 0;
% SCT template path
scttemplatepath = '/data/pt_02098/ZShim_BIDS_0807/derivatives/template/SCT/';
% path to raw data (to get the list of subjects)
rawdatapath     = '/data/pt_02098/ZShim_BIDS_0807/rawdata/';
% path to processed data
processdatapath = '/data/pt_02098/ZShim_BIDS_0807/derivatives/';
% 1 = save the figures as .svg, 0 = do not save
printStatus = 0;
% if printStatus is true, save the results in this folder
% leave empty if you do not want to save figures
figurepath  = '/data/pt_02098/ZShim_BIDS_0807/derivatives/figures/';

%where FSL is downloaded
fslDir  = '/afs/cbs.mpg.de/software/fsl/5.0.11/ubuntu-xenial-amd64/';
%set FSL directory
setenv('FSLDIR', fslDir);
fslDir = getenv('FSLDIR');
fsldirmpath = sprintf('%s/etc/matlab',fslDir);
path(path, fsldirmpath);
setenv('FSLOUTPUTTYPE', 'NIFTI_GZ'); % output type nii gz
%% 3.2.1 EPI-based automation
load(fullfile(datapath, 'signal_templatespace', 'GroupSingle_EPI', 'GM', 'TimeSeries_TSNR', 'Data.mat'))

no = (noLEFT_DORSAL+noLEFT_VENTRAL+noRIGHT_DORSAL+noRIGHT_VENTRAL)./4;
manual = (manualLEFT_DORSAL+manualLEFT_VENTRAL+manualRIGHT_DORSAL+manualRIGHT_VENTRAL)./4;
auto = (autoLEFT_DORSAL+autoLEFT_VENTRAL+autoRIGHT_DORSAL+autoRIGHT_VENTRAL)./4;

notests = 3;
Results_automation{1,1} = 'EPI-based';
Results_automation{1,2} = ZShim_CalculateResults(no,manual,auto,0,notests);

% I. supplementary analysis for mean gray matter signal intensity
clearvars -except Results_automation datapath printStatus figurepath codepath scttemplatepath rawdatapath recalculateResults saveResults processdatapath fslDir
load(fullfile(datapath, 'signal_templatespace', 'GroupSingle_EPI', 'GM', 'SingleVolume_MeanSignal', 'Data.mat'))

no = (noLEFT_DORSAL+noLEFT_VENTRAL+noRIGHT_DORSAL+noRIGHT_VENTRAL)./4;
manual = (manualLEFT_DORSAL+manualLEFT_VENTRAL+manualRIGHT_DORSAL+manualRIGHT_VENTRAL)./4;
auto = (autoLEFT_DORSAL+autoLEFT_VENTRAL+autoRIGHT_DORSAL+autoRIGHT_VENTRAL)./4;

notests = 3;
Results_automation{2,1} = 'EPI-based, Supp I (gray matter signal intensity)';
Results_automation{2,2} = ZShim_CalculateResults(no,manual,auto,0,notests);

%% 3.2.2 FM-based automation
clearvars -except Results_automation datapath printStatus figurepath codepath scttemplatepath rawdatapath recalculateResults saveResults processdatapath fslDir
load(fullfile(datapath, 'signal_templatespace', 'GroupSingle_FM', 'GM', 'TimeSeries_TSNR', 'Data.mat'))

no = (noLEFT_DORSAL+noLEFT_VENTRAL+noRIGHT_DORSAL+noRIGHT_VENTRAL)./4;
manual = (manualLEFT_DORSAL+manualLEFT_VENTRAL+manualRIGHT_DORSAL+manualRIGHT_VENTRAL)./4;
auto = (autoLEFT_DORSAL+autoLEFT_VENTRAL+autoRIGHT_DORSAL+autoRIGHT_VENTRAL)./4;

notests = 3;
Results_automation{3,1} = 'FM-based';
Results_automation{3,2} = ZShim_CalculateResults(no,manual,auto,0,notests);

% I. supplementary analysis for mean gray matter signal intensity
clearvars -except Results_automation datapath printStatus figurepath codepath scttemplatepath rawdatapath recalculateResults saveResults processdatapath
load(fullfile(datapath, 'signal_templatespace', 'GroupSingle_EPI', 'GM', 'SingleVolume_MeanSignal', 'Data.mat'))

no = (noLEFT_DORSAL+noLEFT_VENTRAL+noRIGHT_DORSAL+noRIGHT_VENTRAL)./4;
manual = (manualLEFT_DORSAL+manualLEFT_VENTRAL+manualRIGHT_DORSAL+manualRIGHT_VENTRAL)./4;
auto = (autoLEFT_DORSAL+autoLEFT_VENTRAL+autoRIGHT_DORSAL+autoRIGHT_VENTRAL)./4;

notests = 3;
Results_automation{4,1} = 'FM-based, Supp I (gray matter signal intensity)';
Results_automation{4,2} = ZShim_CalculateResults(no,manual,auto,0,notests);

% ------------------------------------------------------------------------
% Investigation of the suboptimal performance of the FM-based selection
% ------------------------------------------------------------------------
% (explained in detail in Supplementary Material)
% Before reconstructing volumes with modified z-shims based on different
% FM-based procedures (3.2.2.1 - 3.2.2.6) we need a ground truth to compare with.
% Therefore, first step: artificial reconstruction of the volumes 
% based on manual, automated, no z-picks. Then the different modifications
% of the FM-based approach lead to modified z-shims which are used to
% create artifical volumes and we compare their properties to our previous
% implementation.
% If recalculateResults is set to False, no new data is created, but
% instead the saved results from the comparisons against original implementation
% are loaded from extracted_signal/signal_templatespace/GroupWhole/ReconstructedSignal
subGroup = 'epi';
ZShim_Create_ArtificialVolumes(subGroup, scttemplatepath, ...
    rawdatapath,processdatapath,recalculateResults,fslDir);
subGroup = 'fm';
ZShim_Create_ArtificialVolumes(subGroup, scttemplatepath, ...
    rawdatapath,processdatapath,recalculateResults,fslDir);                             

% Run the following steps in order so that necessary files are copied into
% each subject's fmap folder into a subfolder called investigation

% 3.2.2.1 Choice of mask for identifying the spinal cord in the field map phase data
results = ZShim_SuppIII_II_II_I(scttemplatepath,rawdatapath,processdatapath,recalculateResults,fslDir);

Results_automation{5,1} = 'FM-based, Supp- choice of mask';
Results_automation{5,2} = results;

% 3.2.2.2 Choice of parameters employed in the fitting process of the gradient field
clear results
results = ZShim_SuppIII_II_II_II(scttemplatepath,rawdatapath,processdatapath,recalculateResults,fslDir);

Results_automation{6,1} = 'FM-based, Supp- fit using different parameters';
Results_automation{6,2} = results;

% 3.2.2.3 Field gradients in the AP-direction
results = ZShim_SuppIII_II_II_III(scttemplatepath,rawdatapath,processdatapath,recalculateResults,fslDir);
Results_automation{7,1} = 'FM-based, Supp- gradients in AP';
Results_automation{7,2} = results;

% 3.2.2.4 Inhomogeneity-induced mis-localizations between EPIs and field map 
ZShim_SuppIII_II_II_IV(figurepath,rawdatapath,processdatapath,recalculateResults,fslDir)

% 3.2.2.5 Use of different field map
clear results
results = ZShim_SuppIII_II_II_V(scttemplatepath,rawdatapath,processdatapath,recalculateResults,fslDir);
Results_automation{8,1} = 'FM-based, Supp- different field map';
Results_automation{8,2} = results;

% 3.2.2.6 Assessing the reliability of z-shim selection based on FM-based automation
clear results
results = ZShim_SuppIII_II_II_VI(rawdatapath,processdatapath,recalculateResults,fslDir);
Results_automation{9,1} = 'FM-based, Supp- reliability of fm-based picks';
Results_automation{9,2} = results;

%% 3.2.3 Comparing all three approaches
% Compare the three approaches using independent samples t-tests 
clearvars -except Results_automation datapath printStatus figurepath codepath scttemplatepath rawdatapath recalculateResults saveResults processdatapath fslDir

results_3approaches = ZShim_CompareThreeApproaches(rawdatapath,processdatapath, datapath,recalculateResults);

Results_automation{10,1} = 'Comparing all three (manual, EPI-based, FM-based approaches)';
Results_automation{10,2} = results_3approaches;
%%
if saveResults 
    save([datapath filesep 'Results_Automationofzshimming.mat'], 'Results_automation')
end

%% do the plots in Figure 3

ZShim_Figures_FigureIII(printStatus,datapath,figurepath)

