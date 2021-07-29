% In this script, the calculation of the results for the section "3.1 
% Replication and extension of previous findings" are listed for the 
% following manuscript: Kaptan, M., Vannesjo, S. J., Mildner, T., Horn, U., Hartley-Davies, R., 
% Oliva, V., Brooks, J. C. W., Weiskopf, N., Finsterbusch, J., & Eippert, F. (2021). 
% Automated slice-specific z-shimming for fMRI of the human spinal cord. BioRxiv, 
% 2021.07.27.454049. https://doi.org/10.1101/2021.07.27.454049

% The organization of the results were kept consistent with the
% manuscript.

% Merve Kaptan, mkaptan@cbs.mpg.de
% 04.07.2021

% 3.1 Replication and extension of previous findings
clc;clear all;close all
% path to the provided folder with extracted signal data
datapath = '/data/pt_02098/ZShim_BIDS_0807/derivatives/extracted_signal';
% path to analysis code
codepath = '/data/pt_02098/ZShim_BIDS_0807/derivatives/code';
% add code folder & subfolders to path 
addpath(genpath(codepath))
% 1 = save the results as a matrix, 0 = do not save
saveResults = 0;
% 1 = save the figures as .svg, 0 = do not save
printStatus = 0;
% if printStatus is true, save the results in this folder
figurepath  = '/data/pt_02098/ZShim_BIDS_0807/derivatives/figures';
%%
% 3.1.1 Direct replication
load(fullfile(datapath, 'signal_templatespace', 'GroupWhole', 'SingleVolume_MeanSignal', 'Data.mat'))
Results_repandext{1,1} = 'Direct replication';
Results_repandext{1,2} = ZShim_CalculateResults(no,manual,[],1,1);

% supplementary analysis
% I. supplementary analysis for the mean of moco images
clearvars -except Results_repandext datapath printStatus figurepath codepath saveResults
load(fullfile(datapath, 'signal_templatespace', 'GroupWhole', 'SingleVolume_MOCOMeanSignal', 'Data.mat'))
Results_repandext{2,1} = 'Direct replication-supp 1 (mean moco)';
Results_repandext{2,2} = ZShim_CalculateResults(no,manual,[],1,1);

% II. supplementary analysis in the native space
clearvars -except Results_repandext datapath printStatus figurepath codepath saveResults
load(fullfile(datapath, 'signal_nativespace', 'GroupWhole', 'SingleVolume_MeanSignal', 'Data.mat'))
Results_repandext{3,1} = 'Direct replication-supp 2 (native)';
Results_repandext{3,2} = ZShim_CalculateResults(no,manual,[],1,1);

clearvars -except Results_repandext datapath printStatus figurepath codepath saveResults
load(fullfile(datapath, 'signal_nativespace', 'GroupWhole', 'SingleVolume_MOCOMeanSignal', 'Data.mat'))
Results_repandext{4,1} = 'Direct replication-supp 3 (mean moco-native)';
Results_repandext{4,2} = ZShim_CalculateResults(no,manual,[],1,1);

%% 3.1.2 Slice-by-slice characterization of z-shim effects 
clearvars -except Results_repandext datapath printStatus figurepath codepath saveResults
load(fullfile(datapath, 'signal_nativespace', 'GroupWhole', 'SingleVolume_MeanSignal', 'Data_picks.mat'))

notests = 5;

Results_repandext{5,1} = '0 difference';
% for this test, the t-test is not directed (as the differences in the baseline is not expected)
[~,p,~,stats] = ttest(Manual_0, No_0);

results{1,1} = 'mean signal';
results{1,2} = 'p';
results{1,3} = 't-value';
results{1,4} = 'dof';

if notests > 1
    p = p * notests;
end

results{2,2} = p;
results{2,3} = stats.tstat;
results{2,4} = stats.df;

zshimGlobalMeans = mean([No_0 Manual_0]);
meanG = zshimGlobalMeans(1);
meanW = zshimGlobalMeans(2:end);
meanP = 100 * (meanW/meanG);
meanP_Increase = meanP - 100;

results{1,5} = 'percent increase';
results{2,5} = meanP_Increase;

Results_repandext{5,2} = results;
clear results

Results_repandext{6,1} = '1 difference';
Results_repandext{6,2} = ZShim_CalculateResults(No_1,Manual_1,[],0,notests);

Results_repandext{7,1} = '2 difference';
Results_repandext{7,2} = ZShim_CalculateResults(No_2,Manual_2,[],0,notests);

Results_repandext{8,1} = '3 difference';
Results_repandext{8,2} = ZShim_CalculateResults(No_3,Manual_3,[],0,notests);

Results_repandext{9,1} = 'above 3 difference';
Results_repandext{9,2} = ZShim_CalculateResults(No_above3,Manual_above3,[],0,notests);

% supplementary analysis
% I. supplementary analysis: Repated-measures Anova
groups = [repmat({'no'},48,1); repmat({'manual'},48,1)];
t = table(groups, [No_0;Manual_0], [No_1;Manual_1], [No_2;Manual_2], ...
    [No_3;Manual_3], [No_above3;Manual_above3], ...
    'VariableNames', {'group', 'p0', 'p1', 'p2', 'p3', 'pa3'});

Meas = table({'0'; '1'; '2'; '3'; 'm3'},'VariableNames',{'Picks'});
rm = fitrm(t,'p0-pa3~group','WithinDesign',Meas);
Results_repandext{10,1} = 'Supp I - Repeated measures Anova';
Results_repandext{10,2} = ranova(rm, 'WithinModel', 'Picks');

% II. supplementary analysis: repetition of the same analysis for the mean of
% moco images
clearvars -except Results_repandext datapath printStatus figurepath codepath saveResults
load(fullfile(datapath, 'signal_nativespace', 'GroupWhole', 'SingleVolume_MOCOMeanSignal', 'Data_picks.mat'))

notests = 5;

% for this test, the t-test is not directed (as the differences in the baseline is not expected)
[h,p,ci,stats] = ttest(Manual_0, No_0);

results{1,1} = 'mean signal';
results{1,2} = 'p';
results{1,3} = 't-value';
results{1,4} = 'dof';

if notests > 1
    p = p * notests;
end

results{2,2} = p;
results{2,3} = stats.tstat;
results{2,4} = stats.df;

zshimGlobalMeans = mean([No_0 Manual_0]);
meanG = zshimGlobalMeans(1);
meanW = zshimGlobalMeans(2:end);
meanP = 100 * (meanW/meanG);
meanP_Increase = meanP - 100;

results{1,5} = 'percent increase';
results{2,5} = meanP_Increase;

Results_repandext{11,1} = 'Supp II- 0 difference';
Results_repandext{11,2} = results;
clear results

Results_repandext{12,1} = 'Supp II- 1 difference';
Results_repandext{12,2} = ZShim_CalculateResults(No_1,Manual_1,[],0,notests);

Results_repandext{13,1} = 'Supp II- 2 difference';
Results_repandext{13,2} = ZShim_CalculateResults(No_2,Manual_2,[],0,notests);

Results_repandext{14,1} = 'Supp II- 3 difference';
Results_repandext{14,2} = ZShim_CalculateResults(No_3,Manual_3,[],0,notests);

Results_repandext{15,1} = 'Supp II- above 3 difference';
Results_repandext{15,2} = ZShim_CalculateResults(No_above3,Manual_above3,[],0,notests);

groups = [repmat({'no'},48,1); repmat({'manual'},48,1)];
t = table(groups, [No_0;Manual_0], [No_1;Manual_1], [No_2;Manual_2], ...
    [No_3;Manual_3], [No_above3;Manual_above3], ...
    'VariableNames', {'group', 'p0', 'p1', 'p2', 'p3', 'pa3'});

Meas = table({'0'; '1'; '2'; '3'; 'm3'},'VariableNames',{'Picks'});
rm = fitrm(t,'p0-pa3~group','WithinDesign',Meas);
Results_repandext{16,1} = 'Supp II - Repeated measures Anova';
Results_repandext{16,2} = ranova(rm, 'WithinModel', 'Picks');

%% 3.1.3 z-shim effects across different echo times
clearvars -except Results_repandext datapath printStatus figurepath codepath saveResults
load(fullfile(datapath, 'signal_templatespace', 'GroupWhole', 'TE_SingleVolume_MeanSignal', 'Data.mat'))

notests = 2;

Results_repandext{17,1} = 'EchoTime-TE30';
Results_repandext{17,2} = ZShim_CalculateResults(no30,manual30,[],0,notests);

Results_repandext{18,1} = 'EchoTime-TE50';
Results_repandext{18,2} = ZShim_CalculateResults(no50,manual50,[],0,notests);

% supplementary analysis
% I. supplementary analysis for the mean of moco images

clearvars -except Results_repandext datapath printStatus figurepath codepath saveResults
load(fullfile(datapath, 'signal_templatespace', 'GroupWhole', 'TE_MeanVolume_MeanSignal', 'Data.mat'))

notests = 2;

Results_repandext{19,1} = 'EchoTime-TE30-supp 1';
Results_repandext{19,2} = ZShim_CalculateResults(no30,manual30,[],0,notests);

Results_repandext{20,1} = 'EchoTime-TE50-supp1';
Results_repandext{20,2} = ZShim_CalculateResults(no50,manual50,[],0,notests);

%% 3.1.4 z-shim effects in gray matter regions
clearvars -except Results_repandext datapath printStatus figurepath codepath saveResults
load(fullfile(datapath, 'signal_templatespace', 'GroupWhole', 'GM', 'SingleVolume_MeanSignal', 'Data.mat'))

noD = (noLEFT_DORSAL+noRIGHT_DORSAL)./2;
manualD = (manualLEFT_DORSAL+manualRIGHT_DORSAL)./2;
noV = (noLEFT_VENTRAL+noRIGHT_VENTRAL)./2;
manualV = (manualLEFT_VENTRAL+manualRIGHT_VENTRAL)./2;

noDMean = mean(noD,2);
manualDMean = mean(manualD,2);
noVMean = mean(noV,2);
manualVMean = mean(manualV,2);

notests = 3;

Results_repandext{21,1} = 'Gray matter-dorsal';
Results_repandext{21,2} = ZShim_CalculateResults(noD,manualD,[],0,notests);


Results_repandext{22,1} = 'Gray matter-ventral';
Results_repandext{22,2} = ZShim_CalculateResults(noV,manualV,[],0,notests);

Results_repandext{23,1} = 'Gray matter-ventral vs dorsal';
Results_repandext{23,2} = ZShim_CalculateResults([manualD-noD], [manualV-noV],[],0,notests);


% Supplementary material
% I. Repeated-measures ANOVA for mean of signal intensity
groups = [repmat({'dorsal'},48,1); repmat({'ventral'},48,1)];
t = table(groups, [noDMean;noVMean], [manualDMean;manualVMean], ...
    'VariableNames', {'group', 'no', 'manual'});

Meas = table({'No'; 'Manual'},'VariableNames',{'zmode'});
rm = fitrm(t,'no-manual~group','WithinDesign',Meas);

Results_repandext{24,1} = 'Gray matter-supp 1 (mean signal intensity)';
Results_repandext{24,2} = ranova(rm, 'WithinModel', 'zmode');

% II. Repeated-measures ANOVA for variation of signal intensity
zshimSubjectMeans = [noDMean manualDMean noVMean manualVMean];
zshimSubjectStds = [std(noD,[],2) std(manualD,[],2) std(noV,[],2) std(manualV,[],2)];

noDCoeffVar = zshimSubjectStds(:,1)./zshimSubjectMeans(:,1);
manualDCoeffVar = zshimSubjectStds(:,2)./zshimSubjectMeans(:,2);

noVCoeffVar = zshimSubjectStds(:,3)./zshimSubjectMeans(:,3);
manualVCoeffVar = zshimSubjectStds(:,4)./zshimSubjectMeans(:,4);

groups = [repmat({'dorsal'},48,1); repmat({'ventral'},48,1)];
t = table(groups, [noDMean;noVMean], [manualDMean;manualVMean], ...
    'VariableNames', {'group', 'no', 'manual'});

Meas = table({'No'; 'Manual'},'VariableNames',{'zmode'});
rm = fitrm(t,'no-manual~group','WithinDesign',Meas);

Results_repandext{25,1} = 'Gray matter-supp 2 (coeff var of signal intensity)';
Results_repandext{25,2} = ranova(rm, 'WithinModel', 'zmode');

% III. negative control anaysis - Left vs Right
noL = (noLEFT_DORSAL+noLEFT_VENTRAL)./2;
manualL = (manualLEFT_DORSAL+manualLEFT_VENTRAL)./2;
noR = (noRIGHT_DORSAL+noRIGHT_VENTRAL)./2;
manualR = (manualRIGHT_DORSAL+manualRIGHT_VENTRAL)./2;

noLMean = mean(noL,2);
manualLMean = mean(manualL,2);
noRMean = mean(noR,2);
manualRMean = mean(manualR,2);


groups = [repmat({'dorsal'},48,1); repmat({'ventral'},48,1)];
t = table(groups, [noLMean;noRMean], [manualLMean;manualRMean], ...
    'VariableNames', {'group', 'no', 'manual'});

Meas = table({'No'; 'Manual'},'VariableNames',{'zmode'});
rm = fitrm(t,'no-manual~group','WithinDesign',Meas);

Results_repandext{26,1} = 'Gray matter-supp 3 (Control- Left vs Right)';
Results_repandext{26,2} = ranova(rm, 'WithinModel', 'zmode');

% IV. Repetition of the analysis in the mean moco image
clearvars -except Results_repandext datapath printStatus figurepath codepath saveResults
load(fullfile(datapath, 'signal_templatespace', 'GroupWhole', 'GM', 'SingleVolume_MOCOMeanSignal', 'Data.mat'))

noD = (noLEFT_DORSAL+noRIGHT_DORSAL)./2;
manualD = (manualLEFT_DORSAL+manualRIGHT_DORSAL)./2;
noV = (noLEFT_VENTRAL+noRIGHT_VENTRAL)./2;
manualV = (manualLEFT_VENTRAL+manualRIGHT_VENTRAL)./2;

noDMean = mean(noD,2);
manualDMean = mean(manualD,2);
noVMean = mean(noV,2);
manualVMean = mean(manualV,2);

notests = 3;

groups = [repmat({'dorsal'},48,1); repmat({'ventral'},48,1)];
t = table(groups, [noDMean;noVMean], [manualDMean;manualVMean], ...
    'VariableNames', {'group', 'no', 'manual'});

Meas = table({'No'; 'Manual'},'VariableNames',{'zmode'});
rm = fitrm(t,'no-manual~group','WithinDesign',Meas);

Results_repandext{27,1} = 'Gray matter-supp 1';
Results_repandext{27,2} = ranova(rm, 'WithinModel', 'zmode');

Results_repandext{28,1} = 'Gray matter-dorsal';
Results_repandext{28,2} = ZShim_CalculateResults(noD,manualD,[],0,notests);


Results_repandext{29,1} = 'Gray matter-ventral';
Results_repandext{29,2} = ZShim_CalculateResults(noV,manualV,[],0,notests);

Results_repandext{30,1} = 'Gray matter-ventral vs dorsal';
Results_repandext{30,2} = ZShim_CalculateResults([noD-manualD], [noV-manualV],[],0,notests);


%% 3.1.5 z-shim effects on time-series data
clearvars -except Results_repandext datapath printStatus figurepath codepath saveResults
load(fullfile(datapath, 'signal_templatespace', 'GroupWhole', 'GM', 'TimeSeries_TSNR', 'Data.mat'))

notests = 1;

no = (noLEFT_DORSAL+noLEFT_VENTRAL+noRIGHT_DORSAL+noRIGHT_VENTRAL)./4;
manual = (manualLEFT_DORSAL+manualLEFT_VENTRAL+manualRIGHT_DORSAL+manualRIGHT_VENTRAL)./4;

Results_repandext{31,1} = 'Time-series';
Results_repandext{31,2} = ZShim_CalculateResults(no,manual,[],0,notests);

% II. find the most affected slices
clearvars -except Results_repandext datapath printStatus figurepath codepath saveResults
load(fullfile(datapath, 'signal_nativespace', 'GroupWhole', 'TimeSeries_TSNR', 'Data_picks.mat'))

noMean = mean(No_above3,2);
manualMean = mean(Manual_above3,2);

zshimSubjectMeans = [noMean manualMean];

meanGall = zshimSubjectMeans(:,1);
meanWall = zshimSubjectMeans(:,2);
meanPall = 100 * (meanWall./meanGall);
meanP_Increaseall = meanPall - 100;

zshimGlobalMeans = mean(zshimSubjectMeans);
meanG = zshimGlobalMeans(:,1);
meanW = zshimGlobalMeans(:,2);
meanP = 100 * (meanW./meanG);
meanP_Increase = meanP - 100;

Results_repandext{32,1} = 'Time-series/most affected slices (mean, min,max)';
Results_repandext{32,2} = [ meanP_Increase min(meanP_Increaseall) max(meanP_Increaseall)];

% Supplementary material 
% I.tSNR after censoring
clearvars -except Results_repandext datapath printStatus figurepath codepath saveResults
load(fullfile(datapath, 'signal_templatespace', 'GroupWhole', 'GM', 'TimeSeries_TSNR', 'Data_censored.mat'))

notests = 1;

Results_repandext{33,1} = 'Time-series-supp 1 (tSNR after censoring) ';
Results_repandext{33,2} = ZShim_CalculateResults(no,manual,[],0,notests);

% II. find the most affected slices
clearvars -except Results_repandext datapath printStatus figurepath codepath saveResults
load(fullfile(datapath, 'signal_nativespace', 'GroupWhole', 'TimeSeries_TSNR', 'Data_picks_censored.mat'))

noMean = mean(No_above3,2);
manualMean = mean(Manual_above3,2);

zshimSubjectMeans = [noMean manualMean];

meanGall = zshimSubjectMeans(:,1);
meanWall = zshimSubjectMeans(:,2);
meanPall = 100 * (meanWall./meanGall);
meanP_Increaseall = meanPall - 100;

zshimGlobalMeans = mean(zshimSubjectMeans);
meanG = zshimGlobalMeans(:,1);
meanW = zshimGlobalMeans(:,2);
meanP = 100 * (meanW./meanG);
meanP_Increase = meanP - 100;

Results_repandext{34,1} = 'Time-series-supp 2/most affected slices (mean, min,max)';
Results_repandext{34,2} = [ meanP_Increase min(meanP_Increaseall) max(meanP_Increaseall)];
%%
if saveResults 
    save([datapath filesep 'Results_Replication_and_extensionofpreviousfindings.mat'], 'Results_repandext')
end
%% do the plots in Figure 2
ZShim_Figures_FigureII(printStatus,datapath,figurepath)
