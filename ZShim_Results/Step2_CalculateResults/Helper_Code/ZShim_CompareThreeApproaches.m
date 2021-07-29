function results_3approaches = ZShim_CompareThreeApproaches(rawdatapath,processdatapath, datapath, recalculateResults)
% Compare two automated approaches and manual approach
% use data acquired with two automated methods and the reconstructed 
% artificial volumes

% ----------
% Inputs:
% ----------

% scttemplatepath:   fullpath, string, PAM50 template location
% rawdatapath:       fullpath, string, raw data location
% processeddatapath: fullpath, string, processed data location
% recalculateResults: True or False

% Merve Kaptan, mkaptan@cbs.mpg.de

%% part I
%1. Compare the three approaches using independent samples t-tests
load(fullfile(datapath, 'signal_templatespace', 'GroupSingle_EPI', 'GM', 'TimeSeries_TSNR', 'Data.mat'))

no_epi      = (noLEFT_DORSAL+noLEFT_VENTRAL+noRIGHT_DORSAL+noRIGHT_VENTRAL)./4;
manual_epi  = (manualLEFT_DORSAL+manualLEFT_VENTRAL+manualRIGHT_DORSAL+manualRIGHT_VENTRAL)./4;
auto_epi    = (autoLEFT_DORSAL+autoLEFT_VENTRAL+autoRIGHT_DORSAL+autoRIGHT_VENTRAL)./4;

clear no manual auto

load(fullfile(datapath, 'signal_templatespace', 'GroupSingle_FM', 'GM', 'TimeSeries_TSNR', 'Data.mat'))

no_fm      = (noLEFT_DORSAL+noLEFT_VENTRAL+noRIGHT_DORSAL+noRIGHT_VENTRAL)./4;
manual_fm  = (manualLEFT_DORSAL+manualLEFT_VENTRAL+manualRIGHT_DORSAL+manualRIGHT_VENTRAL)./4;
auto_fm    = (autoLEFT_DORSAL+autoLEFT_VENTRAL+autoRIGHT_DORSAL+autoRIGHT_VENTRAL)./4;

no_epi_mean     = mean(no_epi,2);
auto_epi_mean   = mean(auto_epi,2);
manual_epi_mean = mean(manual_epi,2);

no_fm_mean     = mean(no_fm,2);
auto_fm_mean   = mean(auto_fm,2);
manual_fm_mean = mean(manual_fm,2);

nooftests = 4;

[~,p] = ttest2(no_fm_mean, no_epi_mean);
p = p*nooftests;

results{1,1} = 'mean signal intensity';
results{1,2} = 'p value (Bonferroni corrected)';

results{2,1} = 'no vs no (FM vs EPI-based groups)';
results{2,2} = p;

clear p
[~,p] = ttest2(manual_fm_mean-no_fm_mean,manual_epi_mean-no_epi_mean);
p = p*nooftests;

results{3,1} = 'manual-no vs manual-no (FM vs EPI-based groups)';
results{3,2} = p;

clear p
[~,p] = ttest2(auto_fm_mean-no_fm_mean,auto_epi_mean-no_epi_mean);
p = p*nooftests;

results{4,1} = 'auto-no vs auto-no (FM vs EPI-based groups)';
results{4,2} = p;

clear p
[~,p] = ttest2(manual_fm_mean-auto_fm_mean,manual_epi_mean-auto_epi_mean);
p = p*nooftests;

results{5,1} = 'manual-auto vs manual-auto (FM vs EPI-based groups)';
results{5,2} = p;

zshimSubjectMeans_epi = [no_epi_mean auto_epi_mean manual_epi_mean];
zshimSubjectStds_epi = [std(no_epi,[],2) std(auto_epi,[],2) std(manual_epi,[],2)];

no_epi_coeffVar = zshimSubjectStds_epi(:,1)./zshimSubjectMeans_epi(:,1);
auto_epi_coeffVar = zshimSubjectStds_epi(:,2)./zshimSubjectMeans_epi(:,2);
manual_epi_coeffVar = zshimSubjectStds_epi(:,3)./zshimSubjectMeans_epi(:,3);

zshimSubjectMeans_fm = [no_fm_mean auto_fm_mean manual_fm_mean];
zshimSubjectStds_fm = [std(no_fm,[],2) std(auto_fm,[],2) std(manual_fm,[],2)];

no_fm_coeffVar = zshimSubjectStds_fm(:,1)./zshimSubjectMeans_fm(:,1);
auto_fm_coeffVar = zshimSubjectStds_fm(:,2)./zshimSubjectMeans_fm(:,2);
manual_fm_coeffVar = zshimSubjectStds_fm(:,3)./zshimSubjectMeans_fm(:,3);

results{6,1} = 'coefficient of var';

[~,p] = ttest2(no_epi_coeffVar, no_fm_coeffVar);
p = p*nooftests;

results{7,1} = 'no vs no (FM vs EPI-based groups)';
results{7,2} = p;

[~,p] = ttest2(manual_fm_coeffVar-no_fm_coeffVar,manual_epi_coeffVar-no_epi_coeffVar);
p = p*nooftests;

results{8,1} = 'manual-no vs manual-no (FM vs EPI-based groups)';
results{8,2} = p;

[~,p] = ttest2(auto_fm_coeffVar-no_fm_coeffVar,auto_epi_coeffVar-no_epi_coeffVar);
p = p*nooftests;

results{9,1} = 'auto-no vs auto-no (FM vs EPI-based groups)';
results{9,2} = p;

[~,p] = ttest2(manual_fm_coeffVar-auto_fm_coeffVar,manual_epi_coeffVar-auto_epi_coeffVar);
p = p*nooftests;

results{10,1} = 'manual-auto vs manual-auto (FM vs EPI-based groups)';
results{10,2} = p;

results_3approaches{1,1} = 'Comparison of three approaches- ttests';
results_3approaches{1,2} = results;

%% part II
clear results
% 2. compare the similarity of z-shim moments

manual_picks = readtable([processdatapath filesep 'manualzshimPicks_duringScan.csv']);
auto_picks   = readtable([processdatapath filesep 'autozshimPicks_duringScan.csv']);

manual_picks = table2array(manual_picks(:,3:end));
auto_picks   = table2array(auto_picks(:,3:end));

for sub = 1:size(manual_picks,1)
    
    corrs(sub,1) = corr(manual_picks(sub,:)', auto_picks(sub,:)','type', 'Spearman');
    dists(sub,1) = norm(manual_picks(sub,:)'-auto_picks(sub,:)');
    
end

corrs_epi = corrs([3:24 47 48]); %indices for epi group
corrs_fm  = corrs([1 2 25:46]);  %indices for fm group

dist_epi = dists([3:24 47 48]); %indices for epi group
dist_fm  = dists([1 2 25:46]);  %indices for fm group

results{1,1} = 'Correlations EPI;FM mean,min,max';
results{1,2} = [mean(corrs_epi) min(corrs_epi) max(corrs_epi); mean(corrs_fm) min(corrs_fm) max(corrs_fm)];

results{2,1} = 'Comparison EPI vs FM- correlation';
results{2,2} = 'p';
results{2,3} = 't-value';
results{2,4} = 'dof';

[~,p,~,stats] = ttest2(atanh(corrs_epi), atanh(corrs_fm));

results{3,2} = p;
results{3,3} = stats.tstat;
results{3,4} = stats.df;

results_3approaches{2,1} = 'Comparison of three approaches- correlation';
results_3approaches{2,2} = results;

clear results

results{1,1} = 'Euclidian distance EPI;FM mean,min,max';
results{1,2} = [mean(dist_epi) min(dist_epi) max(dist_epi); mean(dist_fm) min(dist_fm) max(dist_fm)];

results{2,1} = 'Comparison EPI vs FM- ED';
results{2,2} = 'p';
results{2,3} = 't-value';
results{2,4} = 'dof';

[~,p,~,stats] = ttest2(atanh(dist_epi), atanh(dist_fm));

results{3,2} = p;
results{3,3} = stats.tstat;
results{3,4} = stats.df;

results_3approaches{3,1} = 'Comparison of three approaches- euclidian distance';
results_3approaches{3,2} = results;

clear results
%% part III

% 3. compare the signal in the reconstructed volumes
% note that with the ZShim_Create_ArtificialVolumes(subGroup, scttemplatepath,
% rawdatapath,processdatapath) the reconstructed volumes are already
% prepared (create reconstructed volumes, normalize, extract signal in gray matter).
% load the signal and compare the signal intensity
if recalculateResults

cd(rawdatapath)
subjects = dir('*sub-*'); % get the list of the subjects

reconsName  = 'ZShim1_Recons';
reconsMode    = {'gre', 'zref', 'manual', 'no'};

recons_signal = ZShim_Load_ReconsSignal(processdatapath,subjects,reconsName,reconsMode);

for r = 1:numel(reconsMode)
    
    eval([reconsMode{r} '1'  '= recons_signal(' num2str(r) ');'])
    
end

reconsName  = 'ZShim2_Recons';
reconsMode    = {'gre', 'zref', 'manual', 'no'};
recons_signal = ZShim_Load_ReconsSignal(processdatapath,subjects,reconsName,reconsMode);

for r = 1:numel(reconsMode)
    
    eval([reconsMode{r} '2'  '= recons_signal(' num2str(r) ');'])
    
end

else
    
    load(fullfile(processdatapath, 'extracted_signal', 'signal_templatespace', 'GroupWhole', ...
        'ReconstructedSignal', 'Compare_ThreeApproaches', 'results.mat'))
   
end

% compare baselines of no z-shimming and benefit of z-shimming over time
[~,p1, ~,stat1] = ttest(nanmean(no1,2), nanmean(no2,2));
p1 = p1 *4;
[~,p2, ~,stat2] = ttest([nanmean(manual1,2)- nanmean(no1,2)],[nanmean(manual2,2)- nanmean(no2,2)]);
p2 = p2 *4;
[~,p3, ~,stat3] = ttest([nanmean(zref1,2)- nanmean(no1,2)],[nanmean(zref2,2)- nanmean(no2,2)]);
p3 = p3 *4;
[~,p4, ~,stat4] = ttest([nanmean(gre1,2)- nanmean(no1,2)],[nanmean(gre2,2)- nanmean(no2,2)]);
p4 = p4 *4;

results1{1,1} = 'no 1 vs no 2';
results1{1,2} = 'p';
results1{1,3} = 't-value';
results1{1,4} = 'dof';

results1{2,2} = p1;
results1{2,3} = stat1.tstat;
results1{2,4} = stat1.df;

results2{1,1} = 'manual-no 1 vs manual-no 2';
results2{1,2} = 'p';
results2{1,3} = 't-value';
results2{1,4} = 'dof';

results2{2,2} = p2;
results2{2,3} = stat2.tstat;
results2{2,4} = stat2.df;

results3{1,1} = 'epi based-no 1 vs epi based -no 2';
results3{1,2} = 'p';
results3{1,3} = 't-value';
results3{1,4} = 'dof';

results3{2,2} = p3;
results3{2,3} = stat3.tstat;
results3{2,4} = stat3.df;

results4{1,1} = 'fm based-no 1 vs fm based -no 2';
results4{1,2} = 'p';
results4{1,3} = 't-value';
results4{1,4} = 'dof';

results4{2,2} = p4;
results4{2,3} = stat4.tstat;
results4{2,4} = stat4.df;

results = {results1;results2;results3;results4};

results_3approaches{4,1} = 'Comparison of three approaches- z-shim effects over time';
results_3approaches{4,2} = results;

clear results results1 results2 results3 results4
% compare benefit of z-shimming - reconstruction of 1st zref
[~,p1, ~,stat1] = ttest(nanmean(manual1,2),nanmean(no1,2));
p1 = p1 *3;
[~,p2, ~,stat2] = ttest(nanmean(zref1,2), nanmean(no1,2));
p2 = p2 *3;
[~,p3, ~,stat3] = ttest(nanmean(gre1,2), nanmean(no1,2));
p3 = p3 *3;

zshimSubjectMeans = [nanmean(no1,2) nanmean(manual1,2) nanmean(zref1,2) nanmean(gre1,2)];
zshimGlobalMeans = mean(zshimSubjectMeans);
meanG = zshimGlobalMeans(1);
meanW = zshimGlobalMeans(2:end);
meanP = 100 * (meanW/meanG);
meanP_Increase = meanP - 100;

results{1,1} = 'Zref1 manual vs no';
results{1,2} = 'p';
results{1,3} = 't-value';
results{1,4} = 'dof';

results{2,2} = p1;
results{2,3} = stat1.tstat;
results{2,4} = stat1.df;

results{3,1} = 'Zref1 epi based vs no';
results{3,2} = 'p';
results{3,3} = 't-value';
results{3,4} = 'dof';

results{4,2} = p2;
results{4,3} = stat2.tstat;
results{4,4} = stat2.df;

results{5,1} = 'Zref1 fm based vs no';
results{5,2} = 'p';
results{5,3} = 't-value';
results{5,4} = 'dof';

results{6,2} = p3;
results{6,3} = stat3.tstat;
results{6,4} = stat3.df;

results{1,5} = 'mean % increase manual (vs no)';
results{2,5} = meanP_Increase(1);

results{3,5} = 'mean % increase epi-based (vs no)';
results{4,5} = meanP_Increase(2);

results{5,5} = 'mean % increase fm-based (vs no)';
results{6,5} = meanP_Increase(3);

results_3approaches{5,1} = 'Comparison of reconstructed volumes- ZRef I';
results_3approaches{5,2} = results;

clear results stat1 stat2 stat3
%  compare benefit of z-shimming - reconstruction of 2nd zref
[~,p1, ~,stat1] = ttest(nanmean(manual2,2),nanmean(no2,2));
p1 = p1 *3;
[~,p2, ~,stat2] = ttest(nanmean(zref2,2), nanmean(no2,2));
p2 = p2 *3;
[~,p3, ~,stat3] = ttest(nanmean(gre2,2), nanmean(no2,2));
p3 = p3 *3;

zshimSubjectMeans = [nanmean(no2,2) nanmean(manual2,2) nanmean(zref2,2) nanmean(gre2,2)];

zshimGlobalMeans = mean(zshimSubjectMeans);
meanG = zshimGlobalMeans(1);
meanW = zshimGlobalMeans(2:end);
meanP = 100 * (meanW/meanG);
meanP_Increase = meanP - 100;

results{1,1} = 'Zref1 manual vs no';
results{1,2} = 'p';
results{1,3} = 't-value';
results{1,4} = 'dof';

results{2,2} = p1;
results{2,3} = stat1.tstat;
results{2,4} = stat1.df;

results{3,1} = 'Zref1 epi based vs no';
results{3,2} = 'p';
results{3,3} = 't-value';
results{3,4} = 'dof';

results{4,2} = p2;
results{4,3} = stat2.tstat;
results{4,4} = stat2.df;

results{5,1} = 'Zref1 fm based vs no';
results{5,2} = 'p';
results{5,3} = 't-value';
results{5,4} = 'dof';

results{6,2} = p3;
results{6,3} = stat3.tstat;
results{6,4} = stat3.df;

results{1,5} = 'mean % increase manual (vs no)';
results{2,5} = meanP_Increase(1);

results{3,5} = 'mean % increase epi-based (vs no)';
results{4,5} = meanP_Increase(2);

results{5,5} = 'mean % increase fm-based (vs no)';
results{6,5} = meanP_Increase(3);

results_3approaches{6,1} = 'Comparison of reconstructed volumes- ZRef II';
results_3approaches{6,2} = results;
