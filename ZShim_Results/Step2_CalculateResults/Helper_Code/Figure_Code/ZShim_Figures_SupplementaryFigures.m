%Supplementary figures
clc;clear all;close all
datapath = '/data/pt_02098/ZShim_BIDS_0807/derivatives/extracted_signal/'; %where the organized folder is to produce the results
codepath = '/data/pt_02098/ZShim_BIDS_0807/derivatives/code/';             %where all the codes are
addpath(genpath(codepath))                                                 %add code folder & subfolders to path 
processdatapath = '/data/pt_02098/ZShim_BIDS_0807/derivatives/';           %where processed data is 
printStatus = 0;                                                           %to print figures as svg
figurepath  = '/data/pt_02098/ZShim_BIDS_0807/derivatives/figures/';       %set the path for saving the figures
%% Supplementary Figure 1
%1. load the data
load([datapath 'signal_nativespace'  filesep 'GroupWhole' filesep 'SingleVolume_MeanSignal' filesep 'Data.mat'])

figYLim      = [30 150];
figPosition  = [0 0 1200 600];
figColor     = [1 1 1];
figYLabel    = 'Slices (inferior to superior)';
figXLabel    =  'Subjects';

%2.plot intensity matrices
figure; hold on; set(gcf, 'color', figColor , 'position', figPosition); box off

subplot(1,2,1)
imagesc(no')
colormap('hot')
caxis(figYLim)
ylabel(figYLabel)
xlabel(figXLabel)

subplot(1,2,2)
imagesc(manual')
colormap('hot')
caxis(figYLim)

if printStatus
   
    cd(path)
    
    print(gcf, '-painters','-dsvg', ...
        [ path filesep 'Figure_Supp_1.svg']);
end
%% Supplementary Figure 2
%1. load data
load([datapath 'signal_templatespace'  filesep 'GroupWhole' filesep 'SingleVolume_MeanSignal' filesep 'Data.mat'])

figPosition  = [801 1 300 700];
figColor     = [1 1 1];
figXLim      = [30 110];
figYLim      = [1 226];
figXLabel    = 'Whole cord signal intensity (a.u.)';
figYLabel    = 'Inferior to Superior';

zNo = no;

% CREATE FIGURE WITH CURVES OF SLICE-BY-SLICE DATA
numSubjects = size(zNo,1);
numSlices   = size(zNo,2);

zNoStats     = [mean(zNo); std(zNo)/sqrt(numSubjects)];

figure; hold on; set(gcf, 'color', figColor , 'position', figPosition); grid on; box off

f1 = fill([zNoStats(1,:)-zNoStats(2,:) fliplr(zNoStats(1,:)+zNoStats(2,:))], ...
    [1:numSlices numSlices:-1:1], 'r');
set(f1, 'facecolor', [1 0 0], 'facealpha', 0.25, 'edgecolor', 'none');
plot(zNoStats(1,:), 1:numSlices, 'color', [1 0 0], 'linewidth', 2);

set(gca,'xticklabel',[],'yticklabel',[],'YColor', [0.5 0.5 0.5],'XColor', [0.5 0.5 0.5],...
    'GridAlpha',0.3, 'xlim', figXLim, 'ylim', figYLim)
xlabel(figXLabel)
ylabel(figYLabel)


if printStatus
   
    cd(path)
    
    print(gcf, '-painters','-dsvg', ...
        [ path filesep 'Figure_Supp_2.svg']);
end


%% Supplementary Figure 5

%I. load data

load([datapath 'signal_templatespace'  filesep 'GroupSingle_EPI' filesep 'GM' filesep 'TimeSeries_TSNR' filesep 'Data.mat'])

zNo  =(noLEFT_DORSAL+noLEFT_VENTRAL+noRIGHT_DORSAL+noRIGHT_VENTRAL)./4;
zManual  =(manualLEFT_DORSAL+manualLEFT_VENTRAL+manualRIGHT_DORSAL+manualRIGHT_VENTRAL)./4;
zAuto = (autoLEFT_DORSAL+autoLEFT_VENTRAL+autoRIGHT_DORSAL+autoRIGHT_VENTRAL)./4;

zshimNo_zref = zNo;
zshimAuto_zref = zAuto;
zshimManual_zref = zManual;

clear zNo zAuto zManual

load([datapath 'signal_templatespace'  filesep 'GroupSingle_FM' filesep 'GM' filesep 'TimeSeries_TSNR' filesep 'Data.mat'])

zNo  =(noLEFT_DORSAL+noLEFT_VENTRAL+noRIGHT_DORSAL+noRIGHT_VENTRAL)./4;
zManual  =(manualLEFT_DORSAL+manualLEFT_VENTRAL+manualRIGHT_DORSAL+manualRIGHT_VENTRAL)./4;
zAuto = (autoLEFT_DORSAL+autoLEFT_VENTRAL+autoRIGHT_DORSAL+autoRIGHT_VENTRAL)./4;

zshimNo_fm = zNo;
zshimManual_fm = zManual;
zshimAuto_fm = zAuto;

noMean_zref = mean(zshimNo_zref,2);
autoMean_zref = mean(zshimAuto_zref,2);
manualMean_zref = mean(zshimManual_zref,2);

noMean_fm = mean(zshimNo_fm,2);
autoMean_fm = mean(zshimAuto_fm,2);
manualMean_fm = mean(zshimManual_fm,2);


figYLim      = [10 20];
figPosition  = [0 0 600 600];
figColor     = [1 1 1];
figXLabel    = [];
figXLim      = [0 4];
XTLabels     = {'z-ref-based group', 'FM-based group'};
figYLabel    =  'mean gray matter tSNR';

figure; hold on; set(gcf, 'color', figColor , 'position', figPosition); box off


positions1 = [0.5 2.5 4.5 6.5]+0.1;
positions2 = [1 3 5 7 ]-0.1 + 0.3;

colors1 = {[0 1 0] [1 0 0 ]};
colors = {[0.5 0.5 0.5]}

% 
data1 = [autoMean_zref, autoMean_fm];
data2 = [manualMean_zref, manualMean_fm ];
datasets = {autoMean_zref manualMean_zref autoMean_fm manualMean_fm }


for c = 1:2
    h=boxplot(data1(:,c), 'positions',positions1(c), 'Colors',colors1{1} ,'Widths', [0.1] ,'Symbol','.r');
    g=boxplot(data2(:,c), 'positions', positions2(c) , 'Colors', colors1{2},'Widths', [0.1] ,'Symbol','.r');
    set(h,'LineWidth',2 )
    set(g,'LineWidth', 2)
 
    clear h g
    

end
xticks([(positions1+positions2)./2])
positions = {positions1(1)+0.07 positions2(1)-0.07 positions1(2)+0.07 positions2(2)-0.07 ...
    positions1(3)+0.07 positions2(3)-0.07 positions1(4)+0.07 positions2(4)-0.07}

counter   = 0;

for p = [1:2:4]
     for j = 1:length(datasets{1})
        line([positions{p} positions{p+1}], [datasets{p}(j) datasets{p+1}(j)], ...
            'color', colors{1}, 'linewidth', 0.5);
    end
end

set(gca,'xticklabel',XTLabels,...
    'GridAlpha',0.5, 'xlim', figXLim, 'ylim', figYLim)
xlabel(figXLabel)
ylabel(figYLabel)

if printStatus
   
    cd(path)
    
    print(gcf, '-painters','-dsvg', ...
        [ path filesep 'Figure_Supp_5.svg']);
end
