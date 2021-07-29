function ZShim_Figures_FigureIII(printStatus,datapath,figurepath)

% Figure 3

% Merve Kaptan
% mkaptan@cbs.mpg.de

if  printStatus
    
    path =  fullfile(figurepath, 'Figure3');
    
    if ~exist(path)
        mkdir(path)
    end
    
end
%% Line plots

load([datapath 'signal_templatespace'  filesep 'GroupSingle_EPI' filesep 'GM' filesep 'TimeSeries_TSNR' filesep 'Data.mat'])

figPosition  = [801 1 300 700];
figColor     = [1 1 1];
figXLim      = [5 20];
figYLim      = [1 226];
figXLabel    = 'Gray matter tSNR';
figYLabel    = 'Inferior to Superior';
numSlices   = 226;

figure; hold on; set(gcf, 'color', figColor , 'position', figPosition); grid on; box off

auto  =(autoLEFT_DORSAL+autoLEFT_VENTRAL+autoRIGHT_DORSAL+autoRIGHT_VENTRAL)./4;
manual  =(manualLEFT_DORSAL+manualLEFT_VENTRAL+manualRIGHT_DORSAL+manualRIGHT_VENTRAL)./4;
no    = (noLEFT_DORSAL+noLEFT_VENTRAL+noRIGHT_DORSAL+noRIGHT_VENTRAL)./4;

numSubjects = 24;

zNoStats     = [mean(no); std(no)/sqrt(numSubjects)];
zManualStats = [mean(manual); std(manual)/sqrt(numSubjects)];
zAutoStats = [mean(auto); std(auto)/sqrt(numSubjects)];

xlim(figXLim); ylim(figYLim);
f1 = fill([zNoStats(1,:)-zNoStats(2,:) fliplr(zNoStats(1,:)+zNoStats(2,:))], ...
    [1:numSlices numSlices:-1:1], 'r');
f2 = fill([zManualStats(1,:)-zManualStats(2,:) fliplr(zManualStats(1,:)+zManualStats(2,:))], ...
    [1:numSlices numSlices:-1:1], 'b');
f3 = fill([zAutoStats(1,:)-zAutoStats(2,:) fliplr(zAutoStats(1,:)+zAutoStats(2,:))], ...
    [1:numSlices numSlices:-1:1], 'g');

set(f1, 'facecolor', [1 0 0], 'facealpha', 0.25, 'edgecolor', 'none');
set(f2, 'facecolor', [0 0 1], 'facealpha', 0.25, 'edgecolor', 'none');
set(f3, 'facecolor', [0 1 0], 'facealpha', 0.25, 'edgecolor', 'none');

p1 = plot(zNoStats(1,:), 1:numSlices, 'color', [1 0 0], 'linewidth', 2); 
p2 = plot(zManualStats(1,:), 1:numSlices, 'color', [0 0 1], 'linewidth', 2);
p3 = plot(zAutoStats(1,:), 1:numSlices, 'color', [0 1 0], 'linewidth', 2);

set(gca, 'YColor', [0.8 0.8 0.8]);
set(gca, 'XColor', [0.8 0.8 0.8]);

xlabel(figXLabel)
ylabel(figYLabel)


if printStatus

    cd(path)   
    print(figure, '-painters','-dsvg', ...
        [ path filesep 'Figure_3_EPI_I.svg']);
end


clearvars -except printStatus datapath figurepath
load(fullfile(datapath, 'signal_templatespace', 'GroupSingle_FM', 'GM', 'TimeSeries_TSNR', 'Data.mat'))

figPosition  = [801 1 300 700];
figColor     = [1 1 1];
figXLim      = [5 20];
figYLim      = [1 226];
figXLabel    = 'Gray matter tSNR';
figYLabel    = 'Inferior to Superior';
numSlices   = 226;


figure; hold on; set(gcf, 'color', figColor , 'position', figPosition); grid on; box off

auto  =(autoLEFT_DORSAL+autoLEFT_VENTRAL+autoRIGHT_DORSAL+autoRIGHT_VENTRAL)./4;
manual  =(manualLEFT_DORSAL+manualLEFT_VENTRAL+manualRIGHT_DORSAL+manualRIGHT_VENTRAL)./4;
no    = (noLEFT_DORSAL+noLEFT_VENTRAL+noRIGHT_DORSAL+noRIGHT_VENTRAL)./4;

numSubjects = 24;

zNoStats     = [mean(no); std(no)/sqrt(numSubjects)];
zManualStats = [mean(manual); std(manual)/sqrt(numSubjects)];
zAutoStats = [mean(auto); std(auto)/sqrt(numSubjects)];

xlim(figXLim); ylim(figYLim);
f1 = fill([zNoStats(1,:)-zNoStats(2,:) fliplr(zNoStats(1,:)+zNoStats(2,:))], ...
    [1:numSlices numSlices:-1:1], 'r');
f2 = fill([zManualStats(1,:)-zManualStats(2,:) fliplr(zManualStats(1,:)+zManualStats(2,:))], ...
    [1:numSlices numSlices:-1:1], 'b');
f3 = fill([zAutoStats(1,:)-zAutoStats(2,:) fliplr(zAutoStats(1,:)+zAutoStats(2,:))], ...
    [1:numSlices numSlices:-1:1], 'g');

set(f1, 'facecolor', [1 0 0], 'facealpha', 0.25, 'edgecolor', 'none');
set(f2, 'facecolor', [0 0 1], 'facealpha', 0.25, 'edgecolor', 'none');
set(f3, 'facecolor', [0 1 0], 'facealpha', 0.25, 'edgecolor', 'none');

p1 = plot(zNoStats(1,:), 1:numSlices, 'color', [1 0 0], 'linewidth', 2); 
p2 = plot(zManualStats(1,:), 1:numSlices, 'color', [0 0 1], 'linewidth', 2);
p3 = plot(zAutoStats(1,:), 1:numSlices, 'color', [0 1 0], 'linewidth', 2);

set(gca, 'YColor', [0.8 0.8 0.8]);
set(gca, 'XColor', [0.8 0.8 0.8]);

xlabel(figXLabel)
ylabel(figYLabel)

if printStatus

    cd(path)   
    print(figure, '-painters','-dsvg', ...
        [ path filesep 'Figure_3_FM_I.svg']);
end

%% Bland-Altmann plots

clearvars -except printStatus datapath figurepath
load([datapath 'signal_templatespace'  filesep 'GroupSingle_EPI' filesep 'GM' filesep 'TimeSeries_TSNR' filesep 'Data.mat'])

figPosition   = [0 0 600 600];
figColor      = [1 1 1];
figXLim1      = [10 20];
figYLim1      = [10 20];
figXLim2      = [10 20];
figYLim2      = [-3 2];
figXLabel1    = 'Gray matter tSNR: manual';
figYLabel1    = 'Gray matter tSNR: auto';
figXLabel2    = 'Gray matter tSNR: manual';
figYLabel2    = 'Gray matter tSNR: auto - manual';
SizeMarker = 200;

auto  =(autoLEFT_DORSAL+autoLEFT_VENTRAL+autoRIGHT_DORSAL+autoRIGHT_VENTRAL)./4;
manual  =(manualLEFT_DORSAL+manualLEFT_VENTRAL+manualRIGHT_DORSAL+manualRIGHT_VENTRAL)./4;

auto_mean = nanmean(auto,2);
manual_mean = nanmean(manual,2);

figure; hold on; set(gcf, 'color', figColor , 'position', figPosition); grid on; box off
subplot(2,1,1); hold on;

scatter(manual_mean, auto_mean, SizeMarker, [0.5 0.5 0.5], 'filled', 'MarkerFaceAlpha',0.5)
xlim(figXLim1)
ylim(figYLim1)
axis square

x = figXLim1(1):0.05:figXLim1(2);
y = figYLim1(1):0.05:figYLim1(2);
plot(x,y, 'Color', [0.5 0.5 0.5])

xlabel(figXLabel1)
ylabel(figYLabel1)

subplot(2,1,2);hold on;
scatter(manual_mean, auto_mean-manual_mean, SizeMarker, [0.5 0.5 0.5], 'filled', 'MarkerFaceAlpha',0.5)
xlim(figXLim2)
ylim(figYLim2)
axis square

yline(mean(auto_mean-manual_mean))
yline(mean(auto_mean-manual_mean)+1.96*std(auto_mean-manual_mean))
yline(mean(auto_mean-manual_mean)-1.96*std(auto_mean-manual_mean))

xlabel(figXLabel2)
ylabel(figYLabel2)

if printStatus

    cd(path)   
    print(figure, '-painters','-dsvg', ...
        [ path filesep 'Figure_3_EPI_II.svg']);
end

clc
clearvars -except printStatus datapath figurepath

load(fullfile(datapath, 'signal_templatespace', 'GroupSingle_FM', 'GM', 'TimeSeries_TSNR', 'Data.mat'))

figPosition   = [0 0 600 600];
figColor      = [1 1 1];
figXLim1      = [10 20];
figYLim1      = [10 20];
figXLim2      = [10 20];
figYLim2      = [-3 2];
figXLabel1    = 'Gray matter tSNR: manual';
figYLabel1    = 'Gray matter tSNR: auto';
figXLabel2    = 'Gray matter tSNR: manual';
figYLabel2    = 'Gray matter tSNR: auto';
SizeMarker = 200;

auto  =(autoLEFT_DORSAL+autoLEFT_VENTRAL+autoRIGHT_DORSAL+autoRIGHT_VENTRAL)./4;
manual  =(manualLEFT_DORSAL+manualLEFT_VENTRAL+manualRIGHT_DORSAL+manualRIGHT_VENTRAL)./4;

auto_mean = nanmean(auto,2);
manual_mean = nanmean(manual,2);

figure; hold on; set(gcf, 'color', figColor , 'position', figPosition); grid on; box off
subplot(2,1,1); hold on;

scatter(manual_mean, auto_mean, SizeMarker, [0.5 0.5 0.5], 'filled', 'MarkerFaceAlpha',0.5)
xlim(figXLim1)
ylim(figYLim1)
axis square

x = figXLim1(1):0.05:figXLim1(2);
y = figYLim1(1):0.05:figYLim1(2);
plot(x,y, 'Color', [0.5 0.5 0.5])

xlabel(figXLabel1)
ylabel(figYLabel1)

subplot(2,1,2);hold on;
scatter(manual_mean, auto_mean-manual_mean, SizeMarker, [0.5 0.5 0.5], 'filled', 'MarkerFaceAlpha',0.5)
xlim(figXLim2)
ylim(figYLim2)
axis square
yline(mean(auto_mean-manual_mean))
yline(mean(auto_mean-manual_mean)+1.96*std(auto_mean-manual_mean))
yline(mean(auto_mean-manual_mean)-1.96*std(auto_mean-manual_mean))

xlabel(figXLabel2)
ylabel(figYLabel2)

if printStatus

    cd(path)   
    print(figure, '-painters','-dsvg', ...
        [ path filesep 'Figure_3_FM_II.svg']);
end

end


