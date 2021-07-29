function ZShim_Figures_FigureII(printStatus,datapath,figurepath)

% Figure 2

% Merve Kaptan
% mkaptan@cbs.mpg.de

if  printStatus
    
    path =  fullfile(figurepath, 'Figure2');
    
    if ~exist(path)
        mkdir(path)
    end
    
end
%% A - line plots 
load(fullfile(datapath, 'signal_templatespace', 'GroupWhole', 'SingleVolume_MeanSignal', 'Data.mat'))

figPosition  = [801 1 300 700];
figColor     = [1 1 1];
figXLim      = [30 110];
figYLim      = [1 226];
figXLabel    = 'Whole cord signal intensity (a.u.)';
figYLabel    = 'Inferior to Superior';

zNo = no;
zManual = manual;

% CREATE FIGURE WITH CURVES OF SLICE-BY-SLICE DATA
numSubjects = size(zNo,1);
numSlices   = size(zNo,2);

zNoStats     = [mean(zNo); std(zNo)/sqrt(numSubjects)];
zManualStats = [mean(zManual); std(zManual)/sqrt(numSubjects)];

figure; hold on; set(gcf, 'color', figColor , 'position', figPosition); grid on; box off

f1 = fill([zNoStats(1,:)-zNoStats(2,:) fliplr(zNoStats(1,:)+zNoStats(2,:))], ...
    [1:numSlices numSlices:-1:1], 'r');
f2 = fill([zManualStats(1,:)-zManualStats(2,:) fliplr(zManualStats(1,:)+zManualStats(2,:))], ...
    [1:numSlices numSlices:-1:1], 'b');
set(f1, 'facecolor', [1 0 0], 'facealpha', 0.25, 'edgecolor', 'none');
set(f2, 'facecolor', [0 0 1], 'facealpha', 0.25, 'edgecolor', 'none');
plot(zNoStats(1,:), 1:numSlices, 'color', [1 0 0], 'linewidth', 2);
plot(zManualStats(1,:), 1:numSlices, 'color', [0 0 1], 'linewidth', 2);

set(gca,'xticklabel',[],'yticklabel',[],'YColor', [0.8 0.8 0.8],'XColor', [0.8 0.8 0.8],...
    'GridAlpha',0.5, 'xlim', figXLim, 'ylim', figYLim)
xlabel(figXLabel)
ylabel(figYLabel)

if printStatus

    cd(path)   
    print(figure, '-painters','-dsvg', ...
        [ path filesep 'Figure_2A.svg']);
end

%% B

clc
clearvars -except printStatus datapath figurepath

load(fullfile(datapath, 'signal_nativespace', 'GroupWhole', 'SingleVolume_MeanSignal', 'Data_picks.mat'))

figXLim      = [0.5 5.5];
figYLim      = [0 125];
figPosition  = [0 0 1000 600];
figColor     = [1 1 1];
figXLabel    = 'Step difference in z-shim indices';
figYLabel    = 'Whole cord signal intensity (a.u.)';

colors =  {[1 0 0],[0 0 1],[1 0 0],[0 0 1],[1 0 0],[0 0 1],[1 0 0],[0 0 1],[1 0 0],[0 0 1]};

figure; hold on; set(gcf, 'color', figColor , 'position', figPosition); grid on; box off

Y = [nanmean([No_0,Manual_0]);nanmean([No_1,Manual_1]);...
    nanmean([No_2,Manual_2]); nanmean([No_3,Manual_3]); ...
    nanmean([No_above3,Manual_above3])];
bar_handle = bar([1;2;3;4;5],Y, 'grouped','linewidth', 2 );
set(bar_handle(1),'FaceColor',[1,1,1],'EdgeColor', [1 0 0]);
set(bar_handle(2),'FaceColor',[1,1,1],'EdgeColor', [0 0 1]);
datasets = {No_0,Manual_0,No_1,Manual_1,...
    No_2,Manual_2,No_3,Manual_3, ...
    No_above3,Manual_above3};
positions = {1,2,3,4,5};
counter   = 0;

for j = 1:numel(datasets)
    
    sem= nanstd(datasets{j})/sqrt(size(datasets{j},1));
    
    if ismember(j,[1 2])
        posit = 1;
    elseif ismember(j,[3 4])
        posit = 2;
    elseif ismember(j,[5 6])
        posit = 3;
    elseif ismember(j,[7 8])
        posit = 4;
    elseif ismember(j,[9 10])
        posit = 5;
    end
    
    if mod(j,2) ==0
        pp = [posit+0.15 posit+0.15];
    else
        pp = [posit-0.15 posit-0.15];
        
    end
    
    e(j) = line(pp, ...
        [nanmean(datasets{j})-sem nanmean(datasets{j})+sem], ...
        'linewidth', 2, 'color', colors{j});
end


for p = [1:2:10]
    counter   = counter+1;
    for j = 1:length(datasets{1})
        line([positions{counter}-0.09 positions{counter}+0.09], [datasets{p}(j) datasets{p+1}(j)], ...
            'color', [0.5 0.5 0.5 0.5], 'linewidth', 2);
    end
end

xticks([1:1:5])
set(gca,'xticklabel',{'0', '1', '2', '3', '>3'},'yticklabel',[],...
    'GridAlpha',0.5, 'xlim', figXLim, 'ylim', figYLim)
xlabel(figXLabel)
ylabel(figYLabel)

if printStatus
    
    cd(path)   
    print(gcf, '-painters','-dsvg', ...
        [ path filesep 'Figure_2B.svg']);
end
%% C

clc
clearvars -except printStatus datapath figurepath

load(fullfile(datapath, 'signal_templatespace', 'GroupWhole', 'GM', 'SingleVolume_MeanSignal', 'Data.mat'))


figXLim      = [0.5 5.5];
figYLim      = [60 130];
figPosition  = [0 0 600 600];
figColor     = [1 1 1];
figXLabel    = [];
figYLabel    = 'Gray matter  signal intensity (a.u.)';

positions1 = [1 4];
positions2 = [2 5];
colors1 = {[1 0 0] [0 0 1]};

figure; hold on; set(gcf, 'color', figColor , 'position', figPosition); box off

No_D = (noLEFT_DORSAL+noRIGHT_DORSAL)./2;
Manual_D = (manualLEFT_DORSAL+manualRIGHT_DORSAL)./2;
No_V = (noLEFT_VENTRAL+noRIGHT_VENTRAL)./2;
Manual_V = (manualLEFT_VENTRAL+manualRIGHT_VENTRAL)./2;

no_D = mean(No_D,2);
manual_D = mean(Manual_D,2);

no_V = mean(No_V,2);
manual_V = mean(Manual_V,2);

data1 = [no_D, no_V];
data2 = [manual_D, manual_V];

for c = 1:2
    
    h=boxplot(data1(:,c), 'positions',positions1(c), 'Colors',colors1{1} ,'Widths', [0.1] ,'Symbol','.r');
    g=boxplot(data2(:,c), 'positions', positions2(c) , 'Colors', colors1{2},'Widths', [0.1] ,'Symbol','.r');
    set(h,'LineWidth',2 )
    set(g,'LineWidth', 2)
    
    clear h g
    
end

plot([positions1(1) positions2(1)], mean(data1), '-','LineWidth',3, 'color',[0.5 0.5 0.5])
plot([positions1(2) positions2(2)], mean(data2), '-','LineWidth',3, 'color',[0.5 0.5 0.5])

distributionPlot(data1,'distWidth',0.1,'showMM',0, 'color',colors1{1},  'widthDiv',[2 1], 'histOri','left', 'xValues',  positions1- 0.2)
distributionPlot(data2,'distWidth',0.1,'showMM',0, 'color', colors1{2},  'widthDiv',[2 2],'histOri','right', 'xValues',  positions2 +0.2)

plotSpread(colors1{1}, data1, 'spreadWidth', 0.2, 'binWidth', 0.2, 'xValues', positions1+ 0.2);
plotSpread(colors1{2}, data2, 'spreadWidth', 0.2, 'binWidth', 0.2, 'xValues',  positions2- 0.2);

xticks([(positions1(1)+positions2(1))/2 (positions1(2)+positions2(2))/2])
set(gca,'xticklabel',{'dorsal', 'ventral'},'yticklabel',[],...
    'GridAlpha',0.5, 'xlim', figXLim, 'ylim', figYLim)
xlabel(figXLabel)
ylabel(figYLabel)

if printStatus
   
    cd(path)
    
    print(gcf, '-painters','-dsvg', ...
        [ path filesep 'Figure_2C.svg']);
end

%% D

clc
clearvars -except printStatus datapath figurepath

load(fullfile(datapath, 'signal_templatespace', 'GroupWhole', 'GM', 'TimeSeries_TSNR', 'Data.mat'))

figXLim      = [0.5 2.5];
figYLim      = [10 20];
figPosition  = [0 0 300 600];
figColor     = [1 1 1];
figXLabel    = [];
figYLabel    = 'Mean gray matter tSNR';

positions1   = 1;
positions2   = 2;
colors = {[1 0 0] [0 0 1] [0.8 0.8 0.8]};

figure; hold on; set(gcf, 'color', figColor , 'position', figPosition); box off

zNo  =(noLEFT_DORSAL+noLEFT_VENTRAL+noRIGHT_DORSAL+noRIGHT_VENTRAL)./4;
zManual  =(manualLEFT_DORSAL+manualLEFT_VENTRAL+manualRIGHT_DORSAL+manualRIGHT_VENTRAL)./4;

zshimNo     = zNo;
zshimManual = zManual;

noMean = mean(zshimNo,2);
manualMean = mean(zshimManual,2);

data1 = [noMean];
data2 = [manualMean];
datasets = {noMean, manualMean};


h=boxplot(data1(:,1), 'positions',positions1(1), 'Colors',colors{1} ,'Widths', [0.1] ,'Symbol','.r');
g=boxplot(data2(:,1), 'positions', positions2(1) , 'Colors', colors{2},'Widths', [0.1] ,'Symbol','.r');
set(h,'LineWidth',2 )
set(g,'LineWidth', 2)

distributionPlot(data1,'distWidth',0.3,'showMM',0, 'color', colors{1}, 'widthDiv',[2 1], 'histOri','left', 'xValues', positions1-0.2)
distributionPlot(data2,'distWidth',0.3,'showMM',0, 'color', colors{2}, 'widthDiv',[2 2], 'histOri','right', 'xValues', positions2 + 0.2)

positions = {positions1(1)+0.07 positions2(1)-0.07};
counter   = 0;

for p = [1]
    for j = 1:length(datasets{1})
        line([positions{p} positions{p+1}], [datasets{p}(j) datasets{p+1}(j)], ...
            'color', colors{3}, 'linewidth', 1);
    end
end


xticks([positions1(1) positions2])
set(gca,'xticklabel',[],...
    'GridAlpha',0.5, 'xlim', figXLim, 'ylim', figYLim)
xlabel(figXLabel)
ylabel(figYLabel)

if printStatus
    
    cd(path)
    
    print(gcf, '-painters','-dsvg', ...
        [ path filesep 'Figure_2D.svg']);
end


