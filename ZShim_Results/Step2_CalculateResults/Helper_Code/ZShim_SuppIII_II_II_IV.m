function ZShim_SuppIII_II_II_IV(figurepath,rawdatapath,processdatapath,recalculateResults)

cd(rawdatapath)
subjects = dir('*sub-*'); % get the list of the subjects

if recalculateResults
%%
for sub = 1:size(subjects,1)
    cd(fullfile(processdatapath,subjects(sub).name,'func','reconstructZRef',filesep));
    
    % load the manual and fm-based picks
    load('ManualSelection.mat');
    load('AutoSelection_FM_gre8.mat');
    
    % absolute difference between picks
    abs_diff(sub,:) = abs(manual - auto_gre);
end
%% find the indices of  >=3 pick differences 
[row,col] = find(abs_diff>=3);
subj_idx = unique(row);
off_subj_dire = subjects(subj_idx,:);
diffs = abs_diff(subj_idx,:);
%%
fitSli = 8; % this is what I input to  Johanna's code
beta = 0;
res = [2.2 1 1]*1e-3;
basis = [1:4];
smoothfactor = '1';

for sub = 1:size(off_subj_dire,1)
    
   cd '/data/pt_02098/ZShim_BIDS_Organized/processed_data'
    cd(['Recons1/' off_subj_dire(s).name '/'])

    if isempty(smoothfactor)
        
        fmFilename = ['phase_swap.nii.gz'];
        
    else
        
        fmFilename = ['phase'  smoothfactor      '.nii.gz'];
    end
    
    maskFilename = ['T2_mask.nii.gz'];
    
    fm = read_avw(fmFilename);
    mask = read_avw(maskFilename);
    mask = logical(mask);
    
    fm = read_avw(fmFilename);
    mask = read_avw(maskFilename);
    mask = logical(mask);
    
    % Rescale field map
    
    dTE = 2.46e-3; %difference in echo time
    HzMax = 1/(2*dTE);
    fmMax = max(abs(fm(:)));
    fm = fm*HzMax/fmMax;

    % Cut the field map 5 slices above and 4 slices below of the coverage of the EPI slice
    % In our acquisition the slice stack was placed in a way that
    % the center of the EPI volume (24 slices) corresponds to the center of the
    % fm acquisition (180 sagittal slices) --> therefore,the middle point of the
    % acquisitions are identical
    fm = fm(:,:,(fmSli/2)-(epiSli/2*5)-4:(fmSli/2)+(epiSli/2*5)+5);
    mask = mask(:,:,(fmSli/2)-(epiSli/2*5)-4:(fmSli/2)+(epiSli/2*5)+5);
    
    % Compute the shims - main code
    [shims,~,~,~,shimSli] = compute_shims(fm,mask,basis,fitSli*res(3),res,beta);
    
    Hz_scale = [1 5e-3 5e-3 5e-3]; % thickness of our EPI slice
    shims = diag(Hz_scale(basis))*shims; % convert shims to Hz    
    offsets(sub,:) = shims(1,:); % get the field offset for each subject 

 
   clear fitArray Zshim_picks fm mask Zshim_picks auto_gre
end

else
    
    load([processdatapath filesep 'extracted_signal' filesep 'signal_templatespace' filesep 'GroupWhole' ...
        filesep 'ReconstructedSignal' filesep 'FM_inhomogeneity' filesep 'results.mat'])
    
    

end

%find the middle slices -- as slices wi
vector = [NaN(1,1) repmat([1 2 3 4 5], 1, 24) NaN(1,1)]';
sli_idx =find(vector==3);

 
 figPosition2 = [0 0 1000 1000];
 figure(2); hold on; set(gcf, 'color', [1 1 1], 'position', figPosition2);  box off
 ylimValues = [-300 900];
 
 for sp = 1:size(offsets,1)
     
     subplot(5,2,sp); hold on;
     yyaxis left
     plot(offsets(sp,:), '-k', 'LineWidth', 2)
     ylim(ylimValues)
     scatter(sli_idx, offsets(sp,sli_idx), 20,'MarkerEdgeColor',[0 0 0],...
         'MarkerFaceColor',[0 0 0],...
         'LineWidth',1.5)
     yticks(-200:200:800)
     if sp ==1
         ylabel('Field offset (Hz)')
         
     end

     yyaxis right
     ylim([-3 9])
     scatter(sli_idx,diffs(sp,:)*2,30,'MarkerEdgeColor',[.5 .5 .5],...
         'MarkerFaceColor',[0.5 .5 .5],...
         'LineWidth',1.5)
     yticks([0:2:8])
     yticklabels({'0', '1', '2', '3', '4'})
     yline(0)


     ax = gca;
     ax.YAxis(1).Color = 'k';
     ax.YAxis(2).Color = 'k';
     
     xlim([2 121])
     xticks([2 121])
     xticklabels({'1', '120'})
     
     if sp ==1
         ylabel('Absolute difference')
         xlabel('Slices (Field-map resolution: 1mm)')
         legend({'Field offset', 'Center of EPI slice', 'Step difference'})
     end     
  
 end
 
 %% save the figure
 
 if ~exist(figurepath)
     mkdir(figurepath)
 end
 
 cd(figurepath)
 print(gcf, '-painters','-dsvg', ...
     ['suppFigure_IV_inhomogenity.svg']);

end

