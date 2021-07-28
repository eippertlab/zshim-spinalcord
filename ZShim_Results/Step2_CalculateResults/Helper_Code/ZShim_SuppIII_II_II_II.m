function results = ZShim_SuppIII_II_II_II(scttemplatepath,rawdatapath,processdatapath,recalculateResults)
% Use different parameters in field map processing and fitting procedure:
%   i)  Set smoothing kernel to 0, 1 or 2 mm (previously 1mm)
%   ii) Set slab thickness to 5, 9 or 13 mm (previously 9mm)
%   iii) Weighting of voxel in fitting procedure equally, or weighted by a
%        raised cosine kernel of FWHM=slab thickness and roll-off beta = 0.5
%        (previously equal weighting)
% Calculate new zshims, create artificial volumes based on these modified
% zshims and extract signal.
% Compare modified implementation signal to previous one and return the result.
% If recalculateResults is set to False, no new data will be created, but
% rather the pre-saved results will be loaded.

if recalculateResults
    cd(rawdatapath)
    subjects = dir('*sub-*');
    
    % prepare the parameter combinations
    % fitSli, beta, totalSli, smoothing
    combs(1,:) = {8, 0, 8, '0'};
    combs(2,:) = {6, 1/2, 8, '0'};
    combs(3,:) = {4, 0, 4, '0'};
    combs(4,:) = {3, 1/2, 4, '0'};
    combs(5,:) = {12, 0, 12, '0'};
    combs(6,:) = {9, 1/2, 12, '0'};
    
    % Smoothing1_NoMaskWeighting
    combs(7,:) = {6, 1/2, 8,  '1'};
    combs(8,:) = {4, 0, 4,  '1'};
    combs(9,:) = {3, 1/2, 4, '1'};
    combs(10,:) = {12, 0, 12 , '1'};
    combs(11,:) = {9, 1/2, 12, '1'};
    
    % Smoothing_2_NoMaskWeighting
    combs(12,:) = {8, 0, 8,  '2'};
    combs(13,:) = {6, 1/2, 8,  '2'};
    combs(14,:) = {4, 0, 4,  '2'};
    combs(15,:) = {3, 1/2, 4 '2'};
    combs(16,:) = {12, 0, 12, '2'};
    combs(17,:) = {9, 1/2, 12, '2'};
    %%
    % smooth the phase map with 2mm kernel
    for sub = 1:size(subjects,1)
        
        cd(fullfile(processdatapath,subjects(sub).name,'fmap','investigation'));
        system('sct_maths -i phase_swap.nii.gz -smooth 2,2,2 -o phase2.nii.gz');
        
    end
    
    %% calculate FM-based picks for different combination of parameters
    for cbs = 1:size(combs,1)
        
        for sub = 1:size(subjects,1)
            
            cd(fullfile(processdatapath,subjects(sub).name,'fmap','investigation'));
            
            fm_type = 'gre';
            ZShims_Calculate_FMbasedpicks(fm_type,[],combs(cbs,:))
            
        end
        
    end
    %% Reconstruct volumes based on these z-shims
    reconsnames = dir('*smooth*.mat*'); % get the list of all
    
    for sub = 1:size(subjects,1)
        
        for r = 1:size(reconsnames,1)
            
            outdir = fullfile(processdatapath,subjects(sub).name,'func','reconstructZRef',filesep);
            cd(fullfile(processdatapath,subjects(sub).name,'fmap','investigation'));
            
            % copy the files to the reconstruction of z-shim values folder
            copyfile(reconsnames(r).name, outdir)
            
            cd(outdir)
            
            ZrefName = 'ZRef1.nii.gz';
            
            ZShim_Fitting_ReconstructZShimSeries(dirData, ZrefName , 24, ...
                [reconsnames(r).name], ['ZShim1_Recons_' reconsnames(r).name(1:end-4)]);
            
            
            system(['sct_apply_transfo '  ...
                ' -i ZShim1_Recons_' reconsnames(r).name(1:end-4) '.nii.gz'  ...
                ' -d ' processdatapath filesep subjects(sub).name filesep 'func' filesep '/moco2volumes_mean.nii.gz ' ...
                ' -w ZRef1Reg' filesep 'warp_ZRef1_MEAN2moco2volumes_mean.nii.gz ' ...
                ' -x spline ' ...
                ' -o ZShim1_Recons_'  reconsnames(r).name(1:end-4) '_MocoReg.nii.gz'])
            
            system(['sct_apply_transfo -i  ZShim1_Recons_'  reconsnames(r).name(1:end-4) '_MocoReg.nii.gz'  ...
                ' -d ' scttemplatepath filesep 'PAM50' filesep 'template' filesep  'PAM50_t2.nii.gz ' ...
                ' -w ' processdatapath filesep subjects(sub).name filesep 'func' filesep 'warps' filesep 'warp_moco2volumes_mean2PAM50_t2.nii.gz' ...
                ' -x spline' ...
                ' -o normalizedVols' filesep 'ZShim1_Recons_' reconsnames(r).name(1:end-4) '_normalized.nii.gz'])
            
            system(['fslmeants -i normalizedVols' filesep 'ZShim1_Recons_' reconsnames(r).name(1:end-4) '_normalized ' ...
                ' -m ' scttemplatepath filesep 'PAM50' filesep 'atlas' filesep  'WHOLE_GM.nii.gz ' ...
                ' -o normalizedVols' filesep 'ZShim1_Recons_'  reconsnames(r).name(1:end-4) '_normalized_GM.txt ' ...
                ' --showall'])
            
        end
        
    end
    
    % extract signal
    reconsName    = 'ZShim1_Recons';
    reconsMode    = dir(['normalizedVols' filesep '*smooth*.txt']); %get all the files
    reconsMode    = {reconsMode.name,'gre'}';
    
    recons_signal = ZShim_Load_ReconsSignal(processdatapath,subjects,reconsName,reconsMode);
    
    recons_signal_mean = mean(recons_signal,3);
    
else
    
    load([processdatapath filesep 'extracted_signal' filesep 'signal_templatespace' filesep 'GroupWhole' ...
        filesep 'ReconstructedSignal' filesep 'FM_parameters' filesep 'results.mat'])
    
    parameterChange_names{18}= 'gre';   
    reconsMode    = parameterChange_names;
    
    resultsDiffParameters(:,18,:) = gre;
    recons_signal = resultsDiffParameters;
    
    recons_signal_mean = mean(resultsDiffParameters,3);
    
end

% 18 is the original implementation
% test all modifications for improvement over original in terms of signal
% mean
for tt = 1:size(recons_signal_mean,2)-1
    
    [~,p(tt), ~,stat(tt)] = ttest(recons_signal_mean(:,tt), recons_signal_mean(:,18), ...
        'Tail', 'right'); 
    
end

idx1 = find(p==min(p)); % index of the smallest p;

% test the improvement over original implementation for best modification
results_meansignal = ZShim_CalculateResults(recons_signal_mean(:,18),recons_signal_mean(:,idx1),[],0,1); 

results_std = std(recons_signal,[],3);

results_coeff = results_std./recons_signal_mean;

clear h p stat idx

% test all modifications for improvement over original in terms of
% coefficients of variation
for tt = 1:size(recons_signal_mean,2)-1
    
    [~,p(tt), ~,stat(tt)] = ttest(results_coeff(:,tt),results_coeff(:,18), ...
        'Tail', 'left');
    
end

idx2 = find(p==min(p)); % index of the smallest p;

% test the improvement over original implementation for best modification
results_coeffvar = ZShim_CalculateResults(results_coeff(:,18),results_coeff(:,idx2),[],0,1); 

results{1,1} = 'results_meansignal';
results{1,2} = results_meansignal;
results{1,3} = ['name of the parameter set: ' reconsMode{idx1} ];

results{2,1} = 'results_coeffvar';
results{2,2} = results_coeffvar;
results{2,3} = ['name of the parameter set: ' reconsMode{idx2} ];


end

