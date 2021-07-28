function results = ZShim_SuppIII_II_II_I(scttemplatepath,rawdatapath,processdatapath,recalculateResults)
% Use a modified mask for voxel selection: 
% either eroded by 1 voxel or dilated by 1 voxel.
% Calculate new zshims, create artificial volumes based on these modified
% zshims and extract signal.
% Compare modified implementation signal to previous one and return the result.
% If recalculateResults is set to False, no new data will be created, but
% rather the pre-saved results will be loaded.

if recalculateResults
    cd(rawdatapath)
    subjects = dir('*sub-*'); % get the list of the subjects
    
    for sub = 1:size(subjects,1)
        
        outdir = fullfile(processdatapath,subjects(sub).name,'fmap','investigation');
        mkdir(outdir)
        cd(outdir)
        
        copyfile(fullfile(rawdatapath, subjects(sub).name, 'anat', '*T2w*.nii.gz'))
        % rename to T2 for further processing
        tmp_dir = dir('*T2w*.nii.gz');
        movefile(tmp_dir(1).name, 'T2.nii.gz')
        clear tmp_dir
        
        % segment T2
        system('sct_propseg -i T2.nii.gz  -c t2 ');
        
        % copy the phase map
        copyfile(fullfile(rawdatapath, subjects(sub).name, 'fmap', '*phasediff*.nii.gz'))
        % rename the phase map for further processing
        tmp_dir = dir('*phasediff*.nii.gz');
        movefile(tmp_dir(1).name, 'phase_ori.nii.gz')
        clear tmp_dir
        
        % swap dimensions of the phase image
        system('fslswapdim phase_ori RL PA IS  phase_swap.nii.gz');
        
        % resample the T2 segmentation to the template space
        system('sct_register_multimodal -i T2.nii.gz -d phase_swap.nii.gz -identity 1 -o RegT2.nii.gz');
        system('sct_apply_transfo -i T2_seg.nii.gz -d phase_swap.nii.gz -w warp_T22phase_swap.nii.gz -o T2_mask.nii.gz -x nn');
        
        % smooth the phase image
        system('sct_maths -i phase_swap.nii.gz -smooth 1,1,1 -o phase1.nii.gz');
        
        % erode the segmentation
        system('sct_maths -i T2_seg.nii.gz -erode 1 -o T2_seg_eroded.nii.gz');
        
        % dilate the segmentation
        system('sct_maths -i T2_seg.nii.gz -dilate 1 -o T2_seg_dilated.nii.gz');
        
        % resample these segmentations to the phase map
        system(['sct_apply_transfo -i T2_seg_dilated.nii.gz -d '...
            ' phase_swap.nii.gz -w warp_T22phase_swap.nii.gz -o T2_mask_dilated.nii.gz -x nn']);
        
        system(['sct_apply_transfo -i T2_seg_eroded.nii.gz -d '...
            ' phase_swap.nii.gz -w warp_T22phase_swap.nii.gz -o T2_mask_eroded.nii.gz -x nn']);
        
        % calculate the new z-shim selections based on these eroded and masks
        % type of the field map, mask, name of the field map
        ZShim_Calculate_FMbasedpicks('gre', 'dilated',{8, 0, 8, '1'}); 
        ZShim_Calculate_FMbasedpicks('gre', 'eroded' ,{8, 0, 8, '1'});
        
        
    end
    %% Register Z-reference EPI image and mean of motion-corrected volumes
    
    for s = 1:size(dire,1)
        
        cd(fullfile(processdatapath,subjects(sub).name,'fmap','investigation'));
        outdir = fullfile(processdatapath,subjects(sub).name,'func','reconstructZRef',filesep);
        
        % copy the files to the reconstruction of z-shim values folder
        copyfile('*dilated.mat', outdir)
        copyfile('*eroded.mat', outdir)
        
        cd(outdir)
        
        ZrefName = 'ZRef1.nii.gz';
        
        % reconstruct the artificial volumes based on eroded and dilated masks
        ZShim_Fitting_ReconstructZShimSeries(outdir, ZrefName , 24, ...
            'AutoSelection_FM_gre8_eroded.mat' , 'ZShim1_Recons_gre_eroded')
        
        ZShim_Fitting_ReconstructZShimSeries(outdir, ZrefName , 24, ...
            'AutoSelection_FM_gre8_dilated.mat' , 'ZShim1_Recons_gre_dilated')
        
        % to normalize these to template space
        % 1. register to moco image (as moco mask will be used)
        system(['sct_apply_transfo ' ...
            ' -i ZShim1_Recons_gre_eroded.nii.gz '   ...
            ' -d ' processdatapath filesep subjects(sub).name filesep 'func' filesep  'moco2volumes_mean.nii.gz ' ...
            ' -w ZRef1Reg' filesep 'warp_ZRef1_MEAN2moco2volumes_mean.nii.gz  ' ...
            ' -x spline ' ...
            ' -o gre_eroded_MocoReg.nii.gz']);
        
        system(['sct_apply_transfo ' ...
            ' -i ZShim1_Recons_gre_dilated.nii.gz '   ...
            ' -d ' processdatapath filesep subjects(sub).name filesep 'func' filesep  'moco2volumes_mean.nii.gz ' ...
            ' -w ZRef1Reg' filesep 'warp_ZRef1_MEAN2moco2volumes_mean.nii.gz  ' ...
            ' -x spline ' ...
            ' -o gre_dilated_MocoReg.nii.gz']);
        
        % 2. apply the warping field (moco volume to PAM50 --> to normalize to the template space)
        system(['sct_apply_transfo -i gre_eroded_MocoReg.nii.gz '  ...
            ' -d ' scttemplatepath filesep 'PAM50' filesep 'template' filesep  'PAM50_t2.nii.gz ' ...
            ' -w ' processdatapath filesep subjects(sub).name filesep 'func' filesep 'warps' filesep 'warp_moco2volumes_mean2PAM50_t2.nii.gz' ...
            ' -x spline ' ...
            ' -o normalizedVols' filesep 'gre_eroded_normalized.nii.gz'])
        
        system(['sct_apply_transfo -i gre_dilated_MocoReg.nii.gz '  ...
            ' -d ' scttemplatepath filesep 'PAM50' filesep 'template' filesep  'PAM50_t2.nii.gz ' ...
            ' -w ' processdatapath filesep subjects(sub).name filesep 'func' filesep 'warps' filesep 'warp_moco2volumes_mean2PAM50_t2.nii.gz' ...
            ' -x spline' ...
            ' -o normalizedVols' filesep 'gre_dilated_normalized.nii.gz'])
        
        % 3. extract the signal in the gray matter
        system(['fslmeants -i normalizedVols' filesep 'gre_dilated_normalized ' ...
            ' -m ' scttemplatepath filesep 'PAM50' filesep 'atlas' filesep  'WHOLE_GM.nii.gz ' ...
            ' -o normalizedVols' filesep 'ZShim1_Recons_gre_dilated_normalized_GM.txt ' ...
            ' --showall']);
        
        system(['fslmeants -i normalizedVols' filesep 'gre_eroded_normalized ' ...
            ' -m ' scttemplatepath filesep 'PAM50' filesep 'atlas' filesep  'WHOLE_GM.nii.gz ' ...
            ' -o normalizedVols' filesep 'ZShim1_Recons_gre_eroded_normalized_GM.txt ' ...
            ' --showall']);
        
    end
    
    % extract signal
    reconsName  = 'ZShim1_Recons';
    reconsMode    = {'gre', 'gre_eroded', 'gre_dilated', 'no'};
    recons_signal = ZShim_Load_ReconsSignal(processdatapath,subjects,reconsName,reconsMode);
    
    for r = 1:numel(reconsMode)
        
        eval([reconsMode{r}  '= recons_signal(' num2str(r) ');'])
        
    end
    
else
    load(fullfile(processdatapath, 'extracted_signal', 'signal_templatespace', 'GroupWhole', ...
        'ReconstructedSignal', 'FM_maskchange', 'results.mat'))
end

% compare eroded and dilated version to the previous implementation
results_eroded  =  ZShim_CalculateResults(no,gre,gre_eroded,0,1);
results_dilated =  ZShim_CalculateResults(no,gre,gre_dilated,0,1);

results{1,1} = 'results_eroded';
results{1,2} = results_eroded;

results{2,1} = 'results_dilated';
results{2,2} = results_dilated;

end