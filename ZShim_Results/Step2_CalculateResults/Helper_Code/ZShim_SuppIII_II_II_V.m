function results = ZShim_SuppIII_II_II_V(scttemplatepath,rawdatapath,processdatapath,recalculateResults,fslDir)
% Use in-house field map instead of vendor-provided field map.
% Calculate new zshims, create artificial volumes based on these modified
% zshims and extract signal.
% Compare modified implementation signal to previous one and return the result.
% If recalculateResults is set to False, no new data will be created, but
% rather the pre-saved results will be loaded.

% ----------
% Inputs:
% ----------

% scttemplatepath:   fullpath, string, PAM50 template location
% rawdatapath:       fullpath, string, raw data location
% processeddatapath: fullpath, string, processed data location
% recalculateResults: True or False

% Merve Kaptan, mkaptan@cbs.mpg.de

if recalculateResults
    
    cd(rawdatapath)
    subjects = dir('*sub-*'); % get the list of the subjects
    
    for sub = 1:size(subjects,1)
        
        cd(fullfile(processdatapath,subjects(sub).name,'fmap','investigation'));
        
        copyfile([rawdatapath filesep subjects(sub).name filesep 'fmap' filesep '*run-01_fieldmap*.nii.gz'])
        tmp_dir = dir('*fieldmap*.nii.gz');
        movefile(tmp_dir(1).name, 'phase_ori_cbsfl_1.nii.gz')
        clear tmp_dir
        
        % swap dimensions of the phase image
        system('sct_image -i phase_ori_cbsfl_1.nii.gz -setorient RPI -o phase_cbsfl_swap.nii.gz');
        
        % smooth the phase image
        system('sct_maths -i phase_ori_cbsfl_1.nii.gz -smooth 1,1,1 -o phase_cbsfl_1_1.nii.gz');
        
        % to use the segmentation of T2, we have to register it
        system('sct_register_multimodal -i T2.nii.gz -d phase_cbsfl_1_1.nii.gz -identity 1 -o RegT2_cbsfl.nii.gz');
        system('sct_apply_transfo -i T2_seg.nii.gz -d phase_cbsfl_1_1.nii.gz -w warp_T22phase_cbsfl_1_1.nii.gz -o T2_mask_cbsfl.nii.gz -x nn');
        
        fm_type = 'cbsfl_1';
        ZShim_Calculate_FMbasedpicks(fm_type,'cbsfl', {8, 0, 8, '1'},fslDir);
        
        outdir = fullfile(processdatapath,subjects(sub).name,'func','reconstructZRef',filesep);
        
        % copy the files to the reconstruction of z-shim values folder
        copyfile('AutoSelection_FM_cbsfl_1.mat', outdir)
        
        cd(outdir)
        ZrefName = 'ZRef1.nii.gz';
        
        % reconstruct the artificial volumes based on in-house (cbsfl) fieldmap
        ZShim_Fitting_ReconstructZShimSeries(outdir, ZrefName , 24, ...
            'AutoSelection_FM_cbsfl_1.mat' , 'ZShim1_Recons_cbsfl_1')
        
        % to normalize these to template space
        % 1. register to moco image (as moco mask will be used)
        system(['sct_apply_transfo ' ...
            ' -i ZShim1_Recons_cbsfl_1.nii.gz '   ...
            ' -d ' processdatapath filesep subjects(sub).name filesep 'func' filesep  'moco2volumes_mean.nii.gz ' ...
            ' -w ZRef1Reg' filesep 'warp_ZRef1_MEAN2moco2volumes_mean.nii.gz  ' ...
            ' -x spline ' ...
            ' -o ZShim1_Recons_cbsfl_1_MocoReg.nii.gz']);
        
        % 2. apply the warping field (moco volume to PAM50 --> to normalize to the template space)
        system(['sct_apply_transfo -i ZShim1_Recons_cbsfl_1_MocoReg.nii.gz '  ...
            ' -d ' scttemplatepath filesep 'PAM50' filesep 'template' filesep  'PAM50_t2.nii.gz ' ...
            ' -w ' processdatapath filesep subjects(sub).name filesep 'func' filesep 'warps' filesep 'warp_moco2volumes_mean2PAM50_t2.nii.gz' ...
            ' -x spline ' ...
            ' -o normalizedVols' filesep 'ZShim1_Recons_cbsfl_1_normalized.nii.gz'])
        
        % 3. extract the signal in the gray matter
        system(['fslmeants -i normalizedVols' filesep 'ZShim1_Recons_cbsfl_1_normalized ' ...
            ' -m ' scttemplatepath filesep 'PAM50' filesep 'atlas' filesep  'WHOLE_GM.nii.gz ' ...
            ' -o normalizedVols' filesep 'ZShim1_Recons_cbsfl_1_normalized_GM.txt ' ...
            ' --showall']);
        
    end
    
    reconsName  = 'ZShim1_Recons';
    reconsMode    = {'gre', 'cbsfl_1', 'manual', 'no'};
    recons_signal = ZShim_Load_ReconsSignal(processdatapath,subjects,reconsName,reconsMode);
    
    for r = 1:numel(reconsMode)
        
        eval([reconsMode{r}  '= recons_signal(' num2str(r) ');'])
        
    end
else
    
    load(fullfile(processdatapath, 'extracted_signal', 'signal_templatespace', 'GroupWhole', ...
        'ReconstructedSignal', 'FM_fmchange', 'results.mat'))
    
end
results_fmcomparison1 =  ZShim_CalculateResults(no,gre,cbsfl_1,0,1);

results{1,1} = 'results_fmcomparison-against no';
results{1,2} = results_fmcomparison1;


results_fmcomparison2 =  ZShim_CalculateResults(no,gre,manual,0,1);

results{2,1} = 'results_fmcomparison-gre vs manual';
results{2,2} = results_fmcomparison2;

results_fmcomparison3 =  ZShim_CalculateResults(no,cbsfl_1,manual,0,1);

results{3,1} = 'results_fmcomparison-cbsfl vs manual';
results{3,2} = results_fmcomparison3;


end