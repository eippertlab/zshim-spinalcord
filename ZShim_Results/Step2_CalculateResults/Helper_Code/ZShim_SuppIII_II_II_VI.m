function results = ZShim_SuppIII_II_II_VI(rawdatapath,processdatapath,recalculateResults,fslDir)
% Compare the z-shims obtained by the in-house field map acquired at the beginning of the scan session
% with the z-shim obtained by the in-house field map acquired at the end of the scan session 
% use Spearman rank-correlation toevaluate how similar these values are.
% If recalculateResults is set to False, no new data will be created, but
% rather the pre-saved results will be loaded.

% ----------
% Inputs:
% ----------

% scttemplatepath:   fullpath, string, PAM50 template location
% rawdatapath:       fullpath, string, raw data location
% processeddatapath: fullpath, string, processed data location
% recalculateResults: True or False

if recalculateResults
    cd(rawdatapath)
    subjects = dir('*sub-*'); % get the list of the subjects
    
    % remember that the first 3 subjects do not have the second field map
    % acquisition!
    subjects = subjects(4:end);
    
    for sub = 1:size(subjects,1)
        
        cd(fullfile(processdatapath,subjects(sub).name,'fmap','investigation'));
        
        % rename the phase image
        copyfile(fullfile(rawdatapath, subjects(sub).name, 'fmap', '*run-02_fieldmap*.nii.gz'))
        tmp_dir = dir('*run-02_fieldmap*.nii.gz');
        movefile(tmp_dir(1).name, 'phase_ori_cbsfl_2.nii.gz')
        clear tmp_dir
        
        % swap dimensions of the phase image
        system('sct_image -i phase_ori_cbsfl_2.nii.gz -setorient RPI -o phase_cbsfl_swap_2.nii.gz');
        
        % smooth the phase image
        system('sct_maths -i phase_ori_cbsfl_2.nii.gz -smooth 1,1,1 -o phase_cbsfl_2_1.nii.gz');
        
        % register T2 and cbsfl field map
        system('sct_register_multimodal -i T2.nii.gz -d phase_cbsfl_2_1.nii.gz -identity 1 -o RegT2_cbsfl.nii.gz');
        system('sct_apply_transfo -i T2_seg.nii.gz -d phase_cbsfl_2_1.nii.gz -w warp_T22phase_cbsfl_2_1.nii.gz -o T2_mask_cbsfl2.nii.gz -x nn');
        
        % as the registered T2_mask is already there, it is not necessary to
        % register it again
        fm_type = 'cbsfl_2';
        ZShim_Calculate_FMbasedpicks(fm_type,'cbsfl2', {8, 0, 8, '1'},fslDir);
        
        load('AutoSelection_FM_cbsfl_1.mat');
        cbsfl1_picks = auto_cbsfl;
        clear auto_cbsfl
        
        load('AutoSelection_FM_cbsfl_2.mat');
        cbsfl2_picks = auto_cbsfl;
        clear auto_cbsfl
        
        subjectCorr(s,1) = corr(cbsfl1_picks,cbsfl2_picks, 'Type', 'Spearman');
        
        cbsfl1_picks = [];
        cbsfl2_picks = [];
        
    end
    
else
    
    load(fullfile(processdatapath, 'extracted_signal', 'signal_templatespace', 'GroupWhole', ...
        'ReconstructedSignal', 'FM_relofpicks', 'results.mat'))
    
end

results{1,1} = 'Reliability of FM-based picks [mean, min, max correlation]';
results{1,2} = [mean(subjectCorr) min(subjectCorr) max(subjectCorr)];


end