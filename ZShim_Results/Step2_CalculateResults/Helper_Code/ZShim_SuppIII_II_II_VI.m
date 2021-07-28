function results = ZShim_SuppIII_II_II_VI(rawdatapath,processdatapath,recalculateResults)

if recalculateResults
    cd(rawdatapath)
    subjects = dir('*sub-*'); % get the list of the subjects
    
    %remember that the first 3 subjects do not have the second field map
    %acquisition!
    subjects = subjects(3:end);
    
    for sub = 1:size(subjects,1)
        
        cd(fullfile(processdatapath,subjects(sub).name,'fmap','investigation'));
        
        %rename the phase image
        copyfile([rawdatapath filesep subjects(sub).name filesep 'fmap' filesep '*run-02_fieldmap*.nii.gz'])
        tmp_dir = dir('*run-02_fieldmap*.nii.gz');
        movefile(tmp_dir(1).name, 'phase_ori_cbsfl_2.nii.gz')
        clear tmp_dir
        
        %swap dimensions of the phase image
        system('sct_image -i phase_ori_cbsfl_2.nii.gz -setorient RPI -o phase_cbsfl_swap_2.nii.gz');
        
        %smooth the phase image
        system('sct_maths -i phase_ori_cbsfl_2.nii.gz -smooth 1,1,1 -o phase_cbsfl_2_1.nii.gz');
        
        % register T2 and cbsfl field map
        system('sct_register_multimodal -i T2.nii.gz -d phase_cbsfl_2_1.nii.gz -identity 1 -o RegT2_cbsfl.nii.gz');
        system('sct_apply_transfo -i T2_seg.nii.gz -d phase_cbsfl_2_1.nii.gz -w warp_T22phase_cbsfl_2_1.nii.gz -o T2_mask_cbsfl2.nii.gz -x nn');
        
        %as the registred T2_mask is already there, it is not necessary to
        %register it again
        fm_type = 'cbsfl_2';
        ZShims_Calculate_FMbasedpicks(fm_type,'cbsfl2', {8, 0, 8, '1'});
        
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
    
    load([processdatapath filesep 'extracted_signal' filesep 'signal_templatespace' filesep 'GroupWhole' ...
        filesep 'ReconstructedSignal' filesep 'FM_relofpicks' filesep 'results.mat'])
    
end

results{1,1} = 'Reliability of FM-based picks [mean, min, max correlation]';
results{1,2} = [mean(subjectCorr) min(subjectCorr) max(subjectCorr)];


end