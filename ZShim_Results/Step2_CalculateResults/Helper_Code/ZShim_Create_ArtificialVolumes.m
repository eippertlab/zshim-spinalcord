function ZShim_Create_ArtificialVolumes(subGroup, scttemplatepath,rawdatapath,processdatapath,recalculateResults)
% Creates artificial volumes for the original approaches to be able to
% later compare the modified FM-based approaches (and the resulting
% artifical volumes) with the original one
% The output is stored in the func folder of each subject in a subfolder
% called 'reconstructZRef'.
% 
% ----------
% Inputs:
% ----------
% subGroup:          string, 'epi' for EPI-based automation, 'fm' for FM-based
%                    automation
% scttemplatepath:   fullpath, string, PAM50 template location
% rawdatapath:       fullpath, string, raw data location
% processeddatapath: fullpath, string, processed data location

if recalculateResults
    cd(rawdatapath)
    subjects = dir('*sub-*'); % get the list of the subjects
    
    for sub = 1:size(subjects,1)
        
        outdir = fullfile(processdatapath,subjects(sub).name,'func','reconstructZRef');
        mkdir(outdir)
        cd(outdir)
        
        % as the middle index is 11 out of 21, this corresponds to no z-shim
        no = (ones(24,1))*11 ;
        save('NoShim.mat','no')
        
        % copy and rename the z-ref EPI measurements
        copyfile(fullfile(rawdatapath, subjects(sub).name, 'func', '*_task-rest_acq-zref_run-01_bold*.nii.gz'))
        tmp_dir = dir('*-zref_run-01_bold*.nii.gz');
        movefile(tmp_dir(1).name, 'ZRef1.nii.gz')
        clear tmp_dir
        
        copyfile(fullfile(rawdatapath, subjects(sub).name,  'func', '*_task-rest_acq-zref_run-02_bold*.nii.gz'))
        tmp_dir = dir('*-zref_run-02_bold*.nii.gz');
        movefile(tmp_dir(1).name, 'ZRef2.nii.gz')
        clear tmp_dir
        
        if isequal(subGroup, 'epi')
            subIdx = [3:24 47 48];
        elseif isequal(subGroup, 'fm')
            subIdx = [1 2 25:46];
        end
        
        % for EPI-based group, calculate FM-based picks
        if isequal(subGroup, 'epi')
            
            ManualPicks = readtable([processdatapath filesep 'manualzshimPicks_duringScan.csv']);
            AutoPicks   = readtable([processdatapath filesep 'autozshimPicks_duringScan.csv']);
            
            ManualPicks = table2array(ManualPicks(subIdx,3:end));
            AutoPicks   = table2array(AutoPicks(subIdx,3:end));
            
            manual = ManualPicks(sub,:)';
            save('ManualSelection.mat', 'manual');
            
            auto_epi = AutoPicks(sub,:)';
            save('AutoSelection_EPI.mat', 'auto_epi');
            
            copyfile([rawdatapath  subjects(sub).name filesep 'anat' filesep '*T2w*.nii.gz'])
            % rename to T2 for further processing
            tmp_dir = dir('*T2w*.nii.gz');
            movefile(tmp_dir(1).name, 'T2.nii.gz')
            clear tmp_dir
            
            % segment T2
            system('sct_propseg -i T2.nii.gz  -c t2 ');
            
            % copy the phase map and rename it
            copyfile([rawdatapath  subjects(sub).name filesep 'fmap' filesep '*phasediff*.nii.gz'])
            % rename the phase map for further processing
            tmp_dir = dir('*phasediff*.nii.gz');
            movefile(tmp_dir(1).name, 'phase_ori.nii.gz')
            clear tmp_dir
            
            % swap dimensions of the phase image
            system('fslswapdim phase_ori RL PA IS  phase_swap.nii.gz');
            
            % register/resample the T2 segmentation to the template space
            system('sct_register_multimodal -i T2.nii.gz -d phase_swap.nii.gz -identity 1 -o RegT2.nii.gz');
            system('sct_apply_transfo -i T2_seg.nii.gz -d phase_swap.nii.gz -w warp_T22phase_swap.nii.gz -o T2_mask.nii.gz -x nn');
            
            % smooth the phase image
            system('sct_maths -i phase_swap.nii.gz -smooth 1,1,1 -o phase1.nii.gz');
            ZShim_Calculate_FMbasedpicks('gre', '', {8, 0, 8, '1'}) % type of the field map, mask, name of the
            
        % for FM-based group, calculate EPI-based picks
        elseif isequal(subGroup, 'fm') 
            
            ManualPicks = readtable([processdatapath filesep 'manualzshimPicks_duringScan.csv']);
            AutoPicks   = readtable([processdatapath filesep 'autozshimPicks_duringScan.csv']);
            
            ManualPicks = table2array(ManualPicks(subIdx,3:end));
            AutoPicks   = table2array(AutoPicks(subIdx,3:end));
            
            manual = ManualPicks(sub,:)';
            save('ManualSelection.mat', 'manual');
            
            auto_gre = AutoPicks(sub,:)';
            auto_gre(auto_gre<=0) = 1;  % fit can result in negative values,
                                        % but the minimum index for the sequence is 1
            save('AutoSelection_FM_gre8.mat', 'auto_gre');
            
            auto_epi = ZShims_Calculate_EPIbasedpicks('ZRef1');
            save('AutoSelection_EPI.mat', 'auto_epi');
            
        end
    end
    
    %%
    for sub = 1:size(subjects,1)
        
        outdir = fullfile(processdatapath,subjects(sub).name,'func','reconstructZRef');
        dirData = outdir;
        cd(outdir)
        
        ZrefName = 'ZRef1.nii.gz';
        
        ZShim_Fitting_ReconstructZShimSeries(dirData, ZrefName , 24, ...
            'AutoSelection_FM_gre8.mat' , 'ZShim1_Recons_gre')
        
        ZShim_Fitting_ReconstructZShimSeries(dirData, ZrefName, 24, ...
            'ManualSelection.mat', 'ZShim1_Recons_manual')
        
        ZShim_Fitting_ReconstructZShimSeries(dirData, ZrefName, 24, ...
            'AutoSelection_EPI.mat', 'ZShim1_Recons_zref')
        
        ZShim_Fitting_ReconstructZShimSeries(dirData, ZrefName, 24, ...
            'NoShim.mat', 'ZShim1_Recons_no')
        
        ZrefName = 'ZRef2.nii.gz';
        
        ZShim_Fitting_ReconstructZShimSeries(dirData, ZrefName , 24, ...
            'AutoSelection_FM_gre8.mat' , 'ZShim2_Recons_gre')
        
        ZShim_Fitting_ReconstructZShimSeries(dirData, ZrefName, 24, ...
            'ManualSelection.mat', 'ZShim2_Recons_manual')
        
        ZShim_Fitting_ReconstructZShimSeries(dirData, ZrefName, 24, ...
            'AutoSelection_EPI.mat', 'ZShim2_Recons_zref')
        
        ZShim_Fitting_ReconstructZShimSeries(dirData, ZrefName, 24, ...
            'NoShim.mat', 'ZShim2_Recons_no')
        
        % motion-correction
        system(['sct_register_multimodal -i ZRef1_MEAN.nii.gz ' ...
            ' -d ' processdatapath filesep subjects(sub).name filesep 'func' filesep 'moco2volumes_mean.nii.gz ' ...
            ' -param step=1,type=im,algo=slicereg,poly=2,deformation=0x0x0,gradStep=1 ' ...
            '  -x spline ' ...
            ' -m ' processdatapath filesep subjects(sub).name filesep 'func' filesep 'moco_mask.nii.gz '...
            ' -ofolder ZRef1Reg' filesep]);
        
        dir_vol = dir('*ZShim1*.nii.gz');
        
        for v = 1:size(dir_vol,1)
            
            system(['sct_apply_transfo ' ...
                ' -i ' dir_vol(v).name  ...
                ' -d ' processdatapath filesep subjects(sub).name filesep 'func' filesep 'moco2volumes_mean.nii.gz ' ...
                ' -w ZRef1Reg' filesep 'warp_ZRef1_MEAN2moco2volumes_mean.nii.gz  ' ...
                ' -x spline ' ...
                ' -o ' dir_vol(v).name(1:end-7)  '_MocoReg.nii.gz']);
        end
        
        
        dir_mocos = dir('*MocoReg*.nii.gz');
        
        mkdir('normalizedVols')
        
        for m = 1:size(dir_mocos,1)
            
            system(['sct_apply_transfo -i ' dir_mocos(m).name ...
                ' -d ' scttemplatepath filesep 'PAM50' filesep 'template' filesep 'PAM50_t2.nii.gz ' ...
                ' -w ' processdatapath subjects(sub).name  filesep 'func' filesep 'warps' filesep 'warp_moco2volumes_mean2PAM50_t2.nii.gz' ...
                ' -x spline' ...
                ' -o normalizedVols' filesep dir_vol(m).name(1:end-7) '_normalized.nii.gz'])
        end
        
        clear dir_vol
        dir_vol = dir(['normalizedVols' filesep '*_normalized*.nii.gz']);
        
        for m = 1:size(dir_vol,1)
            
            dirData = fullfile(processdatapath, subjects(sub).name, 'func', 'reconstructZRef', 'normalizedVols');
            cd(dirData)
            
            system(['fslmeants -i  ' dir_vol(m).name ...
                ' -m ' scttemplatepath filesep 'PAM50' filesep 'atlas' filesep 'WHOLE_GM.nii.gz ' ...
                ' -o '  dir_vol(m).name(1:end-7)  '_GM.txt --showall'])
        end
    end
    
    %% now same for the 2nd Reference measurement
    
    dir_vol = dir('*ZShim2*.nii.gz');
    
    for sub = 1:size(subjects,1)
        
        outdir = fullfile(processdatapath,subjects(sub).name,'func','reconstructZRef');
        cd(outdir)
        
        system('fslmaths ZRef2 -Tmean ZRef2_MEAN')
        
        system(['sct_register_multimodal -i ZRef2_MEAN.nii.gz ' ...
            ' -d ' processdatapath filesep subjects(sub).name filesep 'func' filesep 'moco2volumes_mean.nii.gz ' ...
            ' -param step=1,type=im,algo=slicereg,poly=2,deformation=0x0x0,gradStep=1 ' ...
            '  -x spline ' ...
            ' -m ' processdatapath filesep subjects(sub).name filesep 'func' filesep 'moco_mask.nii.gz '...
            ' -ofolder ZRef2Reg' filesep]);
        
        
        mkdir('ZRef2MocoReg')
        
        for v = 1:size(dir_vol,1)
            
            system(['sct_apply_transfo ' ...
                ' -i ' dir_vol(v).name  ...
                ' -d ' processdatapath filesep subjects(sub).name filesep 'func' filesep 'moco2volumes_mean.nii.gz ' ...
                ' -w ZRef2Reg' filesep 'warp_ZRef2_MEAN2moco2volumes_mean.nii.gz  ' ...
                ' -x spline ' ...
                ' -o ZRef2MocoReg' filesep  dir_vol(v).name(1:end-7)  '_MocoReg.nii.gz']);
        end
        
        dir_mocos = dir(['ZRef2MocoReg' filesep '*MocoReg*.nii.gz']);
        
        mkdir('normalizedVols2')
        
        for m = 1:size(dir_mocos,1)
            
            system(['sct_apply_transfo -i ZRef2MocoReg' filesep dir_mocos(m).name ...
                ' -d ' scttemplatepath filesep 'PAM50' filesep 'template' filesep 'PAM50_t2.nii.gz ' ...
                ' -w ' processdatapath subjects(sub).name  filesep 'func' filesep 'warps' filesep 'warp_moco2volumes_mean2PAM50_t2.nii.gz' ...
                ' -o normalizedVols2' filesep dir_vol(m).name(1:end-7) '_normalized.nii.gz' ...
                ' -x spline'])
            
        end
        
    end
    
    clear dir_vol
    dir_vol = dir(['normalizedVols2' filesep '*_normalized*.nii.gz']);
    
    for sub = 1:size(subjects,1)
        for m = 1:size(dir_vol,1)
            
            dirData = fullfile(processdatapath, subjects(sub).name, 'func', 'reconstructZRef', 'normalizedVols2');
            cd(dirData)
            
            system(['fslmeants -i  ' dir_vol(m).name ...
                ' -m ' scttemplatepath filesep 'PAM50' filesep 'atlas' filesep 'WHOLE_GM.nii.gz ' ...
                ' -o '  dir_vol(m).name(1:end-7)  '_GM.txt --showall'])
        end
    end
    
end
end
