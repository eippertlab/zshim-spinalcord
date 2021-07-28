% Preprocessing- methods section 2.5 

% In this code, the preprocessing steps are listed for the following
% manuscript: WEBLINK
% The organization of the preprocessing steps was kept consistent with the
% manuscript.

% In order to run this code the Spinal Cord Toolbox (version 4.2.2 or higher -
% note that for the manuscript version 4.2.2 was used)
% and FSL should be added to the bash profile of the user to call their
% functions from MATLAB.

% When manual processing was performed, the respective outputs can be found here: 
% https://openneuro.org/datasets/ds003743
% under the "derivatives" parent directory and subject specific
% subdirectories

% Merve Kaptan, mkaptan@cbs.mpg.de; mervekaptan5@gmail.com
% 04.07.2021
%%
% set paths and directories

clc;clear all; close all;


s = system('sct_check_dependencies');
if s == 127
    warning('please make sure that you can call SCT functions from MATLAB')
end

scttemplatepath = '/data/pt_02098/ZShim_BIDS_0807/derivatives/template/SCT/';  %SCT version 4.2.2
rawdatapath     = '/data/pt_02098/ZShim_BIDS_0807/rawdata/';                   %where raw data is
processdatapath = '/data/pt_02098/ZShim_BIDS_0807/derivatives/';               %where processed data is 

mkdir(processdatapath)

cd(rawdatapath)

subjects = dir('*sub-*');  % get the list of the subjects

%% 2.5.1 Motion-correction of fMRI data 

for sub = 1:size(subjects,1)
    
    outdir = fullfile(processdatapath,subjects(sub).name,'func')
    mkdir(outdir)
    cd(fullfile(rawdatapath,subjects(sub).name,'func'))
    
    % Motion-correction for timeseries data with different echo times

    % copy original files to the output directory
    copyfile('*TE*.nii.gz', outdir);
    
    cd(outdir)

    % merge the mean of motion corrected images (from step 1) 
    % and the timeseries with different TEs
    system(['fslmerge -t mergedallTEs_MocoMean.nii.gz moco1volumes_mean.nii.gz ' ...
        '*TE*.nii.gz'])
    % motion correction of TE timeseries using moco1volumes_mean as a target
    % (obtained from 1st step of the motion-correction)
    system(['sct_fmri_moco -i  mergedallTEs_MocoMean.nii.gz -m moco_mask.nii.gz -param iterAvg=0 -x spline'])
    
    % remove unnecessary files
    system('rm -rf *_bold.nii.gz*')
    system('rm -rf *_T0*.nii.gz')
    system('rm -rf *task-rest*zshimTE*.nii.gz')
    
    % split these moco volumes and merge them back for each z-shim condition
    system('fslsplit mergedallTEs_MocoMean_moco.nii.gz ')
    
    % remove the first volume (=mean image)
    % then merge different groups of the 25 volumes as motion-corrected TEs)
    delete('vol0000.nii.gz');
    
    dir_img = dir('*vol0*');
    
    for v =1:25 % first group of 25 volumes were auto z-shim timeseries with TE = 30
        filenameauto1{v} = [dir_img(v).name ];
    end
    
    system(['fslmerge -t auto_moco30.nii.gz '  cell2mat(join(filenameauto1))]);
    clear filenameauto1
    
    for v = 26:50 % 2. group of 25 volumes were auto z-shim timeseries with TE = 40
        filenameauto2{v} = [dir_img(v).name ];
    end
    
    filenameauto2 = filenameauto2(26:50);
    system(['fslmerge -t auto_moco40.nii.gz '  cell2mat(join(filenameauto2))]);
    clear filenameauto2
    
    for v = 51:75 % 3. group of 25 volumes were auto z-shim timeseries with TE = 40
        filenameauto3{v} = [dir_img(v).name ];
    end
    
    filenameauto3 = filenameauto3(51:75);
    system(['fslmerge -t auto_moco50.nii.gz '  cell2mat(join(filenameauto3))]);
    clear filenameauto3
    
    
    for v = 76:100 % 4. group of 25 volumes were manual z-shim timeseries with TE = 30
        filenamemanual1{v} = [dir_img(v).name ];
    end
    
    filenamemanual1 = filenamemanual1(76:100);
    system(['fslmerge -t manual_moco30.nii.gz '  cell2mat(join(filenamemanual1))]);
    clear filenamemanual1
    
    
    for v = 101:125 % 5. group of 25 volumes were manual z-shim timeseries with TE = 40
        filenamemanual2{v} = [dir_img(v).name ];
    end
    filenamemanual2 = filenamemanual2(101:125);
    system(['fslmerge -t manual_moco40.nii.gz '  cell2mat(join(filenamemanual2))]);
    clear filenamemanual2
    
    
    for v = 126:150 % 6. group of 25 volumes were manual z-shim timeseries with TE = 40
        filenamemanual3{v} = [dir_img(v).name ];
    end
    filenamemanual3 = filenamemanual3(126:150);
    system(['fslmerge -t manual_moco50.nii.gz '  cell2mat(join(filenamemanual3))]);
    clear filenamemanual3
    
    for v = 151:175 % 7. group of 25 volumes were no z-shim timeseries with TE = 30
        filenameno1{v} = [dir_img(v).name ];
    end
    
    filenameno1 = filenameno1(151:175);
    system(['fslmerge -t no_moco30.nii.gz '  cell2mat(join(filenameno1))]);
    clear filenameno1
    
    
    for v = 176:200 % 8. group of 25 volumes were no z-shim timeseries with TE = 40
        filenameno2{v} = [dir_img(v).name ];
    end
    
    filenameno2 = filenameno2(176:200);
    system(['fslmerge -t no_moco40.nii.gz '  cell2mat(join(filenameno2))]);
    clear filenameno2
    
    
    for v = 201:225 % 9. group of 25 volumes were no z-shim timeseries with TE = 50
        filenameno3{v} = [dir_img(v).name ];
    end
    
    filenameno3 = filenameno3(201:225);
    system(['fslmerge -t no_moco50.nii.gz '  cell2mat(join(filenameno3))]);
    clear filenameno3
    
    % remove unnecessary output
    system('rm -rf vol0*.nii.gz')
    system('rm -rf *mergedallTEs*.nii.gz')
end
%% 2.5.2 Segmentation
% Segmentation of T2w image

for sub = 1:size(subjects,1)
    
    subid = subjects(sub).name;
    
    % create the output directory
    outdir = fullfile(processdatapath,subid,'anat');
    mkdir(outdir)
    
    % copy the T2-weighted anatomical image there and cd to that directory
    cd(fullfile(rawdatapath,subid))
    copyfile([ 'anat' filesep '*_T2w.nii.gz*'], [outdir filesep]);
    cd(outdir)
    
    % make a tmp directory to keep the folder organization
    mkdir('tmp')
    copyfile('*_T2w*.nii.gz', ['tmp' filesep])
    cd tmp
    
    % first segment the T2 initially
    % then smooth this segmentation with 8mm kernel
    % and use the smoothed image for final segmentation
    system('sct_deepseg_sc -i    *_T2w*  -c t2 -brain 0 ' );
    system('sct_smooth_spinalcord -i *_T2w.nii.gz* -s  *seg.nii.gz* -smooth 0,0,8 ');
    system('sct_deepseg_sc -i *smooth* -c t2 -brain 0 ');
    
    % move the main segmentation to the parent folder
    movefile('*T2w_smooth_seg*', [outdir filesep] )
    
    cd(outdir)
    % remove unnecessary ones
    system('rm -rf tmp')
end

% Segmentation of single volumes under different z-shim conditions:
% the single volumes under different z-shim conditions were averaged
% and the segmentation was done manually using functionalities of fsleyes
for sub = 1:size(subjects,1)
    
    outdir = fullfile(processdatapath,subjects(sub).name,'func')
    cd(fullfile(rawdatapath,subjects(sub).name,'func'))

    % copy original files to the output directory
    copyfile('*onevol*.nii.gz', outdir);
    cd(outdir)
    
    % merge individual volumes
    system(['fslmerge -t mergedall_onevol.nii.gz *onevol*.nii.gz'])
    
    % take the mean of these merged volumes
    system(['fslmaths mergedall_onevol.nii.gz -Tmean onevol_merged_mean.nii.gz'])
    
    % segment manually (name of the file: onevol_merged_mean_manualseg.nii.gz)

    % remove unnecessary files
    delete('mergedall_onevol.nii.gz')
end

% Segmentation of the mean of motion corrected time series under different
% z-shim conditions :
% The segmentation was done manually using functionalities of fsleyes.
% The mean of motion corrected volumes was saved as moco2volumes_mean.nii.gz under
% processdatapath/subjects(sub).name/'func' for each subject. This was used
% for manual segmentation.
% The result can be found in the same folders saved as moco2volumes_mean_manualseg.nii.gz

%% 2.5.3 Registration to template space
% Step a): register T2w image to the template space 

for sub = 1:size(subjects,1)
    
    subid = subjects(sub).name;
    
    % go to the output directory
    outdir = fullfile(processdatapath,subid,'anat');
    cd(outdir)
    
    % label vertebrae automatically
    system(['sct_label_vertebrae -i *T2w.nii.gz* -s ' outdir ...
        filesep '*T2w_smooth_seg.nii* -c t2 ']);
    % remove unnecessary labels
   system (['sct_label_utils -i *T2w_smooth_seg_labeled_discs.nii* '...
        '-keep 2,3,4,5,6,7,8,9 -o ' outdir  filesep 'disc_labels.nii.gz'])
    % note that the  labels were manually corrected if necessary (the name of
    % the file was not changed! The labels used for the normalization can be found 
    % under derivatives/subid/anat/ folder)
    
    % register T2w image to the template
     system(['sct_register_to_template -i *T2w.nii.gz* -s *_smooth_seg.nii* -ldisc disc_labels.nii.gz ' ...
        ' -param step=1,type=seg,algo=slicereg,metric=MeanSquares,' ...
        'iter=10,smooth=2,gradStep=0.5,slicewise=0,smoothWarpXY=2,' ...
        'pca_eigenratio_th=1.6:step=2,type=seg,algo=bsplinesyn,metric=MeanSquares,' ...
        'iter=3,smooth=1,gradStep=0.5,slicewise=0,smoothWarpXY=2,pca_eigenratio_th=1.6 ' ...
        ' -c t2']);
    % rename the normalized T2w image
    system('mv anat2template.nii.gz T2w_normalized.nii.gz')
    
    % remove unnecessary files
    delete('*T2w.nii.gz*')
    delete('*.cache*')
    delete('straight_ref.nii.gz')
    
    % make a directory for warping fields for organizational purposes
    mkdir('warps')
    system(['mv *warp*.nii* warps' filesep])
    % rename the warp fields
    cd warps
    system('mv warp_template2anat.nii.gz warp_PAM502T2.nii.gz')
    system('mv warp_anat2template.nii.gz warp_T22PAM50.nii.gz')

end

% step b) bring onevol functional images into template space
%   step b) i) register mean of onevol images (irrespective of condition)
%           using warping fields of step a) as a starting point
%   step b) ii) normalize onevol images of each condition using warping
%           fields of step b) i)

for sub = 1:size(subjects)
    
    outdir = [processdatapath subjects(sub).name filesep 'func'];
    cd(outdir)
    
    system(['sct_register_multimodal ' ...
        ' -i ' scttemplatepath   filesep 'PAM50' filesep 'template' filesep 'PAM50_t2.nii.gz -d onevol_merged_mean.nii.gz ' ...
        ' -iseg ' scttemplatepath   filesep 'PAM50' filesep 'template' filesep 'PAM50_cord.nii.gz -dseg onevol_merged_mean_manualseg.nii.gz ' ...
        ' -param step=1,type=seg,algo=centermass:step=2,type=seg,algo=bsplinesyn,metric=MeanSquares,smooth=1,slicewise=1,iter=3 ' ...
        ' -initwarp ' processdatapath subjects(sub).name filesep 'anat' filesep 'warps' filesep 'warp_PAM502T2.nii.gz ' ...
        ' -initwarpinv ' processdatapath subjects(sub).name filesep 'anat' filesep 'warps' filesep 'warp_T22PAM50.nii.gz -x spline ']);
    
    system(['sct_apply_transfo -i *task-onevol_acq-nozshim*.nii.gz  -d ' ...
        scttemplatepath filesep 'PAM50' filesep 'template' filesep 'PAM50_t2.nii.gz -w ' ...
        ' warp_onevol_merged_mean2PAM50_t2.nii.gz -o no_onevol_normalized.nii.gz '])
    
    system(['sct_apply_transfo -i *task-onevol_acq-manualzshim*.nii.gz  -d ' ...
        scttemplatepath filesep 'PAM50' filesep 'template' filesep 'PAM50_t2.nii.gz -w ' ...
        'warp_onevol_merged_mean2PAM50_t2.nii.gz -o manual_onevol_normalized.nii.gz '])
    
    system(['sct_apply_transfo -i *task-onevol_acq-autozshim*.nii.gz  -d ' ...
        scttemplatepath filesep 'PAM50' filesep 'template' filesep 'PAM50_t2.nii.gz -w ' ...
        ' warp_onevol_merged_mean2PAM50_t2.nii.gz -o auto_onevol_normalized.nii.gz '])
    
    mkdir('warps')
    system(['mv *warp*.nii* warps' filesep])    
    
end

% step c) bring resting state EPIs and tSNR maps into template space
%   step c) i) register mean moco image (of all time series irrespective of condition) 
%               to the template using warping fields of step a) as a
%               starting point

for sub = 1:size(subjects,1)
    
    cd([processdatapath subjects(sub).name filesep 'func'])
    
    system(['sct_register_multimodal -i '  scttemplatepath   filesep 'PAM50' filesep 'template' filesep 'PAM50_t2.nii.gz ' ...
        ' -d moco2volumes_mean.nii.gz' ...
        ' -iseg ' scttemplatepath   filesep 'PAM50' filesep 'template' filesep 'PAM50_cord.nii.gz ' ...
        ' -dseg moco2volumes_mean_manualseg.nii.gz ' ...
        ' -param step=1,type=seg,algo=centermass:step=2,type=seg,algo=bsplinesyn,metric=MeanSquares,smooth=1,slicewise=1,iter=3 ' ...
        ' -initwarp ' processdatapath subjects(sub).name filesep 'anat' filesep 'warps' filesep 'warp_PAM502T2.nii.gz ' ...
        ' -initwarpinv ' processdatapath subjects(sub).name filesep 'anat' filesep 'warps' filesep 'warp_T22PAM50.nii.gz -x spline ' ...
        ' -ofolder warps' filesep]);
end

%   step c) ii) normalize the mean images and tSNR maps for each condition 
%               using the warping field of step c) i)

zmode = {'no', 'manual', 'auto'};

for sub = 1:size(subjects,1)
    for z = 1:numel(zmode)
        
        cd([processdatapath filesep subjects(sub).name ])
        
        system(['sct_apply_transfo -i func' filesep zmode{z} '_moco_mean.nii.gz ' ...
            ' -d ' scttemplatepath   filesep 'PAM50' filesep 'template' filesep 'PAM50_t2.nii.gz ' ...
            ' -w   func' filesep 'warps' filesep 'warp_moco2volumes_mean2PAM50_t2.nii.gz ' ...
            ' -o   func' filesep zmode{z} '_moco_mean_normalized.nii.gz ' ...
            ' -x spline '])       
        
        system(['sct_apply_transfo -i func' filesep zmode{z} '_moco_tsnr.nii.gz ' ...
            ' -d ' scttemplatepath   filesep 'PAM50' filesep 'template' filesep 'PAM50_t2.nii.gz ' ...
            ' -w   func' filesep 'warps' filesep 'warp_moco2volumes_mean2PAM50_t2.nii.gz ' ...
            ' -o   func' filesep zmode{z} '_moco_tsnr_normalized.nii.gz ' ...
            ' -x spline '])
    end
end

% step d) bring timeseries with different echo times into template space
%   step d) i) take the motion corrected first image of TE timeseries 
%   step d) ii) normalize the TE images using the warping field of step c) i)
% 
zmode = {'no', 'manual', 'auto'};
TEs   = {'30', '40', '50'};

for sub = 1:size(subjects,1)
    cd([processdatapath subjects(sub).name filesep 'func'])
    for z = 1:numel(zmode)
        for t = 1:numel(TEs)
            
            system(['sct_image -i '  zmode{z} '_moco' TEs{t} '.nii.gz '...
                ' -keep-vol  0 '  ...
                ' -o ' zmode{z} '_moco' TEs{t}  '_singlevol.nii.gz ' ]);
             
            system(['sct_apply_transfo -i ' zmode{z} '_moco'  TEs{t} '_singlevol.nii.gz ' ...
                ' -d ' scttemplatepath   filesep 'PAM50' filesep 'template' filesep 'PAM50_t2.nii.gz ' ...
                ' -w warps' filesep 'warp_moco2volumes_mean2PAM50_t2.nii.gz ' ...
                ' -o  '  zmode{z} '_moco'  TEs{t}  '_singlevol_normalized.nii.gz ' ...
                ' -x spline '])
            
        end
    end
end

% step e) calculate MEAN of the moco images under different z-shim
% conditions & TEs, normalize these to the template space using the 
% warping field of step c) i)
   
for sub = 1:size(subjects,1)
    
    cd([processdatapath subjects(sub).name filesep 'func'])
    
    for z = 1:numel(zmode)
        for t = 1:numel(TEs)
            
            system(['fslmaths ' zmode{z} '_moco' TEs{t} ' -Tmean '  zmode{z} '_moco' TEs{t} '_mean'])
            
            system(['sct_apply_transfo -i ' zmode{z} '_moco'  TEs{t} '_mean.nii.gz ' ...
                ' -d ' scttemplatepath   filesep 'PAM50' filesep 'template' filesep 'PAM50_t2.nii.gz ' ...
                ' -w   warps' filesep 'warp_moco2volumes_mean2PAM50_t2.nii.gz ' ...
                ' -o ' zmode{z} '_moco'  TEs{t}  '_mean_normalized.nii.gz ' ...
                ' -x spline '])
                      
        end
    end
end
 
% step f) normalization of field map to the template space

for sub = 1:size(subjects,1)
    
    cd([rawdatapath subjects(sub).name  filesep 'fmap'])
    outdir = [processdatapath subjects(sub).name filesep 'fmap'];
    mkdir(outdir)
    
    copyfile('*phasediff*.nii.gz', outdir)
    copyfile('*magnitude1*.nii.gz', outdir)
    
    cd([processdatapath subjects(sub).name filesep 'fmap'])
    
    system (['sct_register_multimodal -i *phasediff*.nii.gz ' ...
        ' -d ' rawdatapath subjects(sub).name filesep 'anat' filesep '*T2w*.nii.gz '...
        ' -identity 1 ' ...
        ' -o phase_reg2T2.nii.gz ' ...
        ' -x nn']);
      
    system (['sct_apply_transfo -i phase_reg2T2.nii.gz ' ...
        ' -d ' scttemplatepath filesep 'PAM50' filesep 'template' filesep 'PAM50_t2.nii.gz'...
        ' -w ' processdatapath subjects(sub).name filesep 'anat' filesep 'warps' filesep 'warp_T22PAM50.nii.gz ' ...
        ' -o   phase_normalized.nii.gz ' ...
        ' -x nn']);
    
end

%% 2.5.4 EPI signal extraction

% threshold and binarize the gray matter masks

% where SCT atlases are
cd([scttemplatepath filesep 'PAM50' filesep 'atlas' ])
% threshold the gray matter maps
system('fslmaths PAM50_atlas_30.nii.gz -thrp 90 PAM50_atlas_30_thr90.nii.gz')
system('fslmaths PAM50_atlas_31.nii.gz -thrp 90 PAM50_atlas_31_thr90.nii.gz')
system('fslmaths PAM50_atlas_34.nii.gz -thrp 90 PAM50_atlas_34_thr90.nii.gz')
system('fslmaths PAM50_atlas_35.nii.gz -thrp 90 PAM50_atlas_35_thr90.nii.gz')

% binarize them, see info_label.txt for the names
system('fslmaths PAM50_atlas_30_thr90.nii.gz -bin LEFT_VENTRAL ')
system('fslmaths PAM50_atlas_31_thr90.nii.gz -bin RIGHT_VENTRAL ')
system('fslmaths PAM50_atlas_34_thr90.nii.gz -bin LEFT_DORSAL ')
system('fslmaths PAM50_atlas_35_thr90.nii.gz -bin RIGHT_DORSAL ')

% combine them to get the whole cord gray matter map
system('fslmaths LEFT_VENTRAL -add  RIGHT_VENTRAL -add LEFT_DORSAL -add RIGHT_DORSAL WHOLE_GM')

% extract signal
% 1. Extract signal for single volume images

zmode = {'no', 'manual', 'auto'};
masks = {'LEFT_DORSAL', 'RIGHT_DORSAL', 'LEFT_VENTRAL', 'RIGHT_VENTRAL'};

for sub = 1:size(subjects,1)
    
cd([processdatapath subjects(sub).name filesep 'func'])

    for z = 1:numel(zmode)
        
        system(['fslmeants -i ' zmode{z}  '_onevol_normalized.nii.gz ' ...
            ' -m ' scttemplatepath filesep 'PAM50' filesep 'template' filesep 'PAM50_cord.nii.gz ' ...
            ' -o ' zmode{z}  '_onevol_normalized_PAM50_cord.txt --showall ']);

        system(['fslmeants -i ' zmode{z}  '_onevol_normalized.nii.gz ' ...
            ' -m   moco2_volumes_mean_manualseg.nii.gz' ...
            ' -o ' zmode{z}  '_onevol_native.txt --showall ']);
        
        for m = 1:numel(masks)
            
            system(['fslmeants -i ' zmode{z}  '_onevol_normalized.nii.gz ' ...
                ' -m ' scttemplatepath filesep 'PAM50' filesep 'atlas' filesep masks{m} ...
                ' -o   '  zmode{z}  '_onevol_normalized_' masks{m} '.txt --showall ']);            
        end
        
    end
end

% 2. Extract signal for motion-corrected time-series EPI

for sub = 1:size(subjects,1)
    
    cd([processdatapath subjects(sub).name filesep 'func'])
    
    for z = 1:numel(zmode)
        
        system(['fslmeants -i ' zmode{z}  '_moco_mean_normalized.nii.gz ' ...
            ' -m ' scttemplatepath filesep 'PAM50' filesep 'template' filesep 'PAM50_cord.nii.gz ' ...
            ' -o   '  zmode{z}  '_moco_mean_normalized_PAM50_cord.txt --showall ']);

        system(['fslmeants -i ' zmode{z}  '_moco_mean.nii.gz ' ...
            ' -m   moco2_volumes_mean_manualseg.nii.gz' ...
            ' -o ' zmode{z}  '_moco_mean_native.txt --showall ']);
     
        system(['fslmeants -i ' zmode{z}  '_moco_mean.nii.gz ' ...
            ' -m   moco2_volumes_mean_manualseg.nii.gz' ...
            ' -o ' zmode{z}  '_moco_mean_native.txt --showall ']);
        
        for m = 1:numel(masks)
            
            
            system(['fslmeants -i ' zmode{z}  '_moco_mean_normalized.nii.gz ' ...
                ' -m ' scttemplatepath filesep 'PAM50' filesep 'atlas' filesep  masks{m} ...
                ' -o   '  zmode{z}  '_moco_mean_normalized_'  masks{m}  '.txt --showall ']);
            
            system(['fslmeants -i ' zmode{z}  '_moco_tsnr_normalized.nii.gz ' ...
                ' -m ' scttemplatepath filesep 'PAM50' filesep 'atlas' filesep  masks{m}...
                ' -o   '  zmode{z}  '_moco_tsnr_normalized_'  masks{m}  '.txt --showall ']);
            
        end
        
        
    end
end

% 3. Extract signal for motion-corrected time-series with different TEs

zmode = {'no', 'manual', 'auto'};
TEs   = {'30', '40', '50'};

for sub = 1:size(subjects,1)
    
    cd([processdatapath subjects(sub).name filesep 'func'])
    
    for z = 1 :numel(zmode)
        for t = 1:numel(TEs)
            
            system(['fslmeants -i ' zmode{z} '_moco'  TEs{t} '_singlevol_normalized.nii.gz ' ...
                ' -m ' scttemplatepath filesep 'PAM50' filesep 'template' filesep 'PAM50_cord.nii.gz ' ...
                ' -o ' zmode{z} '_moco'  TEs{t}  '_singlevol_normalized_PAM50_cord.txt --showall '])
            
            system(['fslmeants -i ' zmode{z} '_moco'  TEs{t} '_mean_normalized.nii.gz ' ...
                ' -m ' scttemplatepath filesep 'PAM50' filesep 'template' filesep 'PAM50_cord.nii.gz ' ...
                ' -o ' zmode{z} '_moco'  TEs{t}  '_mean_normalized_PAM50_cord.txt --showall '])
            
        end
    end
end

%%
% Please note that the extracted signal is organized to calculate Results
% and the organized data could be found here:  under "derivatives" parent
% directory and "extracted_signal" directory
% for the results in the template space (signal_templatespace/) as well as in the native space (signal_nativespace/)
% Both for native & template space results GroupWhole contains data from
% all participants(N = 48), and GroupSingle_EPI and GroupSingle_FM folders
% for EPI based automation and field map based automation, respectively.
% GM/ folder contains the results for Gray Matter masks

% create directories for each result (an example is shown below/repeated for 
% different results/directories)

if ~exist([processdatapath filesep 'extracted_signal'])
    mkdir([processdatapath filesep 'extracted_signal'])
end

subjGroup = 'GroupWhole';                   % change depending on the group of the participants
saveDire  = 'SingleVolume_MeanSignalS';     % change depending on extracted data type
saveName  = 'Data.mat';                     % for gray matter (GM) saved under GM folder
                                            
zmode     = {'no', 'manual', 'auto'};       % z-shim mode
masks     = {'PAM50_cord'};                 % mask name
fileName  = 'onevol_normalized';

%cut where all subjects have data - for template space only!
ZMin = 710;
ZMax = 935;

for sub = 1:size(subjects,1)
    
    outdir = fullfile([processdatapath subjects(sub).name filesep 'func']);
    cd(outdir)
    
    for z = 1:numel(zmode)
        for m = 1:numel(masks)
            tmp = load([zmode{z} '_' fileName '_' masks{m} '.txt']);
            tmp = tmp(3:4,:)';
            minz = find(tmp(:,1)==ZMin);
            maxz = find(tmp(:,1) == ZMax);
            tmp = tmp(minz(1):maxz(end),:);
            tmp(:,1) = tmp(:,1) + 1;
            slices = unique(tmp(:,1));
            
            for s = 1:size(slices,1)
                signalsli = [];
                signalsli =  nanmean(tmp(find(tmp(:,1) == slices(s)),2));
                results(sub,z,m,s) = signalsli;
            end
            clear tmp;
        end
    end
end

if ~exist([processdatapath filesep 'extracted_signal' filesep saveDire])
    mkdir([processdatapath filesep 'extracted_signal' filesep saveDire])
end

cd([processdatapath filesep 'extracted_signal' filesep saveDire])

if numel(masks)==1
    results = squeeze(results);
    for z = 1:numel(zmode)
        eval([zmode{z}  '= squeeze(results(:,' num2str(z) ',:));'])
    end
    
    if numel(zmode) ==2
        save(saveName, zmode{1}, zmode{2})
    elseif numel(zmode) ==3
        save(saveName, zmode{1}, zmode{2}, zmode{3})
    end    
        
elseif numel(masks1)>1
    for r = 1:numel(zmode)
        for m = 1:numel(masks)
            eval([zmode{r} masks{m}  '= squeeze(results(:,' num2str(z),  'm,:));'])
        end
    end
    
end

%or save the data manually

