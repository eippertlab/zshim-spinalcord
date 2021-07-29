% This is the code to calculate the EPI-based automated z-shim selection
% based on the z-shim reference scan. The output is z-shim indices that maximizes the signal intensity in the cord.

% The indices are written to a .txt file which can then be read by the sequence.

% For more details and description of the method please see the following paper:
% Kaptan, M., Vannesjo, S. J., Mildner, T., Horn, U., Hartley-Davies, R., Oliva, V., Brooks, J. C. W., Weiskopf, N., Finsterbusch, J., & Eippert, F. (2021). 
% Automated slice-specific z-shimming for fMRI of the human spinal cord. BioRxiv, 2021.07.27.454049. 
% https://doi.org/10.1101/2021.07.27.454049

% To use the code please make sure that FSL, SCT and DCM2NIIX are downloaded and that they
% can be called from MATLAB.

% 27.06.2021
% Merve Kaptan
% mkaptan@cbs.mpg.de
% Copyright 2021 Merve Kaptan, MPI for Human Cognitive and Brain Sciences, Leipzig

clc;clear all;close all
%% Settings

% Flags global
% ------------
dcmDirname = '/nobackup/enigma2/prismasys/FE/';                                                 % where dicoms will be sent (from scanner)
oFolderPath =  '/data/pt_02306/fMRI_pilot/ZShims/';                                             % the general output directory
subjectNumber= input('Enter the subject number or name please (needs to be a string input):  ') % name of the subject
oFolderDir = [oFolderPath '/' num2str(subjectNumber)];                                          % create the folder for this subject
mkdir(oFolderDir);                                                                              % make the output directory

% Flags software
% --------------
% add fsl to the path
fslDir = '/afs/cbs.mpg.de/software/fsl/5.0.11/ubuntu-xenial-amd64/';
setenv('FSLDIR', fslDir);
pathFSL = getenv('FSLDIR');
fsldirmpath = sprintf('%s/etc/matlab',pathFSL);
path(path, fsldirmpath);
setenv('FSLOUTPUTTYPE', 'NIFTI_GZ'); % this to tell what the output type would be
% define the path to dicom 2 nifti converter DCM2NIIX
dcm2niixpath    = '/afs/cbs.mpg.de/software/scripts/';
dcm2niixcommand = 'DCM2NIIX --version 20180622 dcm2niix';
% Also: do not forget to 'add' Spinal Cord Toolbox to your bash profile
% before starting MATLAB
s = system('sct_check_dependencies');
if s == 127
    warning('please make sure that you can call SCT functions from MATLAB')
end


% converting DICOM data and putting the output nifti files in the user-defined folder
cmd = [dcm2niixpath dcm2niixcommand ' -b y -z y -o '  oFolderDir ' -f "%t_%p_%s" "'  dcmDirname  '"'  ];
system(cmd);


% Flags scan
% --------------
zrefSeqName = 'ZShimRef'; %the identifier for the z-reference sequence, so that after dcm2niix conversion, code can automatically find it
ZRefDir     = dir([oFolderDir '/*' zrefSeqName '*.nii.gz']);    % the actual name of the z-ref scan after conversion
% make sure that there is only one z-shim reference file in the directory
if size(ZRefDir,1) ~= 1 %there should NOT be more than one zref files so that there is not mistake
    error('Multiple reference scans are detected, please remove the extra ones!')
end

% Flags analysis
% --------------
referenceScan   = ZRefDir(1).name;            % Filename of z-shim series
numSlices       = 24;                         % Number of slice for EPIs
zShimMoments    = 21;                         % Number of z-shim moments

step1ImgMean    = 'ZRef_MEAN.nii.gz';
step1ImgMeanSeg = 'ZRef_MEAN_seg.nii.gz';
step1TxtOutput  = 'IntensitiesEPI.txt';
step1Threshold  = 100;                        % If set to 100, no data are discarded; otherwise data above
% this percentile are discarded considering that it may belong to CSF

%% Calculation part

cd (oFolderDir)
disp(['Working on subject ' subjectNumber]);

% Step 1: create a look-up table of mean signal intensity for each slice and z-shim moment
% ----------------------------------------------------------------------------------------

% Create mean (for better segmentation), segment and extract signals
[status, ~] = unix([pathFSL 'bin/fslmaths ' referenceScan ' -Tmean ' step1ImgMean]);
[status, ~] = unix(['sct_propseg -i ' step1ImgMean ' -c t2s -radius 5']);
[status, ~] = unix([pathFSL 'bin/fslmeants -i ' referenceScan ' -m ' step1ImgMeanSeg ' --showall -o ' step1TxtOutput]);

% Load data (remember that first 3 rows are coordinates)
data     = load(step1TxtOutput);
dataBOLD = data(4:end,:);

if length(unique(data(3,:)))== numSlices
    
    fprintf('Detected %d slices in the mask, no adjustment needed \n', numSlices)
else
    warning ('Not enough slices [ number of slices in the mask =  %d ], the mask will be adjusted \n', length(unique(data(3,:))))
    % In the following part, if the automatic segmentation does not
    % propogate the options for sct_propseg function that affect the z-propogation
    % will be modified.
    
    % if the segmentation does not propogate, first modify the -max-area
    % that affects maximum cross-sectional area (see sct_propseg for details)
    [status, ~] = unix(['sct_propseg -i ' step1ImgMean ' -c t2s -radius 5 -max-area 150']);
    [status, ~] = unix([pathFSL 'bin/fslmeants -i ' referenceScan ' -m ' step1ImgMeanSeg ' --showall -o ' step1TxtOutput]);
    % Load data (remember that first 3 rows are coordinates)
    data     = load(step1TxtOutput);
    dataBOLD = data(4:end,:);
    
    if length(unique(data(3,:)))~= numSlices
        % if the segmentation still does not propogate,  modify the -max-area
        % that affects maximum cross-sectional area and min-contrast (see sct_propseg for
        % details)
        [status, ~] = unix(['sct_propseg -i ' step1ImgMean ' -c t2s -radius 5 -max-area 150 -min-contrast 30']);
        [status, ~] = unix([pathFSL 'bin/fslmeants -i ' referenceScan ' -m ' step1ImgMeanSeg ' --showall -o ' step1TxtOutput]);
        % Load data (remember that first 3 rows are coordinates)
        data     = load(step1TxtOutput);
        dataBOLD = data(4:end,:);
        
        if length(unique(data(3,:)))~= numSlices
            % if the segmentation still does not propogate,  modify the -max-area
            % that affects maximum cross-sectional area, min-contrast and max-deformation (see sct_propseg for
            % details)
            [status, ~] = unix(['sct_propseg -i ' step1ImgMean ' -c t2s -radius 5 -max-area 150 -min-contrast 30 -max-deformation 5']);
            [status, ~] = unix([pathFSL 'bin/fslmeants -i ' referenceScan ' -m ' step1ImgMeanSeg ' --showall -o ' step1TxtOutput]);
            % Load data (remember that first 3 rows are coordinates)
            data     = load(step1TxtOutput);
            dataBOLD = data(4:end,:);
            if length(unique(data(3,:)))~= numSlices
                error ('Not enough slices in the mask, please manually adjust your mask!!')
                % If the image contrast is low or automatic segmentation
                % does not work for some reason, the mask needs to be
                % manually adjusted (fsleyes can be used).
            end
        end
        
    end
    
end

% Find which voxels come from which slice

slices = cell(numSlices,1);
for s = 1:numSlices
    slices{s} = find(data(3,:) == s-1);
end

% Remove CSF-part from data (by assuming prctile X is likely CSF)
dataTmp = dataBOLD(:);
threshold = prctile(dataTmp,step1Threshold);
dataBOLD(dataBOLD > threshold) = NaN;

% Create look-up table containing mean signal for each slice and zshim
lookupTable = NaN*ones(zShimMoments,numSlices);
for s = 1:numSlices
    lookupTable(:,s) = mean(dataBOLD(:,slices{s}),2,'omitnan');
end
csvwrite(['LookupTable_Thr' sprintf('%03d',step1Threshold) '.csv'], lookupTable);

% Clean up
clear status data dataBOLD slices s dataTmp threshold lookupTable;


% Step 2: find best z-shim value based on z-shim reference scan
% -------------------------------------------------------------

% Load look-up table
lookupTable = csvread(['LookupTable_Thr' sprintf('%03d',step1Threshold) '.csv']);

% Determine slice-specific z-shim giving maximum signal
[~,zShimsAuto] = max(lookupTable);

%print these values to the screen
fprintf('%s %s\n','slicesNumber  ',' selectedShims  ')
[(1:numSlices)', zShimsAuto']

%save the z-shim values to a .csv file
csvwrite(['ZShim_AutoRef_Thr' sprintf('%03d',step1Threshold) '.csv'], zShimsAuto');

% now create the file that will be read by the sequence
fid = fopen( 'nin_z-shim_test.shm', 'w' );
fprintf( fid ,'%d \n',numSlices);         %first write the number of slices and go to a new line
fprintf( fid, '%d ',zShimsAuto(1:end-1)); %write the selected shim values (but not the last one so it does not have space afterwards)
fprintf( fid, '%d',zShimsAuto(end));      %now write the last shim value
fclose(fid);                              %close the file


% Clean up
clear lookupTable zShimsAuto;

cd(oFolderDir);

%% remove all the dicoms after calculation

if ~isempty(dir(dcmDirname))
    cd(dcmDirname)
    warning('Removing dicoms sent from scanner...')
    delete *.IMA
    delete *.log
    if ~isempty(dir(dcmDirname))
        warning(['Check ' dcmDirname ' !!'])
    end
end


