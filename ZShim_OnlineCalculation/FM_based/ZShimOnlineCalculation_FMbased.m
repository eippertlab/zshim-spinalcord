% This is the code to calculate the FM-based automated z-shim selection
% based on field map data. The output is z-shim indices that
% compensate the field-gradient in z-direction.
% The indices are written to a .txt file which can then be read by the
% sequence.

% For more details and description of the method please see the following paper:
% Kaptan, M., Vannesjo, S. J., Mildner, T., Horn, U., Hartley-Davies, R., Oliva, V., Brooks, J. C. W., Weiskopf, N., Finsterbusch, J., & Eippert, F. (2021). 
% Automated slice-specific z-shimming for fMRI of the human spinal cord. BioRxiv, 2021.07.27.454049. https://doi.org/10.1101/2021.07.27.454049

% To use the code please make sure that FSL and SCT are downloaded and that they
% can be called from MATLAB.

% 27.06.2021
% Merve Kaptan,Johanna Vannesjo
% mkaptan@cbs.mpg.de
% Copyright 2021 Merve Kaptan, MPI for Human Cognitive and Brain Sciences,
% Leipzig; Johanna Vannesjo, NTNU, Trondheim
clc;clear all;close all
%% Settings

% Flags global
% ------------
dcmDirname = '/data/pt_02098/ZSHIM_SpeedTest/26946.dd_20190628_115833.PRISMA/FILEZ';            % where dicoms will be sent (from scanner)
oFolderPath =  '/data/pt_02306/fMRI_pilot/ZShims/';                                             % the general output directory
subjectNumber= input('Enter the subject number or name please (needs to be a string input):  ') % name of the subject
oFolderDir = fullfile(oFolderPath,num2str(subjectNumber));                                      % create the folder for this subject
mkdir(oFolderDir);                                                                              % make the output directory

% Flags software
% --------------
% add fsl to the path
setenv('FSLDIR', '/afs/cbs.mpg.de/software/fsl/5.0.11/ubuntu-xenial-amd64/');
pathFSL = getenv('FSLDIR');
fsldirmpath = sprintf('%s/etc/matlab',pathFSL);
path(path, fsldirmpath);
setenv('FSLOUTPUTTYPE', 'NIFTI_GZ'); % this to tell what the output type would be
% define the path to dicom 2 nifti converter DCM2NIIX
dcm2niixpath  = '/afs/cbs.mpg.de/software/scripts/';
dcm2niixcommand = 'DCM2NIIX --version 20180622 dcm2niix'
% Also: do not forget to 'add' Spinal Cord Toolbox to your bash profile
% before starting MATLAB
s = system('sct_check_dependencies');
if s == 127
    warning('please make sure that you can call SCT functions from MATLAB')
end

% Flags files
% -------------
% file identifiers that were used in the raw image names
T2SeqName  = 'T2w';                 % identifier for T2-weighted scan
fmSeqName  = 'gre_field_mapping';   % identifier for phase map scan

% Flags analysis
% ---------------
% flags that are used for the slice-wise fit itself
smoothing = 1;                      % if the data will be smoothed or not 
smoothingSigma = '1,1,1';           % enter the smoothing sigma a text and seperated by a comma
fitSli = 8;                         % how many slices (plus& minus the target slice)
                                    % will be used for the slice-wise fit (in the resolution of phase map)!
beta = 0;                           % the weighting of the slices
zshimSteps = 21;                    % how many z-shim steps were acquired
basis = [1:4];                      % basis for the fit (offset, X,Y and Z)
res = [2.2 1 1]*1e-3;               % the resolution of the phase map in meters
% the following parameter is used to make sure for finding the slice in the field-map corresponding to the center of the EPI slice
% depending on the weighting of the slices (beta), the number of slices
% used for fitting (fitSli) can change the total number of slices used for
% fitting
if beta==0 && fitSli==8
    totalSli = 8;
elseif  beta==0 && fitSli==4
    totalSli = 4;
elseif beta==0 && fitSli==2
    totalSli = 2;
elseif beta==0 && fitSli==6
    totalSli = 6;    
end
epiSli = 24;                     % number of axial slices in the EPI acquisition  

%% data preparation
% 1. convert DICOM data and put the output nifti files in the user-defined folder
system([dcm2niixpath dcm2niixcommand ' -b y -z y -o '  ...
    oFolderDir ' -f "%t_%p_%s" "'  dcmDirname  '"'  ]);

cd(oFolderDir)
% 2. prepare the T2-weighted scan
% name of T2 after dcm2niix conversion:
T2Dir  = dir(fullfile(oFolderDir, ['*' T2SeqName '*.nii.gz']));
% if there is more than one T2-weighted image, warn the user
if size(T2Dir,1) ~= 1 
    warning('Multiple T2s are detected, using the first one!')
end
movefile(T2Dir(1).name, 'T2.nii.gz'); % rename T2 image
system('sct_propseg -i  T2.nii.gz -c t2 '); % segment T2 image

% 3. prepare the phase map
% name of phase map after dcm2niix conversion:
fmDir = dir(fullfile(oFolderDir, ['*' fmSeqName '*.nii.gz']));
movefile(fmDir(1).name, 'phase.nii.gz'); % rename phase map
% swap dimensions of the phase image because it was a sagittal acquisition
system([fullfile(pathFSL,'bin','fslswapdim')  ' phase.nii.gz RL PA IS  phase_swapped.nii.gz']);

% 4. register the T2 to phase-map to be able to use the segmentation of T2
system('sct_register_multimodal -i T2.nii.gz -d phase_swapped.nii.gz -identity 1 -o RegT2.nii.gz');
system('sct_apply_transfo -i T2_seg.nii.gz -d phase_swapped.nii.gz -w warp_T22phase_swapped.nii.gz -o phase_seg.nii.gz -x nn');

%% compute shims

if smoothing
    system(['sct_maths -i phase_swapped.nii.gz -smooth ' smoothingSigma ' -o phase_swapped_smoothed.nii.gz']);
    warning('Phase image is smoothed.')
    fmFilename = 'phase_swapped_smoothed.nii.gz';
else 
    fmFilename = 'phase_swapped.nii.gz';
end

fm = read_avw(fmFilename);
mask = read_avw('phase_seg.nii.gz');
mask = logical(mask);

% Rescale field map (please note that this rescaling is only valid for the vendor based field-map)
dTE = 2.46e-3;            % difference in echo time
HzMax = 1/(2*dTE);
fmMax = max(abs(fm(:)));
fm = fm*HzMax/fmMax;
fmSli = size(fm,3);       % number of axial slices (after swapping the dimensions)

% Cut the field map 2 slices above and 4 slices below of the coverage of the EPI slice
% In our acquisition the slice stack was placed in a way that
% the center of the EPI volume (24 slices) corresponds to the center of the
% fm acquisition (180 sagittal slices) --> therefore,the middle point of the
% acquisitions are identical
fm = fm(:,:,(fmSli/2)-(epiSli/2*5)-4:(fmSli/2)+(epiSli/2*5)+2);
mask = mask(:,:,(fmSli/2)-(epiSli/2*5)-4:(fmSli/2)+(epiSli/2*5)+2);

% Compute the shims - main code
[shims,~,~,~,shimSli] = compute_shims(fm,mask,basis,fitSli*res(3),res,beta);

Hz_scale = [1 5e-3 5e-3 5e-3]; % thickness of our EPI slice
shims = diag(Hz_scale(basis))*shims; % convert shims to Hz
Zshim = shims(basis==4,:); % find the z-shims


mP = 0.21e-3;                % [Ts/m] maximum Z-shim gradient moment
mP = mP/floor(zshimSteps/2); % [Ts/m] single Z-shim step gradient moment
gamma = 42.6e6;              % [Hz/T] gyromagnetic ratio
slice = 5e-3;                % [m] slice thickness
mHz = mP*gamma*slice;        % [1]Z-shim frequency range within slice
pHz = mHz;                   % [Hz]
Zshim_picks = Zshim/pHz + ceil(zshimSteps/2);
Zshim_picks = round(Zshim_picks);

%% Create results matrix fitArray
% rescaled, rounded Zshim picks for each slice in the saggital fieldmap
% (will be used as an input to sequence)
fitArray(:,1)= Zshim_picks;
% which slices were used for fitting - output from main code
fitArray(:,2) = shimSli;    

% which 5 slices of the field map correspond to the one EPI slice
if totalSli == 4
    vector = [NaN(1,3) repmat([1 2 3 4 5], 1, 24)]';
elseif totalSli == 6
    vector = [NaN(1,2) repmat([1 2 3 4 5], 1, 24)]';
    vector(end) = [];
elseif totalSli == 8
    vector = [NaN(1,1) repmat([1 2 3 4 5], 1, 24)]';
    vector(end-1:end) = [];
elseif totalSli ==2
    vector = [NaN(1,4) repmat([1 2 3 4 5], 1, 24) NaN(1,1)]';
end
fitArray(:,3) = vector;

Zpicks = fitArray((fitArray(:,3)==3),1);
% calculated moment should not be smaller than 1 although the fit can
% result in negative values
Zpicks(Zpicks<1) = 1; 
sli= [1:24]';

% print values to command window
%print these values to the screen
fprintf('%s %s\n','slicesNumber  ',' selectedShims  ')
[sli, Zpicks]

% save this to a .mat file
save('AutoSelection.mat','Zpicks')

% now create the file that will be read by the sequence
fid = fopen( 'nin_z-shim_test.shm', 'w' );
fprintf( fid ,'%d \n',numSlices);         %first write the number of slices and go to a new line
fprintf( fid, '%d ',Zpicks(1:end-1));     %write the selected shim values (but not the last one so it does not have space afterwards)
fprintf( fid, '%d',Zpicks(end));          %now write the last shim value
fclose(fid);                              %close the file

