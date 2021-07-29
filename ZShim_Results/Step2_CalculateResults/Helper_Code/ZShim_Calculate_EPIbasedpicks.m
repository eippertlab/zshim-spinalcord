function EPI_picks = ZShim_Calculate_EPIbasedpicks(referenceScanName)

% the code to calculate the EPI-based automated z-shim selection
% based on the z-shim reference scan. 
% adapted version for post-hoc calculations

% input
%---------

% referenceScanName: string, the name of the reference measurement

% Merve Kaptan, mkaptan@cbs.mpg.de

numSlices       = 24;                         % Number of slice for EPIs
zShimMoments    = 21;                         % Number of z-shim moments
step0           = 0;                          % Create reconstructed z-shim series for no z-shim and manual z-shim
step1           = 1;                          % Create a look-up table of mean signal intensity for each slice and z-shim moment, so I can later plot the signal traces (involves segmentation etc)

step1ImgMean    = [referenceScanName '_MEAN.nii.gz'];
step1ImgMeanSeg = [referenceScanName '_MEAN_seg.nii.gz'];
step1TxtOutput  = 'IntensitiesEPI.txt';
step1Threshold  = 100;
referenceScan   = [referenceScanName '.nii.gz'];


% Create mean (for better segmentation), segment and extract signals
[status, ~] = unix([' fslmaths ' referenceScan ' -Tmean ' step1ImgMean]);
[status, ~] = unix(['sct_propseg -i ' step1ImgMean ' -c t2s -radius 5']);
[status, ~] = unix([' fslmeants -i ' referenceScan ' -m ' step1ImgMeanSeg ' --showall -o ' step1TxtOutput]);

% Load data (remember that first 3 rows are coordinates)
data     = load(step1TxtOutput);
dataBOLD = data(4:end,:);

if length(unique(data(3,:)))== 24
    fprintf('Detected 24 slices \n')
else
    warning ('Not enough slices, mask will be adjusted \n')
    
end

if length(unique(data(3,:)))~= 24
    [status, ~] = unix(['sct_propseg -i ' step1ImgMean ' -c t2s -radius 5 -max-area 150']);
    [status, ~] = unix([  ' fslmeants -i ' referenceScan ' -m ' step1ImgMeanSeg ' --showall -o ' step1TxtOutput]);
    % Load data (remember that first 3 rows are coordinates)
    data     = load(step1TxtOutput);
    dataBOLD = data(4:end,:);
    
    if length(unique(data(3,:)))~= 24
        [status, ~] = unix(['sct_propseg -i ' step1ImgMean ' -c t2s -radius 5 -max-area 150 -min-contrast 30']);
        [status, ~] = unix([  ' fslmeants -i ' referenceScan ' -m ' step1ImgMeanSeg ' --showall -o ' step1TxtOutput]);
        % Load data (remember that first 3 rows are coordinates)
        data     = load(step1TxtOutput);
        dataBOLD = data(4:end,:);
        if length(unique(data(3,:)))~= 24
            [status, ~] = unix(['sct_propseg -i ' step1ImgMean ' -c t2s -radius 5 -max-area 150 -min-contrast 30 -max-deformation 5']);
            [status, ~] = unix([  ' fslmeants -i ' referenceScan ' -m ' step1ImgMeanSeg ' --showall -o ' step1TxtOutput]);
            % Load data (remember that first 3 rows are coordinates)
            data     = load(step1TxtOutput);
            dataBOLD = data(4:end,:);
        end
        
    end
    
end

if length(unique(data(3,:)))~= 24
    subjMaskProblem{sub,1} = direSub(sub).name;
else
    
    % Find which voxels come from which slice
    slices = cell(numSlices,1);
    for sub = 1:numSlices
        slices{sub} = find(data(3,:) == sub-1);
    end
    
    
    % Create look-up table containing mean signal for each slice and zshim
    lookupTable = NaN*ones(zShimMoments,numSlices);
    for sub = 1:numSlices
        lookupTable(:,sub) = mean(dataBOLD(:,slices{sub}),2,'omitnan');
    end
    csvwrite(['LookupTable_Thr' sprintf('%03d',step1Threshold) '.csv'], lookupTable);
    
    % Load look-up table
    lookupTable = csvread(['LookupTable_Thr' sprintf('%03d',step1Threshold) '.csv']);
    
    % Determine slice-specific z-shim giving maximum signal
    [~,zShimsAuto] = max(lookupTable);
    EPI_picks = zShimsAuto';
end
end