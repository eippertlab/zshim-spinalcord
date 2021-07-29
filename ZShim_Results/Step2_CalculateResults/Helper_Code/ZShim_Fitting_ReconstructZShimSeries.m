function ZShim_Fitting_ReconstructZShimSeries(dirData, imgRefScan, numSlices, zShimDir, imgRecon)
% Reconstruct artificial volumes by splitting the reference scan into 
% volumes (cycled through z-shim moments) and slices and constructing a new image based on provided 
% zshims for each slice.
% ----------------------
% inputs
% ---------------------
% dirData       = location of zshim reference scan
% imgRefScan    = filename of zshim reference scan
% numSlices     = number of slices
% zShimDir      = csv or mat file with zShims to be selected per slice
% imgRecon      = file name for resulting reconstruction image

% Merve Kaptan, mkaptan@cbs.mpg.de

if contains(zShimDir, '.csv')
    
    zShims = load(zShimDir);
   
elseif ~ contains(zShimDir, '.csv') 
    
    tmpZ = load(zShimDir);
    tmpZ = struct2cell(tmpZ);
    zShims = cell2mat(tmpZ);
    
end

zShims(zShims<= 0) = 1;

% Go to data directoy, create tmp directory, copy reference scan and go there
cd(dirData);
[status, ~] = system('mkdir tmp');
copyfile(imgRefScan, 'tmp')

cd('tmp');
zrefname = imgRefScan;

% Split two times: first according to z-shims then according to slices
[status, ~] = system(['fslsplit  ' zrefname  ' ZRef -t']);
[status, ~] = system(['for j in $(ls ZRef0*); do echo $j; ' 'fslsplit $j $( remove_ext $j) -z; done']);

% Create file-name matrix based on z-shim moment and slice number
fileNames = [];
for sliInd = 1:numSlices
    fileNames = [fileNames 'ZRef' sprintf('%04d', zShims(sliInd)-1) sprintf('%04d', sliInd-1) '.nii.gz '];
end

% Create reconstructed z-shim series
[status, ~] = system(['fslmerge -z ' imgRecon ' ' fileNames]);

% Move back to data directory, copy recon img and delete tmp direcoty
cd(dirData);
copyfile(['tmp' filesep imgRecon '.nii.gz'], [dirData filesep])
[status, ~] = system('rm -Rf tmp');