function ZShim_Calculate_FMbasedpicks(fm_type, maskType,fitParameters)

% inputs
% fm_type : type of the field map, 'gre' (vendor-based) or 'cbsfl' (in-house)
% maskType : 'eroded' or 'dilated' or '' (regular) or 'cbsfl'
% fitParameters: fitSli             

% set FSL environment
% Running FSL from matlab
setenv('FSLDIR', '/afs/cbs.mpg.de/software/fsl/5.0.11/ubuntu-xenial-amd64/');
fsldir = getenv('FSLDIR');
fsldirmpath = sprintf('%s/etc/matlab',fsldir);
path(path, fsldirmpath);
setenv('FSLOUTPUTTYPE', 'NIFTI_GZ'); % this to tell what the output type would be

zshimSteps = 21;
fitSli = fitParameters{1};          % this is what I input to Johanna's code
beta   = fitParameters{2};          % slice weighting
totalSli = fitParameters{3};        % total number of slices
smoothfactor = fitParameters{4};    % smoothing kernel

basis = [1:4];
res = [2.2 1 1]*1e-3; % the resolution of the phase map in meter

epiSli = 24;          % number of axial slices in the EPI acquisition
fmSli = 180;          % number of axial slices in the field map acquisition (after swapping the dimensions)  


if  isequal(fm_type,'gre')
    
    if  smoothfactor == '0'
        fmFilename = ['phase_swap.nii.gz'];
    else
        fmFilename = ['phase'  smoothfactor      '.nii.gz'];
    end
    
elseif isequal(fm_type,'cbsfl_1')
    
    fmFilename = ['phase_cbsfl_1_'  smoothfactor      '.nii.gz'];
    
elseif isequal(fm_type,'cbsfl_2')
    
    fmFilename = ['phase_cbsfl_2_'  smoothfactor      '.nii.gz'];
    
end

if ~isempty(maskType)
    maskFilename = ['T2_mask_' maskType '.nii.gz'];
elseif isempty(maskType)
    maskFilename = ['T2_mask' maskType '.nii.gz'];
end

fm = read_avw(fmFilename);
mask = read_avw(maskFilename);
mask = logical(mask);

% Rescale field map
if isequal(fm_type, 'gre') % Siemens field map
    
    dTE = 2.46e-3; % difference in echo time
    HzMax = 1/(2*dTE);
    fmMax = max(abs(fm(:)));
    fm = fm*HzMax/fmMax;
    
elseif contains(fm_type, 'cbsfl') % in-house field map
    
    fm= (fm.*0.0977)-200;
    
end

% Cut the field map 5 slices above and 4 slices below of the coverage of the EPI slice
% In our acquisition the slice stack was placed in a way that
% the center of the EPI volume (24 slices) corresponds to the center of the
% fm acquisition (180 sagittal slices) --> therefore,the middle point of the
% acquisitions are identical
fm = fm(:,:,(fmSli/2)-(epiSli/2*5)-4:(fmSli/2)+(epiSli/2*5)+5);
mask = mask(:,:,(fmSli/2)-(epiSli/2*5)-4:(fmSli/2)+(epiSli/2*5)+5);

% Compute the shims - main code
[shims,~,~,~,shimSli] = compute_shims(fm,mask,basis,fitSli*res(3),res,beta);

Hz_scale = [1 5e-3 5e-3 5e-3]; % thickness of our EPI slice
shims = diag(Hz_scale(basis))*shims; % convert shims to Hz
Zshim = shims(basis==4,:); % find the z-shims


mP = 0.21e-3;                % [Ts/m] maximum Z-shim gradient moment
mP = mP/floor(zshimSteps/2); % [Ts/m] single Z-shim step gradient moment
gamma = 42.6e6;              % [Hz/T]
slice = 5e-3;                % [m]
mHz = mP*gamma*slice;        % [1] Z-shim phase accrual over slice
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
    vector = [NaN(1,3) repmat([1 2 3 4 5], 1, 24) NaN(1,3)]';
elseif totalSli == 8
    vector = [NaN(1,1) repmat([1 2 3 4 5], 1, 24) NaN(1,1)]';
elseif totalSli == 12
    vector = repmat([1 2 3 4 5], 1, 24)';
    vector([1 end]) = [];
end
  
fitArray(:,3) = vector;

if isequal(fm_type, 'gre')
    
    auto_gre = fitArray((fitArray(:,3)==3),1);
    auto_gre(auto_gre<=0) = 1; % fit can result in negative values,but the minimum index for the sequence is 1

    if ~isempty(maskType)
        save(['AutoSelection_FM_gre' num2str(fitSli) '_' maskType   '.mat'], 'auto_gre')
    elseif isempty(maskType)
        if totalSli == 8 && smoothfactor == '1' && beta == 0 && fitSli == 8 % if the parameter set is not standard save all the parameters
            save(['AutoSelection_FM_gre' num2str(fitSli)    '.mat'], 'auto_gre')
        else
            save(['AutoSelection_FM_gre' num2str(fitSli) '_beta' num2str(beta) '_slices' num2str(totalSli) '_smooth' num2str(smoothfactor) '.mat'], 'auto_gre')
            
        end
    end
    
elseif contains(fm_type, 'cbsfl')
    
    auto_cbsfl = fitArray((fitArray(:,3)==3),1);
    auto_cbsfl(auto_cbsfl<=0) = 1; % fit can result in negative values,but the minimum index for the sequence is 1

    
    if isequal(fm_type,'cbsfl_1')
        
        save('AutoSelection_FM_cbsfl_1.mat', 'auto_cbsfl')
    
    elseif isequal(fm_type,'cbsfl_2')
        
        save('AutoSelection_FM_cbsfl_2.mat', 'auto_cbsfl')
    
    end
        

end

end