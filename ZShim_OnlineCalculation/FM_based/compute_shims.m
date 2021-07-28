function [ SOut, fmShimmed, fmVar, shimMaps, shimSli] = compute_shims( fm,mask,basis,fitWidth,res,beta, mag )
% Fits shim fields to field map within mask
% 
% INPUT:
% fm            field map (3D)
% mask          mask (same size as field map)
% basis         shim field basis to be fitted
% fitWidth      width of fit kernel [m]
% res           resolution [m]
% center        center of coordinate system
% angulation    angulation of acquired map
%
% OUTPUT:
% S             calculated shim values
% fmShimmed     theoretical field map after shimming
% fmVar         field variance within mask, before and after shimming

% compute_shims(fm,mask,basis,fitSli*res(3),res);

afftrans = res;           

% Copyright Johanna Vannesjo, FMRIB, Oxford 2015

% Set default parameters
if nargin<2 || isempty(mask) 
    mask = ones(size(fm)); 
end
if nargin<3, basis = 1:4; end
if nargin<4, fitWidth = 5; end
if nargin<5, afftrans = []; end
if length(afftrans) == 3
    nArray = size(fm);
    center = floor(nArray/2)+1;
    if numel(center) < 3
        center(3) = 1;
    end
    afftrans = [...
        afftrans(1) 0 0 -afftrans(1)*center(1); ...
        0 afftrans(2) 0 -afftrans(2)*center(2); ...
        0 0 afftrans(3) -afftrans(3)*center(3)];     
end
if nargin < 7
    mag = ones(size(fm));
end

% Set sizes & indices
nX = size(fm,1);
nY = size(fm,2);
nSli = size(fm,3);
nBasis = length(basis);

% Compute spherical harmonic shim fields within map
fmCoords = compute_map_coords([nX nY nSli],afftrans);
B  = compute_SH_basis( fmCoords,basis );
B_maps = reshape(B,[nX nY nSli nBasis]);

% Calculate weighting kernel
if fitWidth == 0 || ~isfinite(fitWidth)
    fitSli = 1:nSli;
    nFitSli = nSli;
    shimSli = 0;
    w = ones([nX nY nSli]);
else
    dz = afftrans(3,3);
    beta = beta;
    fitSli1 = floor((1+beta)*fitWidth/(2*dz));
    w = sqrt(raised_cosine(dz*(-fitSli1:fitSli1)',1/fitWidth,beta));
    w = repmat(w,[1 nX nY]);
    w = permute(w,[2 3 1]);
    shimSli = fitSli1+1:(nSli-fitSli1);
    nFitSli = 2*fitSli1 + 1;
    fitSli = -fitSli1:fitSli1;
end

% Perform slice-wise shim fitting
fmShimmed = fm;
S = zeros(nBasis,length(shimSli));
for iS = 1:length(shimSli)
    cSli = shimSli(iS);
    % Weight with magnitude image
    wTemp = w.*mag(:,:,cSli+fitSli);
    fm_temp = wTemp.*fm(:,:,cSli+fitSli);
    B_temp = reshape(repmat(wTemp,[1 1 1 nBasis]).*B_maps(:,:,cSli+fitSli,:),[nX*nY*nFitSli nBasis]);
    mask_temp = logical(mask(:,:,cSli+fitSli));
    maskInds = find(mask_temp);
    
    F = fm_temp(mask_temp);
    S(:,iS) = B_temp(maskInds,:)\F; %#ok<FNDSB>
%     S(:,iS) = pinv(B_temp(maskInds,:))*F;
end

% Calculate shimmed field map
shimMaps = zeros(size(fmShimmed));
if length(fitSli) == nSli
    shimSli = 1:nSli;
    for iB = 1:nBasis
        shimMaps = shimMaps + S(iB)*B_maps(:,:,:,iB);
        fmShimmed = fmShimmed - S(iB)*B_maps(:,:,:,iB);
    end
else
    for iS = 1:length(shimSli)
        for iB = 1:nBasis
            shimMaps(:,:,shimSli(iS)) = shimMaps(:,:,shimSli(iS)) + S(iB,iS)*B_maps(:,:,shimSli(iS),iB);
            fmShimmed(:,:,shimSli(iS)) = fmShimmed(:,:,shimSli(iS)) - S(iB,iS)*B_maps(:,:,shimSli(iS),iB);
        end
    end
end

% Calculate slice-wise standard deviation
fmVar.slice_initial = zeros(nSli,1);
fmVar.slice_shimmed = zeros(nSli,1);
for iS = 1:nSli
    fm_temp = fm(:,:,iS);
    fmShimmed_temp = fmShimmed(:,:,iS);
    fmVar.slice_initial(iS) = std(fm_temp(logical(mask(:,:,iS))));
    fmVar.slice_shimmed(iS) = std(fmShimmed_temp(logical(mask(:,:,iS))));
end

% Calculate global standard deviation
fm_cut = fm(:,:,shimSli); 
fmShimmed_cut = fmShimmed(:,:,shimSli);
fmVar.global_initial = std(fm_cut(logical(mask(:,:,shimSli))));
fmVar.global_shimmed = std(fmShimmed_cut(logical(mask(:,:,shimSli))));
fmVar.globalMean_initial = mean(fm_cut(logical(mask(:,:,shimSli))));
fmVar.globalMean_shimmed = mean(fmShimmed_cut(logical(mask(:,:,shimSli))));

if max(basis) < 5 % First order shimming
    SOut = zeros(4,size(S,2));
elseif max(basis) < 10 % Second order shimming
    SOut = zeros(9,size(S,2));
elseif max(basis) < 17 % Third order shimming
    SOut = zeros(16,size(S,2));
end
SOut(basis,:) = S;

end

