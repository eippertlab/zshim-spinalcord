function [ coords ] = compute_map_coords( nArray,afftrans,coordSystem )
% Compute Cartesian coordinates of array of given size
% 
% INPUT:
% nArray    size of array
% afftrans  affine transformation matrix [m]
%
% OUTPUT:
% coords    coordinates of the array
%
% Copyright Johanna Vannesjo, FMRIB, Oxford 2015

% INFO
% Coords in afftrans matrix, read from nifti-data is defined as:
% L --> R, P --> A, F --> H
% LR = -X, PA = Y, FH = -Z

if nargin < 3
    coordSystem = 'xyz';
end
if nargin<2 || isempty(afftrans)
    center = floor(nArray/2)+1;
    afftrans = [...
        1 0 0 -center(1); ...
        0 1 0 -center(2); ...
        0 0 1 -center(3)]; 
end
if length(afftrans) == 3
    center = floor(nArray/2)+1;
    afftrans = [...
        afftrans(1) 0 0 -afftrans(1)*center(1); ...
        0 afftrans(2) 0 -afftrans(2)*center(2); ...
        0 0 afftrans(3) -afftrans(3)*center(3)];     
end
% 
% vRow =0:nArray(1)-1;
% vCol = 0: nArray(2)-1;
% vMat = 0:nArray(3)-1;

vRow =1:nArray(1);
vCol = 1: nArray(2);
vMat = 1:nArray(3);


[Col,Row,Mat] = meshgrid(vCol,vRow,vMat);
coords = [Row(:) Col(:) Mat(:)];
coords = (afftrans*[coords'; ones(size(Col(:)'))])';
% 4th column of afftrans indicates offset of center of voxel (1,1,1)
if strcmp(coordSystem,'xyz')
    coords(:,1) = -coords(:,1); %flip x-axis
    coords(:,3) = -coords(:,3); %flip z-axis
end

end

