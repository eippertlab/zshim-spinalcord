function [ SH_basis ] = compute_SH_basis( coords,basis )
% Compute spherical harmonic basis function values at input coordinates
% 
% INPUT:
% coords    coordinates for SH calculation [X Y Z]
% basis     selected basis functions (up to 3rd order SH)
%
% OUTPUT:
% maps      cell structure with matrices that each contain the values of
%           one SH basis function at the coordinates defined in 'array'
%
% Copyright Johanna Vannesjo, FMRIB, Oxford 2015

X = coords(:,1);
Y = coords(:,2);
Z = coords(:,3);

nCoords = size(coords,1);
nBasis = length(basis);
SH_basis = zeros(nCoords,nBasis);

for iB = 1:nBasis
    switch basis(iB)
        case 1  % Z0
            SH_basis(:,iB) = ones(nCoords,1);  
        case 2  % X
            SH_basis(:,iB) = X;             
        case 3  % Y
            SH_basis(:,iB) = Y;               
        case 4  % Z
            SH_basis(:,iB) = Z;                
        case 5  % S2 (XY)
            SH_basis(:,iB) = X.*Y;              
        case 6  % ZY
            SH_basis(:,iB) = Z.*Y;
        case 7  % Z2
            SH_basis(:,iB) = 2*Z.*Z  - ( X.*X + Y.*Y);
        case 8  % ZX
            SH_basis(:,iB) = Z.*X;
        case 9  % C2 (X2-Y2)
            SH_basis(:,iB) = X.*X - Y.*Y;
        case 10 % Y3
            SH_basis(:,iB) = 3*Y.*X.*X - Y.*Y.*Y;
        case 11 % ZS2
            SH_basis(:,iB) = X.*Y.*Z;
        case 12 % Z2Y
            SH_basis(:,iB) = 5*Y.*Z.*Z - Y.*(X.*X + Y.*Y + Z.*Z);
        case 13 % Z3
            SH_basis(:,iB) = 2*Z.*Z.*Z - 3*Z.*(X.*X + Y.*Y);
        case 14 % Z2X
            SH_basis(:,iB) = 5*X.*Z.*Z - X.*(X.*X + Y.*Y + Z.*Z);
        case 15 % ZC2
            SH_basis(:,iB) = Z.*(X.*X - Y.*Y);
        case 16 % X3
            SH_basis(:,iB) = X.*X.*X - 3*X.*Y.*Y;
    end
end

end

