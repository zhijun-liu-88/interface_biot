
function [B] = BMatNMM2D_1(xPPs, gxy, MeshShape)
% to compute B matrix at a point in a rectangular or triangular element 
% Input
% xPPs: coordinates of mathematical nodes
% gxy: coordinates of the point 
% Output 
% B: strain matrix 

%global MeshShape 
%global DOF_per_PP

if strcmpi(MeshShape, 'Quad') 
    MeshShape = 'Q4';
end

% rectangular element 
if strcmpi(MeshShape, 'Q4') || strcmpi(MeshShape, 'Q9') 
    
    MaxXPP = max(xPPs(:, 1));
    MinXPP = min(xPPs(:, 1));
    MaxYPP = max(xPPs(:, 2));
    MinYPP = min(xPPs(:, 2));
    
    LocalCoord = zeros(1, 2);
    LocalCoord(1) = -1 + ( gxy(1) - MinXPP ) / (MaxXPP - MinXPP) * 2;
    LocalCoord(2) = -1 + ( gxy(2) - MinYPP ) / (MaxYPP - MinYPP) * 2;

    [~, dNdxi] = lagrange_basis(MeshShape, LocalCoord);

    % Jacobian matrix for the transformation in mathematical
    % element
    J = xPPs' * dNdxi;   % 这里用的是物理覆盖的星点

    % derivatives of N with respect to global coordinates 
    dNdx = dNdxi / J;   

    B = zeros(3, size(xPPs, 1));
    for inode = 1 : size(xPPs, 1)
        B(:, 2 * inode - 1 : 2 * inode) = [ ...
            dNdx(inode, 1)   0;
            0                dNdx(inode, 2);
            dNdx(inode, 2)   dNdx(inode, 1)  ];
    end 
elseif strcmpi(MeshShape, 'IRT3_R') || strcmpi(MeshShape, 'ET3')
    x_tri = xPPs(:, 1);
    y_tri = xPPs(:, 2);
    [B]= BMat_2d_Tri3(x_tri, y_tri);
else
    error("Only quadrilateral and triangular elements are implemented.");
end 