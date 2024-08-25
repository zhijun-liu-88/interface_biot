
function [N] = NMatNMM2D_1(xPPs, gxy, MeshShape)
% to compute N matrix at a point in a rectangular or triangular element 
% Input
% xPPs: coordinates of mathematical nodes
% gxy: coordinates of the point 
% Output 
% N: shape function matrix 

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
    
    % shape function corresponding to each node
    [N1D, ~]=lagrange_basis(MeshShape, LocalCoord);
    
    ndim = 2;   % spatical dimension
    % shape function matrix 
    N = kron(N1D.', eye(ndim));
    
elseif strcmpi(MeshShape, 'IRT3_R') || strcmpi(MeshShape, 'ET3') || ...
        strcmpi(MeshShape, 'T3')
    x_tri = xPPs(:, 1);
    y_tri = xPPs(:, 2);
    gx = gxy(1);  gy = gxy(2);
    [N]=NMat_2d_Tri3(x_tri, y_tri, gx, gy);
else
    error("Only quadrilateral and triangular elements are implemented.");
end 