
function [fu, fp] = LoadUnitMulmaterial_1(Sigma0)
% right-hand side of the Biot's consolidation problem caused by source and
% natural boundary terms. 
% Essential boundary terms are enforced for the displacements on all
% boundary. No natual boundary exists for the solid. Pressure on top
% boundary are prescribed. All other boundaries are impervious

global MeshShape 
global NumElem 
global GDOF_U GDOF_P
global ManiElems PhyPatches
global Mats


% E1 = Mats(1).E;  kappa_1 = Mats(1).kappa;  M1 = Mats(1).M;  alpha_1 = Mats(1).alpha;
% E2 = Mats(2).E;  kappa_2 = Mats(2).kappa;  M2 = Mats(2).M;  alpha_2 = Mats(2).alpha;

fu = zeros(GDOF_U, 1);
fp = zeros(GDOF_P, 1); 


if strcmpi(MeshShape,'BiotQ9Q4') || strcmpi(MeshShape,'Q9')
    MeshShape_u = 'Q9';  MeshShape_p = 'Q4';
elseif strcmpi(MeshShape,'BiotQ4Q4')
    MeshShape_u = 'Q4';  MeshShape_p = 'Q4';
elseif strcmpi(MeshShape,'BiotIRT3_RIRT3_R')
    MeshShape_u = 'T3';  MeshShape_p = 'T3';
end

%% source terms 
% t = current_t;

ngauss = 7; % number of gauss points in 1 direction 
% ndim = 2;   % dimension 
GeoTol = 1.e-13;


% wights and local coordinates 
[wgt_1d, lxs_1d] = CoorWeight_GL(ngauss, 1);


%% maximum and minimum coordinates of the model 
xVertex = ManiElems(1).xVertex;
xMin = min(xVertex(:, 1));
xMax = max(xVertex(:, 1));
yMin = min(xVertex(:, 2));
yMax = max(xVertex(:, 2));
for iele = 2 : NumElem 
    xVertex = ManiElems(iele).xVertex;
    xMinElem = min(xVertex(:, 1));
    xMaxElem = max(xVertex(:, 1));
    yMinElem = min(xVertex(:, 2));
    yMaxElem = max(xVertex(:, 2));
    
    xMin = min(xMin, xMinElem);
    xMax = max(xMax, xMaxElem);
    yMin = min(yMin, yMinElem);
    yMax = max(yMax, yMaxElem);
end

%% Uniform pressure on top boundary 
MaxXLoad = 0.8;   MinXLoad = -0.8;
LengthLoad = 0;
Traction = [ 0;  -Sigma0 ];


for iele = 1 : NumElem
    
    ME = ManiElems(iele);
    
    xy_ele = ME.xVertex;
    xy_ele = [xy_ele; xy_ele(1, :) ];
    
    MaxYElem = max(xy_ele(:, 2));
    MaxXElem = max(xy_ele(:, 1));
    MinXElem = min(xy_ele(:, 1));

    if abs(MaxYElem - yMax) > GeoTol
        continue;
    end
    
    if MinXElem > MaxXLoad + GeoTol || MaxXElem < MinXLoad - GeoTol
        continue;
    end
    
    
    PPs = PhyPatches(ME.PP);
    
%     DOFs_U = ME.DOF_u;
%     DOFs_P = ME.DOF_p;
%     
%     % elemental right-hand vector
%     FeU = zeros(length(DOFs_U), 1);
%     FeP = zeros(length(DOFs_P), 1);
    
    xPPs_u = zeros(length(PPs), 2);
    for ipp = 1 : length(PPs)
        xPPs_u(ipp, :) = PPs(ipp).xNode;
    end
    
    Dofs = ME.DOF_u;

    FeU = zeros(length(Dofs), 1);
    
    for inode = 1 : size(xy_ele, 1) - 1
        if abs(xy_ele(inode, 2) - yMax) >= GeoTol  || ...
                abs(xy_ele(inode + 1, 2) - yMax) >= GeoTol
            continue;
        end
        
        MinXEdge = min(xy_ele(inode : inode+1, 1));
        MaxXEdge = max(xy_ele(inode : inode+1, 1));
        
        if MinXEdge > MaxXLoad + GeoTol
            continue;
        end
        if MaxXEdge < MinXLoad - GeoTol
            continue;
        end
        
        MinXSeg = max(MinXLoad, MinXEdge);
        MaxXSeg = min(MaxXLoad, MaxXEdge);
        
        len_seg = MaxXSeg - MinXSeg;
        if len_seg < GeoTol
            continue;
        end
        
        J = len_seg / 2;    % Jacobian for the 1D gauss integration
        
        LengthLoad= LengthLoad + len_seg;
    
    
        for igauss = 1 : length(wgt_1d) 
            
            lx = lxs_1d(igauss);   % local coordinates
            % global coordinates of the integration point
            gxy = [ MinXSeg + (MaxXSeg - MinXSeg) / 2 * (lx + 1), yMax ];

            % shape function matrix for the solid 
            [Nu] = NMatNMM2D_1(xPPs_u, gxy, MeshShape_u); 
            
            FeU = FeU + Nu.' * Traction * J * wgt_1d(igauss);
            
        end  % l : loop of integration points 
    end  % k : loop of quadrilaterals
    
    fu(Dofs, 1) = fu(Dofs, 1) + FeU;
    
end % iele

fprintf("LengthLoad = %10.4e\n", LengthLoad);




    