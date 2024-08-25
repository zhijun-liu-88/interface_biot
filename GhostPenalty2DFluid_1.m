
function [Kpp_ghost, Mpp_ghost] = GhostPenalty2DFluid_1(gamma_g_p, gamma_g_m)
% to compute ghost penalty terms of the solid with the linear approximation

global MeshShape
global DOF_per_PP_P GDOF_P
global nGhostFace GhostFaces 
global ManiElems PhyPatches


if strcmpi(MeshShape,'BiotQ9Q4') || strcmpi(MeshShape,'BiotQ4Q4') || ...
        strcmpi(MeshShape,'Quad') || strcmpi(MeshShape,'Q4')
    MeshShape_p = 'Q4';  
    interp_order = 1;
    ngauss = 2; % number of gauss points in 1 direction 
elseif strcmpi(MeshShape,'BiotIRT3_RIRT3_R') || strcmpi(MeshShape,'BiotIRT3_LIRT3_L') ...
        || strcmpi(MeshShape,'IRT3_R') || ...
        strcmpi(MeshShape,'IRT3_L') || strcmpi(MeshShape,'ET3') 
    MeshShape_p = 'T3';  
    interp_order = 1;
    ngauss = 2; % number of gauss points in 1 direction 
end


% wights and local coordinates 
[wgt_1d, lxs_1d] = CoorWeight_GL(ngauss, 1);

% initialization of the ghost penalty matrix for solid 
Kpp_ghost = zeros(GDOF_P, GDOF_P);
Mpp_ghost = zeros(GDOF_P, GDOF_P);

for iface = 1 : nGhostFace
    Face = GhostFaces(iface);
    % 'L' and 'R' represent the left and right side of the face 
    IndexElem_L = Face.ManiElem(1);
    IndexElem_R = Face.ManiElem(2);
    
    ME_L = ManiElems(IndexElem_L);
    ME_R = ManiElems(IndexElem_R);
    
    PPs_L = PhyPatches(ME_L.PP);
    PPs_R = PhyPatches(ME_R.PP);
    
    Dofs_L = ME_L.DOF_p;
    Dofs_R = ME_R.DOF_p;
    
    xPPs_L_u = zeros(length(PPs_L), 2);
    for ipp = 1 : length(PPs_L)
        xPPs_L_u(ipp, :) = PPs_L(ipp).xNode;
    end
    
    xPPs_R_u = zeros(length(PPs_R), 2);
    for ipp = 1 : length(PPs_R)
        xPPs_R_u(ipp, :) = PPs_R(ipp).xNode;
    end
    
    if strcmpi(MeshShape,'BiotQ9Q4') || strcmpi(MeshShape,'BiotQ4Q4') || ...
            strcmpi(MeshShape,'Q4') || strcmpi(MeshShape,'Quad')
        xPPs_L = xPPs_L_u(1 : 4, :);
        xPPs_R = xPPs_R_u(1 : 4, :);
    elseif strcmpi(MeshShape,'BiotIRT3_RIRT3_R') || ...
            strcmpi(MeshShape,'BiotIRT3_LIRT3_L') || ...
            strcmpi(MeshShape,'IRT3_L') || strcmpi(MeshShape,'IRT3_R') || ...
            strcmpi(MeshShape,'ET3') 
        xPPs_L = xPPs_L_u(1 : 3, :);
        xPPs_R = xPPs_R_u(1 : 3, :);
    else 
        error("Element not implemented.");
    end
    
    if isempty(Face.n) 
        Face.compute_n_2D(xPPs_L);  % to compute the unit normal vector 
    end
    nVector = reshape(Face.n, 1, 2); % unit normal vector 
    nMat = nVector.' * nVector;
    
    % coordiantes of two endpoints 
    xySeg = Face.xVertex;
    len_seg = norm(xySeg(1, :) - xySeg(2, :)); % length 
    J = len_seg / 2;    % Jacobian for the 1D gauss integration 
    
    h = len_seg;
    
    K_LL = zeros(length(Dofs_L), length(Dofs_L));
    K_LR = zeros(length(Dofs_L), length(Dofs_R));
    K_RL = zeros(length(Dofs_R), length(Dofs_L));
    K_RR = zeros(length(Dofs_R), length(Dofs_R));
    for igauss = 1 : length(wgt_1d) 
        lx = lxs_1d(igauss);   % local coordinates
        % global coordinates of the integration point
        gxy = xySeg(1, :) + ( xySeg(2, :)- xySeg(1, :) ) / 2 * (lx + 1);
        
        [Bu_L] = BMatNMM2D_1(xPPs_L, gxy, MeshShape_p);
        [Bu_R] = BMatNMM2D_1(xPPs_R, gxy, MeshShape_p);
        % B matrix for the fluid in 2D 
        Bp_L = [ Bu_L(1, 1 : 2 : end);
                 Bu_L(2, 2 : 2 : end) ];
        Bp_R = [ Bu_R(1, 1 : 2 : end);
                 Bu_R(2, 2 : 2 : end) ];
        
        % formulation in the paper: fictitious domain methods using cut elements:
        % iii. a stabilized nitsche method for stokes’problem
        K_LL = K_LL + Bp_L.' * nMat * Bp_L * J * wgt_1d(igauss);
        K_LR = K_LR - Bp_L.' * nMat * Bp_R * J * wgt_1d(igauss);
        K_RL = K_RL - Bp_R.' * nMat * Bp_L * J * wgt_1d(igauss);
        K_RR = K_RR + Bp_R.' * nMat * Bp_R * J * wgt_1d(igauss);
        
        % another formulation, from the paper: Fictitious domain finite 
        % element methods using cut elements: II. A stabilized Nitsche method
%         K_11 = K_11 + Bp_1.' * Bp_1 * J * wgt_1d(igauss);
%         K_12 = K_12 - Bp_1.' * Bp_2 * J * wgt_1d(igauss);
%         K_21 = K_21 - Bp_2.' * Bp_1 * J * wgt_1d(igauss);
%         K_22 = K_22 + Bp_2.' * Bp_2 * J * wgt_1d(igauss);
    end  % igauss 
    
    Kpp_ghost(Dofs_L, Dofs_L) = Kpp_ghost(Dofs_L, Dofs_L) + K_LL * h;
    Kpp_ghost(Dofs_L, Dofs_R) = Kpp_ghost(Dofs_L, Dofs_R) + K_LR * h;
    Kpp_ghost(Dofs_R, Dofs_L) = Kpp_ghost(Dofs_R, Dofs_L) + K_RL * h;
    Kpp_ghost(Dofs_R, Dofs_R) = Kpp_ghost(Dofs_R, Dofs_R) + K_RR * h;
    
    Mpp_ghost(Dofs_L, Dofs_L) = Mpp_ghost(Dofs_L, Dofs_L) + K_LL * h^3;
    Mpp_ghost(Dofs_L, Dofs_R) = Mpp_ghost(Dofs_L, Dofs_R) + K_LR * h^3;
    Mpp_ghost(Dofs_R, Dofs_L) = Mpp_ghost(Dofs_R, Dofs_L) + K_RL * h^3;
    Mpp_ghost(Dofs_R, Dofs_R) = Mpp_ghost(Dofs_R, Dofs_R) + K_RR * h^3;
    
end  % iface 

Kpp_ghost = Kpp_ghost * gamma_g_p;
Mpp_ghost = Mpp_ghost * gamma_g_m;
    
        
        
        
        
        
        
        