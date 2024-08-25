

function [Kup_I] = KInterface_up_Nitsche_AguilarFootym()
% To compute the submatrix of the coefficient matrix, Kuu, corresponding to
% the interface term 

global MeshShape NumElem 
global GDOF_U GDOF_P
global ManiElems PhyPatches
global nInterface Interfaces
global Mats

% figure

if strcmpi(MeshShape,'BiotQ9Q4')  
    MeshShape_u = 'Q9';  MeshShape_p = 'Q4';
    ngauss = 3;
elseif strcmpi(MeshShape,'BiotQ4Q4') || strcmpi(MeshShape,'Q4') || ... 
        strcmpi(MeshShape,'Quad')
    MeshShape_u = 'Q4';  MeshShape_p = 'Q4';
    ngauss = 2; % number of gauss points in 1 direction 
elseif strcmpi(MeshShape,'BiotIRT3_RIRT3_R') || strcmpi(MeshShape,'BiotIRT3_LIRT3_L') || ...
        strcmpi(MeshShape,'IRT3_R') || strcmpi(MeshShape,'IRT3_L') 
    MeshShape_u = 'T3';  MeshShape_p = 'T3';
    ngauss = 2; % number of gauss points in 1 direction 
end

Kup_I = zeros(GDOF_U, GDOF_P);

% weights and local coordinates 
[wgt_1d, lxs_1d] = CoorWeight_GL(ngauss, 1);

for iface = 1 : nInterface
    
    Face = Interfaces(iface);
    % 'L' and 'R' represent the left and right side of the face 
    IndexElem_L = Face.ManiElem(1);
    IndexElem_R = Face.ManiElem(2);
    
    ME_L = ManiElems(IndexElem_L);     ME_R = ManiElems(IndexElem_R);
    
    mat_L = Mats(ME_L.mat);       mat_R = Mats(ME_R.mat);  
    alpha_L = mat_L.alpha;   alpha_R = mat_R.alpha;
    
    %kappa_L = mat_L.kappa;   kappa_R = mat_R.kappa;
    % weights for fluid velocity 
    %cp_L = kappa_R/(kappa_L + kappa_R) * alpha_L;   
    %cp_R = kappa_L/(kappa_L + kappa_R) * alpha_R;
    cp_L = 0.5 * alpha_L;   cp_R = 0.5 * alpha_R;
    
    PPs_L = PhyPatches(ME_L.PP); 
    PPs_R = PhyPatches(ME_R.PP);
    
    Dofs_U_L = ME_L.DOF_u;   Dofs_P_L = ME_L.DOF_p;
    Dofs_U_R = ME_R.DOF_u;   Dofs_P_R = ME_R.DOF_p;
    
    xPPs_U_L = zeros(length(PPs_L), 2);
    for ipp = 1 : length(PPs_L)
        xPPs_U_L(ipp, :) = PPs_L(ipp).xNode;
    end
    
    xPPs_U_R = zeros(length(PPs_R), 2);
    for ipp = 1 : length(PPs_R)
        xPPs_U_R(ipp, :) = PPs_R(ipp).xNode;
    end
    
    if strcmpi(MeshShape_p,'Q4')
        xPPs_P_L = xPPs_U_L(1:4, :);  xPPs_P_R = xPPs_U_R(1:4, :);
    elseif strcmpi(MeshShape_p,'T3')
        xPPs_P_L = xPPs_U_L(1:3, :);  xPPs_P_R = xPPs_U_R(1:3, :);
    end
    
    xVertex_L = ME_L.xVertex;
    Face.compute_n_2D(xVertex_L);  % to compute the unit normal vector 

    n_vec = reshape(Face.n, 2, 1);
    
    K_LL = zeros(length(Dofs_U_L), length(Dofs_P_L));
    K_LR = zeros(length(Dofs_U_L), length(Dofs_P_R));
    K_RL = zeros(length(Dofs_U_R), length(Dofs_P_L));
    K_RR = zeros(length(Dofs_U_R), length(Dofs_P_R));
    
    xySeg = Face.xVertex;
    len_seg = norm(xySeg(1, :) - xySeg(2, :)); % length 
    J = len_seg / 2;    % Jacobian for the 1D gauss integration 
    
    for igauss = 1 : length(wgt_1d) 
        lx = lxs_1d(igauss);   % local coorinates
        % global coordinates of the integration point
        gxy = xySeg(1, :) + ( xySeg(2, :)- xySeg(1, :) ) / 2 * (lx + 1);
         
        % shape function matrices
        [Nu_L] = NMatNMM2D_1(xPPs_U_L, gxy, MeshShape_u); 
        [Nu_R] = NMatNMM2D_1(xPPs_U_R, gxy, MeshShape_u); 
        
        [Nu_L1] = NMatNMM2D_1(xPPs_P_L, gxy, MeshShape_p); 
        [Nu_R1] = NMatNMM2D_1(xPPs_P_R, gxy, MeshShape_p); 
        
        Np_L = Nu_L1(1, 1 : 2 : end);
        Np_R = Nu_R1(1, 1 : 2 : end);
        
        K_LL = K_LL + Nu_L.' * n_vec * Np_L * cp_L * J * wgt_1d(igauss);
        K_LR = K_LR + Nu_L.' * n_vec * Np_R * cp_R * J * wgt_1d(igauss);
        K_RL = K_RL - Nu_R.' * n_vec * Np_L * cp_L * J * wgt_1d(igauss);
        K_RR = K_RR - Nu_R.' * n_vec * Np_R * cp_R * J * wgt_1d(igauss);
    end  % igauss 
    
    Kup_I(Dofs_U_L, Dofs_P_L) = Kup_I(Dofs_U_L, Dofs_P_L) + K_LL;
    Kup_I(Dofs_U_L, Dofs_P_R) = Kup_I(Dofs_U_L, Dofs_P_R) + K_LR;
    Kup_I(Dofs_U_R, Dofs_P_L) = Kup_I(Dofs_U_R, Dofs_P_L) + K_RL;
    Kup_I(Dofs_U_R, Dofs_P_R) = Kup_I(Dofs_U_R, Dofs_P_R) + K_RR;
end  % iface 


                
                
    

