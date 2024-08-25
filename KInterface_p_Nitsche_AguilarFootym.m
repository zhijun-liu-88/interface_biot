

function [Kpp_I] = KInterface_p_Nitsche_AguilarFootym(gamma_I_p)
% To compute the submatrix of the coefficient matrix, Kpp, corresponding to
% the interface term 

global MeshShape NumElem 
global GDOF_P
global ManiElems PhyPatches
global nInterface Interfaces
global Mats

% figure

if strcmpi(MeshShape,'BiotQ9Q4')  
    MeshShape_p = 'Q4';
    ngauss = 2;
elseif strcmpi(MeshShape,'BiotQ4Q4') || strcmpi(MeshShape,'Q4') || ... 
        strcmpi(MeshShape,'Quad')
    MeshShape_p = 'Q4';
    ngauss = 2; % number of gauss points in 1 direction 
elseif strcmpi(MeshShape,'BiotIRT3_RIRT3_R') || strcmpi(MeshShape,'BiotIRT3_LIRT3_L') || ...
        strcmpi(MeshShape,'IRT3_R') || strcmpi(MeshShape,'IRT3_L') 
    MeshShape_p = 'T3';
    ngauss = 2; % number of gauss points in 1 direction 
end

Kpp_I = zeros(GDOF_P, GDOF_P);

% weights and local coordinates 
[wgt_1d, lxs_1d] = CoorWeight_GL(ngauss, 1);

for iface = 1 : nInterface
    
    Face = Interfaces(iface);
    % 'L' and 'R' represent the left and right side of the face 
    IndexElem_L = Face.ManiElem(1);
    IndexElem_R = Face.ManiElem(2);
    
    ME_L = ManiElems(IndexElem_L);     ME_R = ManiElems(IndexElem_R);
    
    mat_L = Mats(ME_L.mat);       mat_R = Mats(ME_R.mat);  
    
    kappa_L = mat_L.kappa;   kappa_R = mat_R.kappa;
    % weights for fluid velocity 
    %cp_L = kappa_R/(kappa_L + kappa_R);   
    %cp_R = kappa_L/(kappa_L + kappa_R);
    cp_L = 0.5;   cp_R = 0.5;
    
    Dofs_L = ME_L.DOF_p;
    Dofs_R = ME_R.DOF_p;
    
    PPs_L = PhyPatches(ME_L.PP); 
    PPs_R = PhyPatches(ME_R.PP);
    
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
    
    xVertex_L = ME_L.xVertex;
    Face.compute_n_2D(xVertex_L);  % to compute the unit normal vector 
    
    n_vec = reshape(Face.n, 1, 2);
    
    DMatP_L = kappa_L * eye(2);  DMatP_R = kappa_R * eye(2);
    
    nD_L = n_vec * DMatP_L;  nD_R = n_vec * DMatP_R;
    
    K_LL = zeros(length(Dofs_L), length(Dofs_L));
    K_LR = zeros(length(Dofs_L), length(Dofs_R));
    K_RL = zeros(length(Dofs_R), length(Dofs_L));
    K_RR = zeros(length(Dofs_R), length(Dofs_R));
    
    xySeg = Face.xVertex;
    len_seg = norm(xySeg(1, :) - xySeg(2, :)); % length 
    J = len_seg / 2;    % Jacobian for the 1D gauss integration 
    
    h = len_seg;
    
    for igauss = 1 : length(wgt_1d) 
        lx = lxs_1d(igauss);   % local coorinates
        % global coordinates of the integration point
        gxy = xySeg(1, :) + ( xySeg(2, :)- xySeg(1, :) ) / 2 * (lx + 1);
         
        % shape function and strain matrices for displacement 
        [Nu_L] = NMatNMM2D_1(xPPs_L, gxy, MeshShape_p); 
        [Nu_R] = NMatNMM2D_1(xPPs_R, gxy, MeshShape_p); 
        [Bu_L] = BMatNMM2D_1(xPPs_L, gxy, MeshShape_p);
        [Bu_R] = BMatNMM2D_1(xPPs_R, gxy, MeshShape_p);
        
        % shape function matrices for pressure 
        Np_L = Nu_L(1, 1 : 2 : end);
        Np_R = Nu_R(1, 1 : 2 : end);
        
        Bp_L = [ Bu_L(1, 1 : 2 : end);
                 Bu_L(2, 2 : 2 : end) ];
        Bp_R = [ Bu_R(1, 1 : 2 : end);
                 Bu_R(2, 2 : 2 : end) ];
        
        % matrix for the traction 
        wmat_L = cp_L * nD_L * Bp_L;     wmat_R = cp_R * nD_R * Bp_R;
         
        K_LL = K_LL - wmat_L.' * Np_L * J * wgt_1d(igauss);
        K_LR = K_LR + wmat_L.' * Np_R * J * wgt_1d(igauss);
        K_RL = K_RL - wmat_R.' * Np_L * J * wgt_1d(igauss);
        K_RR = K_RR + wmat_R.' * Np_R * J * wgt_1d(igauss);
        
        K_LL = K_LL - Np_L.' * wmat_L * J * wgt_1d(igauss);
        K_LR = K_LR - Np_L.' * wmat_R * J * wgt_1d(igauss);
        K_RL = K_RL + Np_R.' * wmat_L * J * wgt_1d(igauss);
        K_RR = K_RR + Np_R.' * wmat_R * J * wgt_1d(igauss);
        
        K_LL = K_LL + gamma_I_p * (Np_L.' * Np_L) * J * wgt_1d(igauss) / h;
        K_LR = K_LR - gamma_I_p * Np_L.' * Np_R * J * wgt_1d(igauss) / h;
        K_RL = K_RL - gamma_I_p * Np_R.' * Np_L * J * wgt_1d(igauss) / h;
        K_RR = K_RR + gamma_I_p * (Np_R.' * Np_R) * J * wgt_1d(igauss) / h;
    end  % igauss 
    
    Kpp_I(Dofs_L, Dofs_L) = Kpp_I(Dofs_L, Dofs_L) + K_LL;
    Kpp_I(Dofs_L, Dofs_R) = Kpp_I(Dofs_L, Dofs_R) + K_LR;
    Kpp_I(Dofs_R, Dofs_L) = Kpp_I(Dofs_R, Dofs_L) + K_RL;
    Kpp_I(Dofs_R, Dofs_R) = Kpp_I(Dofs_R, Dofs_R) + K_RR;
end  % iface 


                
                
    

