

function [Kuu_I] = KInterface_u_Nitsche_AguilarFootym(AnalysisType, gamma_I_u)
% To compute the submatrix of the coefficient matrix, Kuu, corresponding to
% the interface term 

global MeshShape NumElem 
global GDOF_U
global ManiElems PhyPatches
global nInterface Interfaces
global Mats


% figure

if strcmpi(MeshShape,'BiotQ9Q4') || strcmpi(MeshShape,'Q9')
    MeshShape_u = 'Q9';  
    ngauss = 3;
elseif strcmpi(MeshShape,'BiotQ4Q4')
    MeshShape_u = 'Q4';  
    ngauss = 2; % number of gauss points in 1 direction 
elseif strcmpi(MeshShape,'BiotIRT3_RIRT3_R')
    MeshShape_u = 'T3';  
    ngauss = 2; % number of gauss points in 1 direction 
end

Kuu_I = zeros(GDOF_U, GDOF_U);

% weights and local coordinates 
[wgt_1d, lxs_1d] = CoorWeight_GL(ngauss, 1);

for iface = 1 : nInterface
    
    Face = Interfaces(iface);
    % 'L' and 'R' represent the left and right side of the face 
    IndexElem_L = Face.ManiElem(1);
    IndexElem_R = Face.ManiElem(2);
    
    ME_L = ManiElems(IndexElem_L);     ME_R = ManiElems(IndexElem_R);
    
    mat_L = Mats(ME_L.mat);       mat_R = Mats(ME_R.mat);  
    
    E_L = mat_L.E;   nu_L = mat_L.nu;
    E_R = mat_R.E;   nu_R = mat_R.nu;
%     lambda_L = E_L * nu_L / (1 + nu_L) / (1 - 2 * nu_L);
%     lambda_R = E_R * nu_R / (1 + nu_R) / (1 - 2 * nu_R);
%     G_L = E_L / 2 / (1 + nu_L);
%     G_R = E_R / 2 / (1 + nu_R);
    
    %eta_L = lambda_L + 2 * G_L;  eta_R = lambda_R + 2 * G_R;
    
    %weights for effective stress 
    %cu_L = eta_R/(eta_L + eta_R);   cu_R = eta_L/(eta_L + eta_R);
    cu_L = 0.5;   cu_R = 0.5;
    
%     kappa_L = mat_L.kappa;   kappa_R = mat_R.kappa;
%     % weights for fluid velocity 
%     cp_L = kappa_R/(kappa_L + kappa_R);   
%     cp_R = kappa_L/(kappa_L + kappa_R);
    
    PPs_L = PhyPatches(ME_L.PP); 
    PPs_R = PhyPatches(ME_R.PP);
    
    Dofs_L = ME_L.DOF_u;
    Dofs_R = ME_R.DOF_u;
    
    xPPs_L = zeros(length(PPs_L), 2);
    for ipp = 1 : length(PPs_L)
        xPPs_L(ipp, :) = PPs_L(ipp).xNode;
    end
    
    xPPs_R = zeros(length(PPs_R), 2);
    for ipp = 1 : length(PPs_R)
        xPPs_R(ipp, :) = PPs_R(ipp).xNode;
    end
    
    xVertex_L = ME_L.xVertex;
    Face.compute_n_2D(xVertex_L);  % to compute the unit normal vector 
    
    nx = Face.n(1);  ny = Face.n(2);
    % matrix of outward unit normal for the computing of tractions
    nMat = [ nx   0    ny;
             0    ny   nx ];
         
    [DMatU_L] = ConstitutiveMat(E_L, nu_L, AnalysisType);
    [DMatU_R] = ConstitutiveMat(E_R, nu_R, AnalysisType);
    
    nD_L = nMat * DMatU_L;
    nD_R = nMat * DMatU_R;

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
         
        % shape function and strain matrices
        [N_L] = NMatNMM2D_1(xPPs_L, gxy, MeshShape_u); 
        [N_R] = NMatNMM2D_1(xPPs_R, gxy, MeshShape_u); 
        [B_L] = BMatNMM2D_1(xPPs_L, gxy, MeshShape_u);
        [B_R] = BMatNMM2D_1(xPPs_R, gxy, MeshShape_u);
        
        % matrix for the traction 
        tmat_L = cu_L * nD_L * B_L;     tmat_R = cu_R * nD_R * B_R;
         
        K_LL = K_LL - tmat_L.' * N_L * J * wgt_1d(igauss);
        K_LR = K_LR + tmat_L.' * N_R * J * wgt_1d(igauss);
        K_RL = K_RL - tmat_R.' * N_L * J * wgt_1d(igauss);
        K_RR = K_RR + tmat_R.' * N_R * J * wgt_1d(igauss);
        
        K_LL = K_LL - N_L.' * tmat_L * J * wgt_1d(igauss);
        K_LR = K_LR - N_L.' * tmat_R * J * wgt_1d(igauss);
        K_RL = K_RL + N_R.' * tmat_L * J * wgt_1d(igauss);
        K_RR = K_RR + N_R.' * tmat_R * J * wgt_1d(igauss);
        
        K_LL = K_LL + gamma_I_u * (N_L.' * N_L) * J * wgt_1d(igauss) / h;
        K_LR = K_LR - gamma_I_u * N_L.' * N_R * J * wgt_1d(igauss) / h;
        K_RL = K_RL - gamma_I_u * N_R.' * N_L * J * wgt_1d(igauss) / h;
        K_RR = K_RR + gamma_I_u * (N_R.' * N_R) * J * wgt_1d(igauss) / h;
    end  % igauss 
    
    Kuu_I(Dofs_L, Dofs_L) = Kuu_I(Dofs_L, Dofs_L) + K_LL;
    Kuu_I(Dofs_L, Dofs_R) = Kuu_I(Dofs_L, Dofs_R) + K_LR;
    Kuu_I(Dofs_R, Dofs_L) = Kuu_I(Dofs_R, Dofs_L) + K_RL;
    Kuu_I(Dofs_R, Dofs_R) = Kuu_I(Dofs_R, Dofs_R) + K_RR;
end  % iface 


                
                
    

