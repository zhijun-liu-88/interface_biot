
function [Kuu_ghost] = GhostPenalty2DSolid_Quadratic1(gamma_g_u)
% to compute ghost penalty terms of the solid with the linear/quadratic approximation

global GDOF_U
global MeshShape
global nGhostFace GhostFaces 
global ManiElems PhyPatches


if strcmpi(MeshShape,'BiotQ9Q4') || strcmpi(MeshShape,'Q9')
    MeshShape_u = 'Q9';  
    interp_order = 2;
    ngauss = 4;
elseif strcmpi(MeshShape,'BiotQ4Q4') || strcmpi(MeshShape,'Quad') || strcmpi(MeshShape,'Q4')
    MeshShape_u = 'Q4';  
    interp_order = 1;
    ngauss = 2; % number of gauss points in 1 direction 
elseif strcmpi(MeshShape,'BiotIRT3_RIRT3_R') || strcmpi(MeshShape,'IRT3_R') || ...
        strcmpi(MeshShape,'IRT3_L') || strcmpi(MeshShape,'ET3') 
    MeshShape_u = 'T3';  
    interp_order = 1;
    ngauss = 2; % number of gauss points in 1 direction 
end


% wights and local coordinates 
[wgt_1d, lxs_1d] = CoorWeight_GL(ngauss, 1);

% initialization of the ghost penalty matrix for solid 
Kuu_ghost = zeros(GDOF_U, GDOF_U);
for iface = 1 : nGhostFace
    Face = GhostFaces(iface);
    % 'L' and 'R' represent the left and right side of the face 
    IndexElem_L = Face.ManiElem(1);
    IndexElem_R = Face.ManiElem(2);
    
    ME_L = ManiElems(IndexElem_L);
    ME_R = ManiElems(IndexElem_R);
    
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
    
    if strcmpi(MeshShape_u,'Q9') || strcmpi(MeshShape_u,'Q4')
        xVertex_L = xPPs_L(1 : 4, :); % coordintes of vertices of the left mathematical element 
    elseif strcmpi(MeshShape_u,'T6') || strcmpi(MeshShape_u,'T3')
        xVertex_L = xPPs_L(1 : 3, :);
    end 
    
    Face.compute_n_2D(xVertex_L);  % to compute the unit normal vector 
    nVector = reshape(Face.n, 1, 2); % unit normal vector 
    nMat = nVector.' * nVector;
    
    % coordiantes of two endpoints 
    xySeg = Face.xVertex;
    len_seg = norm(xySeg(1, :) - xySeg(2, :)); % length 
    J = len_seg / 2;    % Jacobian for the 1D gauss integration 
    
    h = len_seg;
    
    for iorder = 1 : interp_order
        h_exp = 2 * iorder - 1; % exponent of h 
        
        
        if iorder == 2 
            nx = nVector(1);  ny = nVector(2);
            nnVector = [ nx^2  ny^2  nx*ny ];
            nMat = nnVector.' * nnVector;
        end
    
        % matrix corresponding to the 1st order interpolation 
        K_LL = zeros(length(Dofs_L), length(Dofs_L));
        K_LR = zeros(length(Dofs_L), length(Dofs_R));
        K_RL = zeros(length(Dofs_R), length(Dofs_L));
        K_RR = zeros(length(Dofs_R), length(Dofs_R));
        for igauss = 1 : length(wgt_1d) 
            lx = lxs_1d(igauss);   % local coordinates
            % global coordinates of the integration point
            gxy = xySeg(1, :) + ( xySeg(2, :)- xySeg(1, :) ) / 2 * (lx + 1);

            if iorder == 1
                [Bu_L] = BMatNMM2D_1(xPPs_L, gxy, MeshShape_u); % strain matrix of 1st element
                [Bu_R] = BMatNMM2D_1(xPPs_R, gxy, MeshShape_u); % strain matrix of 2nd element 

                % matrix corresponding to gradient of u_x of 1st element 
                Bx_L = [ Bu_L(1, :);  Bu_L(2, 2 : end), 0 ];
                By_L = [ 0, Bu_L(1, 1 : end-1);  Bu_L(2, :) ];   
                % matrix corresponding to gradient of u_x of 2st element 
                Bx_R = [ Bu_R(1, :);  Bu_R(2, 2 : end), 0 ];
                By_R = [ 0, Bu_R(1, 1 : end-1);  Bu_R(2, :) ];
            elseif iorder == 2  
                if strcmpi(MeshShape_u,'Q9')
                    % matrix corresponding to 2nd-order gradient of u_x of 1st element 
                    [Bx_L] = d2Ndx2Q9(xPPs_L, gxy);
                    % matrix corresponding to 2nd-order gradient of u_x of 2nd element 
                    [Bx_R] = d2Ndx2Q9(xPPs_R, gxy);
                else
                    error("element not implemented.");
                end
                
                temp_mat =[ 0; 0; 0 ];
                % matrix corresponding to 2nd-order gradient of u_y of 1st element 
                By_L = [ temp_mat, Bx_L(:, 1 : end - 1) ];
                % matrix corresponding to 2nd-order gradient of u_y of 2st element 
                By_R = [ temp_mat, Bx_R(:, 1 : end - 1) ];
            end  

            % iorder-th order derivative terms 
            % formulation in the papers: 1. fictitious domain methods using cut elements:
            % iii. a stabilized nitsche method for stokes’problem 
            % 2. High-order cut finite elements for the elastic wave equation
            K_LL = K_LL + (Bx_L.' * nMat * Bx_L + By_L.' * nMat * By_L) * J * wgt_1d(igauss);
            K_LR = K_LR - (Bx_L.' * nMat * Bx_R + By_L.' * nMat * By_R) * J * wgt_1d(igauss);
            K_RL = K_RL - (Bx_R.' * nMat * Bx_L + By_R.' * nMat * By_L) * J * wgt_1d(igauss);
            K_RR = K_RR + (Bx_R.' * nMat * Bx_R + By_R.' * nMat * By_R) * J * wgt_1d(igauss);

            % 1st-order derivative terms 
            % another formulation, from the paper: Fictitious domain finite 
            % element methods using cut elements: II. A stabilized Nitsche method
    %         K_11 = K_11 + (Bx_1.' * Bx_1 + By_1.' * By_1) * J * wgt_1d(igauss);
    %         K_12 = K_12 - (Bx_1.' * Bx_2 + By_1.' * By_2) * J * wgt_1d(igauss);
    %         K_21 = K_21 - (Bx_2.' * Bx_1 + By_2.' * By_1) * J * wgt_1d(igauss);
    %         K_22 = K_22 + (Bx_2.' * Bx_2 + By_2.' * By_2) * J * wgt_1d(igauss);


        end  % igauss 
    
        Kuu_ghost(Dofs_L, Dofs_L) = Kuu_ghost(Dofs_L, Dofs_L) + K_LL * h^h_exp;
        Kuu_ghost(Dofs_L, Dofs_R) = Kuu_ghost(Dofs_L, Dofs_R) + K_LR * h^h_exp;
        Kuu_ghost(Dofs_R, Dofs_L) = Kuu_ghost(Dofs_R, Dofs_L) + K_RL * h^h_exp;
        Kuu_ghost(Dofs_R, Dofs_R) = Kuu_ghost(Dofs_R, Dofs_R) + K_RR * h^h_exp;
    end  % iorder 
end  % iface 

Kuu_ghost = Kuu_ghost * gamma_g_u;
    
        
        
        
        
        
        
        