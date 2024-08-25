

function [Kuu_eb] = KBC_u_Nitsche_Mulmaterial_ym(AnalysisType, gamma_D_u)
% To apply essential boundary conditions with Nitsche's method 
% Analytical displacements are imposed on all boundaries

global MeshShape NumElem 
global GDOF_U
global ManiElems PhyPatches
global Mats

% figure

if strcmpi(MeshShape,'BiotQ9Q4') || strcmpi(MeshShape,'Q9')
    MeshShape_u = 'Q9';  
    ngauss = 4;
elseif strcmpi(MeshShape,'BiotQ4Q4')
    MeshShape_u = 'Q4';  
    ngauss = 2; % number of gauss points in 1 direction 
elseif strcmpi(MeshShape,'BiotIRT3_RIRT3_R')
    MeshShape_u = 'T3';  
    ngauss = 2; % number of gauss points in 1 direction 
end


Kuu_eb = zeros(GDOF_U, GDOF_U);

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


%% boundary data 
% number of dirichlet boundaries
n_dbc = 3;  

% boundary data; top (yMax), bottom (yMin), left (xMin) and right (xMax) surfaces 
check_x_value = {yMin, xMax, xMin};
check_dof_value = {2, 1, 1};
n_vec_value = { [0; -1], [1; 0], [-1; 0] };
check_x_opt_value = { 1, 2, 1};   % 2 for maximum and 1 for minimum 
constr_dof_value = { [1 2], [1 2], [1 2] };

dbc = struct( 'check_x', check_x_value, 'check_dof', check_dof_value, ...
    'n_vec', n_vec_value, 'check_x_opt', check_x_opt_value, ...
    'constr_dof', constr_dof_value );


%% elements list 
%dofs_bc = [];   % list of DOFs related to the essential boundary 
eles_bc = [];   % elements related to the essential boundary condition

for ibc = 1 : n_dbc
    
    check_x = dbc(ibc).check_x;
    check_dof = dbc(ibc).check_dof;
    
    for iele = 1 : NumElem  
        
        ME = ManiElems(iele);

        xy_ele = ME.xVertex;
        xy_ele = [xy_ele; xy_ele(1, :) ];
    
        bc_flag = 0;
        for inode = 1 : size(xy_ele, 1) - 1
            if abs(xy_ele(inode, check_dof) - check_x) < GeoTol && ...
                abs(xy_ele(inode + 1, check_dof) - check_x) < GeoTol
                bc_flag = 1;
                break;
            end
        end
        
        if bc_flag == 0
            continue;
        end
        
        %Dofs = ME.DOF_u;
        
        eles_bc = union(eles_bc, iele);
        %dofs_bc = union(dofs_bc, Dofs);
    end  % iele 
end  % ibc

% fprintf("finished computing dofs_bc.\n");
% fprintf("length(eles_bc) = %d\n",length(eles_bc));


%% boundary terms in Nitsche's method 
for ibc = 1 : n_dbc
    
    len_bc = 0;
    
    n_vec = dbc(ibc).n_vec;
    nx = n_vec(1);  ny = n_vec(2);
    % matrix of outward unit normal for the computing of tractions
    nMat = [ nx   0    ny;
             0    ny   nx ];
    
    constr_dof = dbc(ibc).constr_dof;
    % outward unit normal matrix corresponding to the constrained DOFs
    nMat_c = nMat(constr_dof, :);  
    
    check_x = dbc(ibc).check_x;
    check_dof = dbc(ibc).check_dof;
    
    for jele = 1 : length(eles_bc)
        
        iele = eles_bc(jele);
        
        ME = ManiElems(iele);
        
        mat = Mats(ME.mat);
        E = mat.E;  nu = mat.nu;  

        [DMatU] = ConstitutiveMat(E, nu, AnalysisType);
        
        PPs = PhyPatches(ME.PP);
        xPPs_u = zeros(length(PPs), 2);
        for ipp = 1 : length(PPs)
            xPPs_u(ipp, :) = PPs(ipp).xNode;
        end
        
        h = max(xPPs_u(:, 1)) - min(xPPs_u(:, 1)); 
        
        Dofs = ME.DOF_u;
        
        xy_ele = ME.xVertex;
        xy_ele = [xy_ele; xy_ele(1, :) ];

        for inode = 1 : size(xy_ele, 1) - 1

            if abs(xy_ele(inode, check_dof) - check_x) > GeoTol || ...
                abs(xy_ele(inode + 1, check_dof) - check_x) > GeoTol
                continue;
            end

            len_seg = norm(xy_ele(inode, :) - xy_ele(inode + 1, :));
            J = len_seg / 2;    % Jacobian for the 1D gauss integration 

            len_bc = len_bc + len_seg;

            for igauss = 1 : length(wgt_1d) 

                lx = lxs_1d(igauss);   % local coorinates
                % global coordinates of the integration point
                gxy = xy_ele(inode, :) + ...
                    ( xy_ele(inode + 1, :) - xy_ele(inode, :) ) / 2 * ...
                    ( lx + 1 );

                % shape function matrix 
                [N] = NMatNMM2D_1(xPPs_u, gxy, MeshShape_u); 
                
                % shape function matrix corresponding to the constrained DOFs
                N_c = N(constr_dof, :);
                
                [B] = BMatNMM2D_1(xPPs_u, gxy, MeshShape_u);
                
                % matrix corresponding to traction on the boundary 
                t_mat = nMat_c * DMatU * B;
                
                Kuu_eb(Dofs, Dofs) = Kuu_eb(Dofs, Dofs) - ...
                    t_mat.' * N_c * J * wgt_1d(igauss);
                Kuu_eb(Dofs, Dofs) = Kuu_eb(Dofs, Dofs) - ...
                    N_c.' * t_mat * J * wgt_1d(igauss);
                Kuu_eb(Dofs, Dofs) = Kuu_eb(Dofs, Dofs) + ...
                    gamma_D_u * (N_c.' * N_c) * J * wgt_1d(igauss) / h;

            end  % iguass 
        end  % inode 
    end  % jele 
    
    fprintf("ibc = %d, len_bc = %16.8e\n", ibc, len_bc);
%     if abs(len_bc - 1) > 1e-6
%         fprintf("wrong solid boundary length, ibc = %d\n", ibc);
%     end
    
end % ibc
                
                
    

