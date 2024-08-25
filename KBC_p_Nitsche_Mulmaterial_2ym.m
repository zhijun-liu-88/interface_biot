

function [Kpp_eb] = KBC_p_Nitsche_Mulmaterial_2ym(gamma_D_p)
% To apply essential boundary conditions with Nitsche's method 
% Analytical pressures are imposed on all boundaries


global MeshShape NumElem GDOF_P
global ManiElems PhyPatches
global Mats


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


% figure


Kpp_eb = zeros(GDOF_P, GDOF_P);
% fp_eb = zeros(GDOF_P, 1);

% ndim = 2;   % dimension 

GeoTol = 1.e-13;
 
% wights and local coordinates 
[wgt_1d, lxs_1d] = CoorWeight_GL(ngauss, 1);
% [wgt_2d, lxs_2d] = CoorWeight_GL(2, 2);


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
n_dbc = 1;  

% boundary data; top (yMax), bottom (yMin), left (xMin) and right (xMax) surfaces 
check_x_value = {yMax};
check_dof_value = {2};
n_vec_value = { [0; 1] };
check_x_opt_value = { 2 };   % 2 for maximum and 1 for minimum 
constr_dof_value = { 1 };

dbc = struct( 'check_x', check_x_value, 'check_dof', check_dof_value, ...
    'n_vec', n_vec_value, 'check_x_opt', check_x_opt_value, ...
    'constr_dof', constr_dof_value );


MaxXLoad =0.8;   MinXLoad = -0.8;
% horizontal coordinate of two segments with prescribed pressure 
xSegs = [ xMin, MinXLoad;  MaxXLoad, xMax ]; 




%% elements list 
%dofs_bc = [];   % list of DOFs related to the essential boundary 
eles_bc = [];   % elements related to the essential boundary condition

for iele = 1 : NumElem 

    ME = ManiElems(iele);
    
    xy_ele = ME.xVertex;
    xy_ele = [xy_ele; xy_ele(1, :) ];

    MaxYElem = max(xy_ele(:, 2));
    if abs(MaxYElem - yMax) > GeoTol 
        continue;
    end

    bc_flag = 0;

    for inode = 1 : size(xy_ele, 1) - 1
        if abs(xy_ele(inode, 2) - yMax) > GeoTol || ...
                abs(xy_ele(inode + 1, 2) - yMax) > GeoTol
            continue;
        end
        
        MinXEdge = min(xy_ele(inode : inode+1, 1));
        MaxXEdge = max(xy_ele(inode : inode+1, 1));

        for iseg = 1 : 2 
            xPrescribeSeg = xSegs(iseg, :);

            if MinXEdge > xPrescribeSeg(2) + GeoTol 
                continue;
            end
            if MaxXEdge < xPrescribeSeg(1) - GeoTol 
                continue;
            end

            MinXSeg = max(xPrescribeSeg(1), MinXEdge);
            MaxXSeg = min(xPrescribeSeg(2), MaxXEdge);
            len_seg = MaxXSeg - MinXSeg;
            if len_seg < GeoTol
                continue;
            end

            bc_flag = 1;
            break;
        end  % iseg 
        if bc_flag == 1 
            break;
        end
    end  % inode 

    if bc_flag == 0
        continue;
    end

%     PPs=PPIncludeElem(iele,:);
%     Dofs = zeros(1, DOF_per_PP * length(PPs));
%     for ipp = 1:length(PPs)
%         Dofs( (ipp - 1) * DOF_per_PP + 1 : ipp * DOF_per_PP ) = ...
%             ( PPs(ipp) - 1 ) * DOF_per_PP + (1 : DOF_per_PP);
%     end 

    eles_bc = union(eles_bc, iele);
    %dofs_bc = union(dofs_bc, Dofs);
end  % iele       找到除去中间的那个边

% fprintf("finished computing dofs_bc.\n");
% fprintf("length(eles_bc) = %d\n",length(eles_bc));

ndim = 2;
 
%% boundary terms in Nitsche's method 

    
len_bc = 0;

n_vec = dbc(1).n_vec;

%     constr_dof = dbc(ibc).constr_dof;

% sub-normal vector corresponding to the constrained DOFs
%     n_vec_c = n_vec(constr_dof);


for jele = 1 : length(eles_bc)
    
    iele = eles_bc(jele);
    
    ME = ManiElems(iele);
    
    mat = Mats(ME.mat);
    kappa = mat.kappa;
    DMatP = kappa * eye(ndim);
    
    PPs = PhyPatches(ME.PP);
    xPPs_u = zeros(length(PPs), 2);
    for ipp = 1 : length(PPs)
        xPPs_u(ipp, :) = PPs(ipp).xNode;
    end
    
    if strcmpi(MeshShape,'BiotQ9Q4') || strcmpi(MeshShape,'BiotQ4Q4') || ...
            strcmpi(MeshShape,'Q4') || strcmpi(MeshShape,'Quad')
        xPPs_p = xPPs_u(1 : 4, :);
    elseif strcmpi(MeshShape,'BiotIRT3_RIRT3_R') || ...
            strcmpi(MeshShape,'BiotIRT3_LIRT3_L') || ...
            strcmpi(MeshShape,'IRT3_L') || strcmpi(MeshShape,'IRT3_R')
        xPPs_p = xPPs_u(1 : 3, :);
    else
        error("Element not implemented.");
    end
    
    h = max(xPPs_u(:, 1)) - min(xPPs_u(:, 1));
    
    Dofs = ME.DOF_p;
    
    xy_ele = ME.xVertex;
    xy_ele = [xy_ele; xy_ele(1, :) ];
    
    for inode = 1 : size(xy_ele, 1) - 1
        
        if abs(xy_ele(inode, 2) - yMax) > GeoTol || ...
                abs(xy_ele(inode + 1, 2) - yMax) > GeoTol
            continue;
        end
        
        
        MinXEdge = min(xy_ele(inode : inode+1, 1));
        MaxXEdge = max(xy_ele(inode : inode+1, 1));
        
        for iseg = 1 : 2
            xPrescribeSeg = xSegs(iseg, :);
            
            if MinXEdge > xPrescribeSeg(2) + GeoTol
                continue;
            end
            if MaxXEdge < xPrescribeSeg(1) - GeoTol
                continue;
            end
            
            MinXSeg = max(xPrescribeSeg(1), MinXEdge);
            MaxXSeg = min(xPrescribeSeg(2), MaxXEdge);
            len_seg = MaxXSeg - MinXSeg;
            if len_seg < GeoTol
                continue;
            end
            
            
            %             len_seg = norm(xy_ele(inode, :) - xy_ele(inode + 1, :));
            
            J = len_seg / 2;    % Jacobian for the 1D gauss integration
            
            len_bc = len_bc + len_seg;
            
            for igauss = 1 : length(wgt_1d)
                
                lx = lxs_1d(igauss);   % local coordinates
                % global coordinates of the integration point
                gxy = [ MinXSeg + (MaxXSeg - MinXSeg) / 2 * (lx + 1), yMax ];
                
                % shape function matrix for the solid
                [Nu_1] = NMatNMM2D_1(xPPs_p, gxy, MeshShape_p);
                % shape function matrix for pressure
                Np = Nu_1(1, 1 : 2 : end);
                
                % B matrix for the solid
                [Bu_1] = BMatNMM2D_1(xPPs_p, gxy, MeshShape_p);
                Bp = [ Bu_1(1, 1 : 2 : end);
                    Bu_1(2, 2 : 2 : end) ];
                
                
                w_mat = n_vec.' * DMatP * Bp;
                
                Kpp_eb(Dofs, Dofs) = Kpp_eb(Dofs, Dofs) - ...
                    w_mat.' * Np * J * wgt_1d(igauss);
                Kpp_eb(Dofs, Dofs) = Kpp_eb(Dofs, Dofs) - ...
                    Np.' * w_mat * J * wgt_1d(igauss);
                Kpp_eb(Dofs, Dofs) = Kpp_eb(Dofs, Dofs) + ...
                    gamma_D_p * (Np.' * Np) * J * wgt_1d(igauss) / h;
                
            end  % iguass
        end  % inode
    end  % jele
    
%     fprintf("ibc = %d, len_bc = %16.8e\n", ibc, len_bc);
    %     if abs(len_bc - 1) > 1e-6
    %         fprintf("wrong fluid boundary length, ibc = %d\n", ibc);
    %     end
    
end % ibc
                
                
    

