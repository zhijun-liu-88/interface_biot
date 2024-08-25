
function [Kup_b]=CoupleMatUBC_AguilarFootym()
% To compute the coupling matrix on the Dirichlet boundary of displacement. 
% Only used when essential boundary condition for displacement is not
% satisfied by the interpolation.  


global MeshShape NumElem 
global GDOF_U GDOF_P
global ManiElems PhyPatches
global Mats


if strcmpi(MeshShape,'BiotQ9Q4')  
    MeshShape_u = 'Q9';  MeshShape_p = 'Q4';
    ngauss = 4;
elseif strcmpi(MeshShape,'BiotQ4Q4') || strcmpi(MeshShape,'Q4') || ... 
        strcmpi(MeshShape,'Quad')
    MeshShape_u = 'Q4';  MeshShape_p = 'Q4';
    ngauss = 2; % number of gauss points in 1 direction 
elseif strcmpi(MeshShape,'BiotIRT3_RIRT3_R') || strcmpi(MeshShape,'BiotIRT3_LIRT3_L') || ...
        strcmpi(MeshShape,'IRT3_R') || strcmpi(MeshShape,'IRT3_L') 
    MeshShape_u = 'T3';  MeshShape_p = 'T3';
    ngauss = 2; % number of gauss points in 1 direction 
end


Kup_b = zeros(GDOF_U, GDOF_P);


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



%% bounary data 
n_dbc = 3;  

% boundary data; top (yMax), bottom (yMin), left (xMin) and right (xMax) surfaces 
check_x_value = {yMin, xMax, xMin};
check_dof_value = {2, 1, 1};
n_vec_value = { [0; -1], [1; 0], [-1; 0] };
check_x_opt_value = { 1, 2, 1};   % 2 for maximum and 1 for minimum 
constr_dof_value = { [1 2], [1, 2], [1, 2] };

dbc = struct( 'check_x', check_x_value, 'check_dof', check_dof_value, ...
    'n_vec', n_vec_value, 'check_x_opt', check_x_opt_value, ...
    'constr_dof', constr_dof_value );
    
GeoTol = 1.e-13;

% wights and local coordinates 
[wgt_1d, lxs_1d] = CoorWeight_GL(ngauss, 1);


%% coupling matrix 
for ibc = 1 : n_dbc
    
    %fprintf("%d-th boundary.\n", ibc);
    
%     len_bc = 0;
    
    n_vec = dbc(ibc).n_vec;
    
    constr_dof = dbc(ibc).constr_dof;

    % sub-normal vector corresponding to the constrained DOFs
    n_vec_c = n_vec(constr_dof);
    
    check_x = dbc(ibc).check_x;
    check_dof = dbc(ibc).check_dof;
    check_x_opt = dbc(ibc).check_x_opt;
    
    for iele = 1 : NumElem
        
        ME = ManiElems(iele);
        mat = Mats(ME.mat);
        alpha = mat.alpha; 
        
        xy_ele = ME.xVertex;
        xy_ele = [xy_ele; xy_ele(1, :) ];
        
        if check_x_opt == 1
            ele_check_x = min(xy_ele(:, check_dof));
        else
            ele_check_x = max(xy_ele(:, check_dof));
        end 

        if abs(ele_check_x - check_x) > GeoTol
            continue;
        end
        
        PPs = PhyPatches(ME.PP);
    
        DOFs_U = ME.DOF_u;
        DOFs_P = ME.DOF_p;

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
        
        KeUP = zeros(length(DOFs_U), length(DOFs_P));   

        for inode = 1 : size(xy_ele, 1) - 1

            if abs(xy_ele(inode, check_dof) - check_x) > GeoTol || ...
                abs(xy_ele(inode + 1, check_dof) - check_x) > GeoTol
                continue;
            end

            len_seg = norm(xy_ele(inode, :) - xy_ele(inode + 1, :));
            J = len_seg / 2;    % Jacobian for the 1D gauss integration 


            for igauss = 1 : length(wgt_1d) 

                lx = lxs_1d(igauss);   % local coordinates
                % global coordinates of the integration point
                gxy = xy_ele(inode, :) + ...
                    ( xy_ele(inode + 1, :) - xy_ele(inode, :) ) / 2 * ...
                    ( lx + 1 );

                % shape function matrix 
                [Nu] = NMatNMM2D_1(xPPs_u, gxy, MeshShape_u); 
                
                [Nu_1] = NMatNMM2D_1(xPPs_p, gxy, MeshShape_p); 
                % shape function matrix for pressure 
                Np = Nu_1(1, 1 : 2 : end);

                % shape function matrix corresponding to the constrained DOFs
                Nu_c = Nu(constr_dof, :);
                
                KeUP = KeUP + Nu_c.' * n_vec_c * Np * J * wgt_1d(igauss);

            end  % iguass 
        end  % inode 
        
        Kup_b(DOFs_U, DOFs_P) = Kup_b(DOFs_U, DOFs_P) + KeUP * alpha;
        
    end  % iele 
%     
%     fprintf("In CoupleMatBoundary2D, ibc = %d, len_bc = %16.8e\n", ibc, len_bc);
%     if abs(len_bc - 1) > 1e-6
%         fprintf("wrong boundary length, ibc = %d\n", ibc);
%     end
    
end % ibc



    

