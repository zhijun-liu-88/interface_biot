
function [Kuu, Kpp, Mpp, Mpp_s, Kup]=CoeffMatBiot2Dym(AnalysisType)
% to compute coefficient matries of fully discrete fixed-stress spliting of
% a 2D consolition problem
% input 
% beta: the stabilziation parameter for the fixed-stress splitting 


global MeshShape 
global NumElem 
global DOF_per_PP_U DOF_per_PP_P GDOF_U GDOF_P
global ManiElems PhyPatches
global Mats


Kup = zeros(GDOF_U, GDOF_P);

Kuu = zeros(GDOF_U);    
Kpp = zeros(GDOF_P);
Mpp = zeros(GDOF_P);
Mpp_s = zeros(GDOF_P);

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

ndim = 2;   % dimension 
% wights and local coordinates 
[wgt, lxs] = CoorWeight_GL(ngauss, ndim);

for iele = 1 : NumElem
    
    ME = ManiElems(iele);
    
    mat = Mats(ME.mat);
    E = mat.E;  nu = mat.nu;  
    kappa = mat.kappa;     M = mat.M;  alpha = mat.alpha; 
    
    % Lame's constants 
    lambda = E * nu / (1 + nu) / (1 - 2 * nu);
    G = E / 2 / (1 + nu);

    % drained bulk modulus 
    Kdr = 2*G/ndim + lambda; 

    % the minimum choice of stabilization parameter in the paper: Convergence 
    % of iterative coupling for coupled flow and geomechanics
    % Stabilization parameter for fixed-stress splitting. 
    beta = alpha^2/Kdr/2;
    
    [DMatU] = ConstitutiveMat(E, nu, AnalysisType);
    DMatP = kappa * eye(ndim);
    
    
    PPs = PhyPatches(ME.PP);
    
    DOFs_U = ME.DOF_u;
    DOFs_P = ME.DOF_p;
    
    % element stiffness and right-hand vector
    KeUU = zeros(length(DOFs_U), length(DOFs_U));   
    KeUP = zeros(length(DOFs_U), length(DOFs_P));   
    KePP = zeros(length(DOFs_P), length(DOFs_P)); 
    MePP = zeros(length(DOFs_P), length(DOFs_P));   

    xPPs_u = zeros(length(PPs), 2);
    for ipp = 1 : length(PPs)
        xPPs_u(ipp, :) = PPs(ipp).xNode;    % 与流形单元有关的四个物理覆盖的星点
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
    
    % length of the mathematical element in x and y direction 
    la = max(xPPs_u(:, 1)) - min(xPPs_u(:, 1)); 
    lb = max(xPPs_u(:, 2)) - min(xPPs_u(:, 2));
    
    xy_poly = ME.xVertex;
    xy_poly = [xy_poly; xy_poly(1, :) ];
    [NumIntegQuad, xyIntegQuadVect] = IntegQuad2d(xy_poly);
    
    % element stiffness 
    for iquad=1:NumIntegQuad     % loop of quadrilaterals
                
        xQuad=xyIntegQuadVect(:,:,iquad);  % 流形单元的坐标

%         AreaQuad=0.5*abs(det([1,xQuad(1,:);1,xQuad(2,:);1,xQuad(3,:)]))+...
%                  0.5*abs(det([1,xQuad(1,:);1,xQuad(3,:);1,xQuad(4,:)]));
% 
%         if AreaQuad < 10^(-20) * la * lb                                       
%             continue;
%         end

        for igauss=1:length(wgt)    % loop of integration points
            
            [Nq, dNqdxi]=lagrange_basis('Q4', lxs(igauss,:));
            
            gxy = Nq.' * xQuad;   % global coordinates of the integration points
            
            % Jacobian matrix of the transformation in integration quadrilateral
            J_i = xQuad' * dNqdxi;
            % Jacobian
            detJ_i = abs(det(J_i));
            
            % B matrix 
            [Bu] = BMatNMM2D_1(xPPs_u, gxy, MeshShape_u);
            [Bu_1] = BMatNMM2D_1(xPPs_p, gxy, MeshShape_p);
            Bp = [ Bu_1(1, 1 : 2 : end);
                   Bu_1(2, 2 : 2 : end) ];
            
            % B matrix for the volumetric strain, only valid for plane
            % strain problem 
            Bv = Bu(1, :) + Bu(2, :);
            
            [Nu_1] = NMatNMM2D_1(xPPs_u, gxy, MeshShape_p); 
            % shape function matrix for pressure 
            Np = Nu_1(1, 1 : 2 : end);
            
            KeUU = KeUU + Bu.' * DMatU * Bu * detJ_i * wgt(igauss);
            KeUP = KeUP + alpha * Bv.' * Np * detJ_i * wgt(igauss);
            
            KePP = KePP + Bp.' * DMatP * Bp * detJ_i * wgt(igauss);
            MePP = MePP + (Np.' * Np) * detJ_i * wgt(igauss);
            
        end  % l : loop of integration points 
    end  % k : loop of quadrilaterals
    
%     MePP = MePP * 1/M;
%     MePP_s = MePP * beta;
    
    Kuu(DOFs_U, DOFs_U) = Kuu(DOFs_U, DOFs_U) + KeUU;
    Kpp(DOFs_P, DOFs_P) = Kpp(DOFs_P, DOFs_P) + KePP;
    Mpp(DOFs_P, DOFs_P) = Mpp(DOFs_P, DOFs_P)+ MePP * 1/M;
    Mpp_s(DOFs_P, DOFs_P) = Mpp_s(DOFs_P, DOFs_P) + MePP * beta;
    
    Kup(DOFs_U, DOFs_P) = Kup(DOFs_U, DOFs_P) + KeUP;
            
end % iele

    
    
    
    
    




    
    