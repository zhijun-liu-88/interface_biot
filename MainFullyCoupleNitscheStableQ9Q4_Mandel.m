
% Non-stabilized 2D FEM with unfitted meshes for the quasi-staic modeling of Biot's 
% consolidation problem with u-p formulation. Nitsche's method is applied
% for both displacements and pore pressure. 
% Mandel problem is considered. The analytical solution is taken from the
% paper "On the optimization of the fixed-stress splitting
% for Biot's equations" by Storvik et al.
% Left and right boundaries are traction free and prescribed with zero pore
% pressure. Top and bottom boundaries are impervious. Analytical
% displacements are prescribed on top and bottom boundaries.

clc;
clear;
clear AnaUPMandel;
% clear all;
close all;

InputFile = "input_Mandel.txt";




global MeshShape NumElem CoorElemNode  ElemNodeFL 
global NumPP CoorPPStar PPIncludeElem 
global DOF_per_PP_U DOF_per_PP_P GDOF_U GDOF_P
global ManiElems PhyPatches
global dof_u_full  dof_p_full
%global M alpha k mu delta_t current_t theta




%%
% geometric parameters
Lx = 100;   % half width 
Ly = 10;    % half height 

% Data taken from the paper "On the optimization of the fixed-stress splitting
% for Biot's equations" by Storvik et al.
% material parameters
E = 5.94e9;        % Young's modulus 
nu = 0.2;       % Poisson's raio
kappa = 1e-10;      % hydraulic conductivity 
alpha = 1.0;     % Biot's coefficient
M = 1.65e10;
F = 6e8;         % applied force 

delta_t = 10;      % time duration of each step 
nsteps = 10;         % number of time steps of the simulation 


% E = 1;        % Young's modulus 
% nu = 0.2;       % Poisson's raio
% kappa = 0.01;      % hydraulic conductivity 
% alpha = 1.0;     % Biot's coefficient
% M = 40;
% F = 100;         % applied force 
% 
% delta_t = 10;      % time duration of each step 
% nsteps = 3;         % number of time steps of the simulation 

AnalysisType = "plane_strain";

% penalty factor
PenaltyFactor=10^4 * E;

% parameter for the integration with respect to time, 0, 0.5, 2/3, and 1.0 
% represent the classical forward, central, Galerkin, and backward Euler
theta = 1.0;  

ndim = 2;           % dimension of the model 

% Lame's constants 
lamda = E * nu / (1 + nu) / (1 - 2 * nu);
G = E / 2 / (1 + nu);

% draied bulk modulus 
Kdr = 2*G/ndim + lamda; 

% the minimum choice of stabilization parameter in the paper: Convergence 
% of iterative coupling for coupled flow and geomechanics
beta = alpha^2/Kdr/2;

iteration_tol = 1e-10;   % tolerance for the iteration 

% solution of characteristic equation
nMandelTerm = 20;
[angle_vec] = SolveCharacteristicEq(E, nu, alpha, M, nMandelTerm);

gamma_D_u = 400 * lamda; % Nitsche method parameter for displacement 
gamma_D_p = 50;   % Nitsche method parameter for pressure

gamma_g_u = 0.05*lamda;    % ghost penalty parameter for the solid 
gamma_g_m = 0.5 / M; % ghost penalty parameter for mass of fluid 
gamma_g_p = 0.4 * kappa * delta_t;   % ghost penalty parameter for the fluid 

%% pre-processing 
% Main_preprocessing_1(InputFile);
Main_preprocessing_cut(InputFile);

% TecFile = "Consolidation1D.dat";
% FileID = fopen(TecFile, "w");
% 
% TempLine = 'variables="x","y","ux","uy","p","exx","eyy","exy",';
% TempLine = strcat(TempLine,'"sxx","syy","sxy"');
% fprintf(FileID,'%s\n',TempLine);
% fclose(FileID);


%% constitutive matrix 
[DMatU] = ConstitutiveMat(E, nu, AnalysisType);
DMatP = kappa * eye(ndim);


% d_p_n: values of DoFs of previous time step for pressure 
% d_p_i: values of DoFs of previous iteration for pressure 
% d_u_n: values of DoFs of previous time step for displacement
% d_u_i: values of DoFs of previous iteration for displacement

%% marching time 
% initial displacements and pore pressure 
d_p_n = zeros(GDOF_P, 1); 
d_u_n = zeros(GDOF_U, 1); 
current_t = 0;      % current time 
for ipp = 1 : NumPP 
    PP = PhyPatches(ipp);
    gxy = PP.xNode;
    [ana_u, ana_p] = AnaUPMandel_2(E, nu, alpha, M, kappa, F, ...
        Lx, current_t, gxy, angle_vec, nMandelTerm);
    d_u_n(PP.DOF_u) = ana_u;
    if PP.patchUP(2) == 0 
        continue;
    end
    d_p_n(PP.DOF_p) = ana_p; 
end

% coefficient matrices that keep constant during the simulation 
[Kuu, Kpp, Mpp, Mpp_s, Kup]=CoeffMatBiot2D(DMatU, DMatP, ....
    M, beta, alpha);


[Kuu_eb] = KBC_u_Nitsche_Mandel_2(DMatU, gamma_D_u);
[Kpp_eb] = KBC_p_Nitsche_Mandel_2(DMatP, gamma_D_p);

% Additional coupling matrix on the essential boundary of displacement. 
% Only used when essential boundary condition for displacement is not
% satisfied. 
[Kup_b]=CoupleMatBoundary2D_Mandel(alpha);

% ghost penalty coefficient matrix for solid
[Kuu_ghost] = GhostPenalty2DSolid_Quadratic1(gamma_g_u);
% ghost penalty coefficient matrix for fluid 
[Kpp_ghost, Mpp_ghost] = GhostPenalty2DFluid_1(gamma_g_p, gamma_g_m);

for istep = 1 : nsteps  % time marching 
% for istep = 1 : 1   % time marching 
    
    fprintf("%d-th time step.\n", istep);
    
    current_t = current_t + delta_t;
    
    % initialization of solution of previous iteration 
    d_u_i = d_u_n;
    d_p_i = d_p_n;
    
    % right hand sides caused by sources and natural boundary terms
    fu = zeros(GDOF_U, 1);
    fp = zeros(GDOF_P, 1); 
    
    % Additional terms of RHS in Nitsche's method for essential
    % boundary condition for pore pressure 
    fp_eb = zeros(GDOF_P, 1);
    % Additional terms of RHS in Nitsche's method for essential
    % boundary condition for displacements
    [fu_eb] = FBC_u_Nitsche_Mandel_2(DMatU, current_t, gamma_D_u, E, ...
        nu, alpha, M, kappa, F, Lx, angle_vec, nMandelTerm);
    % Additional terms of RHS of the fluid on the EBC of solid for the 
    % symmetric Nitsche's method 
    [fp_b]=FpCoupleUBC_Mandel_2(E, nu, alpha, M, kappa, F, Lx, ...
        current_t, angle_vec, nMandelTerm);
    
    
%     % Non-symmetric Nitsche method for fully coupled problem 
%     KMat = zeros(GDOF_U + GDOF_P, GDOF_U + GDOF_P);
%     KMat(dof_u_full, dof_u_full) = Kuu + Kuu_eb;
%     KMat(dof_u_full, dof_p_full) = -(Kup - Kup_b);
%     KMat(dof_p_full, dof_u_full) = Kup.';
%     KMat(dof_p_full, dof_p_full) = delta_t * (Kpp + Kpp_eb) + Mpp;
%     
%     FVec = zeros(GDOF_U + GDOF_P, 1);
%     FVec(dof_u_full) = fu + fu_eb;
%     FVec(dof_p_full) = (fp + fp_eb) * delta_t + Mpp * d_p_n ...
%             + Kup.' * d_u_n;
%     solution = KMat \ FVec;


    % Symmetric Nitsche method for fully coupled problem 
    KMat = zeros(GDOF_U + GDOF_P, GDOF_U + GDOF_P);
    KMat(dof_u_full, dof_u_full) = -(Kuu + Kuu_eb +Kuu_ghost);
    KMat(dof_u_full, dof_p_full) = (Kup - Kup_b);
    KMat(dof_p_full, dof_u_full) = (Kup - Kup_b).';
    KMat(dof_p_full, dof_p_full) = delta_t * (Kpp + Kpp_eb) + Mpp + Kpp_ghost + Mpp_ghost;
    
    FVec = zeros(GDOF_U + GDOF_P, 1);
    FVec(dof_u_full) = -(fu + fu_eb);
    FVec(dof_p_full) = (fp + fp_eb) * delta_t + Mpp * d_p_n ...
            + Kup.' * d_u_n - fp_b;

    solution = KMat \ FVec;
    
    d_u = solution(dof_u_full);
    d_p = solution(dof_p_full);
    
    % updating of solutions of the previous time step 
    d_u_n = d_u;
    d_p_n = d_p;


end  % istep: end of time marching

       
    % displacement contour 
    IndexDisp = 1;
    if strcmpi(MeshShape,'BiotQ9Q4') ||strcmpi(MeshShape,'BiotQ4Q4')
        DispContourNMM_1(IndexDisp, d_u);
    else
        DispContourNMM(IndexDisp, d_u);
    end
    % pressure contour 
    if strcmpi(MeshShape,'BiotQ9Q4') ||strcmpi(MeshShape,'BiotQ4Q4')
        PressureContourNMM_1(d_p);
    else
        PressureContourNMM(d_p);
    end
     
    PlotP1DMandel(E, nu, istep, current_t, d_p, alpha, M, kappa, F, ...
        Lx, angle_vec, nMandelTerm);
    PlotV1DMandel(E, nu, istep, current_t, d_u, alpha, M, kappa, F, ...
        Lx, angle_vec, nMandelTerm);


   