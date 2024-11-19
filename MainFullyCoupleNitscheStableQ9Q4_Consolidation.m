
% Stabilized NMM with Q4P4 elements for the quasi-staic modeling of 
% Biot's consolidation problem with u-p formulation, for the 1D consoliation
% problem, which is computed with a 2D model 
% Analytical solution for this problem is taken from the paper "General Theory
% of Three-Dimensional Consolidation" by Biot. 
% Stabilization technique in the paper "Stabilized low-order finite elements
% for strongly coupled poromechanical  problems" by Wentao Li and Changfu Wei is used.

clc;
clear;
clear all;
close all;

InputFile = "input_2DColumn.txt";




global MeshShape NumElem CoorElemNode  ElemNodeFL 
global NumPP CoorPPStar PPIncludeElem 
global DOF_per_PP_U DOF_per_PP_P GDOF_U GDOF_P
global ManiElems PhyPatches
global dof_u_full  dof_p_full
%global M alpha k mu delta_t current_t theta




%%
% geometric parameters
H = 0.7/ cos(pi/4);   % height 
% H = 1;

% Data taken from the paper "On the optimization of the fixed-stress splitting
% for Biot's equations" by Storvik et al.
% material parameters 
% E = 1e-2;        % Young's modulus 
% nu = 0;
% n = 0.30;       % porosity
% Ks = 1.0e3;    % bulk modulus of the solid phase 
% Kf = 2;       % bulk modulus of fluid 
% kappa_w = 1e-13;      % intrinsic permeability
% muw = 1e-12;      % dynamic viscosity 
% alpha = 1.0;    % Biot's coefficient
% M = 1 / ((alpha - n)/Ks + n / Kf);  % Biot's modulus 
% Q = 1e-5;       % pressure applied on the bottom surface 

E = 10;        % Young's modulus 
nu = 0.2;       % Poisson's raio
n = 0.30;       % porosity
Ks = 1.0e6;    % bulk modulus of the solid phase 
Kf = 2e3;       % bulk modulus of fluid 
kappa_w = 1e-13;      % intrinsic permeability
muw = 1e-9;      % dynamic viscosity 
alpha = 1.0;    % Biot's coefficient
M = 1 / ((alpha - n)/Ks + n / Kf);  % Biot's modulus 
Q = 1e-2;       % pressure applied on the bootom surface

AnalysisType = "plane_strain";

kappa = kappa_w/muw;

delta_t = 2;      % time duration of each step 
nsteps = 2;         % number of time steps of the simulation 

% penalty factor
PenaltyFactor=10^4 * E;

% parameter for the integration with respect to time, 0, 0.5, 2/3, and 1.0 
% represent the classical forward, central, Galerkin, and backward Euler
theta = 1.0;  

tao = 0;

ndim = 2;           % dimension of the model 

% Lame's constants 
lamda = E * nu / (1 + nu) / (1 - 2 * nu);
G = E / 2 / (1 + nu);

iteration_tol = 1e-10;   % tolerance for the iteration 

gamma_D_u = 400 * lamda; % Nitsche method parameter for displacement 
gamma_D_p = 500;   % Nitsche method parameter for pressure

gamma_g_u = 0.05*lamda;    % ghost penalty parameter for the solid 
gamma_g_m = 0.5 / M; % ghost penalty parameter for mass of fluid 
gamma_g_p = 0.4 * kappa * delta_t;   % ghost penalty parameter for the fluid 

%% pre-processing 
% Main_preprocessing_1(InputFile);
Main_preprocessing_cut(InputFile);

[nSegDisDirich_t, SegsDisDirich_t, nSegDisDirich_b, SegsDisDirich_b, ...
    nSegDisDirich_r, SegsDisDirich_r, nSegDisDirich_l, SegsDisDirich_l] = BCSegsConsolidation();

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
    [ana_u, ana_p] = AnaUP1DConsolidation(E, alpha, M, muw, kappa_w, nu, Q, gxy, H, current_t);
    d_u_n(PP.DOF_u) = ana_u;
    d_p_n(PP.DOF_p) = ana_p;
end

% coefficient matrices that keep constant during the simulation 
[Kuu, Kpp, Mpp, Kup, Spp]=CoeffMatBiot2D_Consolidation(DMatU, DMatP, ....
    M, alpha, tao);


[Kuu_eb] = KBC_u_Nitsche_Consolidation_2(DMatU, gamma_D_u, nSegDisDirich_t, SegsDisDirich_t, nSegDisDirich_b, SegsDisDirich_b, ...
    nSegDisDirich_r, SegsDisDirich_r, nSegDisDirich_l, SegsDisDirich_l);
[Kpp_eb] = KBC_p_Nitsche_Consolidation_2(DMatP, gamma_D_p, nSegDisDirich_b, SegsDisDirich_b, nSegDisDirich_t, SegsDisDirich_t,...
    nSegDisDirich_r, SegsDisDirich_r, nSegDisDirich_l, SegsDisDirich_l);


% Additional coupling matrix on the essential boundary of displacement. 
% Only used when essential boundary condition for displacement is not
% satisfied. 
[Kup_b]=CoupleMatBoundary2D_Consolidation(alpha, nSegDisDirich_t, SegsDisDirich_t, nSegDisDirich_b, SegsDisDirich_b, ...
    nSegDisDirich_r, SegsDisDirich_r, nSegDisDirich_l, SegsDisDirich_l);

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
    [fu] = LoadConsolidation94(Q, nSegDisDirich_b, SegsDisDirich_b);
    fp = zeros(GDOF_P, 1); 
    
    % Additional terms of RHS in Nitsche's method for essential
    % boundary condition for pore pressure
    fp_eb = zeros(GDOF_P, 1); 
%     [fp_eb] = FBC_p_Nitsche_Consolidation_2(DMatP, E, alpha, M, muw, kappa_w, nu, Q, ...
%     H, current_t, gamma_D_p, nSegDisDirich_t, SegsDisDirich_t, nSegDisDirich_b, SegsDisDirich_b, ...
%     nSegDisDirich_r, SegsDisDirich_r, nSegDisDirich_l, SegsDisDirich_l);
    % Additional terms of RHS in Nitsche's method for essential
    % boundary condition for displacements
    fu_eb = zeros(GDOF_U, 1);
%     [fu_eb] = FBC_u_Nitsche_Consolidation_2(DMatU, E, alpha, M, muw, kappa_w, nu, Q, ...
%     H, current_t, gamma_D_u, nSegDisDirich_t, SegsDisDirich_t, nSegDisDirich_b, SegsDisDirich_b, ...
%     nSegDisDirich_r, SegsDisDirich_r, nSegDisDirich_l, SegsDisDirich_l);
    % Additional terms of RHS of the fluid on the EBC of solid for the 
    % symmetric Nitsche's method 
    [fp_b] = zeros(GDOF_P, 1);
     
%     KMat = zeros(GDOF_U + GDOF_P, GDOF_U + GDOF_P);
%     KMat(dof_u_full, dof_u_full) = - Kuu;
%     KMat(dof_u_full, dof_p_full) = (Kup - Kup_b);
%     KMat(dof_p_full, dof_u_full) = (Kup - Kup_b).';
%     KMat(dof_p_full, dof_p_full) = theta * delta_t * Kpp + Mpp + Spp;
%     
%     FVec = zeros(GDOF_U + GDOF_P, 1);
%     FVec(dof_u_full) = fu + fu_eb;
%     FVec(dof_p_full) = (fp + fp_eb) * delta_t + Mpp * d_p_n ...
%             + Kup.' * d_u_n;
% 
%     [KNew,FNew] = DBC_FEM_Consolidation_1(KMat, FVec, E, alpha, M, muw, kappa_w, nu, Q, H, current_t);   
%     solution = KNew \ FNew;


    KMat = zeros(GDOF_U + GDOF_P, GDOF_U + GDOF_P);
    KMat(dof_u_full, dof_u_full) = -(Kuu + Kuu_eb +Kuu_ghost);
    KMat(dof_u_full, dof_p_full) = (Kup - Kup_b);
    KMat(dof_p_full, dof_u_full) = (Kup).';
    KMat(dof_p_full, dof_p_full) = theta * delta_t * (Kpp + Kpp_eb) + Mpp + Spp + Kpp_ghost + Mpp_ghost;
    
    FVec = zeros(GDOF_U + GDOF_P, 1);
    FVec(dof_u_full) = -(fu + fu_eb);
    FVec(dof_p_full) = (fp + fp_eb) * delta_t + (Mpp + Spp) * d_p_n ...
            + Kup.' * d_u_n;

    solution = KMat \ FVec;
    
    d_u = solution(dof_u_full);
    d_p = solution(dof_p_full);
    
    % updating of solutions of the previous time step 
    d_u_n = d_u;
    d_p_n = d_p;
    

end  % istep: end of time marching

    % displacement contour 
    IndexDisp = 2;
    if strcmpi(MeshShape,'BiotQ9Q4') ||strcmpi(MeshShape,'BiotQ4Q4')
        DispContourNMM_1(IndexDisp, d_u_n);
    else
        DispContourNMM(IndexDisp, d_u);
    end
    % pressure contour 
    if strcmpi(MeshShape,'BiotQ9Q4') ||strcmpi(MeshShape,'BiotQ4Q4')
        PressureContourNMM_1(d_p_n);
    else
        PressureContourNMM(d_p);
    end
     
    PlotP1DConsolidation_2(E, nu, current_t, d_p, alpha, M, kappa_w, muw, ...
        Q, H);
    PlotV1DConsolidation_2(E, nu, current_t, d_u, alpha, M, kappa_w, muw, ...
        Q, H);


   