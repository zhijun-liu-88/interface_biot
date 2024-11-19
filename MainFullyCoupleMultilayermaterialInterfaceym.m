% The poroelasticity analysis is conducted for a footing with three layers


clc;
clear;
% clear all;
close all;

InputFile = "input_Muitilayer_materialInterface1.txt";




global MeshShape 
global NumPP 
global DOF_per_PP_U DOF_per_PP_P GDOF_U GDOF_P
global ManiElems PhyPatches
global dof_u_full  dof_p_full
global Mats
%global M alpha k mu delta_t current_t theta



%%

% to define materials 
E1 =  3e4;   E2 = E1 * 0.5;  E3 = E1 * 0.2;  nu_1 = 0.2;   nu_2 = 0.2;   nu_3 = 0.2;
% E1 =  9;   E2 = 41;  E3 = 13;  nu_1 = 0.14;   nu_2 = 0.21;   nu_3 = 0.08;
kappa_1 = 1e-4;      kappa_2 = kappa_1 * 2;   kappa_3 = kappa_1 * 5 ;% hydraulic conductivity 
M = 1e20;   alpha = 1;       

[Mats] = DefineMulMaterial(E1, E2,  E3, nu_1, nu_2, nu_3, kappa_1, kappa_2, kappa_3, M, alpha);

Sigma0 = 1e4;              % pressure 

delta_t = 1e-3;               % time duratio4n of each step
nsteps = 10000;
% total_t = 1e-2;
% nsteps = total_t/delta_t;         % number of time steps of the simulation 






% coordinates of the interface
NumInterface = 2;
xyInterface1 = [ -4    -2.4000001   4    -2.4000001];
xyInterface2 = [ -4    2.4000001   4    2.4000001];
yInter = [xyInterface1(1, 2),xyInterface2(1, 2)];








% % xInter = xyInterface(1, 1);



AnalysisType = "plane_strain";


% parameter for the integration with respect to time, 0, 0.5, 2/3, and 1.0 
% represent the classical forward, central, Galerkin, and backward Euler
theta = 1.0;  

ndim = 2;           % dimension of the model 


DOF_per_PP_U = ndim;    % number of DOFs of each physical patch for displacement
DOF_per_PP_P = 1;       % number of DOFs of each physical patch for pressure
% DOF_per_PP = DOF_per_PP_U + DOF_per_PP_P;


tic

%% pre-processing 
Main_materialpreprocessing_Interfaceymcs(InputFile, xyInterface1,xyInterface2,NumInterface);
fprintf("finished preprocessing.\n");
toc



%% ghost penalty parameters 


gamma_g_u = 1 * Mats(1).E;    % ghost penalty parameter for the solid 
gamma_g_m = 0.1 / Mats(1).M; % ghost penalty parameter for mass of fluid 
gamma_g_p = 1 * Mats(1).kappa * delta_t;   % ghost penalty parameter for the fluid 

gamma_D_u = 400 * Mats(1).E; % Nitsche method parameter for displacement, for the bondary 
gamma_D_p = 20 * Mats(1).kappa;   % Nitsche method parameter for pressure, for the bondary  

gamma_I_u = 400 * Mats(1).E; % Nitsche method parameter for displacement, for the interface
gamma_I_p = 500 * Mats(1).kappa;   % Nitsche method parameter for pressure, for the interface 




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
    ana_u = [0; 0];   ana_p = 0;
    d_u_n(PP.DOF_u) = ana_u;
    if PP.patchUP(2) == 0 
        continue;
    end
    d_p_n(PP.DOF_p) = ana_p; 
end

% coefficient matrices that keep constant during the simulation 
if strcmpi(MeshShape,'BiotQ9Q4') || strcmpi(MeshShape,'BiotQ4Q4')
    [Kuu, Kpp, Mpp, Mpp_s, Kup]=CoeffMatBiot2Dym(AnalysisType);
else
    [Kuu, Kpp, Mpp, Mpp_s, Kup]=CoeffMatFixedStress2D(DMatU, DMatP, ....
        M, beta, alpha);
end
fprintf("finished coefficient matrices.\n");
toc



% Additional terms of coefficient matrix in Nitsche's method for essential
% boundary condition for displacements 
[Kuu_eb] = KBC_u_Nitsche_Mulmaterial_ym(AnalysisType, gamma_D_u);
% Additional terms of coefficient matrix in Nitsche's method for essential
% boundary condition for pressure 
[Kpp_eb] = KBC_p_Nitsche_Mulmaterial_2ym(gamma_D_p);
fprintf("finished K_eb.\n");
toc
% Additional coupling matrix on the essential boundary of displacement. 
% Used when Nitsche's method is used to get a symmetric weak form 
[Kup_b]=CoupleMatUBC_AguilarFootym();
fprintf("finished Kup_b.\n");
toc

% submatrices of the coefficient matrix, Kuu, Kpp, and Kup, corresponding to
% the interface term 
[Kuu_I] = KInterface_u_Nitsche_AguilarFootym(AnalysisType, gamma_I_u);
[Kpp_I] = KInterface_p_Nitsche_AguilarFootym(gamma_I_p);
[Kup_I] = KInterface_up_Nitsche_AguilarFootym();


% ghost penalty coefficient matrix for solid
[Kuu_ghost] = GhostPenalty2DSolid_Quadratic1(gamma_g_u);
% ghost penalty coefficient matrix for fluid 
[Kpp_ghost, Mpp_ghost] = GhostPenalty2DFluid_1(gamma_g_p, gamma_g_m);
fprintf("finished ghost penalty matrices.\n");
toc

% 
% Symmetric Nitsche method for fully coupled problem 
KMat = zeros(GDOF_U + GDOF_P, GDOF_U + GDOF_P);
KMat(dof_u_full, dof_u_full) = -(Kuu + Kuu_eb + Kuu_I + Kuu_ghost);
KMat(dof_u_full, dof_p_full) = (Kup - Kup_b - Kup_I);
KMat(dof_p_full, dof_u_full) = (Kup - Kup_b - Kup_I).';
KMat(dof_p_full, dof_p_full) = delta_t * (Kpp + Kpp_eb + Kpp_I) + Kpp_ghost ...
    + Mpp + Mpp_ghost;




% coefficient matrix without ghost penalty
KMat_us = zeros(GDOF_U + GDOF_P, GDOF_U + GDOF_P);
KMat_us(dof_u_full, dof_u_full) = -(Kuu + Kuu_eb + Kuu_I);
KMat_us(dof_u_full, dof_p_full) = (Kup - Kup_b - Kup_I);
KMat_us(dof_p_full, dof_u_full) = (Kup - Kup_b - Kup_I).';
% KMat_us(dof_p_full, dof_p_full) = delta_t * (Kpp + Kpp_eb + Kpp_I) + Mpp + Kpp_fpl;
KMat_us(dof_p_full, dof_p_full) = delta_t * (Kpp + Kpp_eb + Kpp_I) + Mpp;


% % 
% kappa_s = cond(KMat, 2);
% kappa_us = cond(KMat_us, 2);
% 
% 
% fprintf("condition number of coefficient matrix: %14.6e\n", kappa_s);
% fprintf("condition number of coefficient matrix without ghost penalty: %14.6e\n", kappa_us);
% toc 
% % 
% % 
% % 
% FileID = fopen('Couple_interface.txt', 'w');
% fprintf(FileID, "step current_t GDOF Ed Ep:\n");
% 
% 
% 
% 
% 
% 
% right hand sides caused by sources and natural boundary terms
[fu, fp] = LoadUnitMulmaterial_1(Sigma0);


% Additional terms of RHS in Nitsche's method for essential
% boundary condition for pore pressure
fp_eb = zeros(GDOF_P, 1);
% Additional terms of RHS in Nitsche's method for essential
% boundary condition for displacements
fu_eb = zeros(GDOF_U, 1);
% Additional terms of RHS of the fluid on the EBC of solid for the
% symmetric Nitsche's method
fp_b = zeros(GDOF_P, 1);
% 
% Ed_0 = 0;    Ep_0 = 0;    Ee = 0;
% nsteps = 100;

% for istep = 1 : nsteps  % time marching 
for istep = 1 : 1   % time marching 
    
    %fprintf("%d-th time step.\n", istep);
    
    current_t = current_t + delta_t;
    
    % initialization of solution of previous iteration 
    d_u_i = d_u_n;
    d_p_i = d_p_n;
    

    
%     fprintf("finished RHS.\n");
%     toc
    
    
    % RHS vector for symmetric Nitsche method for fully coupled problem 
    FVec = zeros(GDOF_U + GDOF_P, 1);
    FVec(dof_u_full) = -(fu + fu_eb);
    FVec(dof_p_full) = (fp + fp_eb) * delta_t + Mpp * d_p_n ...
            + Kup.' * d_u_n - fp_b;
%     FVec(dof_p_full) = (fp + fp_eb) * delta_t + (Mpp + Kpp_fpl) * d_p_n ...
%             + Kup.' * d_u_n - fp_b;

    solution = KMat \ FVec;




    d_u = solution(dof_u_full);
    d_p = solution(dof_p_full);
    
    

    
    
    % updating of solutions of the previous time step 
    d_u_n = d_u;
    d_p_n = d_p;
    
    
    
    
    if mod(istep, 20) == 0 || istep == 1
        file_name = ['solutionmc_', num2str(istep), '.txt'];
        fid = fopen(file_name, 'w');
        fprintf(fid, "%20.14e\n", solution.');
        fclose(fid);
    end
    
    
    
    
    
    
    
    
    
    % displacement contour 
    if  istep == 1
            IndexDisp = 2;
            if strcmpi(MeshShape,'BiotQ9Q4') || strcmpi(MeshShape,'BiotQ4Q4')
                DispContourNMM_1(IndexDisp, d_u);
            else
                DispContourNMM(IndexDisp, d_u);
            end
            % pressure contour
            if strcmpi(MeshShape,'BiotQ9Q4') || strcmpi(MeshShape,'BiotQ4Q4')
                PressureContourNMM_1(d_p);
            else
                PressureContourNMM(d_p);
            end
    end
    fprintf("finished contours.\n");
    toc
    
  

end  % istep: end of time marching


save ('ManiElems.mat','ManiElems')
save ('PhyPatches.mat','PhyPatches')
save ('dof_u_full.mat','dof_u_full')
save ('dof_p_full.mat','dof_p_full')


