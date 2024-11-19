% The hydro-mechanical coupling within a soil slope containing several rock blocks
% is simulated to demonstrate the effectiveness of the proposed methodology
% for complex geotechical problems with numerous material interfaces.

clc;
clear;
% clear all;
close all;

InputFile = "input_Slope_Interface33.txt";




global MeshShape 
global NumPP 
global DOF_per_PP_U DOF_per_PP_P GDOF_U GDOF_P
global ManiElems PhyPatches
global dof_u_full  dof_p_full
global Mats
%global M alpha k mu delta_t current_t theta



%%
% % Parameters from Wenan Wu's paper on CMAME 
% E1 = 7.2;  % Young's modulus of soil skeleton
% E2=10E1;
% nu_1= 0.4;   % Poisson ratio of soil skeleton 
% nu_2 = 0.4;
% rho_s = 2000;   % mass density of soil grain 
% rho_w = 1000;   % mass density of pore water 
% n = 0.3;    % porosity of soil 
% k1 = 4.5e-7;    % Intrinsic permeability of soil
% k2 = 0.9e-7;
% nu1 = 1e-3;      % Dynamic viscosity of fluid 
% kappa_1 = k1/nu1;   % hydraulic conductivity 
% kappa_2 = k2/nu1;
% M = 2.8514e9;   % Biot's modulus 
% alpha = 1;    % Biot's coefficient



% to define materials 
% E1 = 3e4;     E2 = E1*5;   nu_1 = 0.2;   nu_2 = 0.2;
% % kappa_1 = 1e-1;      kappa_2 = kappa_1 * 10;   % hydraulic conductivity 
% kappa_1 = 1e-4;      kappa_2 = 5e-5; 
% M = 1e20;   alpha = 1; 
E1 = 5;     E2 = E1 * 0.2;   nu_1 = 0.2;   nu_2 = 0.2;
kappa_1 = 1e-1;      kappa_2 = kappa_1 * 10;   % hydraulic conductivity 
M = 1e2;   alpha = 1;       
% Mats(2, 1) = PorousMaterial; 
% Mats(1) = PorousMaterial(E1, nu_1, kappa_1, M, alpha);
% Mats(2) = PorousMaterial(E2, nu_2, kappa_2, M, alpha);
[Mats] = DefineSquareMaterial(E1, E2, nu_1, nu_2, kappa_1, kappa_2, M, alpha);



% % to define materials 
% E1 = 5;     E2 = 1e0;   nu_1 = 0.2;   nu_2 = 0.2;
% kappa_1 = 1e-1;      kappa_2 = 0.1;   % hydraulic conductivity 
% M = 1e2;   alpha = 1;       
% % Mats(2, 1) = PorousMaterial; 
% % Mats(1) = PorousMaterial(E1, nu_1, kappa_1, M, alpha);
% % Mats(2) = PorousMaterial(E2, nu_2, kappa_2, M, alpha);
% [Mats] = DefineSquareMaterial(E1, E2, nu_1, nu_2, kappa_1, kappa_2, M, alpha);

rho_s = 20;   % mass density of soil grain 
rho_w = 10;   % mass density of pore water 
g_gravity = 10;         % gravity acceleration

Rainfall = 0.4;              

delta_t = 0.1;              % time duration of each step 
% total_t = 1;
% nsteps = total_t/delta_t;         % number of time steps of the simulation 
nsteps = 100000;

% % paramter for the footing problem in my paper in Computational Geosciences
% E1 = 1e8;     E2 = E1 * 1;   nu_1 = 0;   nu_2 = 0;
% kappa_1 = 9.56e-9;      kappa_2 = kappa_1 * 1;   % hydraulic conductivity 
% M = 1e2;   alpha = 1;       
% % Mats(2, 1) = PorousMaterial; 
% % Mats(1) = PorousMaterial(E1, nu_1, kappa_1, M, alpha);
% % Mats(2) = PorousMaterial(E2, nu_2, kappa_2, M, alpha);
% [Mats] = DefineSquareMaterial(E1, E2, nu_1, nu_2, kappa_1, kappa_2, M, alpha);
% 
% delta_t = 1e-0;      % time duration of each step 
% total_t = 1e-2;
% nsteps = total_t/delta_t;         % number of time steps of the simulation 



% xyregion = [0    0 ;60     0;60      5 ;50      5;20     35; 0     35];
% xyregion=[xyregion;xyregion(1,:)];
% coordinates of the interface 
Inface = 7;
ModelInterface(Inface,1) = Interface;
% xyInterface1 = [ 2    2;6     2; 6     6; 2     6];
xyInterface7 = [ 17.000001    20.5; 22   24.5;20.8     26.000001; 17.8    25.000001 ];
xyInterface1 = [20.000001   16; 23.5   15; 24      15;  24    20; 22.5   20;20.000001     18.000001;20.000001   17];
xyInterface2 = [ 19    9.5; 20    9.5; 22    11 ; 21    13; 20    13;  18.000001   10];
% xyInterface3 = [ 25    20; 26    21 ; 27     21; 28    20;  28    19 ; 25    17];
xyInterface3 = [25.5    7.5; 32    9.5; 32.5     12.500001 ; 28     14 ;26     12.5];
xyInterface4 = [ 14    29; 16.000001   28 ;18     30.5; 15.5     32; 13.5   30.8 ];
xyInterface5 = [37.000001    12; 36    9.600001 ;38.5    7.5; 41.5     7.8; 42    10.5];
xyInterface6 = [12    27.000001; 9.5    26 ;11     24 ; 15     24 ;  16     26];
ModelInterface(1).xVertex = xyInterface1;
ModelInterface(2).xVertex = xyInterface2;
ModelInterface(3).xVertex = xyInterface3;
ModelInterface(4).xVertex = xyInterface4;
ModelInterface(5).xVertex = xyInterface5;
ModelInterface(6).xVertex = xyInterface6;
ModelInterface(7).xVertex = xyInterface7;
% ModelInterface(8).xVertex = xyInterface8;
% xyInterface11 = [xyInterface1;xyInterface1(1,:)];
% xyInterface21 = [xyInterface2;xyInterface2(1,:)];
% xyInterface31 = [xyInterface3;xyInterface3(1,:)];
% xyInterface41 = [xyInterface4;xyInterface4(1,:)];
% xyInterface51 = [xyInterface5;xyInterface5(1,:)];
% xyInterface61 = [xyInterface6;xyInterface6(1,:)];
% xyInterface71 = [xyInterface6;xyInterface7(1,:)];
% xyInterface81 = [xyInterface6;xyInterface8(1,:)];
% xInter = xyInterface(1, 1);
 




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
Main_preprocessing_SlopeInterface(InputFile,ModelInterface,Inface);
fprintf("finished preprocessing.\n");
toc

% line segments with prescribed fluid velocity  （具有规定流体速度的线段）
[nSegFluidNeumann, SegsFluidNeumann] = BCSegsSlope_FluidNeumann(Rainfall);


% [nBCSeg, BoundarySegs] = BoundarySegsDisc_4(xyregion);


%% ghost penalty parameters 
% gamma_g_u = 0.8 * lambda;    % ghost penalty parameter for the solid 
% gamma_g_m = 2 / M; % ghost penalty parameter for mass of fluid 
% gamma_g_p = 8 * kappa * delta_t;   % ghost penalty parameter for the fluid 
% 
% gamma_D_u = 400 * lambda; % Nitsche method parameter for displacement 
% gamma_D_p = 20 * kappa;   % Nitsche method parameter for pressure 

% gamma_g_u = 1 * lambda;    % ghost penalty parameter for the solid 
% gamma_g_m = 0.5 / M; % ghost penalty parameter for mass of fluid 
% gamma_g_p = 1 * kappa * delta_t;   % ghost penalty parameter for the fluid 
% 
% gamma_D_u = 400 * lambda; % Nitsche method parameter for displacement 
% gamma_D_p = 20 * kappa;   % Nitsche method parameter for pressure 


gamma_g_u = 1 * Mats(1).E;    % ghost penalty parameter for the solid 
gamma_g_m = 0.1 / Mats(1).M; % ghost penalty parameter for mass of fluid 
gamma_g_p = 1 * Mats(1).kappa * delta_t;   % ghost penalty parameter for the fluid 

gamma_D_u = 400 * Mats(1).E; % Nitsche method parameter for displacement, for the bondary 
gamma_D_p = 20 * Mats(1).kappa;   % Nitsche method parameter for pressure, for the bondary  

gamma_I_u = 400 * Mats(1).E; % Nitsche method parameter for displacement, for the interface
gamma_I_p = 500 * Mats(1).kappa;   % Nitsche method parameter for pressure, for the interface 


% lambda = E2 * nu_2 / (1 + nu_2) / (1 - 2 * nu_2);
% gamma_g_u = 0.5 * lambda;    % ghost penalty parameter for the solid 
% gamma_g_m = 0.1 / M; % ghost penalty parameter for mass of fluid 
% gamma_g_p = 0.4 * 0.1 * delta_t;   % ghost penalty parameter for the fluid 
% 
% gamma_D_u = 50 * lambda; % Nitsche method parameter for displacement 
% gamma_D_p = 10 * 0.1;   % Nitsche method parameter for pressure 






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
    d_u_n(PP.DOF_u) = 0;
    if PP.patchUP(2) == 0 
        continue;
    end
    
    if gxy(2) < 5 
        ana_p = 1e1 * (5 - gxy(2));
    else
        ana_p = 0;
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



% % fluid stiffness of inf-sup stabilization by the fluid pressure Laplacian
% % method 
% [Kpp_fpl]=CoeffMatFPL(delta_t);


% Additional terms of coefficient matrix in Nitsche's method for essential
% boundary condition for displacements 
[Kuu_eb] = KBC_u_Nitsche_Slope_ym(AnalysisType, gamma_D_u);
% Additional terms of coefficient matrix in Nitsche's method for essential
% boundary condition for pressure 
[Kpp_eb] = KBC_p_Nitsche_Slope_2ymg(gamma_D_p);
fprintf("finished K_eb.\n");
toc
% Additional coupling matrix on the essential boundary of displacement. 
% Used when Nitsche's method is used to get a symmetric weak form 
[Kup_b]=CoupleMatUBC_Slopeym();
fprintf("finished Kup_b.\n");
toc

% submatrices of the coefficient matrix, Kuu, Kpp, and Kup, corresponding to
% the interface term 
[Kuu_I] = KInterface_u_Nitsche_Slopeym(AnalysisType, gamma_I_u);
[Kpp_I] = KInterface_p_Nitsche_Slopeym(gamma_I_p);
[Kup_I] = KInterface_up_Nitsche_Slopeym();


% ghost penalty coefficient matrix for solid
[Kuu_ghost] = GhostPenalty2DSolid_Slope(gamma_g_u);
% ghost penalty coefficient matrix for fluid 
[Kpp_ghost, Mpp_ghost] = GhostPenalty2DFluidSlope_1(gamma_g_p, gamma_g_m);
fprintf("finished ghost penalty matrices.\n");
toc


% Symmetric Nitsche method for fully coupled problem 
KMat = zeros(GDOF_U + GDOF_P, GDOF_U + GDOF_P);
KMat(dof_u_full, dof_u_full) = -(Kuu + Kuu_eb + Kuu_I + Kuu_ghost);
KMat(dof_u_full, dof_p_full) = (Kup - Kup_b - Kup_I);
KMat(dof_p_full, dof_u_full) = (Kup - Kup_b - Kup_I).';
% KMat(dof_p_full, dof_p_full) = delta_t * (Kpp + Kpp_eb + Kpp_I) + Kpp_ghost ...
%     + Mpp + Mpp_ghost + Kpp_fpl;
KMat(dof_p_full, dof_p_full) = delta_t * (Kpp + Kpp_eb + Kpp_I) + Kpp_ghost ...
    + Mpp + Mpp_ghost;

% KMat = zeros(GDOF_U + GDOF_P, GDOF_U + GDOF_P);
% KMat(dof_u_full, dof_u_full) = -(Kuu + Kuu_eb + Kuu_ghost);
% KMat(dof_u_full, dof_p_full) = (Kup - Kup_b);
% KMat(dof_p_full, dof_u_full) = (Kup - Kup_b).';
% KMat(dof_p_full, dof_p_full) = delta_t * (Kpp + Kpp_eb) + Kpp_ghost ...
%     + Mpp + Mpp_ghost;




% KMat = zeros(GDOF_U + GDOF_P, GDOF_U + GDOF_P);
% KMat(dof_u_full, dof_u_full) = -(Kuu + Kuu_eb);
% KMat(dof_u_full, dof_p_full) = (Kup - Kup_b);
% KMat(dof_p_full, dof_u_full) = (Kup - Kup_b).';
% KMat(dof_p_full, dof_p_full) = delta_t * (Kpp + Kpp_eb) + Mpp ;


% KMat = zeros(GDOF_U + GDOF_P, GDOF_U + GDOF_P);
% KMat(dof_u_full, dof_u_full) = -(Kuu + Kuu_eb + Kuu_I);
% KMat(dof_u_full, dof_p_full) = (Kup - Kup_b - Kup_I);
% KMat(dof_p_full, dof_u_full) = (Kup - Kup_b - Kup_I).';
% KMat(dof_p_full, dof_p_full) = delta_t * (Kpp + Kpp_eb + Kpp_I) + Mpp;





% coefficient matrix without ghost penalty
KMat_us = zeros(GDOF_U + GDOF_P, GDOF_U + GDOF_P);
KMat_us(dof_u_full, dof_u_full) = -(Kuu + Kuu_eb + Kuu_I);
KMat_us(dof_u_full, dof_p_full) = (Kup - Kup_b - Kup_I);
KMat_us(dof_p_full, dof_u_full) = (Kup - Kup_b - Kup_I).';
% KMat_us(dof_p_full, dof_p_full) = delta_t * (Kpp + Kpp_eb + Kpp_I) + Mpp + Kpp_fpl;
KMat_us(dof_p_full, dof_p_full) = delta_t * (Kpp + Kpp_eb + Kpp_I) + Mpp;






% % % KMat = KMat_us;
% % % 
% kappa_s = cond(KMat, 2);
% kappa_us = cond(KMat_us, 2);
% % 
% % %kmat_us_eig = eig(KMat_us);
% % 
% % 
% fprintf("condition number of coefficient matrix: %14.6e\n", kappa_s);
% fprintf("condition number of coefficient matrix without ghost penalty: %14.6e\n", kappa_us);
% toc 
% % 
% % 
% % 
% FileID = fopen('Couple_interface.txt', 'w');
% fprintf(FileID, "step current_t GDOF Ed Ep:\n");






% right hand sides caused by sources and natural boundary terms
[fu, fp] = LoadUnitSlopeSequence_1(rho_s, g_gravity, ...
    nSegFluidNeumann, SegsFluidNeumann);


% Additional terms of RHS in Nitsche's method for essential
% boundary condition for pore pressure
fp_eb = zeros(GDOF_P, 1);
% [fp_eb] = FBC_p_Nitsche_Slope_2(gamma_D_p);
% Additional terms of RHS in Nitsche's method for essential
% boundary condition for displacements
fu_eb = zeros(GDOF_U, 1);
% Additional terms of RHS of the fluid on the EBC of solid for the
% symmetric Nitsche's method
fp_b = zeros(GDOF_P, 1);
% 
% Ed_0 = 0;    Ep_0 = 0;    Ee = 0;
% nsteps = 100;

for istep = 1 : 1 % time marching 
% for istep = 1 : 1   % time marching 
    
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
    


    
%     % collocation for fully coupled problem 
%     KMat = zeros(GDOF_U + GDOF_P, GDOF_U + GDOF_P);
%     KMat(dof_u_full, dof_u_full) = -(Kuu + Kuu_I);
%     KMat(dof_u_full, dof_p_full) = Kup - Kup_I;
%     KMat(dof_p_full, dof_u_full) = (Kup - Kup_I).';
%     KMat(dof_p_full, dof_p_full) = delta_t * (Kpp + Kpp_I) + Mpp;
%     
%     FVec = zeros(GDOF_U + GDOF_P, 1);
%     FVec(dof_u_full) = -fu;
%     FVec(dof_p_full) = fp * delta_t + Mpp * d_p_n ...
%             + Kup.' * d_u_n;
%    
%     if strcmpi(MeshShape,'BiotQ9Q4') 
%         [KNew,FNew]=DBC_FEM_UnitSquare_1(KMat, FVec, current_t);
%     else
%         [KNew,FNew]=DBC_FEM_UnitSquare(KMat, FVec, current_t);
%     end
%     solution = KNew \ FNew;
%     
%     kappa_s = cond(KNew, 2);
%     fprintf("condition number of coefficient matrix: %16.8e\n", kappa_s);




    d_u = solution(dof_u_full);
    d_p = solution(dof_p_full);
    
    
    
%     d_u1 = zeros(GDOF_U, 1);   d_p1 = zeros(GDOF_P, 1);
%     for ipp = 1 : NumPP 
%         PP = PhyPatches(ipp);
%         gxy = PP.xNode;
%         [ana_u, ana_p] = AnaUPUnitSquare(gxy, current_t, xInter, ...
%             E1, E2, kappa_1, kappa_2);
%         d_u1(PP.DOF_u) = ana_u;
%         if PP.patchUP(2) == 0 
%             continue;
%         end
%         d_p1(PP.DOF_p) = ana_p; 
%     end
%     solution(dof_u_full) = d_u1;
%     solution(dof_p_full) = d_p1;
    
    
    % updating of solutions of the previous time step 
    d_u_n = d_u;
    d_p_n = d_p;
    
    
    if mod(istep, 20) == 0 || istep == 1
        file_name = ['solutionspgd1520_', num2str(istep), '.txt'];
        fid = fopen(file_name, 'w');
        fprintf(fid, "%20.14e\n", solution.');
        fclose(fid);
    end
    
    
    
%     % displacement contour 
%         if istep==1  || istep==100 || istep==500 || istep==1000 || istep==2000 || istep==5000|| istep==10000            
%             IndexDisp = 2;
%             if strcmpi(MeshShape,'BiotQ9Q4') || strcmpi(MeshShape,'BiotQ4Q4')
%                 DispContourNMM_1(IndexDisp, d_u);
%             else
%                 DispContourNMM(IndexDisp, d_u);
%             end
%             IndexDisp = 1;
%             if strcmpi(MeshShape,'BiotQ9Q4') || strcmpi(MeshShape,'BiotQ4Q4')
%                 DispContourNMM_1(IndexDisp, d_u);
%             else
%                 DispContourNMM(IndexDisp, d_u);
%             end
%             % pressure contour
%             if strcmpi(MeshShape,'BiotQ9Q4') || strcmpi(MeshShape,'BiotQ4Q4')
%                 PressureContourNMM_1(d_p);
%             else
%                 PressureContourNMM(d_p);
%             end
%         end
%     fprintf("finished contours.\n");
%     toc
    
%     % denominator for ed 
% 	% [Denormed] = DenormEd_UnitSquare_1(current_t);
%     Denormed = 1;
%     % normalized L^2 norm of the error in displacement
%     [ed]=Ed_UnitSquare_2(d_u, current_t, Denormed, xInter);
%     % denorminator for ep 
% 	% [Denormep] = DenormEp_UnitSquare_1(current_t);
%     Denormep = 1;
%     % normlized L^2 norm of the error in pore pressure 
%     [ep]=Ep_UnitSquare_2(d_p, current_t, Denormep, xInter);
%     % denorminator for normalized error in energy norm 
%     %[Denormee] = DenormEe_UnitSquare_1(DMatU, current_t);
%     Denormee = 1;
%     % normalized error in energy norm 
%     [ee]=Ee_UnitSquare_2(d_u, current_t, AnalysisType, Denormee, xInter);
% %     fprintf("finished error norms.\n");
% %     toc
%     
%     %NumericSE = 1 / 2 * d_u.' * Kuu * d_u;
%     
% %     fprintf("GDOF, NumericSE, ed, ep, ee:\n");
% %     %fprintf("%8d %16.8e %16.8e %16.8e %16.8e\n", GDOF_U + GDOF_P, ...
% %     %    NumericSE, ed, ep, ee);
% %     fprintf("%8d %16.8e %16.8e %16.8e\n", GDOF_U + GDOF_P, ...
% %         ed, ep, ee);
% %     [GDOF_U + GDOF_P,ed, ep, ee]
% 
%     fprintf("%5d %12.4e %8d %16.8e %16.8e %16.8e\n", istep, current_t, ...
%         GDOF_U + GDOF_P, ed, ep, ee);
%     fprintf(FileID, "%5d %12.4e %8d %16.8e %16.8e %16.8e\n", istep, current_t, ...
%         GDOF_U + GDOF_P, ed, ep, ee);
%     
% %     [GDOF_U + GDOF_P, kappa_s, 0,  ed, ep, ee]
%     
% % %     [GDOF_U + GDOF_P, NumericSE, kappa_s, kappa_us,  kappa_us_m, Ed, Ep, Ee];
% %     [GDOF_U + GDOF_P, NumericSE, kappa_s, kappa_us,  0, Ed, Ep, Ee]
% % %     [GDOF_U + GDOF_P, NumericSE, kappa_s, 0,  Ed, Ep, Ee]
% 
% %     fprintf(FileID,"%5d %10.4e %8d %16.8e %16.8e %16.8e\n", istep, current_t, ...
% %         GDOF_U + GDOF_P, ed, ep, ee);
% 
% %     Ed_0 = Ed_0 + ed^2 * delta_t;
% %     Ep_0 = Ep_0 + ep^2 * delta_t;
% %     Ee = Ee + ee^2 * delta_t;

end  % istep: end of time marching


save ('ManiElems.mat','ManiElems')
save ('PhyPatches.mat','PhyPatches')
save ('dof_u_full.mat','dof_u_full')
save ('dof_p_full.mat','dof_p_full')




% Ed_0 = sqrt(Ed_0);  Ep_0 = sqrt(Ep_0);  Ee = sqrt(Ee);
% 
% fprintf("GDOF, Ed_0, Ep_0, Ee:\n");
% fprintf("%8d %16.8e %16.8e %16.8e\n", GDOF_U + GDOF_P, Ed_0, Ep_0, Ee);
% 
% fclose(FileID);
% 
% [GDOF_U + GDOF_P, Ed_0, Ep_0, Ee]
