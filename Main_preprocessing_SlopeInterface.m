function [] = Main_preprocessing_SlopeInterface(InputFile,ModelInterface,Inface)

global MeshShape NumElem CoorElemNode  ElemNodeFL 
global NumPP CoorPPStar PPIncludeElem 
global nGhostFace GhostFaces IfCutElem
global ManiElems PhyPatches
global nInterface Interfaces
global DOF_per_PP_U DOF_per_PP_P GDOF_U GDOF_P
global dof_u_full  dof_p_full

ReadInput(InputFile);

% [PreData]=Preprocess_NMM_2d();
% Preprocess_NMM_2d_cut();
Preprocess_NMM_2d_cut_1();


% xInter = xyInterface(1, 1);

GeoTol = 1e-12;

% material of each element 
for iele = 1 : NumElem 
    k=0;
    ME = ManiElems(iele);
    xy_poly = round(ME.xVertex,14);
    for i=1:Inface
        xyInterface = ModelInterface(i).xVertex;
        [in,on] = inpolygon(xy_poly(:,1),xy_poly(:,2),xyInterface(:,1),xyInterface(:,2));
        if sum(in)==size(xy_poly,1)
            ME.mat = 1;
            break
        else
            k=k+1;
        end
    end
    if k == Inface
        ME.mat = 2;
    end
  
end

% material interface 
[nInterface, Interfaces] = FormInterfaceUnitSquareym3(ModelInterface,Inface);





figure 
set(gca, 'color', 'none')

% to plot manifold elements 
if strcmpi(MeshShape,'BiotQ9Q4') || strcmpi(MeshShape,'Q9') || strcmpi(MeshShape,'BiotQ4Q4')
    for iele = 1 : NumElem
        xVertex = ManiElems(iele).xVertex;
        xVertex = [ xVertex; xVertex(1, :) ];
        plot(xVertex(:, 1), xVertex(:, 2), 'k');
        hold on 
    end 
else
    for i=1:NumElem
        plot(CoorElemNode(ElemNodeFL(i,1):ElemNodeFL(i,2)+1,1),...
             CoorElemNode(ElemNodeFL(i,1):ElemNodeFL(i,2)+1,2), 'k');
        hold on 
    end 
end

% to plot mathematical elements 
if strcmpi(MeshShape,'BiotQ9Q4') || strcmpi(MeshShape,'Q9') || strcmpi(MeshShape,'BiotQ4Q4')
    for iele = 1 : NumElem
        
        ME = ManiElems(iele);
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
        
        xPPs = [xPPs_p; xPPs_p(1, :) ];
        plot(xPPs(:, 1), xPPs(:, 2), 'k', 'LineWidth', 0.8);
        hold on 
        
    end 
else
    for i=1:NumElem
%         plot(CoorElemNode(ElemNodeFL(i,1):ElemNodeFL(i,2)+1,1),...
%              CoorElemNode(ElemNodeFL(i,1):ElemNodeFL(i,2)+1,2), 'k');
%         hold on 

        PPs = PPIncludeElem(i, :);
        xPPs = CoorPPStar(PPs, :);
        xPPs =[ xPPs; xPPs(1, :) ];
        plot(xPPs(:, 1), xPPs(:, 2), 'k', 'LineWidth', 0.8);
        hold on 
    end 
end

% plot physical domain 
ReadInput(InputFile);
global JointBoundLine

for iline = 1 : size(JointBoundLine, 1) 
    line_x = [ JointBoundLine(iline, 1);  JointBoundLine(iline, 3) ];
    line_y = [ JointBoundLine(iline, 2);  JointBoundLine(iline, 4) ];
    plot(line_x, line_y, 'k', 'LineWidth', 2);
    hold on
end

% temp_block_x = [ JointBoundLine(:, 1);  JointBoundLine(end, 3) ];
% temp_block_y = [ JointBoundLine(:, 2);  JointBoundLine(end, 4) ];
% plot(temp_block_x, temp_block_y, 'k', 'LineWidth', 2);

clear JointBoundLine NumHolePoint CoorHolePoint
clear SizeETri HSizeRTri VSizeRTri MinXCenterRTriMC MinYCenterRTriMC
clear HSizeQuad VSizeQuad MinXCenterQuadMC MinYCenterQuadMC




% physical patch star 
if strcmpi(MeshShape,'BiotQ9Q4') || strcmpi(MeshShape,'Q9') || strcmpi(MeshShape,'BiotQ4Q4')
    for ipp = 1 : NumPP 
        PP = PhyPatches(ipp);
        plot(PP.xNode(1), PP.xNode(2), 'bo', 'MarkerSize', 3, 'MarkerFaceColor', 'b');
        hold on
    end
else
    for ipp = 1 : size(CoorPPStar, 1)
        plot(CoorPPStar(ipp, 1), CoorPPStar(ipp, 2), 'bo', 'MarkerSize', 3, 'MarkerFaceColor', 'b');
        hold on
    end
end


for intf=1:nInterface
    xinterface = Interfaces(intf).xVertex;
    plot(xinterface(:, 1), xinterface(:, 2), 'r');
    hold on
end









axis equal 

% % label cut element 
% for iele = 1 : NumElem 
%     if IfCutElem(iele) == 0 
%         continue;
%     end
%     
%     xElem = CoorElemNode(ElemNodeFL(iele,1):ElemNodeFL(iele,2),1);
%     yElem = CoorElemNode(ElemNodeFL(iele,1):ElemNodeFL(iele,2),2);
%     patch(xElem, yElem, 'k');
% end 

% label ghost faces 
for iface = 1 : nGhostFace 
    xVertex = GhostFaces(iface).xVertex; 
    plot(xVertex(:, 1), xVertex(:, 2), 'r-', 'LineWidth', 2);
%     plot(xVertex(:, 1), xVertex(:, 2), 'LineWidth', 3);
end


%% DOF lists of each physical patches 
ndim = 2;   % spatial dimension 
DOF_per_PP_U = ndim;    % number of DOFs of each physical patch for displacement
DOF_per_PP_P = 1;       % number of DOFs of each physical patch for pressure

GDOF_U = 0;
GDOF_P = 0;
GDOF = 0;
% 1st colum: index of first displacement DOF of this patch in all displacement
% DOFs    
% 2nd column: index of first pressure DOF of this patch in all pressure DOFs,
% zero if this patch doesn't contain pressure DOF 
% 3rd column: index of first displacement DOF of this patch in all DOFs    
% 4th column: index of first pressure DOF of this patch in all DOFs, zero if 
% this patch doesn't contain pressure DOF  
% PPFirstDOF = zeros(NumPP, 4);

% indices of displacement and pressure dofs in all the dofs of the model 
dof_u_full = zeros(1, 2*NumPP);
dof_p_full = zeros(1, NumPP);
% for inode = 1 : NumPP 
%     dof_u_full(inode * 2 - 1 : inode * 2) = 3 * (inode-1) + [1; 2];
%     dof_p_full(inode) = 3 * inode; 
% end

for ipp = 1 : NumPP 
    PhyPatches(ipp).DOF_u = (1 : ndim) + GDOF_U;
%     PPFirstDOF(ipp, 1) = GDOF_U + 1;
%     PPFirstDOF(ipp, 3) = GDOF + 1;

    dof_u_full(GDOF_U + 1 : GDOF_U + ndim) = (1:ndim) + GDOF;
    GDOF_U = GDOF_U + ndim;
    
    PatchUP = PhyPatches(ipp).patchUP;
    if PatchUP(2) == 0 
        %PhyPatches(ipp).DOF_p = [];
        PhyPatches(ipp).DOF = (1 : ndim) + GDOF;
        GDOF = GDOF + ndim;
    else
        PhyPatches(ipp).DOF_p = GDOF_P + 1;
        dof_p_full(GDOF_P + 1) = GDOF + ndim + 1;
        GDOF_P = GDOF_P + 1;
        PhyPatches(ipp).DOF = (1 : (ndim + 1)) + GDOF;
        GDOF = GDOF + ndim + 1;
    end 
    
%     if PatchUP(2) == 0 
%         continue;
%     end
%     PPFirstDOF(ipp, 2) = GDOF_P + 1;
%     PPFirstDOF(ipp, 4) = GDOF + 1;
%     GDOF_P = GDOF_P + 1;
%     GDOF = GDOF + 1;
end   

dof_p_full = dof_p_full(1 : GDOF_P);


%% DOF lists of each manifold element 
for iele = 1 : NumElem
    PPs = ManiElems(iele).PP;
    
    DOFs_U = zeros(1, DOF_per_PP_U * length(PPs));
    for ipp = 1:length(PPs)
        DOFs_U( (ipp - 1) * DOF_per_PP_U + (1 : DOF_per_PP_U) ) = ...
            PhyPatches(PPs(ipp)).DOF_u;
    end
    
    ManiElems(iele).DOF_u = DOFs_U;
    
    DOFs_P = zeros(1, DOF_per_PP_P * length(PPs));
    nDOF_P = 0;
    for ipp = 1:length(PPs)
        if PhyPatches(PPs(ipp)).patchUP(2) == 0
            continue;
        end
        DOFs_P(nDOF_P + 1) = PhyPatches(PPs(ipp)).DOF_p;
        nDOF_P = nDOF_P + 1;
    end
    
    ManiElems(iele).DOF_p = DOFs_P(1 : nDOF_P);
end


% % physical patch star 
% if strcmpi(MeshShape,'BiotQ9Q4') || strcmpi(MeshShape,'Q9')
%     for ipp = 1 : NumPP 
%         PP = PhyPatches(ipp);
%         plot(PP.xNode(1), PP.xNode(2), 'ko', 'MarkerSize', 3, 'MarkerFaceColor', 'k');
%         hold on
%         
%         if PhyPatches(ipp).patchUP(2) 
%             plot(PP.xNode(1), PP.xNode(2), 'ko', 'MarkerSize', 10, 'MarkerFaceColor', 'none');
%             hold on
%         end
%             
%     end
% else
%     for ipp = 1 : size(CoorPPStar, 1)
%         plot(CoorPPStar(ipp, 1), CoorPPStar(ipp, 2), 'ko', 'MarkerSize', 3, 'MarkerFaceColor', 'k');
%         hold on
%         plot(CoorPPStar(ipp, 1), CoorPPStar(ipp, 2), 'ko', 'MarkerSize', 8, 'MarkerFaceColor', 'none');
%         hold on
%     end
% end
% 
% 
% axis equal 
    