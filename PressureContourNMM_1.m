
function []=PressureContourNMM_1(solution)
% to draw the solution (displacement or pore pressure) contour for a 2D
% Biot's consolidation problem. 

global MeshShape
global NumElem 
global ManiElems PhyPatches
global DOF_per_PP_U DOF_per_PP_P 


if strcmpi(MeshShape,'BiotQ9Q4') || strcmpi(MeshShape,'BiotQ4Q4') || ...
        strcmpi(MeshShape,'Q4') || strcmpi(MeshShape,'Quad')
    MeshShape_p = 'Q4';
elseif strcmpi(MeshShape,'BiotIRT3_RIRT3_R') || ...
        strcmpi(MeshShape,'ET3') || strcmpi(MeshShape,'IRT3_R')
    MeshShape_p = 'T3';
end

figure   

%% to contour the displacement in each element 
for iele = 1 : NumElem
    
    ME = ManiElems(iele);
    PPs = PhyPatches(ME.PP);
    
    DOFs_P = ME.DOF_p;
    
    xPPs_u = zeros(length(PPs), 2);
    for ipp = 1 : length(PPs)
        xPPs_u(ipp, :) = PPs(ipp).xNode;
    end
    if strcmpi(MeshShape,'BiotQ9Q4') || strcmpi(MeshShape,'BiotQ4Q4')
        xPPs_p = xPPs_u(1 : 4, :);
    elseif strcmpi(MeshShape,'BiotIRT3_RIRT3_R')
        xPPs_p = xPPs_u(1 : 3, :);
    else 
        error("Element not implemented.");
    end
    
    xVertex = ME.xVertex;
    
    % pressure at vertices
    VertexP = zeros(size(xVertex, 1), 1);
    ana_pp = zeros(size(xVertex, 1), 1);
    for iVertex = 1 : size(xVertex, 1)
        
        gxy = xVertex(iVertex, :);
        
        [N] = NMatNMM2D_1(xPPs_p, gxy, MeshShape_p);
        N = N(1, 1 : 2 : end);

        VertexP(iVertex) = N * solution(DOFs_P);

    end
    
    patch(xVertex(:, 1), xVertex(:, 2), VertexP);
    hold on
        
end  % iele


%% figure title 
title ('Contour of pore pressure p');

colorbar('vert');

colormap jet;

shading interp
axis equal

axis off 

set(gca, 'FontName', 'Times New Roman');
set(gca,'fontsize',20);