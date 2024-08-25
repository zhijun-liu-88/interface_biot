
function []=DispContourNMM_1(IndexDisp, solution)
% to contour the displacement for the problem solved by 2d conventional NMM

global MeshShape
global NumElem 
global ManiElems PhyPatches


if strcmpi(MeshShape,'BiotQ9Q4') || strcmpi(MeshShape,'Q9')
    MeshShape_u = 'Q9'; 
elseif strcmpi(MeshShape,'BiotQ4Q4') || strcmpi(MeshShape,'Q4') || ...
        strcmpi(MeshShape,'Quad')
    MeshShape_u = 'Q4';  
elseif strcmpi(MeshShape,'BiotIRT3_RIRT3_R') || ...
        strcmpi(MeshShape,'ET3') || strcmpi(MeshShape,'IRT3_R')
    MeshShape_u = 'T3'; 
end

figure   

%% to contour the displacement in each element 
for iele = 1 : NumElem
    
    ME = ManiElems(iele);
    PPs = PhyPatches(ME.PP);
    
    DOFs_U = ME.DOF_u;
   
    xPPs_u = zeros(length(PPs), 2);
    for ipp = 1 : length(PPs)
        xPPs_u(ipp, :) = PPs(ipp).xNode;
    end
    
    xVertex = ME.xVertex;
    
    Disp = zeros(size(xVertex, 1), 1);
    
    for iVertex = 1 : size(xVertex, 1)
        
        gxy = xVertex(iVertex, :);
        
        
        [N] = NMatNMM2D_1(xPPs_u, gxy, MeshShape_u);
        

        Disp(iVertex) = N(IndexDisp, :) * solution(DOFs_U);

    end
    
    patch(xVertex(:, 1), xVertex(:, 2), Disp);
    hold on
        
end  % iele


%% figure title 
if IndexDisp==1
    title ('Contour of displacement u');
elseif IndexDisp==2
    title ('Contour of displacement v');
end

colorbar('vert');

colormap jet;

shading interp
axis equal

axis off 

set(gca, 'FontName', 'Times New Roman');
set(gca,'fontsize',20);