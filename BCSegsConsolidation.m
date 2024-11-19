
function [nSegDisDirich_t, SegsDisDirich_t, nSegDisDirich_b, SegsDisDirich_b, ...
    nSegDisDirich_r, SegsDisDirich_r, nSegDisDirich_l, SegsDisDirich_l] = BCSegsConsolidation()
% to get the line segements with essential boundary conditions for the
% disc (the circular boundary) 
% input
% RainFall: intensity of the rainfall, m/s.

global NumElem 
global ManiElems


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

GeoTol = 1e-7;


% array of boundary segment objects 
SegsDisDirich_t(1) = BoundarySegment;  
nSegDisDirich_t = 0;   % number of boundary segments 
len_segs_ut = 0.0;

SegsDisDirich_b(1) = BoundarySegment;  
nSegDisDirich_b = 0;   % number of boundary segments 
len_segs_ub = 0.0;

SegsDisDirich_r(1) = BoundarySegment;  
nSegDisDirich_r = 0;   % number of boundary segments 
len_segs_ur = 0.0;

SegsDisDirich_l(1) = BoundarySegment;  
nSegDisDirich_l = 0;   % number of boundary segments 
len_segs_ul = 0.0;
%% top surface
% coordinates of two endpoints of the surface 
xP1 = [ 0.735  0.665 ];
xP2 = [ 0.665  0.735 ];
theta = pi/4;  % inclination of the surface 
n_vec = [ cos(theta);  sin(theta) ];

for iele = 1 : NumElem 
    
%     if iele == 85 
%         fprintf("check position.\n");
%     end
    
    ME = ManiElems(iele);
    
    xy_ele = ME.xVertex;
    xy_ele = [xy_ele; xy_ele(1, :) ];
    
    xMaxElem = max(xy_ele(:, 1));
    yMaxElem = max(xy_ele(:, 2));
    
    if xMaxElem < 0.665 - GeoTol || yMaxElem < 0.665 - GeoTol
        continue;
    end

    for inode = 1 : size(xy_ele, 1) - 1
        
        % to check if first point is on P1P2 
        xP0 = xy_ele(inode, :);
        x01 = xP1 - xP0;
        x02 = xP2 - xP0; 
       
        if abs(x01(1)*x02(2) - x01(2)*x02(1)) > 1e-10
            continue;
        end
        
        % to check if first point is on P1P2 
        xP0 = xy_ele(inode + 1, :);
        x01 = xP1 - xP0;
        x02 = xP2 - xP0; 
       
        if abs(x01(1)*x02(2) - x01(2)*x02(1)) > 1e-10
            continue;
        end
        
%         if ~IfPointOnLineSeg(xy_ele(inode, :), xP1, xP2) || ...
%             ~IfPointOnLineSeg(xy_ele(inode + 1, :), xP1, xP2)
%             continue;
%         end
        
        nSegDisDirich_t = nSegDisDirich_t + 1;
        xSeg = xy_ele(inode : inode + 1, :);
        SegsDisDirich_t(nSegDisDirich_t) = BoundarySegment(xSeg, iele);
        len_segs_ut = len_segs_ut + norm(xSeg(2, :) - xSeg(1, :));
        
        SegsDisDirich_t(nSegDisDirich_t).n = n_vec;
        SegsDisDirich_t(nSegDisDirich_t).LocalIndexDOF = 1;
    end  % inode 
end  % iele 

fprintf("len_segs_ut = %20.12e\n", len_segs_ut);

%% bottom surface
% coordinates of two endpoints of the surface 
xP1 = [ 0.035  -0.035 ];
xP2 = [ -0.035  0.035 ];
theta = -3*pi/4;  % inclination of the surface 
n_vec = [ cos(theta);  sin(theta) ];

for iele = 1 : NumElem 
    
%     if iele == 85 
%         fprintf("check position.\n");
%     end
    
    ME = ManiElems(iele);
    
    xy_ele = ME.xVertex;
    xy_ele = [xy_ele; xy_ele(1, :) ];
    
    xMaxElem = max(xy_ele(:, 1));
    yMaxElem = max(xy_ele(:, 2));
    
    if xMaxElem > 0.04 + GeoTol || yMaxElem > 0.04 + GeoTol
        continue;
    end

    for inode = 1 : size(xy_ele, 1) - 1
        
        % to check if first point is on P1P2 
        xP0 = xy_ele(inode, :);
        x01 = xP1 - xP0;
        x02 = xP2 - xP0; 
       
        if abs(x01(1)*x02(2) - x01(2)*x02(1)) > 1e-10
            continue;
        end
        
        % to check if first point is on P1P2 
        xP0 = xy_ele(inode + 1, :);
        x01 = xP1 - xP0;
        x02 = xP2 - xP0; 
       
        if abs(x01(1)*x02(2) - x01(2)*x02(1)) > 1e-10
            continue;
        end
        
%         if ~IfPointOnLineSeg(xy_ele(inode, :), xP1, xP2) || ...
%             ~IfPointOnLineSeg(xy_ele(inode + 1, :), xP1, xP2)
%             continue;
%         end
        
        nSegDisDirich_b = nSegDisDirich_b + 1;
        xSeg = xy_ele(inode : inode + 1, :);
        SegsDisDirich_b(nSegDisDirich_b) = BoundarySegment(xSeg, iele);
        len_segs_ub = len_segs_ub + norm(xSeg(2, :) - xSeg(1, :));
        
        SegsDisDirich_b(nSegDisDirich_b).n = n_vec;
        SegsDisDirich_b(nSegDisDirich_b).LocalIndexDOF = 1;
    end  % inode 
end  % iele 

fprintf("len_segs_ub = %20.12e\n", len_segs_ub);

%% right surface
% coordinates of two endpoints of the surface 
xP1 = [ 0.035  -0.035 ];
xP2 = [ 0.735  0.665 ];
theta = -pi/4;  % inclination of the surface 
n_vec = [ cos(theta);  sin(theta) ];

for iele = 1 : NumElem 
    
    ME = ManiElems(iele);
    
    xy_ele = ME.xVertex;
    xy_ele = [xy_ele; xy_ele(1, :) ];
    
    xMaxElem = max(xy_ele(:, 1));
    yMaxElem = max(xy_ele(:, 2));
    
    if xMaxElem < 0.035 - GeoTol || yMaxElem > 0.68 + GeoTol
        continue;
    end

    for inode = 1 : size(xy_ele, 1) - 1
        
        % to check if first point is on P1P2 
        xP0 = xy_ele(inode, :);
        x01 = xP1 - xP0;
        x02 = xP2 - xP0; 
       
        if abs(x01(1)*x02(2) - x01(2)*x02(1)) > 1e-10
            continue;
        end
        
        % to check if first point is on P1P2 
        xP0 = xy_ele(inode + 1, :);
        x01 = xP1 - xP0;
        x02 = xP2 - xP0; 
       
        if abs(x01(1)*x02(2) - x01(2)*x02(1)) > 1e-10
            continue;
        end
        
%         if ~IfPointOnLineSeg(xy_ele(inode, :), xP1, xP2) || ...
%             ~IfPointOnLineSeg(xy_ele(inode + 1, :), xP1, xP2)
%             continue;
%         end
        
        nSegDisDirich_r = nSegDisDirich_r + 1;
        xSeg = xy_ele(inode : inode + 1, :);
        SegsDisDirich_r(nSegDisDirich_r) = BoundarySegment(xSeg, iele);
        len_segs_ur = len_segs_ur + norm(xSeg(2, :) - xSeg(1, :));
        
        SegsDisDirich_r(nSegDisDirich_r).n = n_vec;
        SegsDisDirich_r(nSegDisDirich_r).LocalIndexDOF = 1;
    end  % inode 
end  % iele 

fprintf("len_segs_ur = %20.12e\n", len_segs_ur);

% figure 
% for iseg = 1 : nSegDisDirich_r 
%     xVertex = SegsDisDirich_r(iseg).xVertex;
%     plot(xVertex(:, 1), xVertex(:, 2), 'k');
%     hold on 
% end
% 
% title("Segments for Dirichlet boundary condition of displacement.");
% 
% axis equal 
% hold off 


%% left surface
% coordinates of two endpoints of the surface 
xP1 = [ -0.035  0.035 ];
xP2 = [ 0.665  0.735 ];
theta = 3*pi/4;  % inclination of the surface 
n_vec = [ cos(theta);  sin(theta) ];

for iele = 1 : NumElem 
    
    ME = ManiElems(iele);
    
    xy_ele = ME.xVertex;
    xy_ele = [xy_ele; xy_ele(1, :) ];
    
    xMaxElem = max(xy_ele(:, 1));
    yMaxElem = max(xy_ele(:, 2));
    
    if xMaxElem > 0.68 + GeoTol || yMaxElem < 0.035 - GeoTol
        continue;
    end

    for inode = 1 : size(xy_ele, 1) - 1
        
        % to check if first point is on P1P2 
        xP0 = xy_ele(inode, :);
        x01 = xP1 - xP0;
        x02 = xP2 - xP0; 
       
        if abs(x01(1)*x02(2) - x01(2)*x02(1)) > 1e-10
            continue;
        end
        
        % to check if first point is on P1P2 
        xP0 = xy_ele(inode + 1, :);
        x01 = xP1 - xP0;
        x02 = xP2 - xP0; 
       
        if abs(x01(1)*x02(2) - x01(2)*x02(1)) > 1e-10
            continue;
        end
        
%         if ~IfPointOnLineSeg(xy_ele(inode, :), xP1, xP2) || ...
%             ~IfPointOnLineSeg(xy_ele(inode + 1, :), xP1, xP2)
%             continue;
%         end
        
        nSegDisDirich_l = nSegDisDirich_l + 1;
        xSeg = xy_ele(inode : inode + 1, :);
        SegsDisDirich_l(nSegDisDirich_l) = BoundarySegment(xSeg, iele);
        len_segs_ul = len_segs_ul + norm(xSeg(2, :) - xSeg(1, :));
        
        SegsDisDirich_l(nSegDisDirich_l).n = n_vec;
        SegsDisDirich_l(nSegDisDirich_l).LocalIndexDOF = 1;
    end  % inode 
end  % iele 

fprintf("len_segs_ul = %20.12e\n", len_segs_ul);

% figure 
% for iseg = 1 : nSegDisDirich_l 
%     xVertex = SegsDisDirich_l(iseg).xVertex;
%     plot(xVertex(:, 1), xVertex(:, 2), 'k');
%     hold on 
% end
% 
% title("Segments for Dirichlet boundary condition of displacement.");
% 
% axis equal 
% hold off 