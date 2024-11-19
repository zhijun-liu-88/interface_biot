
function [nSegFluidNeumann, SegsFluidNeumann] = BCSegsSlope_FluidNeumann(Rainfall)
% to get the line segements with essential boundary conditions for the
% disc (the circular boundary) 
% input
% RainFall: intensity of the rainfall, m/s.

global NumElem 
global ManiElems

velocity = [ 0  -Rainfall ];

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
SegsFluidNeumann(1) = BoundarySegment;  
nSegFluidNeumann = 0;   % number of boundary segments 
len_segs = 0.0;

%% top surface 
check_x = yMax;
check_dof = 2;
n_vec = [0;  1];

normal_velocity = velocity * n_vec;

for iele = 1 : NumElem 
    
    ME = ManiElems(iele);
    
    xy_ele = ME.xVertex;
    xy_ele = [xy_ele; xy_ele(1, :) ];
    
    yMaxElem = max(xy_ele(:, 2));
    
    if abs(yMaxElem - check_x) > GeoTol
        continue;
    end
    
    for inode = 1 : size(xy_ele, 1) - 1
        if abs(xy_ele(inode, check_dof) - check_x) > GeoTol || ...
            abs(xy_ele(inode + 1, check_dof) - check_x) > GeoTol
            continue;
        end
        
        nSegFluidNeumann = nSegFluidNeumann + 1;
        xSeg = xy_ele(inode : inode + 1, :);
        SegsFluidNeumann(nSegFluidNeumann) = BoundarySegment(xSeg, iele);
        len_segs = len_segs + norm(xSeg(2, :) - xSeg(1, :));
        
        SegsFluidNeumann(nSegFluidNeumann).n = n_vec;
        SegsFluidNeumann(nSegFluidNeumann).LocalIndexDOF = 1;
        SegsFluidNeumann(nSegFluidNeumann).BC_value = normal_velocity;
    end  % inode 
end  % iele 
    

%% inclined surface
% coordinates of two endpoints of the surface 
xP1 = [ 50  5 ];
xP2 = [ 20  35 ];
theta = pi/4;  % inclination of the surface 
n_vec = [ sin(theta);  cos(theta) ];
normal_velocity=velocity *n_vec;
for iele = 1 : NumElem 
    
%     if iele == 85 
%         fprintf("check position.\n");
%     end
    
    ME = ManiElems(iele);
    
    xy_ele = ME.xVertex;
    xy_ele = [xy_ele; xy_ele(1, :) ];
    
    xMaxElem = max(xy_ele(:, 1));
    yMaxElem = max(xy_ele(:, 2));
    
    if xMaxElem < 20 - GeoTol || yMaxElem < 5 - GeoTol
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
        
        nSegFluidNeumann = nSegFluidNeumann + 1;
        xSeg = xy_ele(inode : inode + 1, :);
        SegsFluidNeumann(nSegFluidNeumann) = BoundarySegment(xSeg, iele);
        len_segs = len_segs + norm(xSeg(2, :) - xSeg(1, :));
        
        SegsFluidNeumann(nSegFluidNeumann).n = n_vec;
        SegsFluidNeumann(nSegFluidNeumann).LocalIndexDOF = 1;
        SegsFluidNeumann(nSegFluidNeumann).BC_value = normal_velocity;
    end  % inode 
end  % iele 
        
        
fprintf("len_segs = %20.12e\n", len_segs);

figure 
for iseg = 1 : nSegFluidNeumann 
    xVertex = SegsFluidNeumann(iseg).xVertex;
    plot(xVertex(:, 1), xVertex(:, 2), 'k');
    hold on 
end

title("Segments for Neumann boundary condition of fluid.");

axis equal 
hold off 
    


