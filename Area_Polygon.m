% to compute the area of a polygon
% the polygon can be convex or concave 



% inpute 
% NumNodes is the number of nodes of the polygon
% CoorNodes records the coordinates of each node of the polygon, the first node is repeatedly
% recorded at last, CoorNodes is a (NumNodes+1)*2 matrix; these nodes are arranged in anticlockwise
% order


% output 
% Area is the area of the polygon



function [Area]=Area_Polygon(NumNodes,CoorNodes)

Area=0;

for i=1:NumNodes
    Area=Area+0.5*det([1, 0,                0;...
                       1, CoorNodes(i,1),   CoorNodes(i,2);...
                       1, CoorNodes(i+1,1), CoorNodes(i+1,2)]);
end