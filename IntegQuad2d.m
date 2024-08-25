function [NumIntegQuad, xyIntegQuadVect] = IntegQuad2d(xy_poly)
% to form integration quadrilaterals for 2d Gauss integration in a polygon
% (convex or concave)
% input 
% xy_poly: coordinates of vertices of the polygon. The input is a closed
% polyon, which means the last vertex and the first vertex are the same. 
% output
% NumIntegQuad: number of quadrilaterals for integration 
% xyIntegQuadVect(:, :, k): coordinates of 4 vertices ofk-th quadrilateral 
% for integration.  




MaxX_poly=max(xy_poly(:, 1));
MinX_poly=min(xy_poly(:, 1));
MaxY_poly=max(xy_poly(:, 2));
MinY_poly=min(xy_poly(:, 2));

IfRect=1;     % indicator of the polygon being a rectangle

GeoTol = 1.e-14;

for k = 1 : size(xy_poly, 1) - 1
    if      ( abs( xy_poly(k,1)-MaxX_poly )<GeoTol           || abs( xy_poly(k,1)-MinX_poly )<GeoTol ) && ...
            ( abs( xy_poly(k,2)-MaxY_poly )<GeoTol           || abs( xy_poly(k,2)-MinY_poly )<GeoTol ) && ...
            ( abs( xy_poly(k,1)-xy_poly(k+1,1) )<GeoTol || abs( xy_poly(k,2)-xy_poly(k+1,2) )<GeoTol )
        continue;
    else 
        IfRect=0;
        break;
    end
end

if IfRect     % rectangular 

    NumIntegQuad=1;

    xyIntegQuadVect=zeros(4,2,NumIntegQuad);

    xyIntegQuadVect(:,:,1)=[MinX_poly,MinY_poly;...
                            MaxX_poly,MinY_poly;
                            MaxX_poly,MaxY_poly;
                            MinX_poly,MaxY_poly];
else                                                        % not rectangular 

%     [Tris] = polygon_triangulation(xy_poly);
% 
%     xyIntegQuadVect = Tris;
%     NumIntegQuad = size(Tris, 3);

    
    xy_poly = xy_poly(1 : end - 1, :);
    Tris = triangulation(polyshape(xy_poly));
    
    NumIntegQuad = size(Tris.ConnectivityList, 1);
    
    xyIntegQuadVect = zeros(4, 2, NumIntegQuad);  
    for i = 1 : NumIntegQuad
        vertices = Tris.ConnectivityList(i, :);
        vertices = [ vertices, vertices(1) ];
        xyIntegQuadVect(:, :, i) = Tris.Points(vertices, :);
    end
    
end  % if IfRect