classdef Face_Gamma < handle
    % the class of interior faces (common face of mathematical elements) of
    % cut element on which the ghost penalty term are defined. 
    % see the following papers: 1. A Stabilized Nitsche Fictitious Domain 
    % Method for the Stokes Problem, 2. CutFEM: Discretizing geometry and 
    % partial differential equations, 3. Fictitious domain finite element 
    % methods using cut elements: II. A stabilized Nitsche method, 4. other
    % papers of Erik Burman 
    
    properties
        ndim;       % spatial dimension
        % 2*2, coordinates of two endpoints for 2D face; coordinates of
        % boundary vertices in anticlockwise order for 3D face, nVertex*3. 
        xVertex;    
        ManiElem;  % indices of manifold elements sharing the face 
        n;          % outward unit normal to the first element 
    end
    
    methods
        function obj = Face_Gamma(xVertex, ManiElem)
            if nargin > 1
                if size(xVertex, 2) ~= 2 && size(xVertex, 2) ~= 3
                    error("dimensions must be 2 or 3.");
                end
                obj.xVertex = xVertex;
                obj.ManiElem = ManiElem;
                obj.ndim = size(xVertex, 2);
            end
        end
        
        % to get the outward unit normal to the first element 
        function compute_n_2D(obj, xVertex1) 
            % xVertex1: the coordinates of vertices of the first 
            % element, in anticlockwise order 
            p1p2 = obj.xVertex(2, :) - obj.xVertex(1, :);
            nVertex1 = size(xVertex1, 1);
            for ivertex = 1 : nVertex1
                p3 = xVertex1(ivertex, :);  
                if ivertex < nVertex1
                    p4 = xVertex1(ivertex + 1, :);
                else
                    p4 = xVertex1(1, :);
                end
                p3p4 = p4 - p3;
                if norm(p3p4) < 1e-14 
                    continue;
                end

                % cross product of p1p2 and p3p4 
                cross_product = p1p2(1) * p3p4(2) - p1p2(2) * p3p4(1);
                if abs(cross_product) > 1e-15  % not parallel 
                    continue;
                end

                p2p3 = p3 - obj.xVertex(2, :);
                % cross product of p1p2 and p2p3 
                cross_product = p1p2(1) * p2p3(2) - p1p2(2) * p2p3(1);
                if abs(cross_product) > 1e-15 
                    continue;
                end

                obj.n = [ p3p4(2), -p3p4(1) ] / norm(p3p4);
                break;
            end  % ivertex 
        end % compute_n_2D
        
        
    end  % methods 
end  % classdef 

