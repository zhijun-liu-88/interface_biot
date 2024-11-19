classdef BoundarySegment < handle
    % the class of 2D boundary segments 
    
    properties
        ndim;       % spatial dimension
        xVertex;    % 2*2, coordinates of two endpoints 
        ManiElem;  % index of the manifold element containing the segment 
        n;          % outward unit normal
        LocalIndexDOF  % local index of the constrained DOF 
        BC_value    % value of the constrained DOF 
        
    end
    
    methods
        function obj = BoundarySegment(xVertex, ManiElem)
            if nargin > 1
                if size(xVertex, 2) ~= 2 && size(xVertex, 2) ~= 3
                    error("dimensions must be 2 or 3.");
                end
                obj.xVertex = xVertex;
                obj.ManiElem = ManiElem;
                obj.ndim = size(xVertex, 2);
            end
        end
        
        function set.LocalIndexDOF(obj, ConstraintDOF) 
            obj.LocalIndexDOF = ConstraintDOF;
        end
        
        function set.BC_value(obj, BC_value) 
            obj.BC_value = BC_value;
        end
        
        function set.n(obj, n_vec) 
            obj.n = n_vec;
        end
        
        function compute_n(obj)
            segment_vector = obj.xVertex(2, :) - obj.xVertex(1, :);
            normal_vector = [ segment_vector(2)  -segment_vector(1) ];
            obj.n = normal_vector / norm(normal_vector);
        end 
        
    end  % methods 
end  % classdef 

