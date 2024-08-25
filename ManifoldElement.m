classdef ManifoldElement < handle
    % the class of manifold element 
    
    properties
        ndim;           % spatial dimension
        xVertex % coordinates of vertices
        PP     % indices of physical patches containing it
        % indices of edge-based smoothing domains containing part of it
        EdgeSD 
        % lists of displacement and pressure DOFs 
        DOF_u;   DOF_p; 
        % list of DOFs 
        %DOF 
        % numbers of PPs for containing u and p DOFs, respectively
        nPP_u;  nPP_p;
        mat  % index of the material type 
    end
    
    methods
        function obj = ManifoldElement(xVertex)
            if nargin > 0 
                if size(xVertex, 2) ~= 2 && size(xVertex, 2) ~= 3
                    error("dimensions must be 2 or 3.");
                end
                obj.xVertex = xVertex;
                obj.ndim = size(xVertex, 2);
            end 
            obj.EdgeSD = [];
        end
        
        function set.PP(obj, IndexPP) 
            obj.PP = IndexPP;
        end
        
        function AddEdgeSD(obj, SDs)
            obj.EdgeSD = union(obj.EdgeSD, SDs);
        end
        
        function set.DOF_u(obj, DOF_u) 
            obj.DOF_u = DOF_u;
        end
        function set.DOF_p(obj, DOF_p) 
            obj.DOF_p = DOF_p;
        end
        function set.mat(obj, mat_id) 
            obj.mat = mat_id;
        end
        
    end
end

