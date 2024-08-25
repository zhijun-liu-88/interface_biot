classdef MathematicalPatch < handle
    % the class of mathematical patch 
    
    properties
        ndim;       % spatial dimension
        xNode;      % coodinates of node for interpolation 
        xVertex;    % coordinates of vertices on the boundary 
        patchUP     % flag of displacement DOF and pressure DOF 
    end
    
    methods
        function obj = MathematicalPatch(xNode)
            if nargin > 0
                if size(xNode, 2) ~= 2 && size(xNode, 2) ~= 3
                    error("dimensions must be 2 or 3.");
                end
                obj.xNode = xNode;
                obj.ndim = size(xNode, 2);
            end
            obj.patchUP = [0, 0];
        end
        
        function set.xVertex(obj, xVertex) 
            obj.xVertex = xVertex;
        end
        
        function set.patchUP(obj, patchUP) 
            obj.patchUP = patchUP;
        end
    end
end

