classdef MathematicalElement < handle
    % the class of mathematical element 
    
    properties
        Type    % element type 
        xNode  % coordinates of nodes for interpolation 
        nNode  % number of nodes for interpolation 
        MP     % indices of mathematical patches (nodes) containing it 
        xVertex % coordinates of vertices of the mathematical element
    end
    
    methods
        function obj = MathematicalElement(Type, MP)
            if nargin > 0
                obj.Type = Type;
                obj.MP = MP;
            end
        end
        
        function set.xNode(obj, xNode) 
            obj.xNode = xNode;
        end
        
        function set.xVertex(obj, xVertex) 
            obj.xVertex = xVertex;
        end
    end
end

