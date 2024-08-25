classdef PhysicalPatch < MathematicalPatch
    % the class of mathematical patch 
    
    properties
        ManiElem  % indices of manifold elements contained by it
        IndexMP   % index of the mathematical patch containing it 
%         % index of first displacement DOF of this patch in all displacement
%         % DOFs 
%         FDOF_u;     
%         % index of first pressure DOF of this patch in all pressure DOFs,
%         % zero if this patch doesn't contain pressure DOF 
%         FDOF_p;
%         % index of first displacement DOF of this patch in all DOFs 
%         FDOF_u_c;     
%         % index of first pressure DOF of this patch in all DOFs, zero if 
%         % this patch doesn't contain pressure DOF  
%         FDOF_p_c;
        % lists of displacement and pressure DOFs 
        DOF_u;   DOF_p; 
        % list of DOFs 
        DOF 
    end
    
    methods
        function PPobj = PhysicalPatch(xNode, IndexMP)
            if nargin == 0 
                xNode = [0, 0];
            end
            PPobj = PPobj@MathematicalPatch(xNode);
            if nargin >= 2 
                PPobj.IndexMP = IndexMP;
            end
           PPobj.DOF_u = [];
           PPobj.DOF_p = [];
           PPobj.DOF = [];
        end
        
        function set.DOF_u(obj, DOF_u) 
            obj.DOF_u = DOF_u;
        end
        function set.DOF_p(obj, DOF_p) 
            obj.DOF_p = DOF_p;
        end
        function set.DOF(obj, DOF) 
            obj.DOF = DOF;
        end
        
    end
end

