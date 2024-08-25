classdef PorousMaterial < handle
    % the class of porous material
    
    properties
        E      % Young's modulus 
        nu     % Poission's ratio 
        kappa  % hydraulic conductivity 
        M      % Biot's modulus 
        alpha   % Biot--Willis coefficient 
        
    end
    
    methods
        function obj = PorousMaterial(E, nu, kappa, M, alpha)
            if nargin > 0
                obj.E = E;
            end
            if nargin > 1
                obj.nu = nu;
            end
            if nargin > 2
                obj.kappa = kappa;
            end
            if nargin > 3
                obj.M = M;
            end
            if nargin > 4
                obj.alpha = alpha;
            end
        end
        
        function set.E(obj, E) 
            obj.E = E;
        end
        
    end  % methods 
end  % classdef 

