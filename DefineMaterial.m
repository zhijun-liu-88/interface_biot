

function [Mats] = DefineMaterial(E1, E2, nu_1, nu_2, kappa_1, kappa_2, M, alpha)

Mats(2, 1) = PorousMaterial;
% material of the left subdomian 
E = 1e0;          % Young's modulus 
% nu = 0.2;       % Poisson's ratio
nu = 0;
kappa = 0.1;    % hydraulic conductivity 
alpha = 1.0;    % Biot's coefficient
M = 100;
Mats(1) = PorousMaterial(E, nu, kappa, M, alpha);

% material of the right subdomian 
E = 709;          % Young's modulus 
% nu = 0.2;       % Poisson's raios
nu = 0;
kappa = 1e-3;    % hydraulic conductivity 
alpha = 1.0;    % Biot's coefficient
% E = 1e0;          % Young's modulus 
% nu = 0;
% kappa = 0.1;    % hydraulic conductivity 
% alpha = 1.0;    % Biot's coefficient
M = 100;
Mats(2) = PorousMaterial(E, nu, kappa, M, alpha);