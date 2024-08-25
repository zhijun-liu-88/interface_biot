function [DMat] = ConstitutiveMat(E, nu, AnalysisType)



if strcmpi(AnalysisType, "plane_stress")
    DMat = E / (1 - nu.^2) *[1 nu 0; nu 1 0; 0 0 (1-nu)/2];
elseif strcmpi(AnalysisType, "plane_strain")
    DMat= E/(1+nu)/(1-2*nu)*[1-nu, nu,   0;...
                                    nu,   1-nu, 0;...
                                    0,    0,    0.5-nu]; 
elseif strcmpi(AnalysisType, "3d_analysis")
    % Lame's constants 
    lamda = E * nu / (1 + nu) / (1 - 2 * nu);
    mu = E / 2 / (1 + nu);

    DMat = [ lamda + 2 * mu,  lamda,  lamda,  0,   0,   0;
             lamda,  lamda + 2 * mu,  lamda,  0,   0,   0;
             lamda,  lamda,  lamda + 2 * mu,  0,   0,   0;
             0,             0,            0,  mu,  0,   0;
             0,             0,            0,  0,   mu,  0;
             0,             0,            0,  0,   0,   mu  ];
end
    