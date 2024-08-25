

function [Mats] = DefineMulMaterial(E1, E2,  E3, nu_1, nu_2, nu_3, kappa_1, kappa_2, kappa_3, M, alpha)

Mats(3, 1) = PorousMaterial;

Mats(1) = PorousMaterial(E1, nu_1, kappa_1, M, alpha);
Mats(2) = PorousMaterial(E2, nu_2, kappa_2, M, alpha);
Mats(3) = PorousMaterial(E3, nu_3, kappa_3, M, alpha);


