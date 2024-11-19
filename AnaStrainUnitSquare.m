function [ana_e] = AnaStrainUnitSquare(gxy, t, xInter, E1, E2)
% analytical strain for a 2D consolidation





if gxy(1) <= xInter
    exx = gxy(1)^2;
else
    exx = gxy(1)^2 + (E1 - E2)/E2 * xInter^2;
end
exx = exx * exp(-t);
ana_e = [exx;  0;  0];