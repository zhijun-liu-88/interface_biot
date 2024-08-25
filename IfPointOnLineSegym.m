


function [Flag]=IfPointOnLineSegym(Coor0,Coor1,Coor2)

% to judge whether or not a point P0 is on a line segment P1P2

% input 
% Coor0 records the coordinates of point in cartesian coordinate system, 1*2 vector, P0
% Coor1 records the coordinates of the first vertex of the line segment, 1*2 vector, P1
% Coor2 records the coordinates of the second vertex of the line segment, 1*2 vector, P2

% output 
% Flag is the indicator representing whether or not on the line segment, '1' for on and '0' for not 


Flag=0;

Coor01=round(Coor1-Coor0,14);     % vector P0P1
Coor02=round(Coor2-Coor0,14);     % vector P0P2

length01 = sqrt(sum(Coor01.^2));
length02 = sqrt(sum(Coor02.^2));
max_len = max(length01, length02);

if abs(Coor01(1)*Coor02(2)-Coor01(2)*Coor02(1))< 10^(-18) * max_len^2     % P0, P1 and P2 on the same line 
    
    if sum(Coor01.*Coor02) <= 10^(-18) * max_len^2      % P0 is between P1 and P2 
        
        Flag=1;
        
    end
    
end

end
    
    