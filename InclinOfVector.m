% the function to compute the inclination of a vector, >=0, <2*pi


% input 
% CoorVector records the coordinates of the begining node and ending node of the vector, is a 1*4
% vector; CoorVecor(1,1:2) is the begining node, CoorVecotor(1,3:4) is the ending node


% output 
% Inclination is the inclination of the vector, [0,2*pi)


function[Inclination]=InclinOfVector(CoorVector)

if CoorVector(1)==CoorVector(3)
    
    if CoorVector(4)>CoorVector(2)
        Inclination=pi/2;
    else 
        Inclination=3*pi/2;
    end
    
elseif CoorVector(2)==CoorVector(4)
    
    if CoorVector(3)>CoorVector(1)
        Inclination=0;
    else
        Inclination=pi;
    end
    
elseif CoorVector(3)>CoorVector(1) && CoorVector(4)>CoorVector(2)
    
    Inclination=atan((CoorVector(4)-CoorVector(2))/(CoorVector(3)-CoorVector(1)));
    
elseif CoorVector(3)<CoorVector(1) && CoorVector(4)>CoorVector(2)
    
    Inclination=atan((CoorVector(4)-CoorVector(2))/(CoorVector(3)-CoorVector(1)))+pi;
    
elseif CoorVector(3)<CoorVector(1) && CoorVector(4)<CoorVector(2)
    
    Inclination=atan((CoorVector(4)-CoorVector(2))/(CoorVector(3)-CoorVector(1)))+pi;
    
elseif CoorVector(3)>CoorVector(1) && CoorVector(4)<CoorVector(2)
    
    Inclination=atan((CoorVector(4)-CoorVector(2))/(CoorVector(3)-CoorVector(1)))+2*pi;
    
end

end
    
    