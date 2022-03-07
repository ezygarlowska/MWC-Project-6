function [As Sk] = empirical_fun(Ur)
%The following equations are the 7.3 and 7.4 from practicals manual

B = 0.857/(1+exp((-0.471-log10(Ur))/(0.297)));
psi = -90+90*tanh(0.815/(Ur.^0.672)); %It will give the result in degrees

As = B*sind(psi);
Sk = B*cosd(psi);
end
