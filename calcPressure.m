function [ p ] = calcPressure(rho, rho0)
%getPressure 
%   Calculates the pressure of a specific particle
%   p0 = pressure constant
cs2 = 400;
p0 = (cs2*rho0/7);
p = p0 * (((rho/rho0)^7) - 1);
end

