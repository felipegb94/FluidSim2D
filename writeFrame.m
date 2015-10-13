function [] = writeFrame(Pos, Vel, Rho_RhoHalf_dRho, numParticles, currStep)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

totalNumParticles = size(Pos,2);
headers = {'x', 'y','xVel','yVel','Rho', 'RhoHalf', 'dRho', 'typeOfParticle'};
data = zeros(totalNumParticles,5);
filename = strcat('C:\Users\Tim\Desktop\Fall_15\SBEL Lab\FluidSim2D_control\Data\Data.csv.', num2str(currStep));
for i = 1:numParticles
    data(i,1) = Pos(1,i);
    data(i,2) = Pos(2,i);
    data(i,3) = Vel(1,i);
    data(i,4) = Vel(2,i);
    data(i,5) = Rho_RhoHalf_dRho(1,i);
    data(i,6) = Rho_RhoHalf_dRho(2,i);
    data(i,7) = Rho_RhoHalf_dRho(3,i);
    data(i,8) = 0;
end
for i = numParticles+1:totalNumParticles
    data(i,1) = Pos(1,i);
    data(i,2) = Pos(2,i);
    data(i,3) = Vel(1,i);
    data(i,4) = Vel(2,i);
    data(i,5) = Rho_RhoHalf_dRho(1,i);
    data(i,6) = Rho_RhoHalf_dRho(2,i);
    data(i,7) = Rho_RhoHalf_dRho(3,i);
    data(i,8) = 1;
end

csvwrite_with_headers(filename, data, headers);

end

