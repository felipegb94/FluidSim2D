function [] = writeFrame(Pos, Vel, numParticles, currStep)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

totalNumParticles = size(Pos,2);
headers = {'x', 'y','xVel','yVel','typeOfParticle'};
data = zeros(totalNumParticles,5);
filename = strcat('/Users/felipegb94/repos/FluidSim2D/data/state_step.csv.', num2str(currStep));
for i = 1:numParticles
    data(i,1) = Pos(1,i);
    data(i,2) = Pos(2,i);
    data(i,3) = Vel(1,i);
    data(i,4) = Vel(2,i);
    data(i,5) = 0;
end
for i = numParticles+1:totalNumParticles
    data(i,1) = Pos(1,i);
    data(i,2) = Pos(2,i);
    data(i,3) = Vel(1,i);
    data(i,4) = Vel(2,i);
    data(i,5) = 1;
end

csvwrite_with_headers(filename, data, headers);

end

