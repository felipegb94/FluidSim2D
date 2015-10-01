function [Pos, Vel, VelHalf, Acc, Rho_RhoHalf_dRho] = initParticleSystem(params)
%initParticleSystem 
%   This function initialize an array of particles.

disp('Initializing Particle System');


numParticles = params.numParticles;
numParticlesPerRow = params.numParticlesPerRow;

dr = params.initialSeparation;
boxWidth = params.boxWidth;
boxHeight = params.boxHeight;

numParticlesLowBound = (boxWidth/dr) + 7;

numParticlesSideBound = (boxHeight / dr) + 1;


numLowerBoundaryParticles = 3*numParticlesLowBound;
numSideBoundaryParticles = 6*numParticlesSideBound;

totalNumParticles = numParticles + numLowerBoundaryParticles + numSideBoundaryParticles;

Pos = zeros(2, totalNumParticles);
Vel = zeros(2, totalNumParticles);


VelHalf = zeros(2, totalNumParticles);


Acc = zeros(2, totalNumParticles);
Rho_RhoHalf_dRho = zeros(3, totalNumParticles);
Rho_RhoHalf_dRho(1:2,:) = params.rho0;


counter = 1;
% counter = 1;
% Pos(1,counter) = dr;
% Pos(2,counter) = dr;
% counter = counter + 1; 
% 
% Pos(1,counter) = 2*dr;
% Pos(2,counter) = dr;
% counter = counter + 1; 

% Place Particles
for x = dr: dr : dr*(numParticlesPerRow)  
    for y = dr: dr : dr*(numParticlesPerRow)
        Pos(1,counter) = x + 4*dr;
        Pos(2,counter) = y + 4*dr;
        counter = counter + 1; 
    end
end



% Place Lower boundary layers
for x = (-2*dr):dr:(boxWidth+2*dr)
    % Place Lower boundary first layer
    Pos(1,counter) = x;
    Pos(2,counter) = 0;
    counter = counter + 1;
    % Place Lower boundary second layer
    Pos(1,counter) = x;
    Pos(2,counter) = -dr;
    counter = counter + 1;
    % Place Lower boundary third layer
    Pos(1,counter) = x;
    Pos(2,counter) = -2*dr;
    counter = counter + 1;
end

%Place side boundary layers
for y = dr:dr:(boxHeight+2*dr)
    % Place left-most boundary layer (first)
    Pos(1,counter) = -2*dr;
    Pos(2,counter) = y;
    counter = counter + 1;
    
    % Place second boundary layer (second)
    Pos(1,counter) = -dr;
    Pos(2,counter) = y;
    counter = counter + 1;
    
    % Place third boundary layer 
    Pos(1,counter) = 0;
    Pos(2,counter) = y;
    counter = counter + 1;
    
    % Place fourht boundary layer 
    Pos(1,counter) = boxWidth;
    Pos(2,counter) = y;
    counter = counter + 1;
    
    % Place fifth boundary layer 
    Pos(1,counter) = boxWidth + dr;
    Pos(2,counter) = y;
    counter = counter + 1;

    % Place right-most boundary layer (sixth)
    Pos(1,counter) = boxWidth + 2*dr;
    Pos(2,counter) = y;
    counter = counter + 1;
    
end

