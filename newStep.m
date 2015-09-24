function [Pos, Vel, VelHalf, Acc, Rho_RhoHalf_dRho] = newStep(Pos, Vel, VelHalf, Acc, Rho_RhoHalf_dRho, params, currStep)
%step Summary of this function goes here
%   Detailed explanation goes here

    % Assign parameters
    numParticles = params.numParticles;
    dt = params.dt;
    dtHalf = 0.5*dt;
    totalNumParticles = size(Pos,2);
    rho0 = params.rho0;
    p0 = params.pressureConstant;
    g = params.g;
    mu = params.mu;
    h = params.h;
    h2 = h*h;
    epsilon = params.epsilon;
    particleMass = params.particleMass;
    

end

