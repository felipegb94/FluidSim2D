
params.numParticlesPerRow = 10; %Filename
params.numParticles = params.numParticlesPerRow*params.numParticlesPerRow; %Filename
params.numSteps = 1000; % Number of frames
params.dt = 1e-4; % Time step
params.h = 0.04; % Interaction radius
params.initialSeparation = params.h;
params.rho0 = 1000; % Reference density
params.particleMass = (params.initialSeparation^2) * params.rho0;
params.pressureConstant = 0.5;
params.mu = 0.001;
params.g = -9.81; % Gravity acceleration
params.epsilon = 1e-2;
params.boxWidth = 20 * params.h; %Filename
params.boxHeight = 20 * params.h; %Filename

numParticles = params.numParticles;
currStep = 0;
[Pos, Vel, VelHalf, Acc, Rho_RhoHalf_dRho] = initParticleSystem(params);
drawFrame(Pos, numParticles, params.boxWidth, params.boxHeight);

% [ sum ] = checkKernel( Pos, Rho_RhoHalf_dRho, numParticles, 16, params.particleMass, params.h )
% [ sum ] = checkKernel( Pos, Rho_RhoHalf_dRho, numParticles, 15, params.particleMass, params.h )
% [ sum ] = checkKernel( Pos, Rho_RhoHalf_dRho, numParticles, 34, params.particleMass, params.h )
% [ sum ] = checkKernel( Pos, Rho_RhoHalf_dRho, numParticles, 1, params.particleMass, params.h )


writeFrame(Pos, Vel, Rho_RhoHalf_dRho, numParticles, currStep);
currStep = currStep + 1;


[Pos, Vel, VelHalf, Acc, Rho_RhoHalf_dRho] = step(Pos, Vel, VelHalf, Acc, Rho_RhoHalf_dRho, params, 0.5);
currStep = currStep + 1;

for i = currStep:5000
    fprintf(strcat('Curr Step = ', num2str(i), '\n'));
    [Pos, Vel, VelHalf, Acc, Rho_RhoHalf_dRho] = step(Pos, Vel, VelHalf, Acc, Rho_RhoHalf_dRho, params, 1);
    if(mod(i,10) == 0)
        writeFrame(Pos, Vel, Rho_RhoHalf_dRho, numParticles, i);
    end
end













