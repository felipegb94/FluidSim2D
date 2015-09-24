params.boxWidth = 0.8; %Filename
params.boxHeight = 0.8; %Filename
params.numParticlesPerRow = 10; %Filename
params.numParticles = params.numParticlesPerRow*params.numParticlesPerRow; %Filename
params.numSteps = 1000; % Number of frames
params.dt = 1e-4; % Time step
params.h = 0.04; % Interaction radius
params.initialSeparation = params.h;
params.mass = 1; % Mass
params.particleMass = params.mass/params.numParticles; % Mass
params.rho0 = 1000; % Reference density
params.pressureConstant = 0.5;
params.mu = 0.001;
params.g = -9.81; % Gravity acceleration
params.epsilon = 1e-2;
params.firstLayer = (params.boxWidth/params.initialSeparation) + 1 + (2*params.boxHeight / params.initialSeparation);

numParticles = params.numParticles;
currStep = 0;
[Pos, Vel, VelHalf, Acc, Rho_RhoHalf_dRho] = initParticleSystem(params);
drawFrame(Pos, numParticles, params.boxWidth, params.boxHeight)
writeFrame(Pos, Vel, numParticles, currStep);
currStep = currStep + 1;


[Pos, Vel, VelHalf, Acc, Rho_RhoHalf_dRho] = firstStep(Pos, Vel, VelHalf, Acc, Rho_RhoHalf_dRho, params);
currStep = currStep + 1;


for i = 1:500
    fprintf(strcat('Curr Step = ', num2str(i), '\n'));
    [Pos, Vel, VelHalf, Acc, Rho_RhoHalf_dRho] = step(Pos, Vel, VelHalf, Acc, Rho_RhoHalf_dRho, params, currStep);
    if(mod(i,10) == 0)
        writeFrame(Pos, Vel, numParticles, currStep);
    end
    currStep = currStep + 1;
end


%Particles = setInitialDensities(Particles, params);












