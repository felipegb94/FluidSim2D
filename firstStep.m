function [Pos, Vel, VelHalf, Acc, Rho_RhoHalf_dRho] = firstStep(Pos, Vel, VelHalf, Acc, Rho_RhoHalf_dRho, params)
%firstStep Summary of this function goes here
%   Detailed explanation goes here

%% Calculate accelerations first
numParticles = params.numParticles;
dt = params.dt;
totalNumParticles = size(Pos,2);
rho0 = params.rho0;
p0 = params.pressureConstant;
g = params.g;
mu = params.mu;
h = params.h;
epsilon = params.epsilon;
particleMass = params.particleMass;

PosNext = zeros(size(Pos,1),numParticles);
Rho_RhoHalf_dRhoNext = zeros(size(Rho_RhoHalf_dRho,1),numParticles);

for i = 1:numParticles
    xAcc = 0;
    yAcc = 0;
    rho_i = Rho_RhoHalf_dRho(1,i);
    drho_i = 0;
    p_i = calcPressure(rho_i, rho0, p0);
    x_i = Pos(1,i);
    y_i = Pos(2,i);
    vx_i = Vel(1,i);
    vy_i = Vel(2,i); 
    for j = 1:totalNumParticles
        x_j = Pos(1,j);
        y_j = Pos(2,j);
        dx = x_i - x_j;
        dy = y_i - y_j;
        r_ij = [dx; dy];
        norm_r_ij = norm(r_ij);
        q = norm_r_ij/h;
        if ((i ~= j) && (q<=2) && (q>=0))
            
            rho_j = Rho_RhoHalf_dRho(1,j);
            p_j = 10000;
            if j <= numParticles
                p_j = calcPressure(rho_j, rho0, p0);
            end
            
            vx_j = Vel(1,j);
            vy_j = Vel(2,j); 
       
            dvx = vx_i - vx_j;
            dvy = vy_i - vy_j;
       
            v_ij = [dvx; dvy];
       
            dw_ij = dW(q, h); 
            grad_a_wab = (dw_ij/(h*h*q)) .* r_ij;
            
            rho_bar = (rho_i+rho_j)/2;
            % Everything is a scalar except unit_r_ij
            pressureTerm = (particleMass*((p_i/(rho_i*rho_i)) + (p_j/(rho_j*rho_j))))...
                        .* grad_a_wab;
                 
            muNumerator = particleMass * 2*mu .* (r_ij' * grad_a_wab) .* v_ij;
            muDenominator = (rho_bar*rho_bar)*(norm_r_ij*norm_r_ij + h*h*epsilon);
            muTerm = (muNumerator./muDenominator);
                   
            xAcc = xAcc - pressureTerm(1) + muTerm(1);
            yAcc = yAcc - pressureTerm(2) + muTerm(2);
        
            drho_i = drho_i + (rho_i/rho_j)*particleMass*v_ij'*grad_a_wab;

       
       end
       
    end
    Acc(1,i) = xAcc;
    Acc(2,i) = yAcc + g;
    VelHalf(1,i) = VelHalf(1,i) + 0.5*dt*xAcc;
    VelHalf(2,i) = VelHalf(2,i) + 0.5*dt*yAcc;
    % rhoHalf
    Rho_RhoHalf_dRhoNext(2,i) = Rho_RhoHalf_dRho(1,i) + 0.5*dt*drho_i;
    % dRho
    Rho_RhoHalf_dRhoNext(3,i) = drho_i;
    PosNext(1,i) = Pos(1,i) + dt * VelHalf(1,i);
    PosNext(2,i) = Pos(2,i) + dt * VelHalf(1,i);
end

Pos(:,1:numParticles) = PosNext;
Rho_RhoHalf_dRho(:,1:numParticles) = Rho_RhoHalf_dRhoNext;

end

