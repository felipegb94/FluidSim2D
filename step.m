function [Pos, Vel, VelHalf, Acc, Rho_RhoHalf_dRho] = step(Pos, Vel, VelHalf, Acc, Rho_RhoHalf_dRho, params, stepTerm)
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
    

    NL = GetNeighborList(Pos, params.h, params.boxWidth, params.boxHeight, params.numParticles);
   

    % Update velocity and rho at t using acceleration at t-dt and 
    % vel and rho at t-0.5*dt.
    % This is not supposed to happend in the first step, however all of the
    % values are 0 at that point so it does not really matter if we run
    % thiis part of the code. An if statement would be an overkill because
    % it would be evaluated at every integration step.
    for i = 1:numParticles
        Vel(1,i) = VelHalf(1,i) + dtHalf*Acc(1,i);
        Vel(2,i) = VelHalf(2,i) + dtHalf*Acc(2,i);
        Rho_RhoHalf_dRho(1,i) = Rho_RhoHalf_dRho(2,i) + dtHalf*Rho_RhoHalf_dRho(3,i); 
    end
    
    % Create temporary arrays for position and density vector. There is no
    % need to do this for Vel or Acc vector because, the Acc vector is not
    % used in the step loop and VelHalf is the vector that is being update
    % in the step loop and the Vel vector is the one used
    PosNext = Pos;
    Rho_RhoHalf_dRhoNext = Rho_RhoHalf_dRho;
    
    Pressures = zeros(1, totalNumParticles);
    for i = 1:numParticles
        rho_i = Rho_RhoHalf_dRho(1,i);
        Pressures(1, i) = calcPressure(rho_i, rho0);
    end
    Pressures(1, numParticles+1:totalNumParticles) = 1000;
    
    %Loops through all of the fluid particles
    for i = 1:numParticles
        xAcc_i = 0; % Reset values for every particle
        yAcc_i = 0;
        rho_i = Rho_RhoHalf_dRho(1,i);
        drho_i = 0;
        p_i = Pressures(1, i);
        x_i = Pos(1,i);
        y_i = Pos(2,i);
        vx_i = Vel(1,i);
        vy_i = Vel(2,i); 
        
        Neighbor = NL{i,1};
        [t, m] = size(Neighbor);
        for k = 2:t
        j = Neighbor(k);
        
        %Loops through all of the neighbors to the fluid
%          for j = 1:totalNumParticles
            
            x_j = Pos(1,j);
            y_j = Pos(2,j);
            dx = x_i - x_j;
            dy = y_i - y_j;
            r_ij = [dx; dy];
            norm_r_ij = norm(r_ij);
            q = norm_r_ij/h;
%             if ((i ~= j) && (q<=2) && (q>=0))
                
                rho_j = Rho_RhoHalf_dRho(1,j);
                p_j = Pressures(1, j);
            
                vx_j = Vel(1,j);
                vy_j = Vel(2,j); 
       
                dvx = vx_i - vx_j;
                dvy = vy_i - vy_j;
       
                v_ij = [dvx; dvy];
       
                dw_ij = dW(q, h); 
                grad_a_wab = (dw_ij/(h2*q)) * r_ij;
            
                rho_bar = (rho_i+rho_j)/2;
                % Everything is a scalar except unit_r_ij
                pressureTerm = (particleMass*((p_i/(rho_i*rho_i)) + (p_j/(rho_j*rho_j))))...
                            .* grad_a_wab;
                 
                muNumerator = particleMass * 2*mu .* (r_ij' * grad_a_wab) .* v_ij;
                muDenominator = (rho_bar*rho_bar)*(norm_r_ij*norm_r_ij + h2*epsilon);
                muTerm = (muNumerator./muDenominator);
                xAcc_i = xAcc_i - pressureTerm(1) + muTerm(1);
                yAcc_i = yAcc_i - pressureTerm(2) + muTerm(2);
        
                drho_i = drho_i + (rho_i/rho_j)*particleMass*v_ij'*grad_a_wab;

%             end
        end
        % Leap frog stuff: 
        % stepTerm = 0.5 for first step. Otherwise stepTerm==1. 
        % Acceleration
        Acc(1,i) = xAcc_i;
        Acc(2,i) = yAcc_i + g;
        % Velocity at t+0.5*dt
        VelHalf(1,i) = VelHalf(1,i) + stepTerm*dt*Acc(1,i);
        VelHalf(2,i) = VelHalf(2,i) + stepTerm*dt*Acc(2,i);
        % rho at t+0.5*dt
        Rho_RhoHalf_dRhoNext(2,i) = Rho_RhoHalf_dRho(2,i) + stepTerm*dt*drho_i;
        % dRho
        Rho_RhoHalf_dRhoNext(3,i) = drho_i;
        % Temporary Position
        PosNext(1,i) = Pos(1,i) + dt * VelHalf(1,i);
        PosNext(2,i) = Pos(2,i) + dt * VelHalf(2,i);
    end
    
    Pos = PosNext;
    Rho_RhoHalf_dRho = Rho_RhoHalf_dRhoNext;
end

