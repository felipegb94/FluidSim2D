function [ sum ] = checkKernel( Pos, Rho_RhoHalf_dRho, numParticles, ParticleNumber, ParticleMass, h )
%checkKernel Summary of this function goes here
%   Detailed explanation goes here


if(ParticleNumber > numParticles)
    assert(-1);
end

x_i = Pos(1,ParticleNumber);
y_i = Pos(2,ParticleNumber);

sum = 0;
for j = 1:numParticles
    x_j = Pos(1,j);
    y_j = Pos(2,j);  
    dx = x_i - x_j;
    dy = y_i - y_j;
    r_ij = [dx; dy];
    norm_r_ij = norm(r_ij);
    q = norm_r_ij/h;
    Wq = W(q, h);
    rhoj = Rho_RhoHalf_dRho(1,j);
    mj_Rhoj = ParticleMass/rhoj;
    sum = sum + mj_Rhoj*Wq;
end

end

