function T1 = IM(U, T, dt, maxite)
%implicitly solve the diffusion part of the hydrodynamics
%input:  T:  vector, temperature of the cell
%        rho: density at t^{n+1}
%        dt: time step
%

global theta
T1 = LI(U, T, dt, maxite, theta);

