clear all;
addpath('../src/EX','../src/IM');
format long;
global dt dr A V dv CFL GAMMA C_v BMl BMr kappa0 a b theta NL_TOL TAU maxite;

%computational time
t_end = 0.01;
%mesh
Nr = 2e2;
R = 1;
dr = R/Nr;
r = 0 : dr : R;
A = 4*pi*r.^2;
V = 4/3*pi*r.^3;
dv = V(2:end) - V(1:end-1);
rc = 0.5*(r(1:end-1) + r(2:end));%r at cell center
%CFL number
CFL = 5e-1;
%ratio of specific heats
GAMMA = 5/4;
%gas constant
GR = 8.3144598;
%specific heat
C_v = 1/(GAMMA-1);
%boundary mark
BMl = -1; BMr = 1;
%coefficient of heat conduction
kappa0 = 1;
%exponent of heat conduction term
a = 0; b = 2.5;
%coefficient of C-N scheme: theta*t^{n+1} + (1-theta)*t^n
theta = 0.5;
%tol for Newton iteration
NL_TOL = 1e-3;
%linear tol for GMRES
TAU = 1e-3;
%number of Newton iteration
maxite = 2;

%initialization
[U T] = initialization(r, rc);

tic;
dt = 0;
t = 0;
deltaT = zeros(size(T));
while t < t_end;
%for j = 1:1
  dt = cal_dt(U, T, dt, t, deltaT);
  if t + dt > t_end
    dt = t_end - t;
  end
  %EX
  U = EX(U, T, dt);

  %IM
  T1 = IM(U, T, dt, maxite);
  deltaT = T1 - T;
  T = T1;
  %relation between hydrodynamics and heat conduction
  U(3, :) = 0.5*U(2,:).^2./U(1,:) + C_v.*U(1,:).*T;

  t = t + dt;
  fprintf('t:%.4e, dt:%.4e\n', t, dt);
end
toc;

save Smooth_the0 rc U T

