clear all;
addpath('../src/EX','../src/IM');
format long;
global dt dr A V dv CFL GAMMA C_v BMl BMr kappa0 a b theta NL_TOL TAU maxite gmres_ite;

%computational time
%t_end = 0.3;
t_end = 1.0;
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
CFL = 0.5;
%ratio of specific heats
GAMMA = 5/4;
%specific heat
C_v = 1.0/(GAMMA-1);
%boundary mark
BMl = -1; BMr = 1;
%coefficient of heat conduction
kappa0 = 1;
%exponent of heat conduction term
a = 0; b = 6.5;
%coefficient of C-N scheme: theta*t^{n+1} + (1-theta)*t^n
theta = 1;
%tol and iterations for Newton iteration
NL_TOL = 1e-7; maxite = 3;
%linear tol and iterations for GMRES
TAU = 1e-3; gmres_ite = 100;

%initialization
[U T] = initialization(rc);

tic;
dt = 0;
t = 0;
deltaT = zeros(size(T));
flag = 0;
while t < t_end;
  dt = cal_dt(flag, U, T, dt, t, deltaT);
  flag = 1;
  if t + dt > t_end
    dt = t_end - t;
  end
  %EX
  %U = EX(U, T, dt);
  %T = U(3,:)./U(1,:)/C_v;

  %IM
  theta = 1;
  T1 = IM(U, T, dt/2, maxite, theta);
  deltaT = T1 - T;
  T = T1;
  %relation between hydrodynamics and heat conduction
  U(3, :) = C_v.*U(1,:).*T;
  %IM
  theta = 0.5;
  T1 = IM(U, T, dt/2, maxite, theta);
  deltaT = T1 - T;
  T = T1;
  %relation between hydrodynamics and heat conduction
  U(3, :) = C_v.*U(1,:).*T;

  t = t + dt;
  fprintf('t:%.4e, dt:%.4e\n', t, dt);
end
toc;

save(['Barenblatt_the',num2str(theta),'_t',num2str(t_end),'_k',num2str(b),'.mat'], 'rc', 'U', 'T');

%Barenblatt exact
E0 = 10;
Q = E0/U(1,1)/C_v;
KAI0 = kappa0/U(1,1)/C_v;
XI0 = ((3*b+2)/(2^(b-1)*b*pi^b))^(1/(3*b+2))*(gamma(5/2+1/b)/gamma(1+1/b)/gamma(3/2))^(b/(3*b+2));
rf = XI0*(KAI0*Q^b*t_end)^(1/(3*b+2));
T_c = Q*XI0^3/rf^3*(b*XI0^2/2/(3*b+2))^(1/b);
T_exact = T_c*(1-min(rc,rf).^2/rf^2).^(1/b);
save(['Barenblatt_exact_E',num2str(E0),'_t',num2str(t_end),'_k',num2str(b),'.mat'], 'rc', 'T_exact');

