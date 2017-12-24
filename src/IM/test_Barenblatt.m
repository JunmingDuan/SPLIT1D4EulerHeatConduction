clear all;
format long;
global dt dr A V dv C_v kappa0 GAMMA b SCHEME theta tau;
t_end = 0.3; Nr = 2e2;
R = 1;
maxite = 5;
GAMMA = 5/4;
C_v = 1.0/(GAMMA-1);
kappa0 = 1;
CFL = 5e-1;
%theta*t^{n+1} + (1-theta)*t^n
theta = 0;
%tol for Newton iteration
NL_TOL = 1e-7;
%linear tol for GMRES
TAU = 1e-2;

dr = R/Nr;
r = 0 : dr : R;
A = 4*pi*r.^2;
V = 4/3*pi*r.^3;
dv = V(2:end) - V(1:end-1);
rc = 0.5*(r(1:end-1) + r(2:end));%r at cell center

%hydrodynamics variables
E0 = 1e1;
E = 1e-4*ones(size(rc));
rho = ones(size(rc));
E(1) = E0/dv(1);
u = zeros(size(rc));
T = initialization_T(rc, E, rho);

tic;
t = 0;
while t < t_end;
%for j = 1:1
  if t < 1e-18
    dt = dr^2/rho(1)/C_v/(kappa0*T(1)^2.5);
  else
    vf = (norm(deltaT,1)/dt)./(norm(T(2:end)-T(1:end-1),1)/dr);
    dt = CFL*dr/vf;
  end
  if t + dt > t_end
    dt = t_end - t;
  end
  T1 = IM(T, E, rho, u, dt, maxite);
  deltaT = T1 - T;
  T = T1;
  %relation between hydrodynamics and heat conduction
  E = C_v.*rho.*T;
  t = t + dt;
  fprintf('t:%.8e, dt:%.8e\n', t, dt);
end
toc;

hold on;
plot(rc, T, 'bo');
save test rc u

%exact
Q = E0/rho(1)/C_v;
KAI0 = kappa0/rho(1)/C_v;
XI0 = ((3*b+2)/(2^(b-1)*b*pi^b))^(1/(3*b+2))*(gamma(5/2+1/b)/gamma(1+1/b)/gamma(3/2))^(b/(3*b+2));
rf = XI0*(KAI0*Q^b*t_end)^(1/(3*b+2));
T_c = Q*XI0^3/rf^3*(b*XI0^2/2/(3*b+2))^(1/b);
T_exact = T_c*(1-min(rc,rf).^2/rf^2).^(1/b);
plot(rc, T_exact, '-k');
legend('numerical','analytical');
axis([0,1,0,1.4]);
xlabel('r');
ylabel('Temperature');
fprintf('err1: %.6e\n', norm(T-T_exact, 1)/length(T));
fprintf('err2: %.6e\n', norm(T-T_exact, 2)/sqrt(length(T)));
fprintf('errinf %.6e\n', norm(T-T_exact, inf));

