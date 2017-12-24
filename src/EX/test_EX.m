clear all;
format long;
global GAMMA CFL dt dr A V dv BMl BMr;
GAMMA = 5/3;
CFL = 0.8;
BMl = -1;
BMr = -1;
t_end = 1.4;
Nr = 4e2;
R = 1;

dr = R/Nr;
r = 0 : dr : R;
A = 4*pi*r.^2;
A = [A;A;A];
V = 4/3*pi*r.^3;
dv = V(2:end) - V(1:end-1);
dv = [dv;dv;dv];
rc = 0.5*(r(1:end-1) + r(2:end));%r at cell center

u = initialization(rc);
t = 0;

while t < t_end
  dt = cal_dt(u);
  if t + dt > t_end
    dt = t_end - t;
  end
  u = EX(r, u, dt, 0);
  t = t + dt;
  fprintf('t:%.5f, dt:%.5f\n', t, dt);
end

save test rc u
plot(rc, u(1,:));

