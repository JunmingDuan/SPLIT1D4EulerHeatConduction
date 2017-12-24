function [U T] = initialization(rc)
global dv GAMMA C_v

%Barenblatt
%rho = ones(size(rc));
%u = zeros(size(rc));
%E0 = 1e1;
%E = 1e-4*ones(size(rc));
%E(1) = E0/dv(1);
%U = [rho; rho.*u; E];
%T = E/C_v./rho;

%Sod
rhol = 1;
ul = 0;
pl = 1;
rhor = 0.125;
ur = 0;
pr = 0.1;
for i = 1:length(rc)
  if rc(i) < 10
    U(1,i) = rhol;
    U(2,i) = rhol*ul;
    U(3,i) = 0.5*rhol*ul^2 + pl/(GAMMA-1);
  else
    U(1,i) = rhor;
    U(2,i) = rhor*ur;
    U(3,i) = 0.5*rhor*ur^2 + pr/(GAMMA-1);
  end
end
T = U(3,:)./U(1,:)/C_v;

%Nod
%rhol = 1;
%ul = -1;
%pl = 1e-5;
%rhor = 1;
%ur = -1;
%pr = 1e-5;

%for i = 1:length(r)
  %if r(i) < 1
    %u(1,i) = rhol;
    %u(2,i) = rhol*ul;
    %u(3,i) = 0.5*rhol*ul^2 + pl/(GAMMA-1);
  %else
    %u(1,i) = rhor;
    %u(2,i) = rhor*ur;
    %u(3,i) = 0.5*rhor*ur^2 + pr/(GAMMA-1);
  %end
%end

