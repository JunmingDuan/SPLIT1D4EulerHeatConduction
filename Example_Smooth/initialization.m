function [U T] = initialization(r, rc)
global dv GAMMA C_v

%Smooth problem
rho = 1./rc;
u = zeros(size(rc));
E0 = 1e2;
c0 = 1/4;
f = @(x) E0*exp(-x.^2/c0^2)/(c0*sqrt(pi))^3;
E = (E0*(erf(r(2:end)/c0)-erf(r(1:end-1)/c0)) - 2*pi*c0^2*(r(2:end).*f(r(2:end)) - ...
r(1:end-1).*f(r(1:end-1))))./dv;
U = [rho; rho.*u; E];
T = E/C_v./rho;

%Sod
%rhol = 1;
%ul = 0;
%pl = 1;
%rhor = 0.125;
%ur = 0;
%pr = 0.1;
%for i = 1:length(rc)
  %if rc(i) < 10
    %U(1,i) = rhol;
    %U(2,i) = rhol*ul;
    %U(3,i) = 0.5*rhol*ul^2 + pl/(GAMMA-1);
  %else
    %U(1,i) = rhor;
    %U(2,i) = rhor*ur;
    %U(3,i) = 0.5*rhor*ur^2 + pr/(GAMMA-1);
  %end
%end
%T = U(3,:)./U(1,:)/C_v;

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

