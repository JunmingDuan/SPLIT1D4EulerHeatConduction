function u = initialization(r)
global gamma

%Sod
rhol = 1;
ul = 0;
pl = 1;
rhor = 0.125;
ur = 0;
pr = 0.1;

for i = 1:length(r)
  if r(i) < 10
    u(1,i) = rhol;
    u(2,i) = rhol*ul;
    u(3,i) = 0.5*rhol*ul^2 + pl/(gamma-1);
  else
    u(1,i) = rhor;
    u(2,i) = rhor*ur;
    u(3,i) = 0.5*rhor*ur^2 + pr/(gamma-1);
  end
end

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
    %u(3,i) = 0.5*rhol*ul^2 + pl/(gamma-1);
  %else
    %u(1,i) = rhor;
    %u(2,i) = rhor*ur;
    %u(3,i) = 0.5*rhor*ur^2 + pr/(gamma-1);
  %end
%end

