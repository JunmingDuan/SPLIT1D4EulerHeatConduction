function T1 = LI(U, T, dt, maxite, theta)
%linear implicit solver
%input:  K:  Newton iteration times
global kappa0 C_v dv dr A TAU NL_TOL

%kappa = kappa0 * (rho(:,[1,1:end])+rho(:,[1:end,end]))/2 .* (T(:,[1,1:end])+T(:,[1:end,end]))/2;

%linear equation
%n = length(T);
%row = [];
%col = [];
%val = [];
%RHS = zeros(n,1);
%for i = 1:n
  
  %row = [row, i];
  %col = [col, i];
  %tmp = C_v*rho(i) + dt/dv(i)/dr*(A(i+1)*K(i+1)+A(i)*K(i));
  %val = [val, tmp];
  %RHS(i) = C_v*rho(i)*T(i) - rho(i)*u(i)^2/2;
  
  %if i ~= 1
    %row = [row, i];
    %col = [col, i-1];
    %tmp = - dt/dv(i)/dr*(A(i)*K(i));
    %val = [val, tmp];
  %else
    %row = [row, i];
    %col = [col, i];
    %tmp = - dt/dv(i)/dr*(A(i)*K(i));
    %val = [val, tmp];
  %end
  
  %if i ~= n
    %row = [row, i];
    %col = [col, i+1];
    %tmp = - dt/dv(i)/dr*(A(i+1)*K(i+1));
    %val = [val, tmp];
  %else
    %row = [row, i];
    %col = [col, i];
    %tmp = - dt/dv(i)/dr*(A(i+1)*K(i+1));
    %val = [val, tmp];
  %end
%end
%RHS(1) = C_v*rho(1)*T(1);
%RHS(n) = C_v*rho(n)*T(n);
%MAT = sparse(row, col, val, n, n, n*3);
%T1 = MAT\RHS;
%T1 = T1';
%full(MAT)
%RHS

%Newton iteration
for k = 1:maxite
  rhs = -F(U, T, dt, T, 0)';
  %[dT, flag, relres] = gmres(@afun, rhs, 10, TAU*norm(rhs,2), 20, [], [], T', T, u, E, rho, dt, theta);
  [dT, flag, relres] = gmres(@afun, rhs, 10, TAU, 20, [], [], T', U, T, dt, theta);
  %dT = gmres(MAT, RHS, 10, 1e-9, 20, [], [], T');
  T = T + dT';
  relres;
  nl_res = norm(F(U, T, dt, T, 0), 2);
  if nl_res < NL_TOL
    break;
  end
end
T1 = T;

