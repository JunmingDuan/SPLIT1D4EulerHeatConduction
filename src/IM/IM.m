function T1 = IM(U, T, dt, maxite, theta)
%linear implicit solver
%input:  K:  Newton iteration times
global kappa0 C_v dv dr A TAU NL_TOL gmres_ite

%LI
%rhs = -F(U, T, dt, T, 0)';
%MAT = LI_mat(U, T, dt);
%T1 = gmres(MAT, rhs, 10, 1e-9, 20, [], [], T');
T1 = LI_mat(U, T, dt);

%Newton iteration
%{for k = 1:maxite%}
  %rhs = -F(U, T, dt, T, 0)';
  %%[dT, flag, relres] = gmres(@afun, rhs, 10, TAU*norm(rhs,2), 20, [], [], T', T, u, E, rho, dt, theta);
  %%[L, U] = ilu()
  %%[dT, flag, relres] = gmres(@afun, rhs, 10, TAU, 20, [], [], T', U, T, dt, theta);
  %T = T + dT';
  %relres;
  %nl_res = norm(F(U, T, dt, T, 0), 2);
  %if nl_res < NL_TOL
    %break;
  %end
%end
%T1 = T;


