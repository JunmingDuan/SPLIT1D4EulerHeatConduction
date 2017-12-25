function T1 = JFNK(U, T, dt, maxite, theta)
%Jacobian free Newton Krylov solver
%input:  K:  Newton iteration times
global kappa0 C_v dv dr A TAU NL_TOL gmres_ite

%[M1, M2] = ilu(MAT, struct('type', 'ilutp', 'droptol', 1e-5));
%[M1, M2] = ilu(MAT);
%MAT = LI_mat(U, T, dt);
%Newton iteration
res0 = norm(F(U, T, dt, T, 0), 2);
for k = 1:maxite
  rhs = -F(U, T, dt, T, 0)';
  MAT = LI_mat(U, T, dt);
  [dT, flag, relres] = gmres(@afun, rhs, 10, TAU, 20, MAT, [], T', U, T, dt, theta);
  %[dT, flag, relres] = gmres(@afun, rhs, 10, TAU, 20, M1, M2, T', U, T, dt, theta);
  %[dT, flag, relres] = gmres(@afun, rhs, 10, TAU, 20, [], [], T', U, T, dt, theta);
  T = T + dT';
  relres;
  nl_res = norm(F(U, T, dt, T, 0), 2);
  if nl_res < NL_TOL*res0
    break;
  end
end
T1 = T;

