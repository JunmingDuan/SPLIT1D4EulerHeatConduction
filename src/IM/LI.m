function T1 = LI(U, T, dt)
% linearized method solve the diffusion term
global C_v A dv dr kappa0
n = length(T);
rho = U(1, :);
u = U(2, :)./U(1, :);

RHS = C_v*rho.*T - rho.*u.^2/2;
RHS(1) = C_v*rho(1)*T(1);
RHS(n) = C_v*rho(n)*T(n);
MAT = LI_mat(U, T, dt);
T1 = MAT\RHS';
T1 = T1';

end

