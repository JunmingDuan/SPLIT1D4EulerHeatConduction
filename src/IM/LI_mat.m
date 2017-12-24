function T1 = LI_mat(U, T, dt)
% return maxtrix formed by the linearized method
global C_v A dv dr kappa0
n = length(T);
row = [];
col = [];
val = [];
RHS = zeros(n,1);
rho = U(1, :);
u = U(2, :)./U(1, :);
K = kappa0 * (rho(:,[1,1:end])+rho(:,[1:end,end]))/2 .* (T(:,[1,1:end])+T(:,[1:end,end]))/2;

for i = 1:n
  row = [row, i];
  col = [col, i];
  tmp = C_v*rho(i) + dt/dv(i)/dr*(A(i+1)*K(i+1)+A(i)*K(i));
  val = [val, tmp];
  RHS(i) = C_v*rho(i)*T(i) - rho(i)*u(i)^2/2;

  if i ~= 1
    row = [row, i];
    col = [col, i-1];
    tmp = - dt/dv(i)/dr*(A(i)*K(i));
    val = [val, tmp];
  else
    row = [row, i];
    col = [col, i];
    tmp = - dt/dv(i)/dr*(A(i)*K(i));
    val = [val, tmp];
  end

  if i ~= n
    row = [row, i];
    col = [col, i+1];
    tmp = - dt/dv(i)/dr*(A(i+1)*K(i+1));
    val = [val, tmp];
  else
    row = [row, i];
    col = [col, i];
    tmp = - dt/dv(i)/dr*(A(i+1)*K(i+1));
    val = [val, tmp];
  end
end
RHS(1) = C_v*rho(1)*T(1);
RHS(n) = C_v*rho(n)*T(n);
MAT = sparse(row, col, val, n, n, n*3);
RHS = F(U, T, dt, T, 0);
RHS = RHS';
T1 = MAT\RHS;
T1 = T1';

end

