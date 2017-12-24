function MAT = LI_mat(U, T, dt)
% return maxtrix formed by the linearized method
global C_v A dv dr kappa0
n = length(T);
row = [];
col = [];
val = [];
RHS = zeros(n,1);
rho = U(1, :);
u = U(2, :)./U(1, :);
K = kappa(rho, T);

for i = 1:n
  row = [row, i];
  col = [col, i];
  tmp = C_v*rho(i) + dt/dv(i)/dr*(A(i+1)*K(i+1)+A(i)*K(i));
  val = [val, tmp];

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
MAT = sparse(row, col, val, n, n, n*3);

end

