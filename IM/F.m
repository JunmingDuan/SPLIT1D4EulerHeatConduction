function g = F(U, T, dt, Tn, theta)
%return residual
%zero flux boundary, i.e., \pd{T}{r} = 0
global C_v A dv dr

rho = U(1,:);
u = U(2,:)./U(1,:);
E = U(3,:);
if theta == 1
  kn1 = kappa(rho, T);
  g = (C_v.*rho.*T+0.5*rho.*u.^2 - E) - (A(2:end).*kn1(2:end).*([T(2:end),T(end)]-T) ...
  - A(1:end-1).*kn1(1:end-1).*(T-[T(1),T(1:end-1)])).*dt./dv./dr;
  g(1) = (C_v.*rho(1).*T(1)+0.5*rho(1).*u(1).^2 - E(1)) - (A(2).*kn1(2).*(T(2)-T(1))).*dt./dv(1)./dr;
  g(end) = (C_v.*rho(end).*T(end)+0.5*rho(end).*u(end).^2 - E(end)) - (-A(end-1).*kn1(end-1).*(T(end)-T(end-1))).*dt./dv(end)./dr;
elseif theta == 0
  kn = kappa(rho, Tn);
  g = (C_v.*rho.*T+0.5*rho.*u.^2 - E) - (A(2:end).*kn(2:end).*([T(2:end),T(end)]-T) ...
  - A(1:end-1).*kn(1:end-1).*(T-[T(1),T(1:end-1)])).*dt./dv./dr;
  g(1) = (C_v.*rho(1).*T(1)+0.5*rho(1).*u(1).^2 - E(1)) - (A(2).*kn(2).*(T(2)-T(1))).*dt./dv(1)./dr;
  g(end) = (C_v.*rho(end).*T(end)+0.5*rho(end).*u(end).^2 - E(end)) - (-A(end-1).*kn(end-1).*(T(end)-T(end-1))).*dt./dv(end)./dr;
else
  kn1 = kappa(rho, T);
  gn1 = (C_v.*rho.*T+0.5*rho.*u.^2 - E) - (A(2:end).*kn1(2:end).*([T(2:end),T(end)]-T) ...
  - A(1:end-1).*kn1(1:end-1).*(T-[T(1),T(1:end-1)])).*dt./dv./dr;
  gn1(1) = (C_v.*rho(1).*T(1)+0.5*rho(1).*u(1).^2 - E(1)) - (A(2).*kn1(2).*(T(2)-T(1))).*dt./dv(1)./dr;
  gn1(end) = (C_v.*rho(end).*T(end)+0.5*rho(end).*u(end).^2 - E(end)) - (-A(end-1).*kn1(end-1).*(T(end)-T(end-1))).*dt./dv(end)./dr;
  kn = kappa(rho, Tn);
  gn = (C_v.*rho.*T+0.5*rho.*u.^2 - E) - (A(2:end).*kn(2:end).*([T(2:end),T(end)]-T) ...
  - A(1:end-1).*kn(1:end-1).*(T-[T(1),T(1:end-1)])).*dt./dv./dr;
  gn(1) = (C_v.*rho(1).*T(1)+0.5*rho(1).*u(1).^2 - E(1)) - (A(2).*kn(2).*(T(2)-T(1))).*dt./dv(1)./dr;
  gn(end) = (C_v.*rho(end).*T(end)+0.5*rho(end).*u(end).^2 - E(end)) - (-A(end-1).*kn(end-1).*(T(end)-T(end-1))).*dt./dv(end)./dr;
  g = theta*gn1 + (1-theta)*gn;
end

