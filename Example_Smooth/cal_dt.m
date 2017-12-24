function dt = cal_dt(U, T, dt, t, deltaT);
global dr CFL GAMMA C_v kappa0 b;
%give speed at the cell boundaries

if t < 1e-18
  %dt = dr^2/U(1,1)/C_v/(kappa0*T(1)^b);
  dt1 = 1e-10;
else
  vf = (norm(deltaT,1)/dt)./(norm(T(2:end)-T(1:end-1),1)/dr);
  dt1 = CFL*dr/vf;
end

v = U(2,:)./U(1,:);
p = (U(3,:) - 0.5.*U(1,:).*v.^2).*(GAMMA-1);
% ????? here we do not use relation c = RT as in the paper
c = sqrt(GAMMA.*p./U(1,:));

alpha = zeros(size(v));
for i = 1:length(alpha)
  if v(i) < 0
    alpha(i) = c(i) - v(i);
  else
    alpha(i) = v(i) + c(i);
  end
end

dt2 = 0.5*dr/max(abs(alpha));

dt = min(dt1, dt2);

