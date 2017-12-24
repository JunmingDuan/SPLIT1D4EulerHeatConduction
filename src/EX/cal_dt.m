function dt = cal_dt(u);
global gamma CFL dr;
%give speed at the cell boundaries

v = u(2,:)./u(1,:);
p = (u(3,:) - 0.5.*u(1,:).*v.^2).*(gamma-1);
% ????? here we do not use relation c = RT as in the paper
c = sqrt(gamma.*p./u(1,:));

alpha = zeros(size(v));
for i = 1:length(alpha)
  if v(i) < 0
    alpha(i) = c(i) - v(i);
  else
    alpha(i) = v(i) + c(i);
  end
end

dt = CFL*dr/max(abs(alpha));

