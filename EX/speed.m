function alpha = speed(reconl, reconr, BDl, BDr);
global GR GAMMA;
%give speed at the cell boundaries

ul = [BDl, reconl];
ur = [reconr, BDr];
vl = ul(2,:)./ul(1,:);
vr = ur(2,:)./ur(1,:);
pl = (ul(3,:) - 0.5.*ul(1,:).*vl.^2).*(GAMMA-1);
pr = (ur(3,:) - 0.5.*ur(1,:).*vr.^2).*(GAMMA-1);
% ????? here we do not use relation c = RT as in the paper
cl = sqrt(GAMMA.*abs(pl)./ul(1,:));
cr = sqrt(GAMMA.*abs(pr)./ur(1,:));
%cl = sqrt(GR*T);
%cr = sqrt(GR*T);

alpha = zeros(size(vl));
tmpl = zeros(size(vl));
tmpr = zeros(size(vr));
for i = 1:length(alpha)
  if vl(i) < 0
    tmpl(i) = cl(i) - vl(i);
  else
    tmpl(i) = vl(i) + cl(i);
  end
  if vr(i) < 0
    tmpr(i) = cr(i) - vr(i);
  else
    tmpr(i) = vr(i) + cr(i);
  end
end
alpha = max(tmpl, tmpr);


