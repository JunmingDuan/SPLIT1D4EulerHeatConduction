function [flux, G] = LLF(ul, ur, alpha)
global GAMMA
%give LLF flux
%input: ul, ur: matrix, values at the left and right of the cell boundaries
%       alpha:  local characteristic velocity

fl = zeros(size(ul));
fl(1,:) = ul(2,:);
v = ul(2,:)./ul(1,:);
fl(2,:) = v.^2.*ul(1,:);
pl = (ul(3,:) - 0.5.*fl(2,:)).*(GAMMA-1);
fl(3,:) = (ul(3,:) + pl).*v;

fr = zeros(size(ur));
fr(1,:) = ur(2,:);
v = ur(2,:)./ur(1,:);
fr(2,:) = v.^2.*ur(1,:);
pr = (ur(3,:) - 0.5.*fr(2,:)).*(GAMMA-1);
fr(3,:) = (ur(3,:) + pr).*v;

flux = 0.5*((fl + fr) - [alpha;alpha;alpha].*(ur - ul));

G = zeros(size(flux));
G(2, :) = 0.5*(pl + pr);

