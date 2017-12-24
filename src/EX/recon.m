function [vl, vr] = recon(u, BDl, BDr)
%linear reconstruction
%input:  u,  matrix
%output: vl, vr,   left and right limits at cell boundaries

slope = [u(:,1:end),BDr] - [BDl,u(:,1:end)];
slope = minmod(slope(:, 1:end-1), slope(:, 2:end));
vl = u + 0.5*slope;
vr = u - 0.5*slope;

