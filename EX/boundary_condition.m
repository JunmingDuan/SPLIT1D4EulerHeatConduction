function [BDl, BDr] = boundary_condition(ur, ul, BMl, BMr, Dl, Dr)
%give the values of the ghost cells at the left and right boundaries
%input:  ur,ul: reconstruction values at the cell boundaries
%        BMl, BMr:  boundary mark
%        0:   Diriclet,  needs two extra parameters
%        1:   outflow or inflow
%        -1:  reflection

if BMl == 0
  BDl = Dl;
elseif BMl == 1
  BDl = ur(:,1);
elseif BMl == -1
  BDl = ur(:,1);
  BDl(2,:) = -BDl(2,:);
end

if BMr == 0
  BDr = Dr;
elseif BMr == 1
  BDr = ul(:,end);
elseif BMr == -1
  BDr = ul(:,end);
  BDr(2,:) = -BDr(2,:);
end

