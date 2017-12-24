function m = minmod(a, b)
%input:  a,b, matrix
%output: m,   minmod function of a,b's components
for i = 1:size(a,1)
  for j = 1:size(a,2)
    if a(i,j)*b(i,j) <= 0
      m(i,j) = 0;
    elseif abs(a(i,j)) < abs(b(i,j))
      m(i,j) = a(i,j);
    else
      m(i,j) = b(i,j);
    end
  end
end

