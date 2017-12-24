function y = afun(x, U, T, dt, theta)
%return J*x for gmres in the Newton iteration process
%input:  x:  column vector
%        T:  column vector

u = U(2,:)./U(1,:);
x = x';
b = 1e-6;
e = b + b/length(T)*sum(abs(u))/norm(x,2);
y = (F(U, T+e*x, dt, T, theta) - F(U, T, dt, T, theta))/e;
y = y';
%y = ones(size(x));

