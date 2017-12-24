function k = kappa(rho, T)
global kappa0 a b
k = kappa0*(([T(1),T]+[T,T(end)])/2).^b;

