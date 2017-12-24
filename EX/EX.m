function U1 = EX(U, T, dt)
global dr A dv BMl BMr;
%explicitly solve the hydrodynamics part for dt
%input:   U = [rho;m;E]: row vectors denote densiti,momentum,energy
%         alpha: characteristic speed at cell boundaries

%RK2
%stage 1
[BDl, BDr] = boundary_condition(U, U, BMl, BMr, 0, 0);%reflective and outflow BD
[reconl, reconr] = recon(U, BDl, BDr);
[BDl, BDr] = boundary_condition(reconr, reconl, BMl, BMr, 0, 0);%reflective and outflow BD
alpha = speed(reconl, reconr, BDl, BDr);
[flux, G] = LLF([BDl,reconl], [reconr,BDr], alpha);
u1 = U - dt*(A(:,2:end).*flux(:,2:end) - A(:,1:end-1).*flux(:,1:end-1))./dv - dt/dr.*(G(:,2:end) - G(:,1:end-1));

%stage 2
[BDl, BDr] = boundary_condition(u1, u1, BMl, BMr, 0, 0);%reflective and outflow BD
[reconl, reconr] = recon(u1, BDl, BDr);
[BDl, BDr] = boundary_condition(reconr, reconl, BMl, BMr, 0, 0);%reflective and outflow BD
alpha = speed(reconl, reconr, BDl, BDr);
[flux, G] = LLF([BDl,reconl], [reconr,BDr], alpha);
u2 = u1 - dt*(A(:,2:end).*flux(:,2:end) - A(:,1:end-1).*flux(:,1:end-1))./dv - dt/dr.*(G(:,2:end) - G(:,1:end-1));

U1 = 0.5*U + 0.5*u2;

%RK1
%[BDl, BDr] = boundary_condition(U, U, BMl, BMr, 0, 0);%reflective and outflow BD
%alpha = speed(U, U, BDl, BDr);
%[flux, G] = LLF([BDl,U], [U,BDr], alpha);
%U1 = U - dt*(A(:,2:end).*flux(:,2:end) - A(:,1:end-1).*flux(:,1:end-1))./dv - dt/dr.*(G(:,2:end) - G(:,1:end-1));

