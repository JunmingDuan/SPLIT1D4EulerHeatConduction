function plot_Barenblatt(E0, t_end, b, theta, CFL)
%data = load(['Barenblatt_the',num2str(theta),'_t',num2str(t_end),'_k',num2str(b),'_CFL',num2str(CFL),'.mat'], 'rc', 'U', 'T');
%data = load(['Barenblatt_PC_the',num2str(theta),'_t',num2str(t_end),'_k',num2str(b),'_CFL',num2str(CFL),'.mat'], 'rc', 'U', 'T');
data = load(['Barenblatt_LI_the',num2str(theta),'_t',num2str(t_end),'_k',num2str(b),'_CFL',num2str(CFL),'.mat'], 'rc', 'U', 'T');
rc = data.rc;
T = data.T;
%Barenblatt exact
GAMMA = 5/4;
C_v = 1.0/(GAMMA-1);
E0 = 10;
Q = E0/1.0/C_v;
KAI0 = 1/C_v;
XI0 = ((3*b+2)/(2^(b-1)*b*pi^b))^(1/(3*b+2))*(gamma(5/2+1/b)/gamma(1+1/b)/gamma(3/2))^(b/(3*b+2));
rf = XI0*(KAI0*Q^b*t_end)^(1/(3*b+2));
T_c = Q*XI0^3/rf^3*(b*XI0^2/2/(3*b+2))^(1/b);
T_exact = T_c*(1-min(rc,rf).^2/rf^2).^(1/b);
save(['Barenblatt_exact_E',num2str(E0),'_t',num2str(t_end),'_k',num2str(b),'.mat'], 'rc', 'T_exact');

m = 1;
plot(rc(1:m:end), T(1:m:end), '-', rc, T_exact, 'k-');
legend('numerical','analytical');
axis([0,1,0,1.4]);
%axis([0,1,0,1]);
xlabel('r');
ylabel('Temperature');

fprintf('err1: %.6e\n', norm(T-T_exact, 1)/length(T));
fprintf('err2: %.6e\n', norm(T-T_exact, 2)/sqrt(length(T)));
fprintf('errf: %.6e\n', norm(T-T_exact, inf));

end

