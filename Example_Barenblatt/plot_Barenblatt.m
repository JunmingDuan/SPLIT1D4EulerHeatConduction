function plot_Barenblatt(E0, t_end, b, theta)
EX1 = load(['Barenblatt_exact_E',num2str(E0),'_t',num2str(t_end),'_k',num2str(b),'.mat']);
rc1 = EX1.rc;
T_exact = EX1.T_exact;
data = load(['Barenblatt_the',num2str(theta),'_t',num2str(t_end),'_k',num2str(b),'.mat']);
rc = data.rc;
T = data.T;

m = 2;
plot(rc(1:m:end), T(1:m:end), 'ro', rc1, T_exact, '-');
legend('numerical','analytical');
axis([0,1,0,1.4]);
%axis([0,1,0,1]);
xlabel('r');
ylabel('Temperature');

fprintf('err1: %.6e\n', norm(T-T_exact, 1)/length(T));
fprintf('err2: %.6e\n', norm(T-T_exact, 2)/sqrt(length(T)));
fprintf('errf: %.6e\n', norm(T-T_exact, inf));

end

