clear all;
EX1 = load('Barenblatt_exact_E10_t1_k6_5.mat');
%EX1 = load('Barenblatt_exact_E10_t0_3_k2_5.mat');
rc1 = EX1.rc;
T1 = EX1.T_exact;
%data = load('Barenblatt1_the0.mat');
data = load('Barenblatt2_the0.mat');
rc = data.rc;
T = data.T;

m = 2;
plot(rc(1:m:end), T(1:m:end), 'ro', rc1, T1, '-');
legend('numerical','analytical');
%axis([0,1,0,1.4]);
axis([0,1,0,1]);
xlabel('r');
ylabel('Temperature');

%fprintf('err1: %.6e\n', norm(T-T_exact, 1)/length(T));
%fprintf('err2: %.6e\n', norm(T-T_exact, 2)/sqrt(length(T)));
%fprintf('errinf %.6e\n', norm(T-T_exact, inf));
%save Barenblatt_exact_E10_t1_k6_5 rc T_exact

