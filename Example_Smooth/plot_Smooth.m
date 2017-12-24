data = load('Smooth_the0.mat');
rc = data.rc;
U = data.U;
T = data.T;

subplot(2,2,1);
plot(rc, U(1,:));
xlabel('r');
ylabel('Density');
set(gca, 'Xtick', 0:0.2:1);
set(gca, 'Ytick', 0:2:18);
subplot(2,2,2);
plot(rc, R*U(1,:).*T);
xlabel('r');
ylabel('Pressure');
set(gca, 'Xtick', 0:0.2:1);
set(gca, 'Ytick', 0:20:120);
subplot(2,2,3);
plot(rc, U(2,:)./U(1,:));
xlabel('r');
ylabel('Velocity');
set(gca, 'Xtick', 0:0.2:1);
set(gca, 'Ytick', 0:0.5:2);
subplot(2,2,4);
plot(rc, T);
xlabel('r');
ylabel('Temperature');
set(gca, 'Xtick', 0:0.2:1);
set(gca, 'Ytick', 0:2:8);

hold on;
save Smooth_the0 rc U T

