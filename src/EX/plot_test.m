data = load('test.mat');
rc = data.rc;
u = data.u;
data2 = load('test2.mat');
rc2 = data2.rc;
u2 = data2.u;
plot(rc, u(1,:), 'ro', rc2, u2(1,:), '-');
axis([5,16,0.1,1]);
set(gca, 'Xtick', 5:1:16);

