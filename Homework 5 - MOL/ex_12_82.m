m = 50;
h = 1.0/m;
mu = 2.4;
k = mu * h;
A = (-1.5*mu + 0.5*mu*mu) * eye(m);
A(2:m,1:m-1) = A(2:m,1:m-1) + (2*mu - mu*mu) * eye(m-1);
A(3:m,1:m-2) = A(3:m,1:m-2) + (-0.5*mu + 0.5*mu*mu) * eye(m-2);
A(1,m) = 2*mu - mu*mu;
A(1,m-1) = -0.5*mu + 0.5*mu*mu;
A(2,m) = -0.5*mu + 0.5*mu*mu;
lambda = eig(A);
scatter(real(lambda), imag(lambda), 'filled');
hold on
plot(cos(0:2*pi/1000:2*pi)-1, sin(0:2*pi/1000:2*pi));
xlim([-2.5,0.5]);
ylim([-1.5,1.5]);
set(gcf,'Units','centimeters','Position',[6 6 12 12]);
saveas(gcf, 'figures/ex_12_82_mu=2.4.eps');
