mu = 1.6;
m = 50;
h = 1.0/m;
p = 1:1:m;
z = exp(-2*pi*1i*h*p).*((mu*mu-2*mu)*(cos(2*pi*h*p)-1) - 1i*mu*sin(2*pi*h*p));
scatter(real(z), imag(z), "filled");
hold on
plot(cos(0:2*pi/1000:2*pi)-1, sin(0:2*pi/1000:2*pi));
xlim([-2.5,0.5]);
ylim([-1.5,1.5]);
set(gcf,'Units','centimeters','Position',[6 6 12 12]);
saveas(gcf, 'figures/ex_12_82_mu=1.6.eps');