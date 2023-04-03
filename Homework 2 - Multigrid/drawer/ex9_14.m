n = 6;
h = 1.0/n;
X = 0:0.01:1;
k1 = 1.5*n;
k2 = 0.5*n;
Y1 = sin(k1*pi*X);
Y2 = sin(k2*pi*X);
Y3 = -Y2;
plot(X, Y1, 'k')
hold on
plot(X, Y2, ':k')
plot(X, Y3, '-.k')
SX = 0:h:1;
SY1 = sin(k1*pi*SX);
scatter(SX,SY1,'k','filled')
plot(SX,SY1,'k','LineWidth',2)
legend("$sin(\frac{3}{2}n\pi x_j)$","$sin(\frac{1}{2}n\pi x_j)$","$-sin(\frac{1}{2}n\pi x_j)$",'$x_j$','interpreter','latex','FontSize',12)
