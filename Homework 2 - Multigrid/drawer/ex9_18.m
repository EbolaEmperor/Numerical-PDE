w = [1/3, 1/2, 2/3, 1];
lstl = ['k', ':k', '--k', '-.k'];
ind = [1,1,2,3,4,6,7,9];
n = 64;
rho = zeros(1,4);
for wk = 1:4
    cw = w(wk);
    k = 0:0.1:n;
    lam = 1-2*cw*sin(k*pi/(2*n)).^2;
    rho(wk) = lam(11);
    plot(k,lam,lstl(ind(2*wk-1):ind(2*wk)))
    hold on
end
disp(rho);
xlabel('$k$','interpreter','latex');
ylabel('$\lambda_k(T_\omega)$','interpreter','latex');
legend('$\omega=1/3$','$\omega=1/2$','$\omega=2/3$','$\omega=1$','interpreter','latex','FontSize',12);