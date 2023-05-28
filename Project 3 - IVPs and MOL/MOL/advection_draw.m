A = load('result.txt');
f = @(x) exp(-20*((x-2).^2))+exp(-((x-5).^2));
t = 15:0.05:25;
n = size(A,1);
plot(t, A(n,342:542), '.-');
hold on
plot(t, f(t-17), '-.');
ylim([-0.4,1.1]);
saveas(gcf, 'report/figures/advection-2-2.eps');