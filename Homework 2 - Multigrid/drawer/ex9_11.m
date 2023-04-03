n = 8;
h = 1/n;
lx = 0:0.01:1;
x = 0:h:1;
u = @(x)(sin((n-1)*pi*x));
ly = u(lx);
y = u(x);
scatter(x,y,'k','filled')
hold on
plot(x,y,'k')
plot(lx,ly,':k')