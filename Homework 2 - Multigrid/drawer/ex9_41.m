n = 128;
w = 2/3;
v1 = 0;
v2 = 0;
lx1 = 0:0.1:n/2;
lx2 = n/2:0.1:n;
y1 = lx1;
y4 = lx2;
y3 = lx1;
y2 = lx2;
rho = lx1;
for lk = 0:5*n
    ind1 = lk+1;
    ind2 = size(lx2,2)-lk;
    k = lk/10;
    lamk = lambda(n,k,w);
    lamk2 = lambda(n,n-k,w);
    sk = (sin(k*pi/(2*n)))^2;
    ck = (cos(k*pi/(2*n)))^2;
    y1(ind1) = lamk^(v1+v2)*sk;
    y4(ind2) = lamk2^(v1+v2)*ck;
    y3(ind1) = lamk2^v1*lamk^v2*ck;
    y2(ind2) = lamk^v1*lamk2^v2*sk;
    T = [y1(ind1),y2(ind2);y3(ind1),y4(ind2)];
    rho(ind1) = max(abs(eig(T)));
end
disp([min(rho),max(rho)]);

plot(lx1,y1,'k')
hold on
plot(lx2,y2,'-.k')
plot(lx1,y3,':k')
plot(lx2,y4,'k','LineWidth',2)
legend('c1','c2','c3','c4')

function lam = lambda(n,k,w)
    lam = 1-2*w*(sin(k*pi/(2*n)))^2;
end