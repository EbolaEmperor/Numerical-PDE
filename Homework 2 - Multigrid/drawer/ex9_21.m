n = 64;
w = 2/3;
T = eye(n-1)*(1-w);
T(1:n-2,2:n-1) = T(1:n-2,2:n-1) + eye(n-2)*(w/2);
T(2:n-1,1:n-2) = T(2:n-1,1:n-2) + eye(n-2)*(w/2);
x = 1:n-1;
y = x;
for k = 1:n-1
    w = mode(n,k);
    wk = w;
    iter = 0;
    for i = 1:100
        w = T*w;
        iter = iter+1;
        if norm(w)/norm(wk)<0.01
            break;
        end
    end
    y(k) = iter;
end
plot(x,y,'k','LineWidth',2);