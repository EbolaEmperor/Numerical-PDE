u = @(x,y)exp(sin(x));
V = load('result.txt');
realV = u(V(:,1));
err = V(:,2) - realV;

plot(V(:,1), err);

function E = error(A)
    r = max(A);
    l = min(A);
    mid = (l+r)/2;
    n = size(A,1);
    E = A - mid*ones(n,1);
end