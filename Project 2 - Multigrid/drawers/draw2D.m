u = @(x,y)exp(sin(x)+y);
V = load('result.txt');
realV = u(V(:,1),V(:,2));
err = V(:,3) - realV;

n = sqrt(size(V,1));
X = reshape(V(:,1),n,n);
Y = reshape(V(:,2),n,n);
Z = reshape(err,n,n);
surf(X, Y, Z);

function E = error(A)
    r = max(A);
    l = min(A);
    mid = (l+r)/2;
    n = size(A,1);
    E = A - mid*ones(n,1);
end