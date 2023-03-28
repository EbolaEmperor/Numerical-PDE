u = @(x,y)exp(sin(x)+y);
V = load('result.txt');
realV = u(V(:,1),V(:,2));
err = V(:,3) - realV;

tri = delaunay(V(:,1),V(:,2));
trisurf(tri,V(:,1),V(:,2),realV);

function E = error(A)
    r = max(A);
    l = min(A);
    mid = (l+r)/2;
    n = size(A,1);
    E = A - mid*ones(n,1);
end