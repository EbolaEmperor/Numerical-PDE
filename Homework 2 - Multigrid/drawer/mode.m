function w = mode(n,k)
    w = zeros(n-1,1);
    for j = 1:n-1
        w(j,1) = sin(j*k*pi/n);
    end
end

