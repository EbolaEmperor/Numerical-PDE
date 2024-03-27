A = [ 1, -1, 1, 0, 0, 0; 1, -4/3, 1/3, 2/3, 0, 0; 1, -18/11, 9/11, -2/11, 6/11, 0; 1, -48/25, 36/25, -16/25, 3/25, 12/25 ];
for s = 1:4
    b = A(s,:);
    sigma = @(x)(b(2+s)*x.^s);
    X = 0:0.01:2*pi;
    X = exp(1i*X);
    rho = b(1)*(X.^s);
    for i = 2:s+1
        rho = rho + b(i)*(X.^(s+1-i));
    end
    Y = rho./sigma(X);
    plot(real(Y),imag(Y));
    xlim([-4,13]);
    ylim([-7,7]);
    hold on
end