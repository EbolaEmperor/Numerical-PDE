B = [ 1, 0, 0, 0; -1/2, 3/2, 0, 0; 5/12, -4/3, 23/12, 0; -3/8, 37/24, -59/24, 55/24 ];
% When producing the first figure, replace the 4:4 to 1:3.
for s = 4:4
    b = B(s,:);
    rho = @(x)(x.^s-x.^(s-1));
    sigma = @(x)(b(1)+b(2)*x+b(3)*(x.^2)+b(4)*(x.^3));
    X = 0:0.01:2*pi;
    Y = rho(exp(1i*X))./sigma(exp(1i*X));
    plot(real(Y),imag(Y));
    xlim([-1,1]);
    ylim([-1,1]);
    hold on
end
% When producing the first figure, the codes below is not necessary.
b = B(4,:);
HX = zeros(201*201, 2);
cnt = 0;
for X = -1:0.01:1
    for Y = -1:0.01:1
        k = X + Y*1i;
        rho = [1, -1, 0, 0, 0];
        sigma = [0, b(4), b(3), b(2), b(1)];
        p = rho - k * sigma;
        if max(abs(roots(p)))<=1
            cnt = cnt + 1;
            HX(cnt,:) = [X, Y];
        end
    end
end
scatter(HX(:,1), HX(:,2), 2, "filled", 'k');