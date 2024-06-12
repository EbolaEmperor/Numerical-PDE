[theta,omega] = meshgrid(-3*pi:pi/5:3*pi,-3:0.2:3);
u = omega;
v = -sin(theta);
figure
quiver(theta,omega,u,v,'Color', '#FF9900')

hold on

% plotting solutions on the vector field
f = @(t,U) [U(2); -sin(U(1))];

for thetaValues = [-3*pi -2*pi 0 2*pi 3*pi]
    for omegaValues = [0 0.25 0.5 0.75 1 1.25 1.5 1.75 2 2.25 -0.25 -0.5 -0.75 -1 -1.25 -1.5 -1.75 -2]
        [ts,ys] = ode45(f,[0,15],[thetaValues;omegaValues]);
        plot(ys(:,1),ys(:,2),'Color', '#FFCC66' ,'LineWidth',1);
    end
end

A = load("result1.txt");
n = size(A, 2);
m = size(A, 1);
t = 1;
cnt = 0;
for i = 1:t:m
    c = 0.1 + 0.1 * cnt;
    cnt = cnt + 1;
    plot(A(i, 1:2:n), A(i, 2:2:n), 'k', 'linewidth', 1)
    fill(A(i, 1:2:n), A(i, 2:2:n), [c c c])
    hold on
end

A = load("result2.txt");
n = size(A, 2);
m = size(A, 1);
cnt = 0;
for i = 1:t:m
    c = 0.5 + 0.1 * cnt;
    cnt = cnt + 1;
    plot(A(i, 1:2:n), A(i, 2:2:n), 'k', 'linewidth', 1)
    fill(A(i, 1:2:n), A(i, 2:2:n), [c c c])
    hold on
end

% settings for the picture
axis equal
set(gca,'xlim',[-1,2.8.*pi]);
set(gca,'ylim',[-2,pi]);

set(gca,'LooseInset',get(gca,'TightInset'))
set(gca,'looseInset',[0 0 0 0])