
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                   %
%  Equation of motion for simple pendulum with g/l = 1:              %
%    \theta''(t) + \sin(\theta) = 0,                                 %
%  which be reduced to a first-order system:                        %
%    (\theta, \omega)' = (\omega, -\sin(\theta)).                    %
%  Solve this system and plot the vector filed (\theta,\omega)       %
%                                                                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% clear all unwanted variables and graphs

clear all; close all; clc

% plotting the phase portraits

[theta,omega] = meshgrid(-3*pi:pi/5:3*pi,-3:0.2:3);
u = omega;
v = -sin(theta);
figure
quiver(theta,omega,u,v,'r')

hold on

% plotting solutions on the vector field
f = @(t,U) [U(2); -sin(U(1))];

for thetaValues = [-3*pi -2*pi 0 2*pi 3*pi]
    for omegaValues = [0 0.25 0.5 0.75 1 1.25 1.5 1.75 2 2.25 -0.25 -0.5 -0.75 -1 -1.25 -1.5 -1.75 -2]
        [ts,ys] = ode45(f,[0,15],[thetaValues;omegaValues]);
        plot(ys(:,1),ys(:,2),'black','LineWidth',1);
    end
end


% settings for the picture
axis equal
set(gca,'xlim',[-.5.*pi,1.8.*pi]);
set(gca,'ylim',[-3,3]);

set(gca,'LooseInset',get(gca,'TightInset'))
set(gca,'looseInset',[0 0 0 0])



