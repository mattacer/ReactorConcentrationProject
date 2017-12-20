%Solves for the time-dependent position-dependent concentration profiles of
%species in a reactor.

clear all;close all;clc;
constant = [0.05 0.03 0.1 1.5 0.02 0.1 0.05 0.1];
N = 20;
x = [0 1];
dx = (x(2)-x(1))/(N-1);
xAxis = x(1,1):dx:x(1,2);
initialValues = 0;
load('steadystatevalues.mat');

writerObj = VideoWriter('ConcentrationsOverTime.mp4','MPEG-4');
writerObj.FrameRate = 30;
writerObj.Quality = 100;
open(writerObj);

%Discretized ODE's in time
df = {@(y1,yA,y3,yB,yU)constant(3)*(dx^-2)*(y1-2*yA+y3)-(yA.*(yB.^2))-constant(5)*yA ...
    @(y1,yB,y3,yA,yU)constant(3)*(dx^-2)*(y1-2*yB+y3)-2*(yA.*(yB.^2))-constant(1)*yB+constant(2)*yU ...
    @(y1,yU,y3,yA,yB)constant(3)*(dx^-2)*(y1-2*yU+y3)+constant(1)*yB-constant(2)*yU};

odeCount = length(df);
unknownCount = N*odeCount;
DeriVals = zeros(unknownCount,1);
y = zeros(N,3);
t = 0;
h = 1e-2;
while (any(~any(y>=0.99*tempY)'))
    for i = 1:odeCount:unknownCount-odeCount+1
        nodeIndex = ceil(i/odeCount);
        %Calculates the derivative at the next time step
        if(i == 1)
            DeriVals(i) = constant(3)*(dx^-2)*(2*y(nodeIndex+1,1)-2*(1+dx)*y(nodeIndex,1)+2*dx)-(y(nodeIndex,1).*(y(nodeIndex,2).^2))-constant(5)*y(nodeIndex,1);
            DeriVals(i+1) = constant(3)*(dx^-2)*(2*y(nodeIndex+1,2)-2*(1+dx)*y(nodeIndex,2)+2*dx*constant(4))-2*(y(nodeIndex,1).*(y(nodeIndex,2).^2))-constant(1)*y(nodeIndex,2)+constant(2)*y(nodeIndex,3);
            DeriVals(i+2) = 0;
        elseif(i == unknownCount-odeCount+1)
            DeriVals(i) = constant(3)*(dx^-2)*(2*y(nodeIndex-1,1)-2*(1+constant(6)*dx)*y(nodeIndex,1))-(y(nodeIndex,1).*(y(nodeIndex,2).^2))-constant(5)*y(nodeIndex,1);
            DeriVals(i+1) = constant(3)*(dx^-2)*(2*y(nodeIndex-1,2)-2*(1+constant(7)*dx*y(nodeIndex,2))*y(nodeIndex,2))-2*(y(nodeIndex,1).*(y(nodeIndex,2).^2))-constant(1)*y(nodeIndex,2)+constant(2)*y(nodeIndex,2);
            DeriVals(i+2) = constant(3)*(dx^-2)*(2*y(nodeIndex-1,3)-2*(1+constant(8)*dx)*y(nodeIndex,3))+constant(1)*y(nodeIndex,2)-constant(2)*y(nodeIndex,3);
        else
            DeriVals(i) = df{1}(y(nodeIndex-1,1),y(nodeIndex,1),y(nodeIndex+1,1),y(nodeIndex,2),y(nodeIndex,3));
            DeriVals(i+1) = df{2}(y(nodeIndex-1,2),y(nodeIndex,2),y(nodeIndex+1,2),y(nodeIndex,1),y(nodeIndex,3));
            DeriVals(i+2) = df{3}(y(nodeIndex-1,3),y(nodeIndex,3),y(nodeIndex+1,3),y(nodeIndex,1),y(nodeIndex,2));
        end
    end
    %Adds the derivative values multipled by h to the initial values
    y = vectorToMatrix(reshape(y.',[],1) + h*DeriVals,odeCount);
    %Plots the profiles at give times
    if(abs(t-0.1) <= h/2 || abs(t-1) <= h/2 || abs(t-5) <= h/2)
        figure(1); hold on;
        plot(xAxis, y(:,1), 'r.-');
        figure(2); hold on;
        plot(xAxis, y(:,2), 'b.-');
        figure(3); hold on;
        plot(xAxis, y(:,3), 'g.-');
    end
    figure(4); cla; hold on;
    title('Position-Dependent Dimensionless Concentration Profiles Over Time');
    xlabel(sprintf('Position along membrane reactor at t=%3.2f', t));
    ylabel('y');
    plot(xAxis, y(:,1), 'r.-', xAxis, y(:,2), 'b.-', xAxis, y(:,3), 'g.-');
    legend('y_A', 'y_B', 'y_U');
    writeVideo(writerObj, getframe(gcf));
    t = t+h
end
%Plots the steady-state
figure(1);
plot(xAxis, y(:,1), 'r.-'); ylabel('y_A');
title('Position-Dependent Dimensionless Concentration Profile y_A Over Time');
xlabel(sprintf('x\nSteady-state reached at t=%3.2f',t));
figure(2);
plot(xAxis, y(:,2), 'b.-'); ylabel('y_B');
title('Position-Dependent Dimensionless Concentration Profile y_B Over Time');
xlabel(sprintf('x\nSteady-state reached at t=%3.2f',t));
figure(3);
plot(xAxis, y(:,3), 'g.-'); ylabel('y_U');
title('Position-Dependent Dimensionless Concentration Profile y_U Over Time');
xlabel(sprintf('x\nSteady-state reached at t=%3.2f',t));
close(writerObj);
