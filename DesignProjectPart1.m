%Solves for the steady state concentrations along the length a reactor.
%Uses multivariate newtons to create a linear system of equations.
%The algebraic system of linear equations is solved using LU decomposition.

clear all;close all;clc;

%Constants used by the functions defined in f
%The design was tested for three sets of constants, represented by each row
constants = [[0 0 0.1 1.5 0.05 0 0 0];...
    [0.05 0 0.1 1.5 0.02 0.1 0.05 0.1];...
    [0.05 0.03 0.1 1.5 0.02 0.1 0.05 0.1]];
constant = constants(3,:);

%Defines the non-second-order terms of each dimensionless concentration ODE
f = {@(yA,yB,yU)-(yA.*(yB.^2))-constant(5)*yA, ...
    @(yA,yB,yU)-2*(yA.*(yB.^2))-(constant(1)*yB)+(constant(2)*yU), ...
    @(yA,yB,yU)constant(1)*yB-constant(2)*yU};

%Defines the derivatives of each function in f with recpect to each variable
%Each column represents each respective variable(yA,yB,yU)
df = {@(yA,yB,yU)-(yB.^2)-constant(5), @(yA,yB,yU)-2*(yA.*yB), @(yA,yB,yU)0;...
    @(yA,yB,yU)-2*(yB.^2), @(yA,yB,yU)-4*(yA.*yB)-constant(1), @(yA,yB,yU)constant(2);...
    @(yA,yB,yU)0, @(yA,yB,yU)constant(1), @(yA,yB,yU)-constant(2)};

%x defines the one dimensional length of the reactor
x = [0 1];

%N defines the number of points along the length of x to calculate the
%concentration for each profile of yA, yB, yC
N = 20;
initialGuess = 0;
dx = (x(2)-x(1))/(N-1);
odeCount = length(f);
unknownVariableCount = N*odeCount;
%Function of the finite differentiated second-order ODE;
F = @(y1,y2,y3,func,yA,yB,yU)(dx^-2)*constant(3)*(y3-(2*y2)+y1)+func(yA,yB,yU);


%Initialization of vectors needed for newtons method
Fv = zeros(unknownVariableCount,1);
J = zeros(unknownVariableCount,unknownVariableCount);
y = initialGuess*ones(N,3);

%Constructs the constant portion of the Jacobian as J
for i = 1:unknownVariableCount
    if(i <= odeCount)
        if(i ~= 3)
            J(i,i+odeCount) = 2*constant(3)*(dx^-2);
            J(i,i) = -2*constant(3)*((1+dx)*(dx^-2));
        else
            J(i,i) = 1;
        end
    elseif(i > unknownVariableCount-odeCount)
        J(i,i-odeCount) = 2*constant(3)*(dx^-2);
        if(i == unknownVariableCount-odeCount+1)
            J(i,i) = -2*constant(3)*((1+constant(6)*dx)*(dx^-2));
        elseif(i == unknownVariableCount-odeCount+3)
            J(i,i) = -2*constant(3)*((1+constant(8)*dx)*(dx^-2));
        end
    else
        J(i,i-odeCount) = (dx^-2)*constant(3);
        J(i,i+odeCount) = (dx^-2)*constant(3);
        J(i,i) = -2*(dx^-2)*constant(3);
    end
end

difference = Inf;
tolerance = 1.5e-7*ones(unknownVariableCount,1);
%Continues creating linear system of equations until convergene is reached
while(any(difference > tolerance))
    Jacobian = J;
    %Formation of the Jacobian matrix based on the current guess
    for i = 1:odeCount:unknownVariableCount-odeCount+1
        nodeIndex = ceil(i/odeCount);
        currentRow = num2cell(y(nodeIndex,:));
        if(i == 1)
            Fv(1) = constant(3)*(dx^-2)*(2*y(2,1)-2*(1+dx)*y(1,1)+2*dx)-(y(1,1).*(y(1,2).^2))-constant(5)*y(1,1);
            Jacobian(1,1) = Jacobian(1,1) + df{1,1}(currentRow{:});
            Jacobian(1,2) = Jacobian(1,2) + df{1,2}(currentRow{:});
            Jacobian(1,3) = Jacobian(1,3) + df{1,3}(currentRow{:});
            Fv(2) = constant(3)*(dx^-2)*(2*y(2,2)-2*(1+dx)*y(1,2)+2*dx*constant(4))-2*(y(1,1).*(y(1,2).^2))-constant(1)*y(1,2)+constant(2)*y(1,3);
            Jacobian(2,2) = Jacobian(2,2) + df{2,2}(currentRow{:});
            Jacobian(2,1) = Jacobian(2,1) + df{2,1}(currentRow{:});
            Jacobian(2,3) = Jacobian(2,3) + df{2,3}(currentRow{:});
            Fv(3) = 0;
        else
            Jacobian(i,i) = Jacobian(i,i) + df{1,1}(currentRow{:});
            Jacobian(i,i+1) = Jacobian(i,i+1) + df{1,2}(currentRow{:});
            Jacobian(i,i+2) = Jacobian(i,i+2) + df{1,3}(currentRow{:});
            Jacobian(i+1,i+1) = Jacobian(i+1,i+1) + df{2,2}(currentRow{:});
            Jacobian(i+1,i+1-1) = Jacobian(i+1,i+1-1) + df{2,1}(currentRow{:});
            Jacobian(i+1,i+1+1) = Jacobian(i+1,i+1+1) + df{2,3}(currentRow{:});
            Jacobian(i+2,i+2) = Jacobian(i+2,i+2) + df{3,3}(currentRow{:});
            Jacobian(i+2,i+2-1) = Jacobian(i+2,i+2-1) + df{3,2}(currentRow{:});
            Jacobian(i+2,i+2-2) = Jacobian(i+2,i+2-2) + df{3,1}(currentRow{:});
            if(i == unknownVariableCount-odeCount+1)
                Fv(i) = constant(3)*(dx^-2)*(2*y(nodeIndex-1,1)-2*(1+constant(6)*dx)*y(nodeIndex,1))-(y(nodeIndex,1).*(y(nodeIndex,2).^2))-constant(5)*y(nodeIndex,1);
                Fv(i+1) = constant(3)*(dx^-2)*(2*y(nodeIndex-1,2)-2*(1+constant(7)*dx*y(nodeIndex,2))*y(nodeIndex,2))-2*(y(nodeIndex,1).*(y(nodeIndex,2).^2))-constant(1)*y(nodeIndex,2)+constant(2)*y(nodeIndex,3);
                Jacobian(i+1,i+1) = Jacobian(i+1,i+1) + constant(3)*(dx^-2)*(-2*(1+2*constant(7)*dx*y(nodeIndex,2)));
                Fv(i+2) = constant(3)*(dx^-2)*(2*y(nodeIndex-1,3)-2*(1+constant(8)*dx)*y(nodeIndex,3))+constant(1)*y(nodeIndex,2)-constant(2)*y(nodeIndex,3);
            else
                Fv(i) = F(y(nodeIndex-1,1),y(nodeIndex,1),y(nodeIndex+1,1),f{1},currentRow{:});
                Fv(i+1) = F(y(nodeIndex-1,2),y(nodeIndex,2),y(nodeIndex+1,2),f{2},currentRow{:});
                Fv(i+2) = F(y(nodeIndex-1,3),y(nodeIndex,3),y(nodeIndex+1,3),f{3},currentRow{:});
            end
        end
    end
    b = Jacobian*reshape(y.',[],1) - Fv;
    tempY = y;
    %Solves algebraic system
    y = vectorToMatrix(LUSolver(Jacobian,b), odeCount);
    %Calculates difference between solutions, required for convergene test
    difference = abs(y(:)-tempY(:));
end
save('steadystatevalues.mat','tempY');

xAxis = x(1,1):(x(1,2)-x(1,1))/(N-1):x(1,2);
plot(xAxis, y(:,1), 'r.-', xAxis, y(:,2), 'b.-', xAxis, y(:,3), 'g.-');
hold on;
title('Position-Dependent Dimensionless Concentration Profiles at Steady State');
legend('y_A', 'y_B', 'y_U'); ylabel('y'); xlabel('x');
saveas(gcf,'DesignProjectPart1.png');