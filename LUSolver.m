function x = LUSolver(A, b)
%Solves for the x vector of a system of equations with the form Ax=b using
%the LU method.
%A must be an N by N matrix.
%b must be an N length column vector.
%The returned value of the x is an N length column vector.

N = size(A);
%Conditional to catch over/under-determined systems
if(N(1) ~= N(2))
    error('Matrix A must be an N by N matrix.');
end
N = N(1);
%Conditional to catch invalid b vectors
if(size(b) ~= [N 1])
    error('Vector b must be a column vector and have the same number of rows as matrix A');
end

x = zeros(N,1);
y = zeros(N,1);
L = A;
U = eye(N,N);
%Pivots rows
Labs = abs(L);
[val m]=max(Labs(1:N,1)); %Find max value from 1 through n to pivot to
L([1 m],:)=L([m 1],:); %Swap the 1st row with the max row
b([1 m])=b([m 1]);
%For diagonal elements, k, eliminate in the column below
for i=1:N-1
    %Forward elimination
    for j=i+1:N
        U(i,j)=L(i,j)/L(i,i); %Find the elimination factor
        L(:,j)=L(:,j)-L(:,i)*U(i,j); %row operation
    end
end
%Forward substitution for y
for i = 1:N
    y(i) = (b(i) - L(i,:)*y)/L(i,i);
end
%Backward substitution for x
for i = N:-1:1
    x(i) = (y(i) - U(i,:)*x)/U(i,i);
end
end