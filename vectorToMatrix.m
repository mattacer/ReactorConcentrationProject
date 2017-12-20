function matrix = vectorToMatrix(vector, N)
%VECTORTOMATRIX Transforms a M by 1 vector into an M/N by N matrix
%   Every N values in the vector are transformed into a row in output
%   matrix.
if(mod(length(vector), N))
    error('N must be a factor of the length of vector.');
end
rows = length(vector)/N;
matrix = zeros(rows, N);
% tic
for i = 1:N:length(vector)
    matrix(ceil(i/N), :) = vector(i:(i+N-1))';
end
%The code below performs slower due to additional operations
% for i = 1:length(vector)
%      matrix(ceil(i/N), i - matrixLength*(ceil(i/N)-1)) = vector(i);
% end
% toc
end

