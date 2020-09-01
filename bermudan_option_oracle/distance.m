function D = distance(A, B)

% DISTANCE Calculates the distance between two matrix

% Description:
%   DISTANCE calculates the distance matrix D between matrix A and matrix B
%   D = DISTANCE(A, B) returns the distance matrix D.

d = size(A,2);

D = zeros(size(A,1),size(B,1));

for i = 1:size(A,1)
    for j = 1:size(B,1)
        d(i,j) = sum((A(i,:) - B(j,:)).^2);
    end
end
end